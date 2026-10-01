"""Descriptive per-channel statistics from 90 minutes onward; corrected CSV input."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT/'20260916_1um_MBP23_30rpm_CH1_CH3_CH4_corrected.csv'
OUT = ROOT/'post_equilibration_90min'
CHANNELS = (1,3,4)

def describe(values):
    s=pd.Series(values,dtype=float).dropna()
    mean=s.mean(); sd=s.std(ddof=1)
    return dict(n=len(s),mean=mean,sd=sd,cv_percent=100*sd/abs(mean) if mean else np.nan,
                median=s.median(),min=s.min(),q25=s.quantile(.25),q75=s.quantile(.75),max=s.max())

def extract_cycles(f,ch):
    t=f.elapsed_min.to_numpy()*60
    q=f[f's{ch}_flow_ml_min'].to_numpy()
    # Hysteresis prevents transition noise/zeros from generating tiny false cycles.
    threshold=.2*np.median(np.abs(q[np.isfinite(q)&(q!=0)]))
    state=0; starts=[]
    for i,value in enumerate(q):
        if value>=threshold:
            if state==-1: starts.append(i)
            state=1
        elif value<=-threshold: state=-1
    rows=[]
    for a,b in zip(starts[:-1],starts[1:]):
        if t[a]<5400: continue
        v=q[a:b]; pos=v[v>=threshold]; neg=v[v<=-threshold]
        if len(pos)<3 or len(neg)<3: continue
        qp=pos[pos>=np.median(pos)].mean()
        qm=neg[neg<=np.median(neg)].mean()
        upper=v[v>=np.quantile(v,.9)].mean(); lower=v[v<=np.quantile(v,.1)].mean()
        rows.append(dict(channel=f'CH{ch}',cycle=len(rows)+1,start_min=t[a]/60,end_min=t[b]/60,
                         period_s=t[b]-t[a],Q_plus=qp,Q_minus=qm,abs_Q_minus=abs(qm),
                         P_to_P=qp-qm,symmetry=qp/abs(qm),
                         upper_decile_Q_plus=upper,lower_decile_Q_minus=lower,
                         decile_P_to_P=upper-lower,decile_symmetry=upper/abs(lower),
                         raw_peak_to_peak=v.max()-v.min(),cycle_mean_flow=v.mean(),
                         push_volume_ml=np.trapezoid(np.maximum(q[a:b+1],0),t[a:b+1])/60,
                         pull_volume_ml=-np.trapezoid(np.minimum(q[a:b+1],0),t[a:b+1])/60,
                         temperature_C=f[f's{ch}_temp_c'].iloc[a:b].mean()))
    result=pd.DataFrame(rows)
    result['net_volume_ml']=result.push_volume_ml-result.pull_volume_ml
    result['volume_symmetry']=result.push_volume_ml/result.pull_volume_ml
    return result,threshold

def main():
    OUT.mkdir(exist_ok=True)
    f=pd.read_csv(SOURCE)
    assert f.time_s.is_monotonic_increasing and not f.time_s.duplicated().any()
    assert f.notna().all().all(), 'Missing data: review continuity before cycle extraction'
    f['elapsed_min']=(f.time_s-f.time_s.iloc[0])/60
    post=f.loc[f.elapsed_min>=90].copy()
    all_cycles=[]; stats=[]; thresholds={}; compact=[]; bins=[]
    for ch in CHANNELS:
        label=f'CH{ch}'
        cycles,threshold=extract_cycles(f,ch)
        assert (cycles.start_min>=90).all() and (cycles.Q_plus>0).all() and (cycles.Q_minus<0).all()
        all_cycles.append(cycles); thresholds[label]=threshold
        temp=describe(post[f's{ch}_temp_c'])
        stats.append(dict(channel=label,metric='temperature_C',unit='°C',**temp))
        stats.append(dict(channel=label,metric='signed_flow_samples',unit='mL/min',**describe(post[f's{ch}_flow_ml_min'])))
        stats.append(dict(channel=label,metric='absolute_flow_samples',unit='mL/min',**describe(post[f's{ch}_flow_ml_min'].abs())))
        for metric in cycles.columns[4:]:
            unit='s' if metric=='period_s' else ('°C' if metric=='temperature_C' else ('mL' if 'volume_ml' in metric else ('ratio' if 'symmetry' in metric else 'mL/min')))
            stats.append(dict(channel=label,metric='cycle_'+metric,unit=unit,**describe(cycles[metric])))
        row=dict(channel=label,n_samples=len(post),n_cycles=len(cycles),temperature_mean_C=temp['mean'],temperature_sd_C=temp['sd'],temperature_cv_percent=temp['cv_percent'])
        for metric in ['Q_plus','Q_minus','P_to_P','symmetry','period_s']:
            for key in ['mean','sd','cv_percent']: row[f'{metric}_{key}']=describe(cycles[metric])[key]
        compact.append(row)
        for start in np.arange(90,post.elapsed_min.max(),30):
            sf=post.loc[post.elapsed_min.between(start,start+30,inclusive='left')]
            sc=cycles.loc[cycles.start_min.between(start,start+30,inclusive='left')]
            for metric,values in [('temperature_C',sf[f's{ch}_temp_c'])]+[(m,sc[m]) for m in ['Q_plus','Q_minus','P_to_P','symmetry']]:
                bins.append(dict(channel=label,start_min=start,end_min=min(start+30,post.elapsed_min.max()),metric=metric,**describe(values)))
    cycles=pd.concat(all_cycles,ignore_index=True)
    summary=pd.DataFrame(compact); full=pd.DataFrame(stats); drift=pd.DataFrame(bins)
    inter=[]
    for metric in ['temperature_mean_C','Q_plus_mean','Q_minus_mean','P_to_P_mean','symmetry_mean']:
        inter.append(dict(metric=metric,**describe(summary[metric])))
    inter=pd.DataFrame(inter)
    differences=[]
    for a,b in [(1,3),(1,4),(3,4)]:
        differences.append(dict(comparison=f'CH{a} − CH{b}',unit='°C',**describe(post[f's{a}_temp_c']-post[f's{b}_temp_c'])))
    differences=pd.DataFrame(differences)
    tables={'channel_summary':summary,'all_statistics':full,'cycle_descriptors':cycles,
            'inter_channel':inter,'temperature_differences':differences,'30min_windows':drift}
    with pd.ExcelWriter(OUT/'post_equilibration_statistics.xlsx') as writer:
        for name,data in tables.items():
            data.to_csv(OUT/f'{name}.csv',index=False)
            data.to_excel(writer,sheet_name=name,index=False)
    notes=f'''Post-equilibration descriptive statistics
Source: {SOURCE.name}
Window: elapsed time >=90 min, relative to first recorded timestamp.
Actual selected range: {post.elapsed_min.min():.6f}–{post.elapsed_min.max():.6f} min.
Samples per channel: {len(post)}. Channels: CH1, CH3, CH4; CH2 is absent.
Data are already corrected; no additional correction factors applied.
No missing entries; timestamps strictly increasing. Maximum sampling interval: {f.time_s.diff().max():.3f} s.
90 minutes is a user-selected cutoff, not a fitted or verified equilibration time.

Cycle detection: negative-to-positive transitions with hysteresis at 20% of
median nonzero absolute flow (thresholds in mL/min: {thresholds}).
Each cycle extends from one positive threshold crossing to the next.
Only fully bounded cycles starting >=90 min are included. Leading/trailing
partial cycles are excluded. At least 3 samples per polarity are required.
Q+ = mean of upper half (>=median) of positive phase samples above threshold.
Q− = mean of lower half (<=median) of negative phase samples below -threshold.
P–P = Q+ − Q−; symmetry = Q+/|Q−|, calculated per cycle before summarizing.
This explicit plateau estimator is comparable in concept to the S2 upper/lower
half method but uses hysteresis to avoid transition noise.
Upper/lower decile estimates are also supplied, matching the decile concept
in the existing plotting script. Its imported helper is absent, so exact
reproduction of its unpublished cycle algorithm is not claimed.
Raw peak-to-peak = per-cycle max minus min; distinct from plateau P–P.
Volumes integrate clipped positive/negative flow over each cycle with the
trapezoidal rule; pull volume is a positive magnitude.

SD is sample SD (ddof=1); CV=100*SD/abs(mean). Temperature CV uses Celsius,
matching the prior figure; this CV depends on the temperature scale.
Temperature/sample-flow statistics use every selected sample, including
transitions. Flow waveform SD/CV is distinct from cycle-to-cycle plateau CV;
signed-flow CV can be large because its mean is near zero.
Inter-channel statistics use the three channel means (n=3), not pooled samples.
Temperature differences use paired simultaneous samples.
30-minute windows assign cycles by cycle start, with a shorter final window.
No inter-day CV is available from this single recording. Repeated time samples
and cycles are not independent biological/experimental replicates; no naive
confidence intervals or significance tests are reported.

Channel summary:
{summary.to_string(index=False,float_format=lambda v:f'{v:.4f}')}

Temperature differences:
{differences.to_string(index=False,float_format=lambda v:f'{v:.4f}')}
'''
    (OUT/'analysis_notes.txt').write_text(notes)
    print(summary.to_string(index=False)); print('\nInter-channel:\n',inter.to_string(index=False)); print('\nTemperature differences:\n',differences.to_string(index=False))
    print('Saved',OUT)

if __name__=='__main__': main()
