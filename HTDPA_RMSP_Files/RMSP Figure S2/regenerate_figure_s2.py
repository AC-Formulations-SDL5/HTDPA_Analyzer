"""Rebuild Figure S2 using supplied cycle descriptors and recorded temperatures."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DAYS = ['HTDPA_3h_140726', 'HTDPA_3h_150726']
COLORS = ['#78b5d2', '#ee9b51', '#87bd77', '#a397c5']
METRICS = ['Q_plus_ml_min', 'Q_minus_ml_min', 'P_to_P_ml_min', 'symmetry_abs_Qp_over_abs_Qm']
LABELS = [r'$Q^+$', r'$|Q^-|$', 'P–P', r'$Q^+/|Q^-|$']
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9, 'axes.titlesize':10,
                     'axes.labelsize':9, 'pdf.fonttype':42, 'svg.fonttype':'none'})

def cv(s):
    return s.std(ddof=1) / abs(s.mean()) * 100

def style(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='y', color='#e5e7eb', linewidth=.6)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(length=3)

def panel(ax, letter):
    ax.text(-.13, 1.12, letter, transform=ax.transAxes, fontsize=16, weight='bold')

def main():
    frames = {d:pd.read_excel(ROOT/'HTDPA_RMSP_140726_normalized.xlsx', sheet_name=d) for d in DAYS}
    cycles = {d:pd.read_csv(ROOT/f'waveform_descriptors_by_cycle_{d}.csv') for d in DAYS}
    for data in cycles.values():
        data['Q_minus_ml_min'] = data['Q_minus_ml_min'].abs()
    fig = plt.figure(figsize=(15,9), layout='constrained')
    outer = fig.add_gridspec(2,2, height_ratios=[1.65,1], width_ratios=[1.7,1], hspace=.12,wspace=.1)
    left = outer[0,0].subgridspec(2,1,hspace=.22)
    ag = left[0].subgridspec(1,3,wspace=.12)
    aa = [fig.add_subplot(ag[0,i]) for i in range(3)]
    x = np.arange(4)
    ticks = ['CH1','CH2','CH3','CH4']
    summary = cycles[DAYS[0]].groupby('channel')[METRICS].agg(['mean','std'])
    for i,m in enumerate(METRICS[:2]):
        aa[0].bar(x+(i-.5)*.36, summary[m]['mean'], .36, yerr=summary[m]['std'],
                  color=COLORS, hatch='///' if i else None, edgecolor='#435363',linewidth=.6,capsize=2)
    aa[0].legend(handles=[Patch(facecolor='#c2dbe8',label='Push',edgecolor='#435363'),
                          Patch(facecolor='#c2dbe8',label='Pull',hatch='///',edgecolor='#435363')],fontsize=8,frameon=False,loc='upper left',ncol=2)
    for ax,m in zip(aa[1:],METRICS[2:]):
        ax.bar(x,summary[m]['mean'],.65,yerr=summary[m]['std'],color=COLORS,edgecolor='#435363',linewidth=.6,capsize=2)
    aa[2].axhline(1,color='#444',ls='--',lw=1)
    for ax,title,ylabel in zip(aa,['Plateau flow','Peak-to-peak amplitude','Push–pull symmetry'],
                               ['Flow (mL/min)','P–P (mL/min)',LABELS[-1]]):
        ax.set(title=title,ylabel=ylabel,xticks=x,xticklabels=ticks)
        ax.set_ylim(0,ax.get_ylim()[1]*1.12)
        style(ax)
    panel(aa[0],'A')
    bg = left[1].subgridspec(1,2,wspace=.13)
    ba = [fig.add_subplot(bg[0,i]) for i in range(2)]
    within = cycles[DAYS[0]].groupby('channel')[METRICS].agg(cv)
    daily = pd.concat([cycles[d].groupby('channel')[METRICS].mean().assign(day=d) for d in DAYS]).reset_index()
    inter = daily.groupby('channel')[METRICS].agg(cv)
    for ax,data,title in zip(ba,[within,inter],['Within-day CV (140726)','Inter-day CV (140726 vs 150726)']):
        for i,m in enumerate(METRICS):
            ax.bar(x+(i-1.5)*.19,data[m],.19,color=COLORS[i],label=LABELS[i])
        ax.set(title=title,ylabel='CV (%)',xticks=x,xticklabels=ticks)
        style(ax)
    ymax=max(within.max().max(),inter.max().max())*1.35
    for ax in ba: ax.set_ylim(0,ymax)
    ba[0].legend(fontsize=7,ncol=2,frameon=False)
    panel(ba[0],'B')
    eg=outer[0,1].subgridspec(4,1,hspace=.05)
    ea=[]
    for i in range(4):
        ax=fig.add_subplot(eg[i]); ea.append(ax)
        f=frames[DAYS[1]]
        f=f.loc[f.time_s.between(45*60,46*60)]
        assert len(f)>0, 'Missing zoom window'
        ax.plot(f.time_s/60,f[f's{i+1}_flow_norm_ml_min'],color='#287ab0',lw=.7)
        tx=ax.twinx()
        tx.plot(f.time_s/60,f[f's{i+1}_temp_c'],color='#d84b51',lw=.85)
        ax.set_ylim(-60,60); tx.set_ylim(27,33)
        ax.set_xlim(45,46); ax.set_ylabel(f'CH{i+1}\nFlow (mL/min)',fontsize=8,color='#287ab0')
        tx.set_ylabel('Temp. (°C)',fontsize=8,color='#d84b51')
        ax.tick_params(labelsize=7); tx.tick_params(labelsize=7)
        if i<3: ax.tick_params(labelbottom=False)
        else: ax.set_xlabel('Time (min)')
    ea[0].set_title('SFL3X signal · 30 RPM · 150726',pad=12)
    panel(ea[0],'E')
    bottom=outer[1,:].subgridspec(1,2,wspace=.15)
    temp_rows=[]
    for j,d in enumerate(DAYS):
        ax=fig.add_subplot(bottom[j]); tx=ax.twinx()
        series=[frames[d][f's{i}_temp_c'].dropna() for i in range(1,5)]
        means=np.array([s.mean() for s in series]); sds=np.array([s.std(ddof=1) for s in series])
        cvs=np.array([cv(s) for s in series])
        mus=np.r_[means,means.mean()]; sd=np.r_[sds,means.std(ddof=1)]
        cvs=np.r_[cvs,means.std(ddof=1)/means.mean()*100]
        xx=np.arange(5)
        ax.bar(xx-.18,mus,.36,yerr=sd,color=COLORS+['#89c8ce'],capsize=3,edgecolor='#435363',linewidth=.6)
        tx.bar(xx+.18,cvs,.36,color='#f5d4dd',edgecolor='#ce5271',linewidth=.6)
        for k in range(5):
            tx.text(xx[k],max(mus[k]+sd[k],20+cvs[k]*4)+.4,f'{mus[k]:.2f} ± {sd[k]:.2f}',transform=ax.transData,ha='center',fontsize=8)
            temp_rows.append(dict(day=d,channel=ticks[k] if k<4 else 'Inter-channel',mean_C=mus[k],sd_C=sd[k],cv_percent=cvs[k]))
        ax.set(ylim=(20,40),xticks=xx,xticklabels=ticks+['Inter'],ylabel='Temperature (°C)',
               title=('Before' if j==0 else 'After')+' fluidic adjustment ('+d[-6:]+')')
        tx.set(ylim=(0,5),ylabel='CV (%)')
        ax.legend(handles=[Patch(facecolor='#89b6d2',label='Temperature: mean ± SD'),Patch(facecolor='#f5d4dd',edgecolor='#ce5271',label='CV')],frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(0,1.01),ncol=2)
        style(ax); panel(ax,'C' if j==0 else 'D')
    fig.savefig(ROOT/'Figure_S2_regenerated.png',dpi=400,facecolor='white')
    fig.savefig(ROOT/'Figure_S2_regenerated.pdf',facecolor='white')
    fig.savefig(ROOT/'Figure_S2_regenerated.svg',facecolor='white')
    fig.savefig(ROOT/'Figure_S2_preview.png',dpi=120,facecolor='white')
    summary.to_csv(ROOT/'Figure_S2_panel_A_summary.csv')
    within.to_csv(ROOT/'Figure_S2_panel_B_within_day_CV.csv')
    inter.to_csv(ROOT/'Figure_S2_panel_B_inter_day_CV.csv')
    pd.DataFrame(temp_rows).to_csv(ROOT/'Figure_S2_temperature_summary.csv',index=False)
    plt.close(fig)
    print('Saved Figure S2 (PNG, PDF, SVG), preview, and summary tables.')

if __name__=='__main__':
    main()
