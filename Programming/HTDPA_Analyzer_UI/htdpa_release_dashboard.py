"""HTD/PA release analysis dashboard - comparative 6-panel version.

Run with:
    streamlit run htdpa_release_dashboard.py

Key analysis rules
------------------
* Only non-negative cumulative-release values (>= 0 %) are retained.
* Quantitative analysis starts at the first positive-release time point where a strict majority (>50%) of replicates are at/above the assay LOQ.
* CV uses only quantifiable replicates (>=LOQ) within majority-valid time points.
* Linearity is assessed before Kp/model analysis.
* If the selected interval does not meet the target R^2, terminal time points
  are removed one at a time until the target is met or the minimum number of
  points remains. Internal points are never removed.
* The accepted common analysis window is then used for Zero, Higuchi and
  First-order fits. KP and Peppas-Sahlin additionally use 0 < F <= 60 % and
  may use a later user-selected start time.
* KP and Peppas-Sahlin use recorded time since the experimental release origin.
  KP is fitted nonlinearly on the original release-fraction scale.
* The six-panel output follows the manuscript-summary style:
  A) release profiles, B) median CV + IQR, C) linear R^2,
  D) model-fit heatmap, E) apparent Kp, F) KP exponent n.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import curve_fit

try:
    import streamlit as st
except ImportError:
    st = None


ASSET_DIR = Path(__file__).resolve().parent


@dataclass
class ReleaseData:
    time_min: np.ndarray
    replicates_pct: np.ndarray
    replicate_names: List[str]
    source_note: str
    detected_receptor_volumes_ml: Optional[List[float]] = None
    detected_donor_dose_ug: Optional[float] = None
    measured_concentration_ug_ml: Optional[np.ndarray] = None


@dataclass
class ConditionResult:
    label: str
    formulation: str
    probe: str
    rpm: float
    temperature_c: float
    donor_dose_ug: float
    donor_concentration_ug_ml: float
    membrane_area_cm2: float
    data: ReleaseData
    mean_pct: np.ndarray
    sd_pct: np.ndarray
    fit_mask: np.ndarray
    fit_start_min: float
    fit_end_min: float
    n_fit_points: int
    trimmed_terminal_points: int
    analysis_window_mode: str
    requested_start_min: float
    requested_end_min: float
    linear_r2: float
    slope_ug_h: float
    kp_cm_h: float
    median_cv_pct: float
    cv_q1_pct: float
    cv_q3_pct: float
    models: Dict[str, float]
    kp_exponent_n: float
    kp_start_min: float
    kp_end_min: float
    kp_n_points: int
    valid_timepoint_mask: np.ndarray
    quantifiable_counts: np.ndarray
    quantifiable_fraction: np.ndarray
    loq_ug_ml: float
    lod_ug_ml: float


def _r2(y: np.ndarray, yhat: np.ndarray) -> float:
    mask = np.isfinite(y) & np.isfinite(yhat)
    y = np.asarray(y)[mask]
    yhat = np.asarray(yhat)[mask]
    if y.size < 2:
        return np.nan
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    if ss_tot <= 0:
        return np.nan
    return float(1.0 - np.sum((y - yhat) ** 2) / ss_tot)


def _numeric_array(values) -> np.ndarray:
    return pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)


def cumulative_from_concentration(time_min, concentrations_ug_ml, receptor_volume_ml, sample_volume_ml, donor_dose_ug):
    out = np.full_like(concentrations_ug_ml, np.nan, dtype=float)
    for r in range(concentrations_ug_ml.shape[0]):
        c = concentrations_ug_ml[r]
        for j in range(c.size):
            if not np.isfinite(c[j]):
                continue
            previous = c[:j]
            previous = previous[np.isfinite(previous)]
            q_ug = receptor_volume_ml * c[j] + sample_volume_ml * np.sum(previous)
            out[r, j] = 100.0 * q_ug / donor_dose_ug
    return out


def clean_release_data(data: ReleaseData) -> ReleaseData:
    """Keep valid release values >= 0%; negative values become NaN."""
    reps = np.asarray(data.replicates_pct, dtype=float).copy()
    reps[~np.isfinite(reps)] = np.nan
    reps[reps < 0] = np.nan
    return ReleaseData(
        time_min=np.asarray(data.time_min, dtype=float),
        replicates_pct=reps,
        replicate_names=list(data.replicate_names),
        source_note=data.source_note + "; negative release excluded",
        detected_receptor_volumes_ml=data.detected_receptor_volumes_ml,
        detected_donor_dose_ug=data.detected_donor_dose_ug,
        measured_concentration_ug_ml=None if data.measured_concentration_ug_ml is None else np.asarray(data.measured_concentration_ug_ml, dtype=float).copy(),
    )


def parse_htdpa_cumulative_sheet(file_bytes: bytes) -> Optional[ReleaseData]:
    try:
        xls = pd.ExcelFile(io.BytesIO(file_bytes), engine="openpyxl")
    except Exception:
        return None
    if "Cumulative Release" not in xls.sheet_names:
        return None

    raw = pd.read_excel(io.BytesIO(file_bytes), sheet_name="Cumulative Release", header=None, engine="openpyxl")
    target_row = None
    for i in range(len(raw)):
        row_text = " ".join(str(v) for v in raw.iloc[i].dropna().tolist())
        if "Sampling-corrected cumulative release (%" in row_text:
            target_row = i
            break
    if target_row is None:
        return None

    time_row = target_row + 1
    time_vals = _numeric_array(raw.iloc[time_row, 2:].tolist())
    valid_t = np.isfinite(time_vals)
    if valid_t.sum() < 3:
        return None
    last = np.where(valid_t)[0][-1] + 1
    time_min = time_vals[:last]

    reps, names = [], []
    r = time_row + 1
    while r < len(raw):
        label = raw.iloc[r, 1] if raw.shape[1] > 1 else None
        label_s = "" if pd.isna(label) else str(label).strip()
        if label_s.lower() in {"average (%)", "s.d. (%)", "c.v. (%)"}:
            break
        vals = _numeric_array(raw.iloc[r, 2:2 + last].tolist())
        if np.isfinite(vals).sum() >= 3:
            reps.append(vals)
            names.append(label_s or f"Rep {len(reps)+1}")
        elif reps:
            break
        r += 1
    if not reps:
        return None

    receptor_vols, donor_dose = [], None
    for i in range(min(15, len(raw))):
        left = "" if pd.isna(raw.iloc[i, 0]) else str(raw.iloc[i, 0]).strip().lower()
        if left.startswith("receptor volume"):
            j = i
            while j < min(i + len(reps), len(raw)):
                v = pd.to_numeric(pd.Series([raw.iloc[j, 1]]), errors="coerce").iloc[0]
                if np.isfinite(v):
                    receptor_vols.append(float(v))
                j += 1
        if left.startswith("initial") and "donor" in left:
            v = pd.to_numeric(pd.Series([raw.iloc[i, 1]]), errors="coerce").iloc[0]
            if np.isfinite(v):
                donor_dose = float(v)

    # Also parse the measured receptor concentration block, when present, so LOD/LOQ
    # rules can be applied to the analytical measurement rather than to cumulative %.
    measured_conc = None
    conc_row = None
    for i in range(len(raw)):
        row_text = " ".join(str(v) for v in raw.iloc[i].dropna().tolist()).lower()
        if "measured concentration" in row_text and "ug/ml" in row_text.replace("µ", "u"):
            conc_row = i
            break
        if "measured concentration" in row_text and "µg/ml" in row_text:
            conc_row = i
            break
    if conc_row is not None:
        tr = conc_row + 1
        conc_time = _numeric_array(raw.iloc[tr, 2:].tolist())
        valid_ct = np.isfinite(conc_time)
        if valid_ct.sum() >= 3:
            last_c = np.where(valid_ct)[0][-1] + 1
            conc_reps = []
            rr = tr + 1
            while rr < len(raw) and len(conc_reps) < len(reps):
                vals = _numeric_array(raw.iloc[rr, 2:2 + last_c].tolist())
                if np.isfinite(vals).sum() >= 3:
                    conc_reps.append(vals)
                elif conc_reps:
                    break
                rr += 1
            if len(conc_reps) == len(reps) and last_c == len(time_min) and np.allclose(conc_time[:last_c], time_min, equal_nan=True):
                measured_conc = np.asarray(conc_reps, dtype=float)

    return clean_release_data(ReleaseData(
        time_min=np.asarray(time_min, dtype=float),
        replicates_pct=np.asarray(reps, dtype=float),
        replicate_names=names,
        source_note="Auto-detected HTD/PA 'Cumulative Release' sheet",
        detected_receptor_volumes_ml=receptor_vols or None,
        detected_donor_dose_ug=donor_dose,
        measured_concentration_ug_ml=measured_conc,
    ))


def parse_simple_table(file_bytes, filename, sheet_name, input_kind, receptor_volume_ml, sample_volume_ml, donor_dose_ug):
    if filename.lower().endswith(".csv"):
        df = pd.read_csv(io.BytesIO(file_bytes))
    else:
        df = pd.read_excel(io.BytesIO(file_bytes), sheet_name=sheet_name or 0, engine="openpyxl")

    numeric = df.apply(pd.to_numeric, errors="coerce")
    numeric_cols = [c for c in numeric.columns if numeric[c].notna().sum() >= 3]
    if len(numeric_cols) < 2:
        raise ValueError("Could not identify a time column plus at least one replicate column.")
    time_col, rep_cols = numeric_cols[0], numeric_cols[1:]
    valid = numeric[time_col].notna()
    time_min = numeric.loc[valid, time_col].to_numpy(dtype=float)
    reps = numeric.loc[valid, rep_cols].to_numpy(dtype=float).T
    order = np.argsort(time_min)
    time_min, reps = time_min[order], reps[:, order]

    if input_kind == "Receptor concentration (ug/mL)":
        reps = cumulative_from_concentration(time_min, reps, receptor_volume_ml, sample_volume_ml, donor_dose_ug)
    elif input_kind == "Cumulative amount (ug)":
        reps = 100.0 * reps / donor_dose_ug

    measured_conc = None
    if input_kind == "Receptor concentration (ug/mL)":
        measured_conc = numeric.loc[valid, rep_cols].to_numpy(dtype=float).T[:, order]

    return clean_release_data(ReleaseData(
        time_min=time_min,
        replicates_pct=reps,
        replicate_names=[str(c) for c in rep_cols],
        source_note=f"Simple table: time='{time_col}', replicates={len(rep_cols)}",
        measured_concentration_ug_ml=measured_conc,
    ))



def analytical_validity_masks(data: ReleaseData, loq_ug_ml: float, receptor_volume_ml: float, donor_dose_ug: float):
    """Return replicate-level quantifiable mask and timepoint majority-valid mask.

    Preferred basis is measured receptor concentration when available. If only cumulative
    release percentages are available, LOQ is converted approximately to a first-sample
    release threshold: LOQ * receptor volume / donor dose * 100.
    A time point is valid when a strict majority (>50%) of replicates are quantifiable
    and the mean cumulative release is >0.
    """
    reps = np.asarray(data.replicates_pct, dtype=float)
    n_rep = reps.shape[0]
    needed = n_rep // 2 + 1
    if loq_ug_ml <= 0:
        quant = np.isfinite(reps) & (reps >= 0)
    elif data.measured_concentration_ug_ml is not None:
        conc = np.asarray(data.measured_concentration_ug_ml, dtype=float)
        quant = np.isfinite(conc) & (conc >= loq_ug_ml)
    else:
        threshold_pct = 100.0 * loq_ug_ml * receptor_volume_ml / donor_dose_ug if donor_dose_ug > 0 else np.inf
        quant = np.isfinite(reps) & (reps >= threshold_pct)
    counts = np.sum(quant, axis=0)
    frac = counts / max(n_rep, 1)
    finite_counts = np.sum(np.isfinite(reps), axis=0)
    sums = np.nansum(reps, axis=0)
    mean = np.divide(sums, finite_counts, out=np.full(reps.shape[1], np.nan, dtype=float), where=finite_counts > 0)
    valid_tp = (counts >= needed) & np.isfinite(mean) & (mean > 0)
    return quant, valid_tp, counts, frac

def optimize_linear_window(time_min, mean_release_pct, valid_timepoint_mask=None, target_r2=0.98, min_points=5,
                           start_min=None, end_min=None, auto_trim=True):
    """Build the quantitative fitting window.

    Automatic mode starts at the first LOQ-valid positive-release point and can
    sequentially trim terminal points until the target R^2 is reached. Manual mode
    applies user-selected start/end boundaries exactly and does not automatically
    shorten the endpoint. LOQ/positivity QC is always enforced.
    """
    t = np.asarray(time_min, dtype=float)
    y = np.asarray(mean_release_pct, dtype=float)
    if valid_timepoint_mask is None:
        valid_timepoint_mask = np.isfinite(t) & np.isfinite(y) & (y > 0)
    else:
        valid_timepoint_mask = np.asarray(valid_timepoint_mask, dtype=bool) & np.isfinite(t) & np.isfinite(y) & (y > 0)
    if not np.any(valid_timepoint_mask):
        return np.zeros_like(valid_timepoint_mask, dtype=bool), np.nan, 0

    if start_min is None:
        first_idx = np.where(valid_timepoint_mask)[0][0]
        start_min = float(t[first_idx])
    if end_min is None:
        end_min = float(np.nanmax(t[valid_timepoint_mask]))

    valid = valid_timepoint_mask & (t >= float(start_min)) & (t <= float(end_min))
    idx = np.where(valid)[0]
    if idx.size < 2:
        return np.zeros_like(valid, dtype=bool), np.nan, 0

    min_points = max(2, int(min_points))
    work = idx.copy()
    trimmed = 0

    def calc_r2(indices):
        if indices.size < 2:
            return np.nan
        p = np.polyfit(t[indices], y[indices], 1)
        return _r2(y[indices], np.polyval(p, t[indices]))

    last_r2 = calc_r2(work)
    if auto_trim:
        while work.size >= min_points:
            last_r2 = calc_r2(work)
            if np.isfinite(last_r2) and last_r2 >= target_r2:
                break
            if work.size == min_points:
                break
            work = work[:-1]
            trimmed += 1

    mask = np.zeros_like(valid, dtype=bool)
    mask[work] = True
    return mask, float(last_r2) if np.isfinite(last_r2) else np.nan, trimmed


def fit_models_with_n(time_min, y_pct, fit_mask, kp_start_min=None, kp_end_min=None):
    recorded_min_all = np.asarray(time_min, dtype=float)
    t_h_all = recorded_min_all / 60.0
    y_all = np.asarray(y_pct, dtype=float)
    base = np.asarray(fit_mask, dtype=bool) & np.isfinite(t_h_all) & np.isfinite(y_all) & (y_all >= 0)
    recorded_min = recorded_min_all[base]
    t, f = t_h_all[base], y_all[base]
    results = {"Zero": np.nan, "Higuchi": np.nan, "First-O": np.nan, "KP": np.nan, "Pep-Sah": np.nan}
    kp_n = np.nan
    if t.size < 3:
        return results, kp_n, 0

    # Match the publication scripts: transformed models use elapsed time,
    # assigning 17 min to the first point in the accepted linear window.
    elapsed_h = (recorded_min - recorded_min[0] + 17.0) / 60.0

    try:
        p = np.polyfit(t, f, 1)
        results["Zero"] = _r2(f, np.polyval(p, t))
    except Exception:
        pass
    try:
        x = np.sqrt(np.clip(elapsed_h, 0, None))
        p = np.polyfit(x, f, 1)
        results["Higuchi"] = _r2(f, np.polyval(p, x))
    except Exception:
        pass

    try:
        if np.all(f < 100):
            transformed = np.log(100.0 - f)
            p = np.polyfit(t, transformed, 1)
            results["First-O"] = _r2(transformed, np.polyval(p, t))
    except Exception:
        pass

    early = (f > 0) & (f <= 60) & (t > 0)
    if kp_start_min is not None and np.isfinite(kp_start_min):
        early &= recorded_min >= float(kp_start_min)
    if kp_end_min is not None and np.isfinite(kp_end_min):
        early &= recorded_min <= float(kp_end_min)
    te, fe = t[early], f[early]

    def kp(x, k, n):
        return k * np.power(x, n)
    if te.size >= 3:
        try:
            fraction = fe / 100.0
            popt, _ = curve_fit(
                kp, te, fraction,
                p0=[max(1e-6, fraction[0]), 0.7],
                bounds=([0, 0], [np.inf, 3]), maxfev=100000,
            )
            results["KP"] = _r2(fraction, kp(te, *popt))
            kp_n = float(popt[1])
        except Exception:
            pass

    def peppas_sahlin(x, k1, k2, m):
        return k1 * np.power(x, m) + k2 * np.power(x, 2 * m)
    if te.size >= 5:
        try:
            popt, _ = curve_fit(
                peppas_sahlin, te, fe / 100.0,
                p0=[max(1e-6, fe[0] / 100.0), 0.0005, 0.5],
                bounds=([-1e4, -1e4, 0.05], [1e4, 1e4, 2.5]),
                maxfev=100000,
            )
            results["Pep-Sah"] = _r2(fe / 100.0, peppas_sahlin(te, *popt))
        except Exception:
            pass
    return results, kp_n, int(te.size)


def permeability_from_mask(time_min, mean_release_pct, fit_mask, donor_dose_ug, donor_concentration_ug_ml, membrane_area_cm2):
    t_h = np.asarray(time_min, dtype=float)[fit_mask] / 60.0
    q_ug = np.asarray(mean_release_pct, dtype=float)[fit_mask] / 100.0 * donor_dose_ug
    if t_h.size < 2:
        return np.nan, np.nan, np.nan
    p = np.polyfit(t_h, q_ug, 1)
    yhat = np.polyval(p, t_h)
    slope_ug_h = float(p[0])
    r2 = _r2(q_ug, yhat)
    kp = np.nan
    if donor_concentration_ug_ml > 0 and membrane_area_cm2 > 0:
        kp = slope_ug_h / membrane_area_cm2 / donor_concentration_ug_ml
    return float(kp), float(r2), slope_ug_h


def cv_summary(reps, fit_mask, quantifiable_mask):
    """Median time-point CV and IQR using only replicates at/above LOQ.

    A time point must already satisfy the strict-majority LOQ rule through fit_mask.
    Within each accepted time point, BLQ replicates are excluded from mean/SD/CV.
    """
    reps = np.asarray(reps, dtype=float)
    qmask = np.asarray(quantifiable_mask, dtype=bool)
    cvs = []
    for j in np.where(np.asarray(fit_mask, dtype=bool))[0]:
        vals = reps[:, j]
        keep = qmask[:, j] & np.isfinite(vals) & (vals >= 0)
        vals = vals[keep]
        if vals.size < 2:
            continue
        mean = float(np.mean(vals))
        if mean <= 0:
            continue
        sd = float(np.std(vals, ddof=1))
        cvs.append(abs(sd / mean) * 100.0)
    if not cvs:
        return np.nan, np.nan, np.nan
    cvs = np.asarray(cvs, dtype=float)
    return float(np.nanmedian(cvs)), float(np.nanpercentile(cvs, 25)), float(np.nanpercentile(cvs, 75))


def analyze_condition(data, formulation, probe, rpm, temperature_c, donor_dose_ug, donor_concentration_ug_ml,
                      membrane_area_cm2, target_r2, min_fit_points, receptor_volume_ml, loq_ug_ml=0.0, lod_ug_ml=0.0, label=None,
                      analysis_start_min=None, analysis_end_min=None, manual_window=False,
                      kp_start_min=None, kp_end_min=None):
    data = clean_release_data(data)
    reps = np.asarray(data.replicates_pct, dtype=float)
    finite_counts = np.sum(np.isfinite(reps), axis=0)
    sums = np.nansum(reps, axis=0)
    mean = np.divide(sums, finite_counts, out=np.full(reps.shape[1], np.nan, dtype=float), where=finite_counts > 0)
    if reps.shape[0] > 1:
        sd = np.full(reps.shape[1], np.nan, dtype=float)
        for j in range(reps.shape[1]):
            vals = reps[:, j]
            vals = vals[np.isfinite(vals)]
            if vals.size >= 2:
                sd[j] = np.std(vals, ddof=1)
    else:
        sd = np.zeros_like(mean)
    quant_mask, valid_tp, counts, frac = analytical_validity_masks(data, loq_ug_ml, receptor_volume_ml, donor_dose_ug)
    fit_mask, _, trimmed = optimize_linear_window(
        data.time_min, mean, valid_tp, target_r2, min_fit_points,
        start_min=analysis_start_min, end_min=analysis_end_min, auto_trim=not manual_window
    )
    kp, linear_r2, slope = permeability_from_mask(
        data.time_min, mean, fit_mask, donor_dose_ug, donor_concentration_ug_ml, membrane_area_cm2
    )
    models, kp_n, kp_n_points = fit_models_with_n(
        data.time_min, mean, fit_mask,
        kp_start_min=kp_start_min, kp_end_min=kp_end_min,
    )
    med_cv, q1_cv, q3_cv = cv_summary(data.replicates_pct, fit_mask, quant_mask)
    fit_times = data.time_min[fit_mask]
    fit_start = float(np.nanmin(fit_times)) if fit_times.size else np.nan
    fit_end = float(np.nanmax(fit_times)) if fit_times.size else np.nan
    display_label = label or f"{formulation}\n{rpm:g} RPM"
    return ConditionResult(
        label=display_label, formulation=formulation, probe=probe, rpm=rpm, temperature_c=temperature_c,
        donor_dose_ug=donor_dose_ug, donor_concentration_ug_ml=donor_concentration_ug_ml,
        membrane_area_cm2=membrane_area_cm2, data=data, mean_pct=mean, sd_pct=sd,
        fit_mask=fit_mask, fit_start_min=fit_start, fit_end_min=fit_end,
        n_fit_points=int(fit_mask.sum()), trimmed_terminal_points=int(trimmed),
        analysis_window_mode=("Manual" if manual_window else "Automatic"),
        requested_start_min=float(analysis_start_min) if analysis_start_min is not None else np.nan,
        requested_end_min=float(analysis_end_min) if analysis_end_min is not None else np.nan,
        linear_r2=linear_r2,
        slope_ug_h=slope, kp_cm_h=kp, median_cv_pct=med_cv, cv_q1_pct=q1_cv, cv_q3_pct=q3_cv,
        models=models, kp_exponent_n=kp_n,
        kp_start_min=(float(kp_start_min) if kp_start_min is not None else fit_start),
        kp_end_min=(float(kp_end_min) if kp_end_min is not None else fit_end),
        kp_n_points=kp_n_points, valid_timepoint_mask=valid_tp,
        quantifiable_counts=counts, quantifiable_fraction=frac, loq_ug_ml=float(loq_ug_ml), lod_ug_ml=float(lod_ug_ml),
    )


def make_comprehensive_figure(results: List[ConditionResult], y_max: Optional[float] = None, title: str = ""):
    if not results:
        raise ValueError("No conditions to plot")

    n = len(results)
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), dpi=150)
    axA, axB, axC, axD, axE, axF = axes.ravel()
    cmap = plt.get_cmap("viridis")
    colors = [cmap(x) for x in np.linspace(0.20, 0.82, max(n, 2))[:n]]

    # A) release profiles
    max_x, max_y = 0, 0
    for i, r in enumerate(results):
        t = r.data.time_min
        m = r.mean_pct / 100.0 * r.donor_dose_ug
        sd = r.sd_pct / 100.0 * r.donor_dose_ug
        axA.plot(t, m, marker="o", markersize=3.6, linewidth=1.7, color=colors[i], label=r.label.replace("\n", " — "))
        if r.data.replicates_pct.shape[0] > 1:
            axA.fill_between(t, np.maximum(0, m - sd), m + sd, color=colors[i], alpha=0.15, linewidth=0)
        max_x = max(max_x, np.nanmax(t))
        max_y = max(max_y, np.nanmax(m + np.nan_to_num(sd)))
    axA.set_title("A)", loc="left", fontweight="bold")
    axA.set_xlabel("Time (min)", fontweight="bold")
    axA.set_ylabel("Mean cumulative release (µg)", fontweight="bold")
    axA.set_xlim(0, max_x * 1.06 if max_x else 1)
    axA.set_ylim(0, y_max if y_max else max(1, max_y * 1.10))
    axA.grid(True, alpha=0.22, linestyle=":")
    axA.legend(fontsize=8, loc="upper left", frameon=True)

    # B) median CV with IQR error bars
    x = np.arange(n)
    med = np.array([r.median_cv_pct for r in results], dtype=float)
    q1 = np.array([r.cv_q1_pct for r in results], dtype=float)
    q3 = np.array([r.cv_q3_pct for r in results], dtype=float)
    low = np.where(np.isfinite(med-q1), med-q1, 0)
    high = np.where(np.isfinite(q3-med), q3-med, 0)
    axB.bar(x, np.nan_to_num(med), color=colors, edgecolor="black", linewidth=1)
    axB.errorbar(x, med, yerr=[low, high], fmt="none", ecolor="black", capsize=4, linewidth=1.2)
    axB.set_title("B)", loc="left", fontweight="bold")
    axB.set_ylabel("Median CV with IQR (%)", fontweight="bold")
    axB.set_xticks(x, [r.label for r in results], rotation=45, ha="right", fontweight="bold")
    axB.set_ylim(0, 50)
    axB.set_yticks(np.arange(0, 51, 10))
    axB.grid(True, axis="y", alpha=0.22)

    # C) linear R2
    r2s = [r.linear_r2 for r in results]
    axC.bar(x, np.nan_to_num(r2s), color=colors, edgecolor="black", linewidth=1)
    axC.axhline(0.99, linestyle="--", linewidth=1.2, color="gray")
    axC.set_title("C)", loc="left", fontweight="bold")
    axC.set_ylabel("R²", fontweight="bold")
    axC.set_xticks(x, [r.label for r in results], rotation=45, ha="right", fontweight="bold")
    finite_r2 = [v for v in r2s if np.isfinite(v)]
    lo = min(0.90, max(0, (min(finite_r2) - 0.01) if finite_r2 else 0.90))
    axC.set_ylim(lo, 1.001)
    axC.grid(True, axis="y", alpha=0.22)

    # D) model-fit heatmap
    names = ["Zero", "Higuchi", "First-O", "KP", "Pep-Sah"]
    arr = np.array([[r.models[k] for k in names] for r in results], dtype=float)
    shown = np.where(np.isfinite(arr), arr, 0.80)
    im = axD.imshow(shown, aspect="auto", vmin=0.80, vmax=1.00, cmap="RdYlGn")
    axD.set_title("D)", loc="left", fontweight="bold")
    axD.set_xticks(range(len(names)), names, rotation=43, ha="right", fontweight="bold")
    axD.set_yticks(range(n), [r.label for r in results], fontweight="bold")
    for i in range(n):
        for j in range(len(names)):
            v = arr[i, j]
            axD.text(j, i, "NA" if not np.isfinite(v) else f"{v:.3f}", ha="center", va="center", fontsize=8.5, fontweight="bold")
    cb = fig.colorbar(im, ax=axD, fraction=0.046, pad=0.02)
    cb.set_label("R²", fontweight="bold")

    # E) Kp
    kps = [r.kp_cm_h for r in results]
    axE.bar(x, np.nan_to_num(kps), color=colors, edgecolor="black", linewidth=1)
    axE.set_title("E)", loc="left", fontweight="bold")
    axE.set_ylabel("Kp (cm/h)", fontweight="bold")
    axE.set_xticks(x, [r.label for r in results], rotation=45, ha="right", fontweight="bold")
    axE.grid(True, axis="y", alpha=0.22)
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_powerlimits((-2, -2))
    axE.yaxis.set_major_formatter(formatter)
    axE.ticklabel_format(axis="y", style="sci", scilimits=(-2, -2))
    axE.set_ylim(bottom=0)

    # F) KP exponent n. Reference lines are shown without mechanism labels because
    # interpretation of n thresholds depends on geometry/system.
    ns = [r.kp_exponent_n for r in results]
    axF.axhspan(0.0, 0.45, alpha=0.10, color="tab:blue")
    axF.axhspan(0.45, 0.89, alpha=0.10, color="tab:purple")
    axF.axhspan(0.89, 1.60, alpha=0.10, color="tab:orange")
    axF.axhline(0.45, linestyle="--", linewidth=1.2, color="gray")
    axF.axhline(0.89, linestyle="--", linewidth=1.2, color="gray")
    axF.axhline(1.00, linestyle="--", linewidth=1.2, color="gray")
    axF.scatter(x, ns, s=115, c=colors, edgecolors="black", linewidths=1.1, zorder=3)
    axF.set_title("F)", loc="left", fontweight="bold")
    axF.set_ylabel("Apparent KP exponent, n", fontweight="bold")
    axF.set_xticks(x, [r.label for r in results], rotation=45, ha="right", fontweight="bold")
    ymax_n = max([1.2] + [v for v in ns if np.isfinite(v)])
    axF.set_ylim(0, max(1.6, ymax_n * 1.12))
    axF.grid(True, axis="y", alpha=0.18)

    for ax in axes.ravel():
        ax.tick_params(labelsize=9)
        for spine in ax.spines.values():
            spine.set_linewidth(0.9)

    if title:
        fig.suptitle(title, fontsize=17, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97] if title else None, h_pad=2.3, w_pad=2.2)
    return fig


def result_table(results: List[ConditionResult]) -> pd.DataFrame:
    rows = []
    for r in results:
        rows.append({
            "Condition": r.label.replace("\n", " "),
            "Window mode": r.analysis_window_mode,
            "Requested start (min)": r.requested_start_min,
            "Requested end (min)": r.requested_end_min,
            "Fit start (min)": r.fit_start_min,
            "LOQ (ug/mL)": r.loq_ug_ml,
            "LOD (ug/mL)": r.lod_ug_ml,
            "Valid time points by >50% LOQ rule": int(np.sum(r.valid_timepoint_mask)),
            "Fit end (min)": r.fit_end_min,
            "Fit points": r.n_fit_points,
            "KP/Pep-Sah start (min)": r.kp_start_min,
            "KP/Pep-Sah end (min)": r.kp_end_min,
            "KP/Pep-Sah points": r.kp_n_points,
            "Terminal points excluded": r.trimmed_terminal_points,
            "Linear R2": r.linear_r2,
            "Median CV (%)": r.median_cv_pct,
            "CV Q1 (%)": r.cv_q1_pct,
            "CV Q3 (%)": r.cv_q3_pct,
            "Kp (cm/h)": r.kp_cm_h,
            "Slope (ug/h)": r.slope_ug_h,
            "KP n": r.kp_exponent_n,
            **{f"R2 {k}": v for k, v in r.models.items()},
        })
    return pd.DataFrame(rows)


def main():
    if st is None:
        raise RuntimeError("Streamlit is not installed. Install requirements and run with streamlit.")

    st.set_page_config(page_title="HTD/PA Release Analyzer", layout="wide")
    logo_column, _ = st.columns([3, 5])
    with logo_column:
        st.image(str(ASSET_DIR / "AC logo.png"), width=440)
    st.title("HTD/PA Release Analysis Dashboard")
    st.caption("Comparative six-panel analysis for one or multiple XLSX/CSV release datasets.")

    uploads = st.file_uploader("Upload one or more XLSX/CSV datasets", type=["xlsx", "xlsm", "csv"], accept_multiple_files=True)

    with st.sidebar:
        _, logo_column, _ = st.columns([1, 3, 1])
        with logo_column:
            st.image(str(ASSET_DIR / "SDL5 logo.png"), width=120)
        st.header("Global analysis rules")
        target_r2 = st.number_input("Target linear R²", value=0.980, min_value=0.0, max_value=1.0, step=0.005, format="%.3f")
        min_fit_points = st.number_input("Minimum fit points", value=5, min_value=3, step=1)
        membrane_area_cm2 = st.number_input("Membrane area (cm²)", value=0.9, min_value=0.001, step=0.1)
        sample_volume_ml = st.number_input("Sample volume replaced (mL)", value=0.25, min_value=0.0, step=0.05)
        y_max = st.number_input("Panel A Y-axis max (µg) — 0 = auto", value=0.0, min_value=0.0, step=5.0)
        figure_title = st.text_input("Optional figure title", value="")

    st.info("Rules: negative cumulative-release values are excluded; a time point is quantitatively valid only when a strict majority (>50%) of replicates are at/above LOQ and mean release is positive; CV uses quantifiable replicates only. In Automatic mode, regression starts at the first valid point and terminal points may be trimmed to meet target R². In Manual mode, you choose the first and last analysis time points and those boundaries are kept fixed; LOQ/positivity QC still applies within the selected range.")

    if not uploads:
        st.write("Upload one or more files to begin. For a comparison figure like your example, upload each condition as a separate file.")
        return

    results = []
    window_times_by_idx = {}
    for idx, uploaded in enumerate(uploads):
        file_bytes = uploaded.getvalue()
        auto_data = parse_htdpa_cumulative_sheet(file_bytes)

        with st.expander(f"Condition {idx+1}: {uploaded.name}", expanded=(idx == 0)):
            c1, c2, c3, c4 = st.columns(4)
            formulation = c1.text_input("Formulation / condition", value="P23", key=f"form_{idx}")
            probe = c2.text_input("Probe / API", value="MB", key=f"probe_{idx}")
            temperature_c = c3.number_input("Temperature (°C)", value=28.0, step=1.0, key=f"temp_{idx}")
            rpm = c4.number_input("Pump speed (RPM)", value=15.0, step=1.0, key=f"rpm_{idx}")

            c5, c6, c7, c8 = st.columns(4)
            receptor_volume_ml = c5.number_input("Receptor volume (mL)", value=20.0, min_value=0.001, step=0.5, key=f"rv_{idx}")
            donor_volume_ml = c6.number_input("Donor volume (mL)", value=0.5, min_value=0.001, step=0.1, format="%.3f", key=f"dv_{idx}")
            donor_conc = c7.number_input("Donor concentration (µg/mL)", value=500.0, min_value=0.001, step=10.0, key=f"dc_{idx}")
            uploaded_stem = uploaded.name.rsplit(".", 1)[0]
            label = c8.text_input("Figure label", value=uploaded_stem, key=f"label_{idx}")
            donor_dose = donor_volume_ml * donor_conc

            q1, q2, q3 = st.columns(3)
            loq_ug_ml = q1.number_input("Assay LOQ (µg/mL; 0 = disabled)", value=0.0, min_value=0.0, step=0.01, format="%.4f", key=f"loq_{idx}")
            lod_ug_ml = q2.number_input("Assay LOD (µg/mL; optional)", value=0.0, min_value=0.0, step=0.01, format="%.4f", key=f"lod_{idx}")
            q3.markdown("**Validity rule**  \nStrict majority: **>50%** of replicates must be ≥ LOQ.")

            data = auto_data
            if data is None:
                input_kind = st.selectbox(
                    "Replicate-column values", ["Cumulative release (%)", "Receptor concentration (ug/mL)", "Cumulative amount (ug)"], key=f"kind_{idx}"
                )
                sheet_name = None
                if not uploaded.name.lower().endswith(".csv"):
                    xls = pd.ExcelFile(io.BytesIO(file_bytes), engine="openpyxl")
                    sheet_name = st.selectbox("Worksheet", xls.sheet_names, key=f"sheet_{idx}")
                try:
                    data = parse_simple_table(file_bytes, uploaded.name, sheet_name, input_kind, receptor_volume_ml, sample_volume_ml, donor_dose)
                except Exception as e:
                    st.error(f"Could not parse {uploaded.name}: {e}")
                    continue

            st.caption(data.source_note)
            window_times_by_idx[idx] = [
                float(v) for v in np.asarray(data.time_min, dtype=float) if np.isfinite(v)
            ]
            if data.measured_concentration_ug_ml is not None:
                st.caption("Measured receptor concentrations detected: LOQ rule will be applied directly to concentration values.")
            elif loq_ug_ml > 0:
                approx_pct = 100.0 * loq_ug_ml * receptor_volume_ml / donor_dose if donor_dose > 0 else np.nan
                st.warning(f"Measured concentration values were not found. LOQ will be approximated as a cumulative-release threshold (~{approx_pct:.3f}% for the current volume/dose). For rigorous LOQ filtering, upload measured concentrations or an HTD/PA workbook containing them.")
            if data.detected_receptor_volumes_ml:
                st.caption("Detected receptor volume(s): " + ", ".join(f"{v:g}" for v in data.detected_receptor_volumes_ml) + " mL")
            if data.detected_donor_dose_ug is not None:
                st.caption(f"Detected donor dose: {data.detected_donor_dose_ug:g} µg; current entered dose: {donor_dose:g} µg")

            st.markdown("**Analysis time window**")
            window_mode = st.radio(
                "Window selection",
                ["Automatic (LOQ-valid start + terminal R² trimming)", "Manual start and end time"],
                horizontal=True, key=f"window_mode_{idx}"
            )
            manual_window = window_mode.startswith("Manual")
            analysis_start_min = None
            analysis_end_min = None
            if manual_window:
                available_times = [float(v) for v in np.asarray(data.time_min, dtype=float) if np.isfinite(v)]
                if len(available_times) < 2:
                    st.error("At least two finite time points are required for manual window selection.")
                    continue
                w1, w2 = st.columns(2)
                start_choice = w1.selectbox(
                    "Starting time point (min)", available_times, index=0,
                    format_func=lambda x: f"{x:g}", key=f"analysis_start_{idx}"
                )
                valid_end_times = [v for v in available_times if v >= start_choice]
                end_choice = w2.selectbox(
                    "Last time point (min)", valid_end_times, index=len(valid_end_times)-1,
                    format_func=lambda x: f"{x:g}", key=f"analysis_end_{idx}"
                )
                analysis_start_min = float(start_choice)
                analysis_end_min = float(end_choice)
                st.caption("Manual boundaries are fixed. Points inside the range must still pass the positive-release and >50% LOQ validity rules. Automatic terminal trimming is disabled.")

            separate_kp_window = st.checkbox(
                "Use a separate KP / Peppas–Sahlin time window",
                value=False, key=f"separate_kp_window_{idx}",
                help="Use only for a predefined model-domain or LOQ reason. KP/Pep-Sah still require 0 < F ≤ 0.60.",
            )
            kp_start_min = None
            kp_end_min = None
            if separate_kp_window:
                kp_times = [float(v) for v in np.asarray(data.time_min, dtype=float) if np.isfinite(v)]
                k1, k2 = st.columns(2)
                default_kp_start = 85.0 if 85.0 in kp_times and "MB" in probe.upper() else kp_times[0]
                kp_start_index = kp_times.index(default_kp_start)
                kp_start_min = float(k1.selectbox(
                    "KP / Peppas–Sahlin start (min)", kp_times, index=kp_start_index,
                    format_func=lambda x: f"{x:g}", key=f"kp_start_{idx}",
                ))
                kp_end_choices = [v for v in kp_times if v >= kp_start_min]
                default_end = 391.0 if 391.0 in kp_end_choices else kp_end_choices[-1]
                kp_end_min = float(k2.selectbox(
                    "KP / Peppas–Sahlin end (min)", kp_end_choices,
                    index=kp_end_choices.index(default_end),
                    format_func=lambda x: f"{x:g}", key=f"kp_end_{idx}",
                ))

            r = analyze_condition(
                data, formulation, probe, rpm, temperature_c, donor_dose, donor_conc,
                membrane_area_cm2, target_r2, int(min_fit_points), receptor_volume_ml,
                loq_ug_ml=loq_ug_ml, lod_ug_ml=lod_ug_ml, label=label,
                analysis_start_min=analysis_start_min, analysis_end_min=analysis_end_min, manual_window=manual_window,
                kp_start_min=kp_start_min, kp_end_min=kp_end_min,
            )
            results.append(r)
            st.write(
                f"Analysis mode: **{r.analysis_window_mode}**; accepted window: **{r.fit_start_min:g}–{r.fit_end_min:g} min**, "
                f"n={r.n_fit_points}; terminal points excluded={r.trimmed_terminal_points}; "
                f"R²={r.linear_r2:.3f}; Kp={r.kp_cm_h:.5f} cm/h; "
                f"KP/Pep-Sah window={r.kp_start_min:g}–{r.kp_end_min:g} min (n={r.kp_n_points}); "
                f"KP n={r.kp_exponent_n:.3f}; "
                f"LOQ-valid time points={int(np.sum(r.valid_timepoint_mask))}"
            )

    if len(uploads) > 1:
        def apply_manual_window_to_all():
            source_mode = st.session_state.get("window_mode_0", "")
            start_value = st.session_state.get("analysis_start_0")
            end_value = st.session_state.get("analysis_end_0")
            if not source_mode.startswith("Manual") or start_value is None or end_value is None:
                st.session_state["window_apply_message"] = "Select Manual start and end time for Condition 1 first."
                return
            incompatible = [
                str(i + 1) for i, times in window_times_by_idx.items()
                if start_value not in times or end_value not in times
            ]
            if incompatible:
                st.session_state["window_apply_message"] = (
                    "The selected range is not available in Condition "
                    + ", ".join(incompatible)
                    + ". Choose shared time points before applying."
                )
                return
            for i in range(len(uploads)):
                st.session_state[f"window_mode_{i}"] = "Manual start and end time"
                st.session_state[f"analysis_start_{i}"] = start_value
                st.session_state[f"analysis_end_{i}"] = end_value
            st.session_state["window_apply_message"] = (
                f"Applied {start_value:g}–{end_value:g} min to all {len(uploads)} conditions."
            )

        st.button(
            "Apply Condition 1 manual time range to all conditions",
            on_click=apply_manual_window_to_all,
            help="Copies the selected Condition 1 manual start and end times to every uploaded condition.",
        )
        if st.session_state.get("window_apply_message"):
            st.info(st.session_state["window_apply_message"])

    if not results:
        return

    fig = make_comprehensive_figure(results, y_max=(y_max if y_max > 0 else None), title=figure_title)
    st.pyplot(fig, width="stretch")

    table = result_table(results)
    st.subheader("Analysis summary")
    st.dataframe(table, width="stretch")

    png = io.BytesIO()
    fig.savefig(png, format="png", dpi=1200, bbox_inches="tight")
    png.seek(0)
    st.download_button("Download comprehensive figure PNG (1200 dpi)", png, file_name="HTDPA_comprehensive_summary.png", mime="image/png")
    st.download_button("Download analysis summary CSV", table.to_csv(index=False).encode("utf-8"), file_name="HTDPA_analysis_summary.csv", mime="text/csv")
    plt.close(fig)


if __name__ == "__main__":
    main()
