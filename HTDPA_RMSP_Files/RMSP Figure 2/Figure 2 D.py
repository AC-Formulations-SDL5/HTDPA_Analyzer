#!/usr/bin/env python3
"""Generate SLF3X_peak_to_peak_vs_RPM.png at 1200 DPI from six source CSV files.

Dependencies: numpy, pandas, matplotlib.
Run: python3 SLF3X_peak_to_peak_vs_RPM.py
Options: --data-dir PATH --output PATH
See SLF3X_peak_to_peak_vs_RPM.md for inputs and calculation details.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
INPUT_FILES = {
    5: "HTDPA_5RPM_30oC_1H.csv",
    10: "HTDPA_10RPM_30oC_1H.csv",
    15: "HTDPA_15RPM_30oC_1H.csv",
    20: "HTDPA_20RPM_30oC_1H.csv",
    25: "HTDPA_25rpm_30oC_1H.csv",
    30: "HTDPA_30rpm_30oC_1.5H.csv",
}
RPM_COLORS = ["#4E79A7", "#59A14F", "#F28E2B", "#E15759", "#B07AA1", "#76B7B2"]
OLD_EMBEDDED_FACTORS = {1: 4.4, 2: 4.2, 3: 2.0, 4: 3.0}
GRAVIMETRIC_FACTORS = {1: 4.508343, 2: 4.474666, 3: 2.320488, 4: 3.340179}

def contiguous_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    edges = np.flatnonzero(np.diff(np.r_[False, mask, False]))
    return list(zip(edges[::2], edges[1::2]))

def lobe_plateaus(flow: pd.Series, polarity: int) -> np.ndarray:
    values = flow.to_numpy(dtype=float)
    finite = np.isfinite(values)
    robust_amplitude = np.nanpercentile(np.abs(values[finite]), 90)
    deadband = 0.15 * robust_amplitude
    mask = values > deadband if polarity > 0 else values < -deadband

    plateaus: list[float] = []
    for start, stop in contiguous_runs(mask):
        lobe = values[start:stop]
        if len(lobe) < 3:
            continue
        midpoint = np.nanmedian(lobe)
        outer_half = (
            lobe[lobe >= midpoint] if polarity > 0 else lobe[lobe <= midpoint]
        )
        plateaus.append(float(np.nanmean(outer_half)))
    return np.asarray(plateaus)

def analyze_flow(data_dir: Path) -> pd.DataFrame:
    rows = []
    for rpm, filename in INPUT_FILES.items():
        frame = pd.read_csv(data_dir / filename)
        for channel in range(1, 5):
            if rpm == 30 and channel == 2:
                continue
            flow = (frame[f"s{channel}_flow_ml_min"]
                    / OLD_EMBEDDED_FACTORS[channel] * GRAVIMETRIC_FACTORS[channel])
            if not np.isfinite(flow.to_numpy(dtype=float)).any():
                raise ValueError(f"No finite flow: {filename}, channel {channel}")
            positive = lobe_plateaus(flow, 1)
            negative = lobe_plateaus(flow, -1)
            if not len(positive) or not len(negative):
                raise ValueError(f"Missing flow lobes: {filename}, channel {channel}")
            q_pos, q_neg = float(positive.mean()), float(negative.mean())
            rows.append(dict(rpm=rpm, channel=channel, q_pos_ml_min=q_pos,
                             q_neg_ml_min=q_neg, peak_to_peak_ml_min=q_pos-q_neg))
    return summarize_peak_to_peak(pd.DataFrame(rows))

def summarize_peak_to_peak(results: pd.DataFrame) -> pd.DataFrame:
    """Calculate the mean, sample SD, and CV across included sensors."""
    rows = []
    for rpm, group in results.groupby("rpm", sort=True):
        values = group["peak_to_peak_ml_min"]
        mean = float(values.mean())
        sd = float(values.std(ddof=1))
        rows.append({
            "rpm": rpm,
            "n_sensors": len(group),
            "peak_to_peak_sensor_mean_ml_min": mean,
            "peak_to_peak_sensor_sd_ml_min": sd,
            "peak_to_peak_sensor_cv_pct": 100.0 * sd / abs(mean) if mean else np.nan,
        })
    return pd.DataFrame(rows)


def render_peak_to_peak_plot(average_results: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(16, 9))
    rpm = average_results["rpm"].to_numpy()
    values = average_results["peak_to_peak_sensor_mean_ml_min"].to_numpy()
    errors = average_results["peak_to_peak_sensor_sd_ml_min"].to_numpy()
    bars = ax.bar(
        rpm,
        values,
        yerr=errors,
        width=3.3,
        capsize=6,
        color=RPM_COLORS[:len(rpm)],
        edgecolor="#334E68",
        linewidth=2.2,
        error_kw={"elinewidth": 2.5, "capthick": 2.5},
    )
    for bar, value, error, cv in zip(
        bars,
        values,
        errors,
        average_results["peak_to_peak_sensor_cv_pct"],
    ):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + error + 1.2,
            f"{value:.1f}\nCV {cv:.1f}%",
            ha="center",
            va="bottom",
            fontsize=18,
            fontweight="bold",
            color="#173F5F",
        )
    ax.set_xlabel("Pump speed (RPM)", fontsize=30, fontweight="bold")
    ax.set_ylabel(
        "Peak-to-peak flow, Q+ − Q− (mL min⁻¹)",
        fontsize=30,
        fontweight="bold",
        labelpad=18,
    )
    ax.yaxis.set_label_coords(-0.075, 0.46)
    ax.set_xticks(rpm)
    ax.set_ylim(0, max(values + errors) * 1.22)
    ax.grid(axis="y", linestyle="--", linewidth=0.7, alpha=0.35)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=24, width=2.5, length=8)
    for spine in ax.spines.values():
        spine.set_linewidth(3.0)
    for label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        label.set_fontweight("bold")
    for text_item in fig.findobj(match=matplotlib.text.Text):
        text_item.set_fontweight("bold")
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.16, top=0.96)
    fig.savefig(output, format="png", dpi=1200, facecolor="white")
    plt.close(fig)

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=ROOT,
                        help="Folder containing the six source CSV files")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "SLF3X_peak_to_peak_vs_RPM.png",
                        help="Output PNG path (existing file will be overwritten)")
    args = parser.parse_args()
    averages = analyze_flow(args.data_dir)
    render_peak_to_peak_plot(averages, args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
