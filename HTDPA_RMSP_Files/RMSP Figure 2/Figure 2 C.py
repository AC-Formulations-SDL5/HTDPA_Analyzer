#!/usr/bin/env python3
"""Reproduce SLF3X_flow_average_and_peak_to_peak_vs_RPM.png.

Install dependencies: python3 -m pip install numpy pandas matplotlib pillow
Run: python3 SLF3X_flow_average_and_peak_to_peak_vs_RPM.py
Optional: --data-dir PATH --output PATH

Inputs: the six CSV files listed in INPUT_FILES below (no Excel reader needed).
Flow is recalibrated using GRAVIMETRIC_FACTORS / OLD_EMBEDDED_FACTORS.
Lobes use a 15% amplitude deadband, at least three samples, and an outer-half
plateau mean. Points and bars are means across sensors; error bars are sample
SD across sensors. Channel 2 at 30 RPM is excluded, matching the original plot.
Peak-to-peak flow is Q+ minus Q-. The inset is |mean Q+| / |mean Q-|,
with SD across individual sensor ratios. CV = 100 * SD / abs(mean).
The script is standalone and does not import other project scripts.
"""

from __future__ import annotations

import argparse
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parent
# The 1200 dpi panels exceed Pillow's default pixel-count safety limit.
Image.MAX_IMAGE_PIXELS = None
INPUT_FILES = {
    5: "HTDPA_5RPM_30oC_1H.csv",
    10: "HTDPA_10RPM_30oC_1H.csv",
    15: "HTDPA_15RPM_30oC_1H.csv",
    20: "HTDPA_20RPM_30oC_1H.csv",
    25: "HTDPA_25rpm_30oC_1H.csv",
    30: "HTDPA_30rpm_30oC_1.5H.csv",
}
RPM_COLORS = ["#4E79A7", "#59A14F", "#F28E2B", "#E15759", "#B07AA1", "#76B7B2"]
DPI = 1200
OLD_EMBEDDED_FACTORS = {1: 4.4, 2: 4.2, 3: 2.0, 4: 3.0}
GRAVIMETRIC_FACTORS = {1: 4.508343, 2: 4.474666, 3: 2.320488, 4: 3.340179}


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
    return analyze_sensor_average_flow(pd.DataFrame(rows))



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


def analyze_sensor_average_flow(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rpm, group in results.groupby("rpm", sort=True):
        if int(rpm) == 30:
            group = group.loc[group["channel"] != 2]
        q_pos_mean = float(group["q_pos_ml_min"].mean())
        q_pos_sd = float(group["q_pos_ml_min"].std(ddof=1))
        q_neg_mean = float(group["q_neg_ml_min"].mean())
        q_neg_sd = float(group["q_neg_ml_min"].std(ddof=1))
        peak_to_peak = group["peak_to_peak_ml_min"]
        peak_to_peak_mean = float(peak_to_peak.mean())
        peak_to_peak_sd = float(peak_to_peak.std(ddof=1))
        channel_q_ratios = group["q_pos_ml_min"].abs() / group["q_neg_ml_min"].abs()
        q_ratio_sd = float(channel_q_ratios.std(ddof=1))
        rows.append(
            {
                "rpm": rpm,
                "n_sensors": len(group),
                "q_pos_sensor_mean_ml_min": q_pos_mean,
                "q_pos_sensor_sd_ml_min": q_pos_sd,
                "q_pos_sensor_cv_pct": 100.0 * q_pos_sd / abs(q_pos_mean),
                "q_neg_sensor_mean_ml_min": q_neg_mean,
                "q_neg_sensor_sd_ml_min": q_neg_sd,
                "q_neg_sensor_cv_pct": 100.0 * q_neg_sd / abs(q_neg_mean),
                "peak_to_peak_sensor_mean_ml_min": peak_to_peak_mean,
                "peak_to_peak_sensor_sd_ml_min": peak_to_peak_sd,
                "peak_to_peak_sensor_cv_pct": (
                    100.0 * peak_to_peak_sd / abs(peak_to_peak_mean)
                ),
                "sensor_average_q_ratio": abs(q_pos_mean) / abs(q_neg_mean),
                "sensor_q_ratio_sd": q_ratio_sd,
                "sensor_q_ratio_cv_pct": (
                    100.0 * q_ratio_sd / abs(float(channel_q_ratios.mean()))
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("rpm")


def render_sensor_average_flow_plot(average_results: pd.DataFrame, output: BytesIO) -> None:
    fig, ax = plt.subplots(figsize=(16, 9))
    rpm = average_results["rpm"].to_numpy()

    ax.plot(
        rpm,
        average_results["q_pos_sensor_mean_ml_min"],
        color="#52606D",
        linewidth=2.4,
        zorder=1,
    )
    ax.plot(
        rpm,
        average_results["q_neg_sensor_mean_ml_min"],
        color="#52606D",
        linewidth=2.4,
        linestyle="--",
        zorder=1,
    )
    for row, condition_color in zip(
        average_results.itertuples(index=False), RPM_COLORS
    ):
        ax.errorbar(
            row.rpm,
            row.q_pos_sensor_mean_ml_min,
            yerr=row.q_pos_sensor_sd_ml_min,
            marker="o",
            markersize=11,
            capsize=6,
            linestyle="none",
            color=condition_color,
            markeredgecolor="#334E68",
            markeredgewidth=1.1,
            zorder=3,
        )
        ax.errorbar(
            row.rpm,
            row.q_neg_sensor_mean_ml_min,
            yerr=row.q_neg_sensor_sd_ml_min,
            marker="s",
            markersize=11,
            capsize=6,
            linestyle="none",
            color=condition_color,
            markeredgecolor="#334E68",
            markeredgewidth=1.1,
            zorder=3,
        )
        ax.annotate(
            f"CV {row.q_pos_sensor_cv_pct:.1f}%",
            (row.rpm, row.q_pos_sensor_mean_ml_min + row.q_pos_sensor_sd_ml_min),
            xytext=(0, 2 if int(row.rpm) == 30 else 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=18,
            fontweight="bold",
            color=condition_color,
        )
        ax.annotate(
            f"CV {row.q_neg_sensor_cv_pct:.1f}%",
            (row.rpm, row.q_neg_sensor_mean_ml_min + row.q_neg_sensor_sd_ml_min),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=18,
            fontweight="bold",
            color=condition_color,
        )

    ax.axhline(0, color="#52606D", linewidth=1.2)
    inset = ax.inset_axes([0.52, 0.32, 0.44, 0.43])
    ratios = average_results["sensor_average_q_ratio"].to_numpy()
    ratio_errors = average_results["sensor_q_ratio_sd"].to_numpy()
    inset_x = np.arange(len(rpm))
    inset.bar(
        inset_x,
        ratios,
        yerr=ratio_errors,
        width=0.68,
        capsize=5,
        color=RPM_COLORS[:len(rpm)],
        edgecolor="#006B4F",
        linewidth=0.8,
        error_kw={"elinewidth": 1.8, "capthick": 1.8, "ecolor": "#243B53"},
    )
    inset.axhline(1.0, color="#334E68", linewidth=1.1, linestyle=":")
    inset.set_title("Push–pull symmetry |Q+|/|Q−|", fontsize=21, fontweight="bold")
    inset.set_ylim(0, 1.05)
    inset.set_yticks(np.arange(0, 1.01, 0.2))
    inset_labels = [f"{value:g}" for value in rpm]
    inset.set_xticks(inset_x, inset_labels)
    inset.text(
        inset_x[-1] + 0.38,
        -0.055,
        "RPM",
        transform=inset.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        color="#243B53",
        clip_on=False,
    )
    inset.tick_params(labelsize=18)
    inset.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    inset.set_facecolor("white")
    for tick_label in [*inset.get_xticklabels(), *inset.get_yticklabels()]:
        tick_label.set_fontweight("bold")
    for spine in inset.spines.values():
        spine.set_color("#7B8794")
        spine.set_linewidth(2.0)
    ax.set_xlabel("Pump speed (RPM)", fontsize=30, fontweight="bold")
    ax.set_ylabel(
        "Mean plateau flow (mL min⁻¹)",
        fontsize=30,
        fontweight="bold",
        labelpad=18,
    )
    ax.yaxis.set_label_coords(-0.075, 0.46)
    ax.set_xticks(rpm)
    ax.grid(True, linestyle="--", linewidth=0.7, alpha=0.35)
    ax.tick_params(labelsize=24)
    for spine in ax.spines.values():
        spine.set_linewidth(2.5)
    for tick_label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        tick_label.set_fontweight("bold")
    q_legend_handles = [
        Line2D([0], [0], color="#52606D", marker="o", linewidth=2.4,
               markersize=10, label="Mean Q+"),
        Line2D([0], [0], color="#52606D", marker="s", linewidth=2.4,
               linestyle="--", markersize=10, label="Mean Q−"),
    ]
    ax.legend(
        handles=q_legend_handles,
        frameon=False,
        loc="upper left",
        ncol=2,
        prop={"size": 20, "weight": "bold"},
    )
    for text_item in fig.findobj(match=matplotlib.text.Text):
        text_item.set_fontweight("bold")
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.16, top=0.96)
    fig.savefig(output, format="png", dpi=DPI, facecolor="white")
    plt.close(fig)


def render_peak_to_peak_plot(average_results: pd.DataFrame, output: BytesIO) -> None:
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
    fig.savefig(output, format="png", dpi=DPI, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, default=ROOT / "SLF3X_flow_average_and_peak_to_peak_vs_RPM.png")
    args = parser.parse_args()
    averages = analyze_flow(args.data_dir)
    with BytesIO() as left_buffer, BytesIO() as right_buffer:
        render_sensor_average_flow_plot(averages, left_buffer)
        render_peak_to_peak_plot(averages, right_buffer)
        left_buffer.seek(0)
        right_buffer.seek(0)
        with Image.open(left_buffer) as left, Image.open(right_buffer) as right:
            height = max(left.height, right.height)
            canvas = Image.new("RGB", (left.width + right.width, height), "white")
            canvas.paste(left, (0, (height-left.height)//2))
            canvas.paste(right, (left.width, (height-right.height)//2))
            canvas.save(args.output, format="PNG", dpi=(DPI, DPI), optimize=True)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
