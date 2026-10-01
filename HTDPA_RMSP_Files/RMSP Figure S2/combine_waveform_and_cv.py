from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parent
WORKBOOK = ROOT / "HTDPA_RMSP_140726_normalized.xlsx"
SHEETS = ("HTDPA_3h_140726", "HTDPA_3h_150726")
OUT_PNG = ROOT / "waveform_and_cv_combined_vertical.png"
DPI = 1200
CHANNELS = (1, 2, 3, 4)
COLORS = ("#4E79A7", "#F28E2B", "#59A14F", "#E15759")
PARAMETERS = ("Q+", "|Q-|", "P-P", "Q+/|Q-|")


def cycle_descriptors(frame: pd.DataFrame, channel: int) -> pd.DataFrame:
    values = pd.to_numeric(frame[f"s{channel}_flow_norm_ml_min"], errors="coerce").to_numpy()
    times = pd.to_numeric(frame["time_s"], errors="coerce").to_numpy()
    finite = np.isfinite(values) & np.isfinite(times)
    values, times = values[finite], times[finite]
    edges = np.flatnonzero(np.diff(np.sign(values)) != 0) + 1
    starts = np.r_[0, edges]
    stops = np.r_[edges, len(values)]
    rows = []
    for index in range(1, len(starts) - 1):
        neg = values[starts[index - 1]:stops[index - 1]]
        pos = values[starts[index]:stops[index]]
        if len(neg) < 3 or len(pos) < 3 or not (np.all(neg < 0) and np.all(pos > 0)):
            continue
        q_minus = float(np.mean(neg[neg <= np.median(neg)]))
        q_plus = float(np.mean(pos[pos >= np.median(pos)]))
        signal = np.r_[neg, pos]
        plateau = np.r_[np.full(len(neg), q_minus), np.full(len(pos), q_plus)]
        amplitude = float(np.mean(np.abs(plateau)))
        rows.append(
            {
                "cycle": len(rows) + 1,
                "time_s": float(np.mean(times[starts[index - 1]:stops[index]])),
                "q_plus": q_plus,
                "q_minus": q_minus,
                "peak_to_peak": q_plus - q_minus,
                "symmetry": abs(q_plus) / abs(q_minus),
                "cv_signal_raw": float(np.std(signal, ddof=1) / amplitude * 100),
                "cv_residual": float(np.std(signal - plateau, ddof=1) / amplitude * 100),
            }
        )
    return pd.DataFrame(rows)


def load_data() -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    if not WORKBOOK.exists():
        raise FileNotFoundError(f"Missing workbook: {WORKBOOK}")
    frames = {}
    descriptors = {}
    for sheet in SHEETS:
        frame = pd.read_excel(WORKBOOK, sheet_name=sheet)
        frames[sheet] = frame
        parts = []
        for channel in CHANNELS:
            values = cycle_descriptors(frame, channel)
            values["channel"] = f"CH{channel}"
            parts.append(values)
        descriptors[sheet] = pd.concat(parts, ignore_index=True)
    return frames, descriptors


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def style_axis(ax: plt.Axes) -> None:
    ax.tick_params(labelsize=7, width=0.8, length=3)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.grid(axis="y", linestyle="--", alpha=0.25, linewidth=0.6)


def draw_panel_a(axes: list[plt.Axes], descriptors: pd.DataFrame) -> None:
    summaries = descriptors.groupby("channel").agg(
        q_plus=("q_plus", "mean"), q_plus_sd=("q_plus", "std"),
        q_minus=("q_minus", "mean"), q_minus_sd=("q_minus", "std"),
        peak_to_peak=("peak_to_peak", "mean"), peak_to_peak_sd=("peak_to_peak", "std"),
        symmetry=("symmetry", "mean"), symmetry_sd=("symmetry", "std"),
    )
    labels = [f"CH{i}" for i in CHANNELS]
    x = np.arange(len(labels))
    width = 0.34
    axes[0].bar(x - width / 2, summaries.q_plus, width, yerr=summaries.q_plus_sd, capsize=3, color=COLORS, edgecolor="#243B53")
    axes[0].bar(x + width / 2, summaries.q_minus.abs(), width, yerr=summaries.q_minus_sd, capsize=3, color=COLORS, edgecolor="#243B53", hatch="///")
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Flow (mL min-1)", fontsize=8)
    axes[0].set_title("Plateau flow rate", fontsize=9, fontweight="bold")
    axes[0].legend(handles=[Patch(facecolor="#AFC7E3", edgecolor="#243B53", label="Push Q+"), Patch(facecolor="#AFC7E3", edgecolor="#243B53", hatch="///", label="Pull |Q-|")], fontsize=6, frameon=False)
    axes[1].bar(x, summaries.peak_to_peak, yerr=summaries.peak_to_peak_sd, capsize=3, color=COLORS, edgecolor="#243B53")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylabel("P-P (mL min-1)", fontsize=8)
    axes[1].set_title("Average peak-to-peak amplitude", fontsize=9, fontweight="bold")
    axes[2].bar(x, summaries.symmetry, yerr=summaries.symmetry_sd, capsize=3, color=COLORS, edgecolor="#243B53")
    axes[2].axhline(1.0, color="#243B53", linestyle="--", linewidth=1)
    axes[2].set_xticks(x, labels)
    axes[2].set_ylabel("Q+/|Q-|", fontsize=8)
    axes[2].set_title("Push-pull symmetry", fontsize=9, fontweight="bold")
    for ax in axes:
        style_axis(ax)
    add_panel_label(axes[0], "A")


def draw_panel_b(ax: plt.Axes, descriptors: dict[str, pd.DataFrame]) -> None:
    metrics = {"Q+": "q_plus", "|Q-|": "q_minus", "P-P": "peak_to_peak", "Q+/|Q-|": "symmetry"}
    x = np.arange(len(PARAMETERS))
    width = 0.24
    for offset, (sheet, label, color) in enumerate(zip(SHEETS, ("140726", "150726"), ("#4E79A7", "#F28E2B"))):
        cvs = []
        data = descriptors[sheet]
        for parameter in PARAMETERS:
            series = data[metrics[parameter]].abs()
            cvs.append(float(series.std(ddof=1) / abs(series.mean()) * 100))
        ax.bar(x + (offset - 0.5) * width, cvs, width, color=color, edgecolor="#243B53", label=label)
    combined = pd.concat([descriptors[SHEETS[0]], descriptors[SHEETS[1]]])
    interday = [float(combined[metrics[p]].abs().std(ddof=1) / abs(combined[metrics[p]].abs().mean()) * 100) for p in PARAMETERS]
    ax.bar(x + width / 2, interday, width, color="#A1A1A1", edgecolor="#243B53", label="Inter-day")
    ax.set_xticks(x, PARAMETERS)
    ax.set_ylabel("Coefficient of variation (%)", fontsize=8)
    ax.set_title("Within-day and inter-day variation", fontsize=9, fontweight="bold")
    ax.legend(fontsize=6, frameon=False, ncol=3)
    style_axis(ax)
    add_panel_label(ax, "B")


def draw_temperature_panel(ax: plt.Axes, frame: pd.DataFrame, label: str, title: str) -> None:
    means = []
    sds = []
    for channel in CHANNELS:
        series = pd.to_numeric(frame[f"s{channel}_temp_c"], errors="coerce").dropna()
        means.append(float(series.mean()))
        sds.append(float(series.std(ddof=1)))
    x = np.arange(len(CHANNELS))
    ax.bar(x, means, yerr=sds, capsize=3, color=COLORS, edgecolor="#243B53")
    inter_mean = float(np.mean(means))
    inter_cv = float(np.std(means, ddof=1) / abs(inter_mean) * 100)
    ax.axhline(inter_mean, color="#243B53", linestyle="--", linewidth=0.8)
    ax.text(0.98, 0.94, f"CVtemp inter-channel = {inter_cv:.2f}%", transform=ax.transAxes, ha="right", va="top", fontsize=6, fontweight="bold")
    ax.set_xticks(x, [f"CH{i}" for i in CHANNELS])
    ax.set_ylabel("Temperature (degC)", fontsize=8)
    ax.set_title(title, fontsize=9, fontweight="bold")
    style_axis(ax)
    add_panel_label(ax, label)


def draw_panel_e(ax: plt.Axes, frame: pd.DataFrame) -> None:
    subset = frame.loc[pd.to_numeric(frame["time_s"], errors="coerce").between(0, 60)].copy()
    time = pd.to_numeric(subset["time_s"], errors="coerce")
    temperature_axis = ax.twinx()
    for channel, color in zip(CHANNELS, COLORS):
        ax.plot(time, subset[f"s{channel}_flow_norm_ml_min"], color=color, linewidth=0.7, label=f"CH{channel} flow")
        temperature_axis.plot(time, subset[f"s{channel}_temp_c"], color=color, linewidth=0.7, linestyle="--", alpha=0.75, label=f"CH{channel} temp")
    ax.set_xlabel("Time (s)", fontsize=8)
    ax.set_ylabel("Flow (mL min-1)", fontsize=8)
    temperature_axis.set_ylabel("Temperature (degC)", fontsize=8)
    ax.set_title("Reciprocating pump signal (RPM = 30)", fontsize=9, fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    handles_2, labels_2 = temperature_axis.get_legend_handles_labels()
    ax.legend(handles + handles_2, labels + labels_2, fontsize=5, ncol=4, frameon=False, loc="upper right")
    style_axis(ax)
    temperature_axis.tick_params(labelsize=7)
    add_panel_label(ax, "E")


def main() -> None:
    frames, descriptors = load_data()
    figure = plt.figure(figsize=(18, 10), constrained_layout=True)
    grid = figure.add_gridspec(3, 6, height_ratios=(1.25, 0.9, 1.0))
    panel_a = [figure.add_subplot(grid[0, start:start + 2]) for start in (0, 2, 4)]
    panel_b = figure.add_subplot(grid[1, :])
    panel_c = figure.add_subplot(grid[2, 0:2])
    panel_d = figure.add_subplot(grid[2, 2:4])
    panel_e = figure.add_subplot(grid[2, 4:6])
    draw_panel_a(panel_a, descriptors[SHEETS[0]])
    draw_panel_b(panel_b, descriptors)
    draw_temperature_panel(panel_c, frames[SHEETS[0]], "C", "Before fluidic adjustment (140726)")
    draw_temperature_panel(panel_d, frames[SHEETS[1]], "D", "After fluidic adjustment (150726)")
    draw_panel_e(panel_e, frames[SHEETS[1]])
    figure.savefig(OUT_PNG, dpi=DPI, facecolor="white")
    plt.close(figure)
    print(f"Saved: {OUT_PNG}")


if __name__ == "__main__":
    main()
