from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKBOOK = Path("HTDPA_RMSP_140726_normalized.xlsx")
SHEET = "HTDPA_3h_150726"
OUT_PNG = Path("temperature_bar_sd_cv_HTDPA_3h_150726_full_period.png")
OUT_PNG_60MIN = Path("temperature_bar_sd_cv_HTDPA_3h_150726_first_60min.png")
STYLE_LINEWIDTH = 1.2
TEMP_COLORS = ["#9ecae1", "#fdae6b", "#a1d99b", "#bcbddc", "#9edae5"]
TEMP_EDGE_COLORS = ["#1f4e79", "#a64b00", "#2f6b2f", "#654a91", "#277b83"]


def plot_temperature_bar_sd_cv(df: pd.DataFrame, title: str, out_png: Path) -> None:
    channels = [1, 2, 3, 4]
    labels = [f"CH{ch}" for ch in channels] + ["Inter"]

    means = []
    sds = []
    cvs = []
    for ch in channels:
        col = f"s{ch}_temp_c"
        if col not in df.columns:
            raise KeyError(f"Missing required column: {col}")
        s = pd.to_numeric(df[col], errors="coerce").dropna()
        mean_val = float(s.mean())
        sd_val = float(s.std(ddof=1)) if len(s) > 1 else 0.0
        cv_val = (sd_val / abs(mean_val) * 100.0) if mean_val != 0 else np.nan
        means.append(mean_val)
        sds.append(sd_val)
        cvs.append(cv_val)

    inter_mean = float(np.mean(means))
    inter_sd = float(np.std(means, ddof=1)) if len(means) > 1 else 0.0
    inter_cv = (inter_sd / abs(inter_mean) * 100.0) if inter_mean != 0 else np.nan

    means_plot = means + [inter_mean]
    sds_plot = sds + [inter_sd]
    cvs_plot = cvs + [inter_cv]

    x = np.arange(len(labels), dtype=float)
    width = 0.38

    fig, ax = plt.subplots(figsize=(12.0, 5.8), constrained_layout=True)
    bars_temp = ax.bar(
        x - width / 2,
        means_plot,
        yerr=sds_plot,
        capsize=5,
        width=width,
        color=TEMP_COLORS,
        edgecolor=TEMP_EDGE_COLORS,
        linewidth=STYLE_LINEWIDTH,
        error_kw={"elinewidth": STYLE_LINEWIDTH, "capthick": STYLE_LINEWIDTH},
        label="Temperature Mean ± SD",
    )

    ax2 = ax.twinx()
    bars_cv = ax2.bar(
        x + width / 2,
        cvs_plot,
        width=width,
        color="#f8d7da",
        edgecolor="#9d174d",
        linewidth=STYLE_LINEWIDTH,
        label="CV (%)",
    )

    ax.set_xlabel("Channel", fontsize=15, fontweight="bold")
    ax.set_ylabel("Temperature (degC)", fontsize=15, fontweight="bold")
    ax2.set_ylabel("CV (%)", fontsize=15, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=14)
    ax.tick_params(axis="x", labelsize=14, width=STYLE_LINEWIDTH)
    ax.tick_params(axis="y", labelsize=13, width=STYLE_LINEWIDTH)
    ax2.tick_params(axis="y", labelsize=13, width=STYLE_LINEWIDTH)
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels() + ax2.get_yticklabels():
        tick_label.set_fontweight("bold")
    for axis in (ax, ax2):
        for spine in axis.spines.values():
            spine.set_linewidth(STYLE_LINEWIDTH)
    ax.set_ylim(20, 40)
    cv_max = float(np.nanmax(cvs_plot)) if len(cvs_plot) else 0.0
    ax2.set_ylim(0, max(5.0, cv_max * 1.35))
    ax.grid(axis="y", alpha=0.3)

    temp_low, temp_high = ax.get_ylim()
    cv_low, cv_high = ax2.get_ylim()
    for bar, mean_val, sd_val, cv_val, text_color in zip(
        bars_temp, means_plot, sds_plot, cvs_plot, TEMP_EDGE_COLORS
    ):
        cv_top_temp_scale = temp_low + (
            (cv_val - cv_low) / (cv_high - cv_low) * (temp_high - temp_low)
        )
        label_y = max(mean_val + sd_val, cv_top_temp_scale) + 0.35
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            label_y,
            f"{mean_val:.2f} ± {sd_val:.2f} °C",
            transform=ax.transData,
            ha="center",
            va="bottom",
            fontsize=11,
            color=text_color,
            fontweight="bold",
        )

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(
        h1 + h2,
        l1 + l2,
        loc="upper left",
        fontsize=13,
        prop={"size": 13, "weight": "bold"},
    )

    fig.savefig(out_png, dpi=1200)
    plt.close(fig)
    print(f"Saved: {out_png}")


def main() -> None:
    if not WORKBOOK.exists():
        raise FileNotFoundError(f"Missing workbook: {WORKBOOK}")

    df = pd.read_excel(WORKBOOK, sheet_name=SHEET)
    plot_temperature_bar_sd_cv(
        df,
        "",
        OUT_PNG,
    )

    if "time_s" not in df.columns:
        raise KeyError("Missing required column: time_s")
    time_min = pd.to_numeric(df["time_s"], errors="coerce") / 60.0
    df_60min = df.loc[time_min.between(0, 60, inclusive="both")]
    plot_temperature_bar_sd_cv(
        df_60min,
        "",
        OUT_PNG_60MIN,
    )


if __name__ == "__main__":
    main()
