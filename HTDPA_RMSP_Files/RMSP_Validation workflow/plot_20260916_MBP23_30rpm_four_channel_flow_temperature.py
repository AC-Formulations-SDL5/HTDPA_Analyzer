#!/usr/bin/env python3
"""Plot four-channel corrected flow, cycle deciles, and temperature for
20260915_1um_dicloP23_30RPM, styled like HTDPA_0043_20_07082026_four_channel_flow_temperature.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from plot_corrected_flow_top_bottom10 import top_bottom_masks

SOURCE = Path("20260916_1um_MBP23_30rpm.csv")
OUTPUT = Path("20260916_1um_MBP23_30rpm_four_channel_flow_temperature.png")
CORRECTION_FACTORS = {1: 2.2, 2: 1.2, 3: 2.0, 4: 2.5}
CHANNELS = [1, 3, 4]


def main() -> None:
    frame = pd.read_csv(SOURCE)
    time_min = (
        pd.to_numeric(frame["time_s"], errors="coerce")
        - pd.to_numeric(frame["time_s"], errors="coerce").iloc[0]
    ) / 60.0

    fig, axes = plt.subplots(
        len(CHANNELS),
        1,
        figsize=(16, 3 * len(CHANNELS)),
        sharex=True,
        gridspec_kw={"hspace": 0.08},
    )

    cycle_counts = []
    for channel, ax in zip(CHANNELS, axes):
        flow = (
            pd.to_numeric(frame[f"s{channel}_flow_ml_min"], errors="coerce")
            * CORRECTION_FACTORS[channel]
        )
        temperature = pd.to_numeric(
            frame[f"s{channel}_temp_c"], errors="coerce"
        )
        top, bottom, cycles = top_bottom_masks(flow)
        cycle_counts.append(cycles)

        ax.plot(
            time_min,
            flow,
            color="#667784",
            linewidth=0.65,
            alpha=0.78,
            zorder=1,
        )
        ax.scatter(
            time_min[top],
            flow[top],
            s=7,
            color="#009E73",
            linewidths=0,
            zorder=3,
        )
        ax.scatter(
            time_min[bottom],
            flow[bottom],
            s=7,
            color="#D55E00",
            linewidths=0,
            zorder=3,
        )
        ax.axhline(0, color="#52606D", linewidth=1.1)
        ax.set_ylabel(
            f"CH{channel} flow\n(mL min⁻¹)",
            fontsize=18,
            fontweight="bold",
            color="#173F5F",
            labelpad=10,
        )
        ax.tick_params(
            axis="y",
            labelcolor="#173F5F",
            labelsize=15,
            width=2.0,
            length=6,
        )
        ax.grid(True, linestyle="--", linewidth=0.65, alpha=0.30)

        ax2 = ax.twinx()
        ax2.plot(
            time_min,
            temperature,
            color="#E15759",
            linewidth=1.15,
            alpha=0.9,
            zorder=4,
        )
        ax2.set_ylim(20, 40)
        ax2.set_yticks([20, 25, 30, 35, 40])
        ax2.set_ylabel(
            "Temperature (°C)",
            fontsize=18,
            fontweight="bold",
            color="#E15759",
            labelpad=10,
        )
        ax2.tick_params(
            axis="y",
            labelcolor="#E15759",
            labelsize=15,
            width=2.0,
            length=6,
        )

        for axis in (ax, ax2):
            for spine in axis.spines.values():
                spine.set_linewidth(2.0)
            for label in axis.get_yticklabels():
                label.set_fontweight("bold")

    axes[-1].set_xlabel("Time (min)", fontsize=21, fontweight="bold", labelpad=10)
    axes[-1].tick_params(axis="x", labelsize=16, width=2.0, length=7)
    for label in axes[-1].get_xticklabels():
        label.set_fontweight("bold")

    for text_item in fig.findobj(match=matplotlib.text.Text):
        text_item.set_fontweight("bold")

    fig.subplots_adjust(left=0.11, right=0.88, bottom=0.08, top=0.99)
    fig.savefig(OUTPUT, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Cycle counts CH1-CH4: {cycle_counts}")
    print(f"Saved: {OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
