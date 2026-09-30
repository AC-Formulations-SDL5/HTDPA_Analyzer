"""
Build visualizations from CFD_Simulation/reciprocating_flow_results.npz
(produced by warp_reciprocating_flow.py):

  1. A static "snapshot grid" figure (velocity / pressure / vorticity rows x
     4 key cycle-phase columns), cropped to the rectangle interior only.
  2. A sequence of per-frame PNGs (one per saved timestep) used to build a
     self-contained interactive HTML scrubber.
"""
import base64
import io
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)
RESULTS_PATH = os.path.join(HERE, "reciprocating_flow_results.npz")
SWEEP_RESULTS_PATH = os.path.join(HERE, "vorticity_sweep_results.npz")


def load_results():
    return np.load(RESULTS_PATH)


def load_sweep_results():
    return np.load(SWEEP_RESULTS_PATH)


def crop_to_rectangle(d):
    """Crop the full (rectangle + pipe stubs) grids down to just the rectangle interior."""
    x_min = float(d["x_min"])
    dx, dy = float(d["dx"]), float(d["dy"])
    rect_w, rect_h = float(d["rect_w"]), float(d["rect_h"])

    # velocity grid uses half-spacing (Q2 nodes)
    i0_u = int(round((0.0 - x_min) / (dx / 2.0)))
    i1_u = int(round((rect_w - x_min) / (dx / 2.0)))
    frames_u = d["frames_u"][:, :, i0_u:i1_u + 1, :]

    i0_p = int(round((0.0 - x_min) / dx))
    i1_p = int(round((rect_w - x_min) / dx))
    frames_p = d["frames_p"][:, :, i0_p:i1_p + 1]
    frames_vort = d["frames_vort"][:, :, i0_p:i1_p + 1]

    extent = (0.0, rect_w, 0.0, rect_h)
    return frames_u, frames_p, frames_vort, extent


def _percentile_clim(arr, lo=2, hi=98):
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return -1.0, 1.0
    a, b = np.percentile(finite, [lo, hi])
    if a == b:
        a, b = a - 1.0, b + 1.0
    return float(a), float(b)


def draw_velocity_panel(ax, u_frame, extent, clim=None, quiver_stride=10):
    umag = np.sqrt(np.nansum(u_frame ** 2, axis=-1))
    if clim is None:
        clim = (0.0, np.nanmax(umag))
    im = ax.imshow(umag, origin="lower", extent=extent, aspect="equal", cmap="viridis", vmin=clim[0], vmax=clim[1])

    ny, nx = umag.shape
    ys = np.linspace(extent[2], extent[3], ny)
    xs = np.linspace(extent[0], extent[1], nx)
    X, Y = np.meshgrid(xs, ys)
    sl = (slice(None, None, quiver_stride), slice(None, None, max(1, quiver_stride // 3)))
    ax.quiver(
        X[sl], Y[sl], u_frame[..., 0][sl], u_frame[..., 1][sl],
        color="white", alpha=0.6, scale=800, width=0.004,
    )
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    return im, clim


def draw_scalar_panel(ax, field_frame, extent, cmap, clim=None, symmetric=False):
    if clim is None:
        clim = _percentile_clim(field_frame)
        if symmetric:
            m = max(abs(clim[0]), abs(clim[1]))
            clim = (-m, m)
    im = ax.imshow(field_frame, origin="lower", extent=extent, aspect="equal", cmap=cmap, vmin=clim[0], vmax=clim[1])
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    return im, clim


def build_snapshot_grid(d, out_path):
    frames_u, frames_p, frames_vort, extent = crop_to_rectangle(d)
    t = d["frames_t"]
    period = float(d["period"])
    t_start = t.min()

    targets = [t_start + f * period for f in (0.02, 0.35, 0.52, 0.85)]
    labels = ["just after reversal\n(inertial transient)", "steady +flow\n(left→right)",
              "just after reversal\n(inertial transient)", "steady -flow\n(right→left)"]
    idx = [int(np.argmin(np.abs(t - tt))) for tt in targets]

    u_sel = frames_u[idx]
    p_sel = frames_p[idx]
    vort_sel = frames_vort[idx]

    u_clim = (0.0, np.nanpercentile(np.sqrt(np.nansum(u_sel ** 2, axis=-1)), 99))
    p_clim = _percentile_clim(p_sel)
    pm = max(abs(p_clim[0]), abs(p_clim[1]))
    p_clim = (-pm, pm)
    vort_clim = _percentile_clim(vort_sel)
    vm = max(abs(vort_clim[0]), abs(vort_clim[1]))
    vort_clim = (-vm, vm)

    ncols = len(idx)
    fig, axes = plt.subplots(3, ncols, figsize=(2.0 * ncols + 1.2, 3 * 6.2), squeeze=False)

    for c, i in enumerate(idx):
        im0, _ = draw_velocity_panel(axes[0][c], u_sel[c], extent, clim=u_clim)
        axes[0][c].set_title(f"t={t[i]:.2f}s\n{labels[c]}", fontsize=9)

        im1, _ = draw_scalar_panel(axes[1][c], p_sel[c], extent, cmap="coolwarm", clim=p_clim)
        im2, _ = draw_scalar_panel(axes[2][c], vort_sel[c], extent, cmap="PuOr", clim=vort_clim)

        for r in range(3):
            axes[r][c].set_xticks([])
            if c > 0:
                axes[r][c].set_yticks([])

    axes[0][0].set_ylabel("velocity |u|\n(y, mm)", fontsize=9)
    axes[1][0].set_ylabel("pressure p/ρ\n(y, mm)", fontsize=9)
    axes[2][0].set_ylabel("vorticity\n(y, mm)", fontsize=9)

    fig.colorbar(im0, ax=axes[0].tolist(), shrink=0.85, label="mm/s")
    fig.colorbar(im1, ax=axes[1].tolist(), shrink=0.85, label="mm²/s²")
    fig.colorbar(im2, ax=axes[2].tolist(), shrink=0.85, label="1/s")

    fig.suptitle(
        f"Reciprocating flow in rectangle (Q={float(d['q_amp_ml_min']):.0f} mL/min, T={period:.0f}s) "
        f"— key cycle phases",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def _fig_to_base64_png(fig, dpi=85):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def build_frame_images(d, max_frames=None):
    frames_u, frames_p, frames_vort, extent = crop_to_rectangle(d)
    t = d["frames_t"]

    n = len(t) if max_frames is None else min(max_frames, len(t))
    stride = max(1, len(t) // n)
    sel = list(range(0, len(t), stride))[:n]

    images = []
    for k in sel:
        fig, axes = plt.subplots(1, 3, figsize=(7.5, 5.6))
        draw_velocity_panel(axes[0], frames_u[k], extent)
        axes[0].set_title("velocity |u| (mm/s)", fontsize=9)
        draw_scalar_panel(axes[1], frames_p[k], extent, cmap="coolwarm", symmetric=True)
        axes[1].set_title("pressure p/ρ (mm²/s²)", fontsize=9)
        draw_scalar_panel(axes[2], frames_vort[k], extent, cmap="PuOr", symmetric=True)
        axes[2].set_title("vorticity (1/s)", fontsize=9)
        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])
        fig.suptitle(f"t = {t[k]:.2f} s", fontsize=11)
        fig.tight_layout()

        images.append(_fig_to_base64_png(fig))

    return images, t[sel].tolist()


def build_vorticity_vs_flowrate_figure(s, out_path):
    """Quantitative vorticity-vs-flow-rate figure from a run_sweep() results npz.

    Left panel: peak & RMS steady-state vorticity (inside the rectangle) vs Q,
    with a linear fit (slope, intercept, R^2) annotated.
    Right panel: peak vorticity normalized by Q -- a flat line means the
    vorticity response stays perfectly linear with flow rate; a trend reveals
    inertial (Reynolds-number-dependent) departure from linearity.
    """
    q = s["q_ml_min"]
    vort_peak = s["vort_peak"]
    vort_rms = s["vort_rms"]
    re = s["re"]

    slope, intercept = np.polyfit(q, vort_peak, 1)
    fit = slope * q + intercept
    ss_res = np.sum((vort_peak - fit) ** 2)
    ss_tot = np.sum((vort_peak - np.mean(vort_peak)) ** 2)
    r2 = 1.0 - ss_res / ss_tot

    slope_rms, intercept_rms = np.polyfit(q, vort_rms, 1)

    fig, axes = plt.subplots(1, 2, figsize=(12, 7.5))

    ax = axes[0]
    ax.plot(q, vort_peak, "o-", color="#E37400", label="Peak |Vorticity|")
    #ax.plot(q, slope * q + intercept, "--", color="#c0392b", alpha=0.5,
    #        label=f"fit: {slope:.3f}·Q + {intercept:.2f}  (R²={r2:.4f})")
    ax.plot(q, vort_rms, "s-", color="#4285F4", label="RMS |Vorticity|")
    #ax.plot(q, slope_rms * q + intercept_rms, "--", color="#2f6fed", alpha=0.5,
    #        label=f"fit: {slope_rms:.3f}·Q + {intercept_rms:.2f}")
    ax.set_xlabel("Flow-rate Amplitude Q (mL/min)", fontsize=14)
    ax.set_ylabel("Vorticity (1/s)", fontsize=14)
    #ax.set_title("Vorticity vs. flow rate")
    ax.grid(False)
    ax.legend(fontsize=12, loc="upper left", frameon=False)
    ax.grid(alpha=0.3)

    ax2 = axes[1]
    ax2.plot(q, vort_peak / q, "o-", color="#0D652D")
    #ax2.axhline(slope, color="gray", linestyle=":", label=f"linear-fit slope = {slope:.3f}")
    ax2.set_xlabel("Flow-rate Amplitude Q (mL/min)", fontsize=18)
    ax2.set_ylabel("Vorticity / Q ", fontsize=18)
    #ax2.set_title("Linearity check (flat = purely linear response)", fontsize=18)
    ax2_top = ax2.twiny()
    ax2_top.set_xlim(ax2.get_xlim())
    re_ticks = np.interp(ax2.get_xticks(), q, re)
    ax2_top.set_xticks(ax2.get_xticks())
    ax2.set_xticklabels([f"{v:.0f}" for v in ax2.get_xticks()], fontsize=14)
    ax2.set_yticklabels([f"{v:.1f}" for v in ax2.get_yticks()], fontsize=14)
    ax2_top.set_xticklabels([f"{v:.0f}" for v in re_ticks], fontsize=14)
    ax2_top.set_xlabel("Reynolds number", fontsize=18)
    #ax2.legend(fontsize=8)
    ax2.grid(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")

    return {"slope": slope, "intercept": intercept, "r2": r2,
            "slope_rms": slope_rms, "intercept_rms": intercept_rms}


if __name__ == "__main__":
    d = load_results()
    build_snapshot_grid(d, os.path.join(HERE, "snapshot_grid.png"))
    s = load_sweep_results()
    build_vorticity_vs_flowrate_figure(s, os.path.join(HERE, "vorticity_vs_flowrate.png"))
