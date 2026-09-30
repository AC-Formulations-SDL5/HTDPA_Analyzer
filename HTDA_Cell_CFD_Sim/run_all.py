"""
One-command entry point for regenerating all HTDA cell CFD data and figures.

Usage (from this folder):
    python run_all.py                  # everything: sim + sweep + figures + viewer
    python run_all.py --stage sim      # single-Q transient run  -> reciprocating_flow_results.npz
    python run_all.py --stage sweep    # flow-rate sweep         -> vorticity_sweep_results.npz
    python run_all.py --stage figures  # snapshot_grid.png + vorticity_vs_flowrate.png
    python run_all.py --stage viewer   # reciprocating_flow_viewer.html

Existing .npz results are reused unless --force is given.
"""
import argparse
import os

import numpy as np

import warp_reciprocating_flow as sim_mod
import plot_reciprocating_flow as pf

HERE = os.path.dirname(os.path.abspath(__file__))


def stage_sim(args, sim=None):
    if os.path.exists(sim_mod.RESULTS_PATH) and not args.force:
        print(f"[sim] using cached {sim_mod.RESULTS_PATH} (use --force to rerun)")
        return sim
    return sim_mod.run(q_amp_ml_min=args.q, n_periods=args.periods, frames_per_period=60, sim=sim)


def stage_sweep(args, sim=None):
    if os.path.exists(sim_mod.SWEEP_RESULTS_PATH) and not args.force:
        print(f"[sweep] using cached {sim_mod.SWEEP_RESULTS_PATH} (use --force to rerun)")
        return sim
    q_values = np.linspace(args.q_min, args.q_max, args.n_q)
    sim_mod.run_sweep(q_values, sim=sim, n_periods=args.periods, steady_margin_frac=0.15)
    return sim


def stage_figures(args):
    pf.build_snapshot_grid(pf.load_results(), os.path.join(HERE, "snapshot_grid.png"))
    pf.build_vorticity_vs_flowrate_figure(pf.load_sweep_results(), os.path.join(HERE, "vorticity_vs_flowrate.png"))


def stage_viewer(args):
    import build_artifact
    build_artifact.main()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stage", choices=["all", "sim", "sweep", "figures", "viewer"], default="all")
    parser.add_argument("--force", action="store_true", help="rerun simulations even if cached .npz files exist")
    parser.add_argument("--q", type=float, default=sim_mod.Q_AMP_ML_MIN, help="flow-rate amplitude for the single run [mL/min]")
    parser.add_argument("--periods", type=int, default=sim_mod.N_PERIODS, help="number of reciprocation periods to simulate")
    parser.add_argument("--q-min", type=float, default=5.0, help="sweep: lowest flow-rate amplitude [mL/min]")
    parser.add_argument("--q-max", type=float, default=50.0, help="sweep: highest flow-rate amplitude [mL/min]")
    parser.add_argument("--n-q", type=int, default=10, help="sweep: number of flow-rate points")
    args = parser.parse_args()

    sim = None
    if args.stage in ("all", "sim"):
        sim = stage_sim(args, sim)
    if args.stage in ("all", "sweep"):
        sim = stage_sweep(args, sim)
    if args.stage in ("all", "figures"):
        stage_figures(args)
    if args.stage in ("all", "viewer"):
        stage_viewer(args)


if __name__ == "__main__":
    main()
