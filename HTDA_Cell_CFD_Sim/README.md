# HTDA Cell CFD Simulation

Supporting CFD files for the HTD/PA (high-throughput diffusion / permeation apparatus) receptor-cell design.

The folder has two parts:

1. **Transient 2D reciprocating-flow model (Python / NVIDIA Warp).** A single receptor cell (a rectangular chamber with inlet and outlet ports) is driven by the bipolar square-wave flow of the multichannel reciprocating pump. The model gives velocity, pressure and vorticity fields over one reciprocation cycle, plus peak and RMS vorticity as a function of flow-rate amplitude.
2. **Steady-state 3D design study (Ansys Workbench / Fluent 2025 R2).** This is the project used to compare candidate cell geometries (velocity contours, streamlines, pressure).

For the scientific background and how the results were used, see [Explanation.md](Explanation.md).

---

## Contents

| File / folder | Description |
|---|---|
| `warp_reciprocating_flow.py` | Solver: geometry, physics parameters, Taylor–Hood Q2–Q1 Navier–Stokes, `run()` (single flow rate) and `run_sweep()` (flow-rate sweep) |
| `plot_reciprocating_flow.py` | Post-processing: snapshot-grid figure and vorticity-vs-flow-rate figure |
| `build_artifact.py` | Builds the standalone interactive HTML frame viewer |
| `run_all.py` | **One-command runner** for all stages (simulation → sweep → figures → viewer) |
| `reciprocating_flow_CFD.ipynb` | Notebook version: explanation, run, and plots |
| `environment.yml` | Conda environment with the required libraries |
| `reciprocating_flow_results.npz` | *Output:* ~60 frames over one period at Q = 40 mL/min (u, p, vorticity) |
| `vorticity_sweep_results.npz` | *Output:* sweep of Q = 5–50 mL/min (peak/RMS vorticity, peak velocity, pressure, Re) |
| `snapshot_grid.png` | *Output:* velocity / pressure / vorticity at 4 key phases of the cycle |
| `vorticity_vs_flowrate.png` | *Output:* vorticity and vorticity/Q vs. flow-rate amplitude |
| `reciprocating_flow_viewer.html` | *Output:* self-contained interactive viewer (open in any browser) |
| `Workbench/` | Ansys Workbench project (`HTDA_Cell_Design.wbpj` + `_files/`: geometry, mesh, Fluent case/data, exported images) |

The precomputed outputs are included, so the figures and viewer can be inspected without rerunning anything.

---

## Requirements

- **Python 3.10 or 3.11** (tested with 3.11 and warp-lang 1.13.0)
- **NVIDIA GPU with a recent CUDA driver** (recommended). Warp also runs on the CPU if no GPU is found, but it is much slower.
- Libraries: `numpy`, `matplotlib`, `warp-lang`, and `jupyterlab`/`ipykernel` for the notebook (all listed in [environment.yml](environment.yml))
- *Optional:* Ansys Workbench 2025 R2 or newer to open the `Workbench/` project

## Installation

With conda (recommended):

```bash
conda env create -f environment.yml
conda activate htda-cfd
```

With pip only:

```bash
pip install numpy matplotlib warp-lang jupyterlab ipykernel
```

To check the install (prints the CUDA device Warp found):

```bash
python -c "import warp as wp; wp.init()"
```

---

## Running

Run every command from inside this folder, because outputs are written next to the scripts.

### Option A: one command

```bash
python run_all.py                    # reuse cached .npz files, regenerate figures + viewer
python run_all.py --force            # rerun both simulations from scratch
```

Individual stages:

```bash
python run_all.py --stage sim        --force --q 40     # single transient run at Q = 40 mL/min
python run_all.py --stage sweep      --force --q-min 5 --q-max 50 --n-q 10
python run_all.py --stage figures                       # PNG figures from the .npz files
python run_all.py --stage viewer                        # interactive HTML viewer
```

Approximate runtimes on a laptop RTX GPU: single run ~2 min; sweep (10 flow rates) ~15–20 min.

### Option B: notebook

Open `reciprocating_flow_CFD.ipynb` in Jupyter or VS Code, select the `htda-cfd` kernel, and run all cells. The notebook uses the cached `.npz` files if they exist. Delete them to force a new simulation.

### Option C: individual scripts

```bash
python warp_reciprocating_flow.py    # single run -> reciprocating_flow_results.npz
python plot_reciprocating_flow.py    # figures (needs both .npz files)
python build_artifact.py             # interactive HTML viewer
```

---

## Model summary

| Parameter | Value | Variable in `warp_reciprocating_flow.py` |
|---|---|---|
| Chamber (W × H) | 7.75 × 27 mm | `RECT_W`, `RECT_H` |
| Port diameter | 4 mm, centred 7 mm above the chamber floor | `PIPE_D`, `Y_C` |
| Inlet stub length | 12 mm (3 × D, keeps the boundary condition away from the junction) | `STUB_LEN` |
| Through-plane depth (2D → flow-rate conversion) | 4 mm | `DEPTH` |
| Fluid | water, ν = 1.0 mm²/s | `NU` |
| Flow forcing | bipolar square wave, amplitude 40 mL/min, period 9 s | `Q_AMP_ML_MIN`, `PERIOD` |
| Grid / time step | 130 × 112 cells (masked `Grid2D`), Δt = 0.009 s | `NX`, `NY`, `DT` |
| Simulated duration | 2 periods; statistics and frames are taken from the 2nd period | `N_PERIODS` |

- **Units:** mm and s throughout. Pressure is kinematic (p/ρ, mm²/s²). Multiply by ρ = 1000 kg/m³ and by 10⁻⁶ to get Pa.
- **Method:** Taylor–Hood Q2/Q1 elements, implicit viscosity, semi-Lagrangian advection, hard velocity-Dirichlet BCs (parabolic profile at the two port tips, no-slip elsewhere), preconditioned CG saddle-point solve.
- **Sweep statistics:** the sweep excludes samples within 15 % of each flow reversal, so the reported vorticity reflects the fully developed flow and not the reversal transient.
- **Changing the device:** edit the constants at the top of `warp_reciprocating_flow.py`, then run `python run_all.py --force`. Velocity, pressure and vorticity magnitudes scale roughly as 1/`DEPTH`.

## Limitations

This is a simplified 2D planar idealization, and the results are design-supporting predictions. They have not been validated experimentally. The square-wave forcing is instantaneous, so it produces a real, sharp inertial pressure transient at each reversal; this is not a numerical artifact.

## Ansys Workbench project

Open `Workbench/HTDA_Cell_Design.wbpj` in Ansys Workbench 2025 R2 or newer, and keep the `HTDA_Cell_Design_files/` folder next to it. Exported contour, streamline and vector images are in `Workbench/HTDA_Cell_Design_files/user_files/`.
