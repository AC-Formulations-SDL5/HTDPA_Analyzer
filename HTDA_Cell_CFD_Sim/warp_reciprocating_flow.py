"""
Transient 2D CFD of reciprocating flow in a rectangular chamber fed by two pipe
stubs, solved with NVIDIA Warp's FEM module (warp.fem): masked Grid2D geometry,
Taylor-Hood (Q2-Q1) incompressible Navier-Stokes, semi-Lagrangian advection,
hard velocity-Dirichlet boundary conditions.

Geometry (mm):
    Rectangle: width RECT_W, height RECT_H.
    Two horizontal pipe stubs (interior "diameter" PIPE_D) attached to the left
    and right sides of the rectangle, centered at height Y_C from the bottom.

Flow (mm, s):
    A bipolar square-wave flow rate Q(t) (amplitude Q_AMP, period PERIOD, zero
    mean) drives the flow: for the first half-period fluid enters from the left
    pipe and exits the right pipe; for the second half-period the direction
    reverses. The imposed profile is the 2D (plane Poiseuille) parabolic
    profile matching Q(t) given an assumed through-plane depth DEPTH.

Units: millimeters and seconds throughout, so that water's kinematic
viscosity is exactly 1.0 mm^2/s. Pressure is "kinematic" (p/rho, mm^2/s^2);
multiply by density (1000 kg/m^3 for water, in SI) and convert mm->m to get
physical pressure in Pa.
"""
import os
import time

import numpy as np
import warp as wp
import warp.fem as fem
import warp.examples.fem.utils as fem_example_utils
from warp.fem.linalg import array_axpy
from warp.sparse import bsr_copy, bsr_mv
from warp.utils import array_cast

# ---------------------------------------------------------------- geometry --
RECT_W = 7.75
RECT_H = 27.0
PIPE_D = 4.0
PIPE_R = PIPE_D / 2.0
Y_C = 7.0                      # pipe centerline height, measured from the bottom
Y_LO = Y_C - PIPE_R
Y_HI = Y_C + PIPE_R
STUB_LEN = 12.0                 # 3x pipe diameter buffer; keeps the imposed BC away from the junction

X_MIN = -STUB_LEN
X_MAX = RECT_W + STUB_LEN
Y_MIN = 0.0
Y_MAX = RECT_H

# ---------------------------------------------------------------- physics --
NU = 1.0                        # mm^2/s  (water kinematic viscosity at ~20C)
DEPTH = 4.0                     # mm, assumed through-plane depth of the 2D slice
Q_AMP_ML_MIN = 40.0
Q_AMP = Q_AMP_ML_MIN * 1000.0 / 60.0         # mm^3/s
U_MEAN_AMP = Q_AMP / (PIPE_D * DEPTH)        # mm/s
U_MAX_AMP = 1.5 * U_MEAN_AMP                  # 2D parabolic (plane Poiseuille) peak
PERIOD = 9.0

# ---------------------------------------------------------------- grid -----
NX, NY = 130, 112
DX = (X_MAX - X_MIN) / NX
DY = (Y_MAX - Y_MIN) / NY

DT = 0.009
N_PERIODS = 2
TOL = 1.0e-3

RESULTS_PATH = os.path.join(os.path.dirname(__file__), "reciprocating_flow_results.npz")
SWEEP_RESULTS_PATH = os.path.join(os.path.dirname(__file__), "vorticity_sweep_results.npz")


def q_ml_min_to_u_max(q_ml_min):
    """Peak 2D (plane-Poiseuille) boundary velocity [mm/s] for a given flow-rate amplitude [mL/min]."""
    q_amp = q_ml_min * 1000.0 / 60.0             # mm^3/s
    u_mean = q_amp / (PIPE_D * DEPTH)            # mm/s
    return 1.5 * u_mean


@fem.integrand
def cell_activity(
    s: fem.Sample,
    domain: fem.Domain,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    rect_w: float,
    y_lo: float,
    y_hi: float,
):
    pos = domain(s)
    x = pos[0]
    y = pos[1]

    if x >= 0.0 and x <= rect_w and y >= y_min and y <= y_max:
        return 1.0

    if y >= y_lo and y <= y_hi:
        if x >= x_min and x <= 0.0:
            return 1.0
        if x >= rect_w and x <= x_max:
            return 1.0

    return 0.0


@wp.func
def parabolic_shape(y: float, y_c: float, r: float):
    d = (y - y_c) / r
    return wp.max(0.0, 1.0 - d * d)


@wp.func
def u_boundary_value(
    x: wp.vec2,
    x_min: float,
    x_max: float,
    y_c: float,
    r: float,
    tol: float,
    u_max_signed: float,
):
    if x[0] <= x_min + tol:
        return wp.vec2(u_max_signed * parabolic_shape(x[1], y_c, r), 0.0)
    if x[0] >= x_max - tol:
        return wp.vec2(u_max_signed * parabolic_shape(x[1], y_c, r), 0.0)
    return wp.vec2(0.0, 0.0)


@fem.integrand
def mass_form(s: fem.Sample, u: fem.Field, v: fem.Field):
    return wp.dot(u(s), v(s))


@fem.integrand
def inertia_form(s: fem.Sample, u: fem.Field, v: fem.Field, dt: float):
    return mass_form(s, u, v) / dt


@fem.integrand
def viscosity_form(s: fem.Sample, u: fem.Field, v: fem.Field, nu: float):
    return 2.0 * nu * wp.ddot(fem.D(u, s), fem.D(v, s))


@fem.integrand
def viscosity_and_inertia_form(s: fem.Sample, u: fem.Field, v: fem.Field, dt: float, nu: float):
    return inertia_form(s, u, v, dt) + viscosity_form(s, u, v, nu)


@fem.integrand
def transported_inertia_form(s: fem.Sample, domain: fem.Domain, u: fem.Field, v: fem.Field, dt: float):
    pos = domain(s)
    vel = u(s)

    conv_pos = pos - 0.5 * vel * dt
    conv_s = fem.lookup(domain, conv_pos, s)
    conv_vel = u(conv_s)

    conv_pos = conv_pos - 0.5 * conv_vel * dt
    conv_vel = u(fem.lookup(domain, conv_pos, conv_s))

    return wp.dot(conv_vel, v(s)) / dt


@fem.integrand
def div_form(s: fem.Sample, u: fem.Field, q: fem.Field):
    return -q(s) * fem.div(u, s)


@fem.integrand
def vorticity_expr(s: fem.Sample, u: fem.Field):
    return fem.curl(u, s)


def q_signal(t, u_max_amp):
    """Bipolar square wave: +u_max_amp for the first half-period, -u_max_amp for the second."""
    phase = (t % PERIOD) / PERIOD
    return u_max_amp if phase < 0.5 else -u_max_amp


def node_grid_indices(pos, x_min, y_min, dx_node, dy_node):
    i = np.round((pos[:, 0] - x_min) / dx_node).astype(int)
    j = np.round((pos[:, 1] - y_min) / dy_node).astype(int)
    return i, j


class Simulation:
    def __init__(self):
        self.geo = fem.Grid2D(res=wp.vec2i(NX, NY), bounds_lo=wp.vec2(X_MIN, Y_MIN), bounds_hi=wp.vec2(X_MAX, Y_MAX))

        cell_space = fem.make_polynomial_space(self.geo, degree=0)
        activity = cell_space.make_field()
        fem.interpolate(
            cell_activity,
            dest=activity,
            values={
                "x_min": X_MIN, "x_max": X_MAX, "y_min": Y_MIN, "y_max": Y_MAX,
                "rect_w": RECT_W, "y_lo": Y_LO, "y_hi": Y_HI,
            },
        )
        mask = wp.array(activity.dof_values.numpy(), dtype=int)
        self.active_partition = fem.ExplicitGeometryPartition(self.geo, mask)

        self.u_space = fem.make_polynomial_space(self.geo, degree=2, dtype=wp.vec2)
        self.p_space = fem.make_polynomial_space(self.geo, degree=1)

        self.u_space_partition = fem.make_space_partition(
            space_topology=self.u_space.topology, geometry_partition=self.active_partition
        )
        self.p_space_partition = fem.make_space_partition(
            space_topology=self.p_space.topology, geometry_partition=self.active_partition
        )

        domain = fem.Cells(geometry=self.active_partition)
        self.boundary = fem.BoundarySides(self.active_partition)

        u_test = fem.make_test(space=self.u_space, space_partition=self.u_space_partition, domain=domain)
        u_trial = fem.make_trial(space=self.u_space, space_partition=self.u_space_partition, domain=domain)
        p_test = fem.make_test(space=self.p_space, space_partition=self.p_space_partition, domain=domain)

        self.u_test = u_test
        self.u_bd_test = fem.make_test(space=self.u_space, space_partition=self.u_space_partition, domain=self.boundary)
        u_bd_trial = fem.make_trial(space=self.u_space, space_partition=self.u_space_partition, domain=self.boundary)

        # ---- one-time assembly --------------------------------------------
        self.u_matrix_orig = fem.integrate(
            viscosity_and_inertia_form, fields={"u": u_trial, "v": u_test}, values={"dt": DT, "nu": NU},
            output_dtype=wp.float64,
        )
        u_matrix_lhs = bsr_copy(self.u_matrix_orig)

        self.div_matrix_orig = fem.integrate(div_form, fields={"u": u_trial, "q": p_test}, output_dtype=wp.float64)
        div_matrix_lhs = bsr_copy(self.div_matrix_orig)

        # Raw (never-normalized) boundary mass matrix, kept pristine: each time the BC
        # value changes we take a fresh copy and normalize it together with that value
        # in a single call (normalizing an already-idempotent projector a second time
        # with a new value does not rescale the value correctly).
        self.u_bd_projector_raw = fem.integrate(
            mass_form, fields={"u": u_bd_trial, "v": self.u_bd_test}, assembly="nodal", output_dtype=wp.float64,
        )

        u_bd_projector_struct = bsr_copy(self.u_bd_projector_raw)
        fem.normalize_dirichlet_projector(u_bd_projector_struct)
        fem.project_system_matrix(u_matrix_lhs, u_bd_projector_struct)
        div_matrix_lhs -= self.div_matrix_orig @ u_bd_projector_struct

        self.saddle_system = fem_example_utils.SaddleSystem(u_matrix_lhs, div_matrix_lhs)

        self.u_field = self.u_space.make_field(space_partition=self.u_space_partition)
        self.p_field = self.p_space.make_field(space_partition=self.p_space_partition)
        self.vort_field = self.p_space.make_field(space_partition=self.p_space_partition)

        # Cache node positions (grid (i,j) index + physical position) for fast
        # scatter into a regular array during post-processing / plotting.
        u_node_idx = self.u_space_partition.space_node_indices().numpy()
        self.u_pos = self.u_space.node_positions().numpy()[u_node_idx]
        p_node_idx = self.p_space_partition.space_node_indices().numpy()
        self.p_pos = self.p_space.node_positions().numpy()[p_node_idx]

        # Boolean masks selecting nodes strictly inside the rectangle (excluding the
        # pipe stubs), used for computing "inside the main rectangle" metrics directly
        # from dof_values without going through the full NaN-padded plotting grid.
        self.u_in_rect = (self.u_pos[:, 0] >= 0.0) & (self.u_pos[:, 0] <= RECT_W)
        self.p_in_rect = (self.p_pos[:, 0] >= 0.0) & (self.p_pos[:, 0] <= RECT_W)

    def reset_fields(self):
        self.u_field.dof_values.zero_()
        self.p_field.dof_values.zero_()
        self.vort_field.dof_values.zero_()

    def rectangle_metrics(self):
        """Summary statistics of the current field state, restricted to the rectangle interior."""
        umag = np.linalg.norm(self.u_field.dof_values.numpy()[self.u_in_rect], axis=-1)
        p = self.p_field.dof_values.numpy()[self.p_in_rect]
        vort = self.vort_field.dof_values.numpy()[self.p_in_rect]
        return {
            "u_peak": float(np.max(umag)),
            "p_peak": float(np.max(np.abs(p))),
            "p_rms": float(np.sqrt(np.mean(p ** 2))),
            "vort_peak": float(np.max(np.abs(vort))),
            "vort_rms": float(np.sqrt(np.mean(vort ** 2))),
        }

    def bc_rhs(self, t, u_max_amp):
        u_signed = q_signal(t, u_max_amp)

        u_bd_field = fem.ImplicitField(
            domain=self.boundary, func=u_boundary_value,
            values={"x_min": X_MIN, "x_max": X_MAX, "y_c": Y_C, "r": PIPE_R, "tol": TOL, "u_max_signed": u_signed},
        )
        u_bd_value = fem.integrate(
            mass_form, fields={"u": u_bd_field, "v": self.u_bd_test}, assembly="nodal", output_dtype=wp.vec2d,
        )

        u_bd_projector = bsr_copy(self.u_bd_projector_raw)
        fem.normalize_dirichlet_projector(u_bd_projector, u_bd_value)

        u_bd_rhs = wp.zeros_like(u_bd_value)
        fem.project_system_rhs(self.u_matrix_orig, u_bd_rhs, u_bd_projector, u_bd_value)

        div_bd_rhs = -(self.div_matrix_orig @ u_bd_value)

        return u_bd_projector, u_bd_rhs, div_bd_rhs

    def step(self, t, u_max_amp):
        u_bd_projector, u_bd_rhs, div_bd_rhs = self.bc_rhs(t, u_max_amp)

        u_rhs = fem.integrate(
            transported_inertia_form, fields={"u": self.u_field, "v": self.u_test}, values={"dt": DT},
            output_dtype=wp.vec2d,
        )
        bsr_mv(u_bd_projector, x=u_rhs, y=u_rhs, alpha=-1.0, beta=1.0)
        array_axpy(x=u_bd_rhs, y=u_rhs, alpha=1.0, beta=1.0)

        p_rhs = div_bd_rhs

        x_u = wp.empty_like(u_rhs)
        x_p = wp.empty_like(p_rhs)
        array_cast(out_array=x_u, in_array=self.u_field.dof_values)
        array_cast(out_array=x_p, in_array=self.p_field.dof_values)

        err, niter = fem_example_utils.bsr_solve_saddle(
            saddle_system=self.saddle_system, tol=1.0e-6, x_u=x_u, x_p=x_p, b_u=u_rhs, b_p=p_rhs, quiet=True,
        )

        array_cast(in_array=x_u, out_array=self.u_field.dof_values)
        array_cast(in_array=x_p, out_array=self.p_field.dof_values)
        return err, niter

    def compute_vorticity(self):
        fem.interpolate(vorticity_expr, dest=self.vort_field, fields={"u": self.u_field})

    def to_grid_u(self, values):
        i, j = node_grid_indices(self.u_pos, X_MIN, Y_MIN, DX / 2.0, DY / 2.0)
        grid = np.full((2 * NY + 1, 2 * NX + 1, 2), np.nan, dtype=np.float32)
        grid[j, i] = values
        return grid

    def to_grid_p(self, values):
        i, j = node_grid_indices(self.p_pos, X_MIN, Y_MIN, DX, DY)
        grid = np.full((NY + 1, NX + 1), np.nan, dtype=np.float32)
        grid[j, i] = values
        return grid


def run(q_amp_ml_min=Q_AMP_ML_MIN, n_periods=N_PERIODS, frames_per_period=60, quiet=False, sim=None):
    wp.init()
    if sim is None:
        sim = Simulation()
    sim.reset_fields()

    u_max_amp = q_ml_min_to_u_max(q_amp_ml_min)

    n_steps = int(round(n_periods * PERIOD / DT))
    save_stride = max(1, int(round((PERIOD / DT) / frames_per_period)))

    frames_t = []
    frames_u = []
    frames_p = []
    frames_vort = []

    t = 0.0
    t0 = time.perf_counter()
    for step in range(n_steps):
        sim.step(t, u_max_amp)
        t += DT

        # Only record frames during the LAST period (once the startup transient has damped out)
        if step >= n_steps - int(round(PERIOD / DT)) and step % save_stride == 0:
            sim.compute_vorticity()
            frames_t.append(t)
            frames_u.append(sim.to_grid_u(sim.u_field.dof_values.numpy()))
            frames_p.append(sim.to_grid_p(sim.p_field.dof_values.numpy()))
            frames_vort.append(sim.to_grid_p(sim.vort_field.dof_values.numpy()))

        if not quiet and (step % 100 == 0 or step == n_steps - 1):
            umax = np.max(np.linalg.norm(np.nan_to_num(sim.to_grid_u(sim.u_field.dof_values.numpy())), axis=-1))
            print(f"step {step}/{n_steps}  t={t:.3f}s  |u|max={umax:.2f} mm/s  "
                  f"elapsed={time.perf_counter()-t0:.1f}s")

    np.savez(
        RESULTS_PATH,
        frames_t=np.array(frames_t),
        frames_u=np.array(frames_u),
        frames_p=np.array(frames_p),
        frames_vort=np.array(frames_vort),
        x_min=X_MIN, x_max=X_MAX, y_min=Y_MIN, y_max=Y_MAX,
        dx=DX, dy=DY, nx=NX, ny=NY,
        rect_w=RECT_W, rect_h=RECT_H, pipe_d=PIPE_D, y_c=Y_C,
        period=PERIOD, q_amp_ml_min=q_amp_ml_min, u_max_amp=u_max_amp, u_mean_amp=u_max_amp / 1.5,
        nu=NU, depth=DEPTH,
    )
    print(f"saved results to {RESULTS_PATH}")
    print(f"total runtime: {time.perf_counter()-t0:.1f}s")
    return sim


def run_sweep(q_values_ml_min, sim=None, n_periods=N_PERIODS, steady_margin_frac=0.15, quiet=False):
    """Run the transient sim once per flow-rate amplitude in q_values_ml_min (reusing the
    same mesh/matrices across all of them) and record steady-state (reversal-transient
    excluded) vorticity/velocity/pressure summary statistics inside the rectangle.
    """
    wp.init()
    if sim is None:
        sim = Simulation()

    steps_per_period = int(round(PERIOD / DT))
    n_steps = int(round(n_periods * PERIOD / DT))
    last_period_start_step = n_steps - steps_per_period

    q_out, vort_peak, vort_rms, u_peak, p_peak, p_rms, re_out = [], [], [], [], [], [], []

    t0 = time.perf_counter()
    for qi, q_ml_min in enumerate(q_values_ml_min):
        sim.reset_fields()
        u_max_amp = q_ml_min_to_u_max(q_ml_min)
        u_mean_amp = u_max_amp / 1.5

        steady_metrics = []
        t = 0.0
        for step in range(n_steps):
            sim.step(t, u_max_amp)
            t += DT

            if step >= last_period_start_step:
                rel = (step - last_period_start_step) / steps_per_period
                is_steady = not (rel < steady_margin_frac or abs(rel - 0.5) < steady_margin_frac)
                if is_steady:
                    sim.compute_vorticity()
                    steady_metrics.append(sim.rectangle_metrics())

        q_out.append(q_ml_min)
        vort_peak.append(max(m["vort_peak"] for m in steady_metrics))
        vort_rms.append(float(np.mean([m["vort_rms"] for m in steady_metrics])))
        u_peak.append(max(m["u_peak"] for m in steady_metrics))
        p_peak.append(max(m["p_peak"] for m in steady_metrics))
        p_rms.append(float(np.mean([m["p_rms"] for m in steady_metrics])))
        re_out.append(u_mean_amp * PIPE_D / NU)

        if not quiet:
            print(f"[{qi+1}/{len(q_values_ml_min)}] Q={q_ml_min:.1f} mL/min  Re={re_out[-1]:.0f}  "
                  f"vort_peak={vort_peak[-1]:.2f} 1/s  vort_rms={vort_rms[-1]:.2f} 1/s  "
                  f"({len(steady_metrics)} steady samples)  elapsed={time.perf_counter()-t0:.1f}s")

    results = dict(
        q_ml_min=np.array(q_out), vort_peak=np.array(vort_peak), vort_rms=np.array(vort_rms),
        u_peak=np.array(u_peak), p_peak=np.array(p_peak), p_rms=np.array(p_rms), re=np.array(re_out),
        pipe_d=PIPE_D, depth=DEPTH, nu=NU, period=PERIOD,
    )
    np.savez(SWEEP_RESULTS_PATH, **results)
    print(f"saved sweep results to {SWEEP_RESULTS_PATH}")
    print(f"total sweep runtime: {time.perf_counter()-t0:.1f}s")
    return results


if __name__ == "__main__":
    run()
