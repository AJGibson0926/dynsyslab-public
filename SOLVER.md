# Standalone solver

The installed native Solver/CLI is the canonical Python-free runtime. Its
controller workflow exposes separate state-space (`lqr`, `sdre`, `pid`,
`mpc`) and Koopman (`klqr`, `kmpc`) families. In particular, KLQR retains
continuous adaptive-RHS evaluation of its state-dependent `B_Psi(x)` and
Riccati gain. The Python API described below is the retained reference and
optional convenience layer.

Use `solver.solve(SolveRequest(...))`. A completed solve is always saved under
`archives/solutions/<system-name>/`; no data are written before completion by default.

```python
from dynsyslib import load_system
from solver import solve, SolveRequest, SolverOptions, HistoryOptions

system = load_system("duffing", frequency_hz=1.0, amplitude=0.3)
result = solve(SolveRequest(
    system=system, kind="ode", time_span=(0.0, 500.0),
    options=SolverOptions(method="DOP853", rtol=1e-10),
    history=HistoryOptions(track_during_solve=True, checkpoint_interval=10.0),
    label="duffing-reference",
))
```

ODE methods: `RK45`, `DOP853`, `Radau`, `BDF`, `LSODA`, `Euler`, and `RK4`.
The fixed-step methods require `step_size`.  `Radau`, `BDF`, and `LSODA` are
the available stiff pathways.

`SolverOptions.execution_backend` controls the numerical implementation:

- `auto` selects an admitted compatible native kernel when the optional wheel
  is installed and otherwise uses the Python/SciPy reference;
- `python` always uses the reference pathway;
- `native` requires native support and fails clearly if the method, model,
  driver, control/history policy, or installed extension is incompatible.

The admitted native pathway supports complete Euler/RK4 loops for dimensional
Duffing/van der Pol and dimensionless Rayleigh--Plesset, with zero or
sinusoidal forcing, including one-call hierarchy cluster batches. RPE keeps
its model coordinates (`tau`, `r`, `dr/dtau`) inside the kernel, evaluates
sinusoidal driver phase in physical seconds, and returns the normal public
result in dimensional coordinates. The optional official SUNDIALS adapter provides explicit
CVODE Adams/BDF execution on Python 3.12+, but is not auto-admitted yet.
The 2-D Ising model has an explicit native kernel whose independent PRNG seed
and engine are recorded; `auto` retains the established NumPy random stream.
DDE and generic PDE systems stay on reference implementations until a
model-specific kernel can cross the admission gates. Every result records the
requested and selected backend capabilities in metadata.

By default, the solver preserves each method's native time samples and does no
interpolation. To resample a continuous ODE, DDE, or PDE solution after it has
finished, pass `InterpolationOptions`. PCHIP is the default: an efficient,
shape-preserving cubic method. `linear` and `nearest` are also available.

```python
from solver import InterpolationOptions

# Uniform output at 0.01 time-unit spacing. The raw and final trajectories are
# both saved in the same solution archive as raw_time/raw_values and time/values.
interpolation=InterpolationOptions(enabled=True, time_step=0.01, method="pchip")

# An explicit grid (or the legacy SolveRequest(sample_times=...)) always requests resampling.
interpolation=InterpolationOptions(times=np.linspace(0.0, 100.0, 10001), method="linear")
```

The archive metadata records `interpolated`, `interpolation`,
`uniform_time_steps` for the final series, and (when relevant)
`raw_uniform_time_steps`. Discrete lattice trajectories cannot be interpolated.

`control` may be a callable `f(t, state) -> vector` or an object with
`value(t, state) -> vector`; it is added to the system's RHS. The returned
vector must have one element per state variable. The system's built-in driver,
such as RPE acoustic pressure or Duffing forcing, remains part of its parameters.

Controller evaluation is an explicit solver policy. Use
`SolverOptions(control_evaluation="sampled", control_sample_time=...)` for an
integer-tick zero-order hold. Use `control_evaluation="continuous_rhs"` only
for a stateless RHS-safe law; it recomputes the estimator, state-dependent input
map, and gain at every internal solver evaluation, including rejected adaptive
trial states. PID, MPC, and stateful Luenberger/LQE/KLQE observers are rejected
in this mode. `"automatic"` is the compatibility default: a supplied sample
time selects sampled control and its absence selects continuous RHS evaluation.
Scientific runs should declare the policy rather than rely on automatic mode.

A trajectory that fails (leaves the system's validity domain, underflows the
step size, or meets a controller error) aborts the whole native solve by
default. With `"solver": {"on_failure": "truncate"}` the solve instead keeps
each failed trajectory's computed samples, archives only those rows, records
`integration_terminated` (`completed_samples`, `time`, `reason`) in the
archive metadata, and continues with the next trajectory. Use it to show where
a controller or model breaks down; the Viewer draws the computed part.

The solution-aligned control history has accepted solver times. In continuous
RHS mode it is evaluated on accepted states; `control_call_time` and
`control_call_values` separately preserve the internal evaluation tape. In
sampled mode the call log contains only actual controller updates, not every
reuse of the held command.

PDE systems must provide `initial_field(mesh)` and `rhs(t, field, mesh)`; the
PDE pathway uses method of lines, supports adaptive SciPy methods plus
fixed-step Euler/RK4, and saves a time sequence of fields. DDE systems
provide their own history-aware integrator. Lattice systems provide `run(steps)`;
Ising and Potts use this pathway and write sweep-indexed observables.

With `track_during_solve=True`, checkpoint arrays and a manifest containing the
system parameters, mesh/domain, time span, control-independent solver settings,
and method are written in `archives/history/<system>/`. This survives a solver
failure. Without it, only the final completed solution archive is written.

Pass a callback in `postprocessors` to analyze a completed archive without
creating a solver dependency on post-processing. For example:

```python
from postprocessing.koopman import dmd_processor
# add postprocessors=(dmd_processor(rank=10),) to SolveRequest
```
