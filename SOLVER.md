# Solver

The installed native Solver/CLI is the canonical Python-free runtime. A solve is
a versioned `dynsyslab.compute-job/v1` document with `"kind": "solve"`, planned
and run by the compiled workflow:

```powershell
dynsyslab-cli validate examples/native/duffing-free-solve.json
dynsyslab-cli plan examples/native/duffing-free-solve.json
dynsyslab-cli run examples/native/duffing-free-solve.json --archive-root C:\data\dynsyslab-archives
```

The native Solver panel edits the same document. A completed solve is saved
under `archives/solutions/<system-name>/`; no data are written before
completion by default.

## Methods

`solver.method` selects `Euler`, `RK4`, adaptive Dormand--Prince `RK45`
(alias `dopri5`), or SUNDIALS/CVODE `BDF` (alias `ode15s`) for stiff systems.
The fixed-step methods require `step_size`; adaptive methods take `rtol`,
`atol`, `max_step`, and `max_internal_steps`. `time.output_step` sets the
archived sample spacing. Every result records the method and settings in its
archive metadata.

## Control

The controller workflow exposes separate state-space (`lqr`, `sdre`, `pid`,
`mpc`) and Koopman (`klqr`, `kmpc`) families, configured under the job's
`controller` section. The system's built-in driver, such as RPE acoustic
pressure or Duffing forcing, remains part of its parameters.

Controller evaluation is an explicit solver policy. Set
`solver.control_evaluation` to `"sampled"` with a `control_sample_time` for an
integer-tick zero-order hold. Use `"continuous_rhs"` only for a stateless
RHS-safe law; it recomputes the estimator, state-dependent input map, and gain
at every internal solver evaluation, including rejected adaptive trial states.
KLQR uses this mode to recompute its state-dependent `B_Psi(x)` and Riccati
gain. PID, MPC, and stateful observers are rejected in this mode.
`"automatic"` is the compatibility default: a supplied sample time selects
sampled control and its absence selects continuous RHS evaluation. Scientific
runs should declare the policy rather than rely on automatic mode.

The solution-aligned control history has accepted solver times. In continuous
RHS mode it is evaluated on accepted states; `control_call_time` and
`control_call_values` separately preserve the internal evaluation tape. In
sampled mode the call log contains only actual controller updates, not every
reuse of the held command.

## Failures

A trajectory that fails (leaves the system's validity domain, underflows the
step size, or meets a controller error) aborts the whole native solve by
default. With `"solver": {"on_failure": "truncate"}` the solve instead keeps
each failed trajectory's computed samples, archives only those rows, records
`integration_terminated` (`completed_samples`, `time`, `reason`) in the
archive metadata, and continues with the next trajectory. Use it to show where
a controller or model breaks down; the Viewer draws the computed part.

## Other system families

PDE systems are solved by the method of lines and save a time sequence of
fields. Delay systems use their own history-aware integrator. Lattice systems
(Ising, Potts) advance by sweeps and write sweep-indexed observables.

Post-processing is a separate job (`"kind": "postprocess"` or an experiment
stage) that reads a completed archive; the solver has no dependency on it. See
[the Solver and Post-processor guide](docs/SOLVER_POSTPROCESSOR.md).
