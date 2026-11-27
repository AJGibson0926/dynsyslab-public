# Native solver

The native ODE, delay, PDE and lattice solves: integration methods, the `solve` job settings,
designed inputs and initial conditions, guards, resampling, checkpoints, pause, resume and
continuation, dimensional integration, preset metadata and the extension entry points. Every
piece lives in the C++ core (`native/`) or the job pipeline (`workflow/native/`); Python only
submits jobs.

## Integration methods

`core::integrateControlled` (`native/include/dynsyslab/core/Integrator.hpp`) is one driver for
every method. It owns the output clock, the zero-order-hold control clock, the guards, the step
budget, checkpoints and the failure policy. A *stepper* only advances a state across one interval
(`native/src/IntegratorBackend.hpp`).

| `solver.method` | Method | Source |
|---|---|---|
| `Euler`, `RK4` | one fixed step per output (or control-clock) interval | classical |
| `RK23` (`bs23`, `ode23`) | embedded Runge-Kutta 3(2), FSAL | Bogacki & Shampine, *Appl. Math. Lett.* 2 (1989) 321-325 |
| `RK45` (`dopri5`, `ode45`) | embedded Runge-Kutta 5(4) | Dormand & Prince, *J. Comput. Appl. Math.* 6 (1980) 19-26 |
| `DOP853` | Runge-Kutta 8(5,3) with the DOP853 error estimator | Hairer, Norsett & Wanner, *Solving ODEs I*, 2nd ed. (1993), Sec. II.10 |
| `Radau` (`radau5`) | Radau IIA order 5, simplified Newton on the transformed real/complex systems | Hairer & Wanner, *Solving ODEs II*, 2nd ed. (1996), Sec. IV.8 |
| `BDF` (`ode15s`) | SUNDIALS CVODE variable-order BDF | Hindmarsh et al., *ACM TOMS* 31 (2005) 363-396 |
| `Adams` (`ode113`) | SUNDIALS CVODE variable-order Adams-Moulton | Hindmarsh et al. 2005 |
| `LSODA` | automatic Adams/BDF switching | after Petzold, *SIAM J. Sci. Stat. Comput.* 4 (1983) 136-148 |

The RK23 and DOP853 tableaus (`native/src/IntegratorTableaus.hpp`) are the published constants.
RK45 keeps the historical max-norm step controller, so earlier archives reproduce bit for bit;
RK23 and DOP853 use the RMS-norm controller of Hairer, Norsett & Wanner (1993), Sec. II.4.
Radau's Jacobian is a forward difference of the controlled right-hand side, refreshed when Newton
converges slowly.

**LSODA.** At each interval end the stepper compares the last step h with the explicit stability
limit 1/||J||_inf. Adams switches to BDF after three intervals whose step is stability-limited
(h ||J|| > 0.5), and BDF returns to Adams once h ||J|| < 0.5. Switches are counted in
`IntegrationStatistics::methodSwitches`. The test is simpler than ODEPACK's step-size comparison,
so step counts differ from SciPy's LSODA, but the solution agrees within tolerance.

CVODE stops exactly at each interval end (`CVodeSetStopTime`), so a held input never leaks across
a sampled-control tick. BDF, Adams and LSODA need a build with `DYNSYSLAB_ENABLE_SUNDIALS`.

## `solve` job settings

| Key | Meaning |
|---|---|
| `time.start`, `time.stop`, `time.output_step` | uniform output clock |
| `time.output_times` | arbitrary, strictly increasing output clock, replacing `output_step`. The integrator stops exactly on these times, with no interpolation |
| `time.coordinate_system` | `automatic` (the system's own), or `dimensional` for a dimensionless model that declares physical scales |
| `time.sample_times` or `interpolation {method, times \| time_step}` | post-solve resampling (below) |
| `solver.method`, `rtol`, `atol`, `max_step`, `step_size` | method and tolerances |
| `solver.max_internal_steps` | step budget per trajectory (accepted plus rejected steps) |
| `solver.state_bounds` | one `[low, high]` (or `null`) per state, checked after every accepted step |
| `solver.max_state_magnitude` | infinity-norm guard |
| `solver.on_failure` | `error` (abort), `truncate` (archive up to the failure, `integration_terminated`), `reject` (drop the trajectory; see the outcome record) |
| `history.checkpoint_interval` | model time between checkpoint (partial) archives |
| `history.pause_file` | pause when this file exists at a checkpoint |
| `history.checkpoint_wall_seconds` | wall-clock seconds between checkpoint archives |
| `history.deadline_unix_seconds` | checkpoint and pause at this UTC time (set by the queue deadline policy, JOB_QUEUE.md) |
| `history.controller_call_log` | archive every controller evaluation per trajectory: `controller_call_time`, `_state`, `_input`, `_kind` (0 right-hand side, 1 sampled tick, 2 output sample), `_status` (1 accepted step, 0 rejected trial step, -1 unknown, e.g. CVODE); `history.controller_call_log_limit` caps the entries (`core/ControllerCallLog.hpp`) |
| `history.keep_partial_history_on_failure` | with `on_failure: error`, the solve still fails but first writes each started trajectory's computed samples to `<archive>.partial.npz` (`kept_on_failure: true`, `termination_reason`) |
| `system.driver` | `none` (default), `preset` (the preset's recommended driver), or `{amplitude, frequency_hz, phase_rad}` |

Every archive records `integration_method`, `integration_failure_policy` and
`control_input_dimension`. When a guard rejects a trajectory, or a pause stops the solve, the
solve writes `<manifest>.outcome.json` (`dynsyslab.solve-outcome/v1`). It holds the accepted-member
mask, each rejected trajectory's reason and, after a pause, `resume_jobs`.

### Post-solve resampling

`time.sample_times` (or `interpolation.times`, or `interpolation.time_step`) resamples every
output-grid array after the solve. `interpolation.method` is `pchip` (the default), `linear`,
`cubic` or `nearest`. The computed arrays are kept as `raw_<name>` plus `raw_time`, and the
archive records `interpolation {method, raw_samples, output_samples}`. A truncated trajectory gives
NaN after its last computed time.

The interpolation lives in `core::resample` (`TimeSeriesOperations.hpp`):

- **pchip** is the monotone piecewise-cubic Hermite interpolant of Fritsch & Carlson (*SIAM J.
  Numer. Anal.* 17 (1980) 238), with the weighted-harmonic-mean slopes of Fritsch & Butland (*SIAM
  J. Sci. Stat. Comput.* 5 (1984) 300) and shape-preserving end slopes (Moler, *Numerical Computing
  with MATLAB*, 2004, Sec. 3.4). It never overshoots the data.
- **nearest** takes the nearest sample; a tie goes to the earlier one.
- **cubic** uses Catmull-Rom slopes.

The `coordinate_transform` kind `resample` (`workflow/native/src/CoordinateResample.cpp`)
applies the same operation to any post-processing source, a physical record or earlier
coordinates. It takes `time_step` (with optional `start` and `stop`) or `times`, and `method`.
It publishes a time-series record with a `resampling {method, raw_samples, output_samples}`
block, and that record names its source. A resampled physical record stays a physical record.
This replaces the Python in-memory `postprocessing/resampling.py`; the Python fixed-lag methods
now point non-uniform sources to it.

### Checkpoints, pause, resume and continuation

Archives are immutable, so a continuation is always a new child archive.

- **Checkpoint.** The running trajectory's samples so far are written atomically to
  `<archive>.partial.npz` (`dynsyslab.solution-checkpoint/v1`). These partial archives are removed
  once the solve is published.
- **Pause.** If `history.pause_file` exists at a checkpoint, the integrator stops after that sample
  (`NativeSolution::paused`). The paused trajectory is archived with `integration_paused`; later
  trajectories are not started. The outcome record's `resume_jobs` holds a continuation job for the
  paused archive and a solve of the members that had not started, with their exact initial states.
  Running those jobs resumes the solve.
- **Continuation and extension.** `initial_conditions {mode: "continue_from", archive | archives}`
  starts each trajectory at its parent's last computed state; `time.start` must equal the parent's
  last time. The child archive records `continuation {parent_archive, parent_archive_sha256,
  parent_time, parent_completed_samples}`. `detail::extensionJob(parent, stop)`
  (`SolveContinuation.hpp`) builds the whole job from the parent's archived system, controller,
  solver and time settings.

Delay (Mackey-Glass) and lattice solves take `history.checkpoint_interval` too, in model time or
in sweeps. The samples so far go to `<archive>.partial.npz`, and for a lattice the current sites
go there as well (`lattice_current`). As in the Python solver, these runs have checkpoints but no
pause (`DelayIntegrationOptions::checkpoint`, `LatticeCheckpoint`).

A restart begins a new integration from the checkpointed state. It is exact for fixed-step
methods and memoryless controllers. An adaptive method restarts its step-size history, and a
stateful controller (an estimator, an integrator state) restarts from its initial state.

## Designed inputs (`controller.kind`)

| Kind | Settings | Archived `input_design` |
|---|---|---|
| `tabulated` | `times`, `values` (one shared row or one row per trajectory), `interpolation` `linear` or `previous`, `time_units`/`input_units` `model` or `physical`. End values are held | - |
| `command_tape` | `interval`, `start`, `values` or `file` (+ `array`) of an `.npz`, `layout` `trajectory_major` or `transition_major` (MATLAB `Ubig`), `scale` (for example 0.5 for the Korda-Mezic source convention) | `tape_sha256`, `scale`, `interval`, `transitions`, `command_peak` |
| `telegraph_multisine` | `seed`, `interval`, `transition_count`, `amplitudes`, `hold_steps`, `levels` (default ±0.5, ±1), `multisine_fraction`, `angular_frequencies` | `tape_sha256`, `hold_steps`, `amplitude`, `command_peak` |
| `symmetric_pulse` | `magnitude`. Trajectory 2k receives -magnitude and 2k+1 receives +magnitude | `pair_id`, `polarity`, `command` |
| `pulse_branches` | `amplitudes` (B branches), `widths` (one value or one per branch, model time), `start`, `input_units`. Trajectory k is branch k mod B of anchor k / B and holds its amplitude for t - start < width, then zero. Pair it with `initial_conditions.repeat: B` and an anchor mode such as `archive_states`, and solve with `solver.control_evaluation: "sampled"` so each pulse is held over whole samples. Branches +a and -a of equal width form a symmetric pair | `anchor`, `branch`, `command`, `width`, `branches`, and for paired branches `pair_id`, `polarity` |

All four accept `active_time_window`.

The telegraph-plus-multisine tape is a balanced multilevel random telegraph with per-member
amplitude and hold, blended with a normalised random-phase multisine, then made zero mean and
peak clipped (Ljung, *System Identification*, 2nd ed. 1999, Sec. 13.3; Pintelon & Schoukens,
*System Identification*, 2nd ed. 2012, Ch. 2). Each member draws from its own `std::mt19937_64`
stream seeded by `(seed, member)`, so a tape does not depend on how many trajectories the job has.
The numbers differ from the NumPy generator the Python code used.

`operations/native_controlled_training.py` writes the standard jobs: the Korda-Mezic Duffing
command tape, the RPE persistent excitation and the RPE symmetric pulses. It runs them through
`dynsyslab_client` and reads the archives back as a `ControlledDataset`. The Python generators
it replaces are gone.

Paired ±u first responses give data-driven input-map samples; the `pair_id` and `polarity` records
are what the phase-E `input_map` estimator consumes. Command tapes for the Duffing case follow Korda
& Mezic, *Automatica* 93 (2018) 149-160.

**Actuation.** `settings.actuation: "additive_state"` with an n x 1 `input_matrix` (or the explicit
LQR `B`) adds B u straight to the state derivative instead of driving the system's input channel.
It works with the prescribed kinds, `pid`, and `lqr` with explicit `A` and `B`. In the core,
`ControlPolicy` carries a vector law, its `inputDimension` and an n x m `additiveMatrix`; the
solution stores m control columns.

## Initial conditions

The modes are `default`, `explicit`, `uniform`, `normal`, `sphere`, `ellipsoid` and the following:

- `stratified`: a Latin hypercube. Each state's range is cut into `count` strata, each stratum is
  used once, the strata are permuted independently per state, and each point is jittered by
  U(`jitter`) inside its stratum (McKay, Beckman & Conover, *Technometrics* 21 (1979) 239-245).
- `target_tube`: points along a multisine target, with x = offset + sum a_k sin(w_k tau + phi_k)
  and its derivative, plus Gaussian `jitter` per state.
- `concatenate`: `parts`, each a full IC spec with its own `count`.
- `continue_from`: see continuation above.
- `archive_states`: `archive` (a trajectory archive of the system's states) and `times`, one per distinct initial state. Each state is the archived sample nearest its time, which must lie within half a sample interval. These are anchors for restarts such as pulse branches.

Two options apply to every mode at the top level:

- `repeat` uses each generated state for that many consecutive trajectories (2 gives the two
  members of a symmetric pulse pair).
- `clip` gives one `[low, high]` or `null` per state.

State rejection gates are `solver.state_bounds` with `on_failure: "reject"`.

## Dimensional integration

`core::systemInCoordinates(model, "dimensional")` (`DimensionalSystem.hpp`) integrates a
dimensionless model in physical units by the exact change of variables
dX/dt = S f(t/T, X/S, U/P) / T. Here T is `secondsPerModelTime`, S is
`physicalStatePerModelState` and P is `physicalInputPerModelInput`. RPE and Marmottant declare
S = (R0, R0/T). Default initial conditions are scaled to match, and archives record
`integration_form: "dimensional"`, with time in seconds and the state in metres and metres per
second. The bubble scales follow Gibson, Calvisi & Yee, *SIAM J. Appl. Dyn. Syst.* 21 (2022),
doi:10.1137/21M1456078.

## Systems and presets

- **Marmottant `shell_law`** (parameter): `0` is linear (Marmottant et al., *J. Acoust. Soc. Am.*
  118 (2005) 3499-3505), `1` is a cubic sigmoid and `2` is Gompertz. The two smoothings come from
  the MATLAB `Bubble.m` reference, which records no source paper. Between the rupture and break
  radii the smoothed laws keep that reference's zero-tension branch.
- **Preset metadata.** `PresetDescriptor` carries the recommended `driver`, the `reference` the
  preset reproduces, and `processingJson`, the per-method settings and provenance
  (`native/src/PresetMetadata.cpp`, transcribed from the Python preset tables). Solves archive
  `preset_reference` and `preset_processing`. The driver applies only when the job asks for it
  with `system.driver: "preset"`.

- **`rayleigh-plesset-dimensional-baseline`.** The frozen pre-nondimensionalisation RPE of the
  retired Python `rpe_dimensional_baseline` (SI units, ambient 101325 Pa, R0 = 10 um, air-water
  20 C, the drive entering as +p_d(t)). It is the `rpe` model with `pressure_input_sign = +1`
  under the dimensional change of variables, not a second copy of the physics
  (`SystemsRPEBaseline.cpp`); initial states are in metres and m/s.

## Partitioned PDE solves

`system.pde_partitioning` (`enabled`, `partition_shape`, `ghost_width`, `local_workers`) splits a
structured mesh into tensor-product subdomains (`core/PDEParallel.hpp`).

- `communication`: `shared_memory` (default; each subdomain reads its ghosted window),
  `halo_exchange` (each subdomain owns its interior and receives face halos rank to rank, axis by
  axis so corners are exact: `core/PDEHaloExchange.hpp`), or `mpi`. MPI is an optional build,
  `-DDYNSYSLAB_WITH_MPI=ON` (`find_package(MPI)`), whose `mpiExchangeRankField` moves the same
  packed halos with `MPI_Isend`/`MPI_Irecv` in an `mpiexec`-launched rank driver; the default build
  and an in-process solve use `halo_exchange`, its thread-parallel equivalent.
- `adaptive` (`enabled`, or `load_balance: "adaptive"`; `rebalance_imbalance`,
  `rebalance_interval_calls`, `checkpoint_before_repartition`, `refinement_threshold`,
  `coarsening_threshold`): partition wall times are measured and, when max/mean exceeds the
  threshold, each axis is recut to equal cumulative cost (Berger & Bokhari 1987; Hendrickson &
  Devine 2000). Repartitioning leaves the numbers unchanged. Archives record
  `pde_communication`, `pde_adaptive_policy` and `pde_repartition_records`
  (`dynsyslab.pde-repartition/v1`, `core/PDEAdaptivePartition.hpp`).

## Delay rollout, adaptive method of lines, extensions

- **Delay rollout.** `core::rolloutDelayModel` (`DelayRollout.hpp`) advances
  x' = F(x(t), x(t - tau_1), ...) with RK4 on its own history: the method of steps (Bellen &
  Zennaro, *Numerical Methods for Delay Differential Equations*, 2003, Ch. 3). Fractional lags and
  half-step stages use linear interpolation; magnitude guards are per component.
- **Adaptive method of lines.** A PDE `solve` whose `solver.method` is not `Euler` or `RK4` runs
  its semi-discrete field through the ODE integrator with that method and its error control
  (Schiesser, *The Numerical Method of Lines*, 1991). Its `pde_method` is `MOL-<method>`.
- **Extension entry points.** `extension_api.h` adds optional delay, field (PDE) and lattice
  tables, read only when `struct_size` covers them, so ODE-only v1 extensions still load.
  `ExtensionRegistry::models()` exposes them. A `solve` whose `system.name` is one of them runs the
  delay model by method of steps (constant prehistory equal to the initial state), the field model
  by the method of lines, or the lattice model sweep by sweep (`ExtensionModelSolve.cpp`).

## Tests

- `native/tests/integrator_methods_tests.cpp`: methods, stiffness, output times, guards, step
  budget, pause, additive actuation, dimensional RPE, PCHIP and nearest, delay rollout, adaptive
  MOL.
- `workflow/native/tests/solve_features_tests.cpp`: the job features end to end, including pause,
  resume and continuation.
- `native/tests/extension_tests.cpp`: the delay, field and lattice entry points.
- `native/tests/pde_parallel_tests.cpp`: halo exchange equals the shared-memory read (edge,
  periodic and zero boundaries, corners), MPI availability, rebalance records, adaptive parity.
- `workflow/native/tests/solve_history_tests.cpp`: controller call log, partial history on
  failure, the dimensional baseline against the frozen Python form, a partitioned PDE job.
