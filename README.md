# dynsyslab

`dynsyslab` provides a plug-in dynamical-system library, single and composite numerical solves, durable
archives, postprocessing, read-only viewing, and separate JobManager and
Campaigner planning/execution layers.

**New here? Read the [DynSysLab overview](docs/OVERVIEW.md)** — what the suite does, its modules and pipelines, and how to analyse and control a dynamical system, with architecture diagrams.

## Using DynSysLab with an AI agent

Ask Claude Code, Codex or another agent: *"Install and run DynSysLab at
https://github.com/AJGibson0926/dynsyslab-public/releases/latest"*. Installing runs
`dynsyslab-cli quickstart`, which prints [the overview](docs/AGENT_WELCOME.md), runs the
bundled examples and reports their results and figures, without opening a window; running
DynSysLab later means the same command. Opening this repository in an agent gives the same
overview on the first session; source conventions are in [AGENTS.md](AGENTS.md).

## Start here

**Install on Windows x64 from [the latest release](https://github.com/AJGibson0926/dynsyslab-public/releases/latest).**
The repository name to give another Codex user is **AJGibson0926/dynsyslab-public**.
Download and run [install.ps1](install.ps1); it verifies the release checksum,
installs both applications, and preserves existing scientific archives. No
Python, compiler, or administrator account is required. See the
[installation and Codex setup guide](docs/INSTALL.md) for the two-command
download, offline installation, and a first calculation.

The primary human front end is the compiled `DynSysLab.exe` Solver and
Post-processor. Solver edits a proposed Job → Batch → System → Cluster →
Trajectory hierarchy; Post-processor and Viewer share the system-centric
archive tree. Developers building from source use:

```powershell
.\tools\Install-DynSysLabNative.ps1 -PinToTaskbar
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\DynSysLab.exe"
```

The optional, dependency-free Python adapter provides
`python -m dynsyslab_client --help` for scripting the compiled installation,
but it is not an application runtime. See
[Start here](docs/START_HERE.md) for the five-minute
human workflow, the AI-agent checklist, and source-tree/installed equivalents.
It passes `validate`, `plan`, `run`, `archive`, and `queue` operations to the
compiled CLI and launches the two native applications directly.
Completed parity checks and case-specific reruns are retained as reports and
immutable artifacts, not permanent commands; see the
[script lifecycle policy](docs/architecture/SCRIPT_LIFECYCLE_POLICY.md).

## Available systems

`rpe`, `marmottant`, `reddy_szeri`, `duffing`, `van_der_pol`,
`nonlinear_pendulum`, `double_pendulum`,
`spring_coupled_pendulums`, `rossler`, `lorenz`, `chua`, `henon_heiles`,
`hindmarsh_rose`, `hodgkin_huxley`, `fput`, and `kuramoto` are selected with the same launcher. Every dynsyslib object
supplies `rhs`, `jacobian`, `initial_state`, `period`, and `invalid`; the
campaign machinery never branches on equations.

Mackey--Glass is in `DELAY_LIBRARY` because a delay equation requires a history
buffer and a different integrator. The 2-D/3-D Ising and Potts models are in `STATISTICAL_LIBRARY`
with Metropolis dynamics; it produces temperature/field phase maps rather than a
deterministic bifurcation diagram or Lyapunov exponent.

Campaign axes are encoded in the common compute-job contract and executed by
the compiled workflow, either directly or through its durable SQLite queue.
`dynsyslab-jobs submit|list|status|cancel|run-next|recover` is an optional thin
client for that queue. The historical `python -m job_manager` scheduler remains
an explicit migration reference for scheduling policies not yet admitted to
the native contract. Numerical data remains exclusively in the archive. A
positive Lyapunov screen remains a candidate, not by itself proof of a chaotic
attractor.

## Adding a system

Add one class in `dynsyslib/lib/` or a separate system module there, inheriting the
matching `ODESystem`, `DelaySystem`, or `LatticeSystem` contract. An ODE plug-in
supplies `rhs`, `jacobian`, `initial_state`, `period`, and `invalid`, then joins
`ODE_LIBRARY`. Attach a `SinusoidalDriver` (or a new driver implementing
`value(t)` and `period`) rather than changing the integrator.

## Standalone solver

The system-agnostic solver is in `solver/`; see `SOLVER.md`. It dispatches ODE,
PDE, delay-equation, and lattice requests and saves every completed result to a
per-system archive. History checkpoints are opt-in.

Numerical execution is owned by the native C++ engine. Euler, RK4, adaptive
Dormand-Prince RK45, and SUNDIALS/CVODE BDF are available through the same
workflow API. Optional Python bindings call this engine; they do not own a
second solver implementation. See [native build and installation](docs/BUILD_NATIVE.md)
and `ADAPTIVE_BACKENDS.md`.

## Independent archives

`archives/` is a separate package and durable data store. Completed solutions
go to `archives/solutions/`; optional in-progress histories go to
`archives/history/`. It is intentionally independent of both the equation
dynsyslib and numerical solver implementation.

## Installation outside Codex

The compiled installer above is the product installation. To add the optional
Python client, create a normal virtual environment and install the wheel or an
editable checkout:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e .
python -m dynsyslab_client doctor
```

The wheel contains only `dynsyslab_client` and has no NumPy or SciPy runtime
dependency. Historical Python research packages remain available only from a
source checkout after installing `requirements/reference-parity.txt` while
parity and retirement work proceeds; they are not distributed as an alternate
scientific engine or advertised as optional product features.

## Post-processing CLI

The optional Python command is a thin client over the compiled workflow ABI.
It accepts the same versioned job document as the native desktop workbench;
planning never performs a decomposition:

```powershell
# Save a job from the native workbench, or start from examples/native.
dynsyslab plan postprocess.json
dynsyslab run postprocess.json --archive-root C:\data\dynsyslab-archives
```

These commands are pass-through conveniences for `dynsyslab-cli plan` and
`dynsyslab-cli run`. See
`examples/native/duffing-hdmd-sindy-postprocess.json` for the standard chained
H-DMD-to-SINDy job. `python -m postprocessing` is retained only as the explicit
Python reference CLI used by migration parity tests.

H-DMD's standard default remains the first scalar state component (radius for
RPE and Marmottant, position for Duffing), but the same workflow accepts any
real scalar entry in the archived system's observable catalog:

Set the H-DMD stage's `settings.observable` to a state-component index or a
named observable such as `far_field_acoustic_echo` or
`inertial_kinetic_energy`.

Use separate immutable jobs with identical matrix, window, rank, sampling,
sorting, and linear-algebra settings to compare spectra from several scalar
observables in Viewer.

The retained Python reference API additionally accepts an `Observable` instance
or an already sampled scalar `TimeSeries` through `observable_series=`. This is
kept for parity experiments involving history-dependent signals while the
compiled workflow remains canonical.

An expensive run requires the explicit `--execute` confirmation flag. Stage
settings expose the SVD driver/thread policy, rank and tolerance, SINDy storage
policy, and method-specific options without adding another execution engine.
Run output is versioned JSON containing compact archive paths; numerical arrays
remain in derived archives.

Native solution archives store the coordinates actually advanced by the
solver. For dimensionless RPE records, H-DMD and SINDy therefore consume `tau`,
`r`, and `dr/dtau`. Derived records are immutable by default and SINDy archives
carry the exact upstream
H-DMD archive lineage needed by Viewer. The standard SINDy archive stores its
polynomial coefficients, dictionary, state transform, and training states; its
large eigenfunction sample matrix is reconstructed exactly and column-wise on
demand. This coefficient-first representation is storage policy, not modal
truncation: all H-DMD/SINDy eigenvalues remain present. Standard execution also
defers construction of that matrix in memory, avoiding roughly 763 MiB for the
25,001-by-2,000 result while retaining callable polynomial eigenfunctions.

## Desktop GUI and application API

The versioned compute-job contract and compiled workflow service are the stable
front-end API. The native queue schedules durable job manifests and exposes
lifecycle state without importing Python package internals. The primary human front end is the native Solver and
Post-processor described below. The older Tk interface remains a compact
compatibility surface; inspect it with `python -m gui --help` and launch it
explicitly with `python -m gui launch` or `dynsyslab-gui launch`.

Sampled solves should set `SolverOptions(control_evaluation="sampled",
control_sample_time=...)`. This creates an integer-tick zero-order-hold clock
instead of changing a PID/observer state at internal RHS evaluations. Stateless
LQR/SDRE/KLQR laws may instead request `control_evaluation="continuous_rhs"`
to recompute their state-dependent gain at every internal solver evaluation.
Archive compression defaults to lossless LZMA; set
`DYNSYSLAB_COMPRESSION=deflate` when write speed matters more than disk space.

### Native Solver and Post-processor

The C++17 workbench is the primary interactive application. Solver and
Post-processor are separate panels over a shared, versioned native boundary;
post-processing an archived solve does not rerun the dynamics.

Solver job names apply to the entire proposed tree. Each System node preserves
its own model, parameters, ICs, controller, target, campaign, and numerical
settings while the user navigates between nodes. The same `solve_set` contract
is accepted by the compiled planner, runner, and durable queue.

```powershell
dynsyslab-compute launch
# Source checkout:
python -m compute_workbench.cli launch
```

Use `template`, `validate`, `plan`, and explicitly gated `run --execute`
subcommands for automation. See the
[user guide](docs/SOLVER_POSTPROCESSOR.md) and
[architecture](docs/architecture/SOLVER_POSTPROCESSOR_CPP_WORKBENCH.md).

### Native Viewer workbench

The comprehensive Viewer workbench is a modular C++17 desktop application for
traversing linked trajectories, training/control datasets, KROMs,
eigenfunctions, control-map regressions and errors, control attempts,
observability analyses, KLQE results, campaigns, and saved reports. It uses the
compiled archive/payload services as its scientific boundary and renders bounded
plot payloads natively, so archive meaning is not duplicated in the UI.

```powershell
dynsyslab-viewer build
dynsyslab-viewer launch
dynsyslab-viewer launch --source <a run .npz from your archive folder>
```

Use `python -m viewer.ui_cli` for the source-tree equivalent. `catalog` and
`payload` subcommands expose the same versioned boundary for tooling and tests.
See [the native Viewer architecture](docs/architecture/VIEWER_CPP_WORKBENCH.md)
for the information model, UI modules, supported views, performance policy, and
first-pass limitations.

## Control and estimation

Canonical dimensionless RPE state targets are exposed through
`control.RPE_CONTROL_TARGETS` and `control.rpe_control_target(name)`. Each
harmonic target returns the analytic pair `(r, dr/dtau)` from angular-frequency
parameters, avoiding numerical differentiation.

Controlled solves archive three distinct records: the solution-aligned actual
control `TimeSeries`, the solution-aligned target `TimeSeries`, and the exact
controller-evaluation log. Continuous-RHS logs include repeated and rejected
adaptive-solver evaluations; sampled logs contain only actual clock updates,
not reuse of held values. Viewer keeps
these roles separate: the target is a dashed red overlay on state and phase
plots, while actual controller output is rendered in its own control
time-series panel(s).

## KROM extraction

`KROM.extract(selection)` and `rom.extract_krom(source, selection)` create a
new KROM from selected nonredundant KROM indices or exact source eigenvalue /
eigenfunction objects. Existing conjugate partners are included by default;
pass `include_conjugate_duals=False` only when a deliberately unpaired complex
coordinate is required. Eigenvalues, eigenfunctions, sampled values,
coefficients, reconstruction modes, and compatible linear/bilinear control
matrices are copied into an independent model. Extraction provenance records
the source spectrum and all source-spectrum indices; the source KROM is not
mutated.

`rom.reserve_krom_archive(...)` publishes that subset as an immutable derived
archive, including coefficient vectors, sampled eigenfunction values when
available, dictionary/state preprocessing, coordinate-system metadata, and
source lineage. `rom.load_krom_archive(...)` reconstructs an executable
polynomial-eigenfunction KROM from either a full or reserved archive.

`rom.select_observable_modes(...)` applies grouped greedy state-observability
selection to polynomial KROMs. Complex conjugate coordinates remain one
indivisible real/imaginary group, control authority is an optional tie-breaker,
and the report records every selected full-spectrum index and cumulative state
`R^2`. The implementation works in the small polynomial-feature space rather
than loading the full eigenfunction-value matrix. Use
`python -m rom select-observable <archive>` to publish an immutable selected
KROM with the complete policy report embedded in its provenance.

KROM archives expose an axis-aligned training domain when training states are
available. Koopman LQR/MPC defaults to `target_domain_policy="error"`; an
out-of-domain state or target must be explicitly changed to `"warn"` or
`"ignore"` for a labeled extrapolation study.

State reconstruction is a separate validated capability. Use
`python -m rom fit-observation <archive>` to fit an affine real/imaginary-split
decoder with held-out quality gates; failed maps are not published unless the
diagnostic override is explicit. `dynsyslab-krom-training` plans a
system-agnostic state-domain expansion for any compatible ODE KROM and writes
no campaign until `--publish` is supplied. Publication requires a matching
`--screen`; the screen uses the non-persistent canonical ODE kernel and checks
candidate ICs plus any supplied state probes. RPE's canonical target envelopes
are now a preset adapter over this same generic workflow. See
[`docs/KROM_OBSERVATION_AND_TRAINING_DOMAIN.md`](docs/KROM_OBSERVATION_AND_TRAINING_DOMAIN.md).

`control/` is an optional solver plug-in package with multi-tone waveform,
continuous LQR/SDRE-style LQR, finite-horizon linear MPC, Koopman-ROM LQR/MPC,
and PID laws. Targets may be constants, functions, or `TimeSeries` objects.
`estimation/` supplies identity, static-map, and linear Luenberger estimators.
Pass either or both to `SolveRequest(controller=..., estimator=...)`; the
solver gives the controller the estimator's current state at each RHS
evaluation. The default estimator is `IdentityEstimator`, so unspecified
estimation means exact full-state feedback. Every controller accepts
`active_time_window=(start, end)` on its own time coordinate; outside that
interval it has no effect. No controller or estimator changes an ordinary solve. A controller
with `actuation="system_input"` calls a model's explicit input channel (RPE
supports this); `additive_state` controllers return a state-derivative vector.
ODE plug-ins expose their input channels as `rhs_with_input(t, x, u)` and an
`input_matrix(x)` accessor. For most models this returns a constant matrix
`B`; RPE returns its state-dependent `B(x)`. RPE input is external acoustic
pressure: positive pressure compresses the bubble. New systems record
`pressure_input_convention="external-acoustic-pressure/v2"`. Archives that
predate this field reconstruct with `legacy-positive-acceleration/v1` so old
solutions and analytic KROM control experiments remain reproducible.
When a controller is used, its actual solve-time output is stored as an aligned
`control_values` time series (and `raw_control_values` when interpolation is
requested); load it with `archives.load_control_time_series`.
The native workflow binds that history to the model's declared coordinates.
For physically nondimensionalized bubble systems, prescribed pressure and
clock settings are accepted in Pa, Hz, and seconds and converted to
`p_control/p0` and model time. Archives retain both the model units and the
scale back to physical units; explicit `model_*` controller settings are
available when an experiment is already specified in nondimensional
coordinates.

`rom.KROM` is an explicit reduced-order model built from a `KoopmanSpectrum`
with callable eigenfunctions. A spectrum always identifies its source dynamical
system; the KROM control mode determines how its lifted input matrix is formed.
KROM capability metadata distinguishes lifted evolution, state reconstruction,
and analytic control. Delay-coordinate H-DMD modes cannot be used as physical
state modes, while SINDy KROMs expose eigenfunction lifting/control without
claiming reconstruction modes. Dimensionless controller adapters accept physical
solver state/time and preserve the system input's physical units.
`KoopmanLQRController(..., target=<state target>)` lifts a canonical model-state
target for feedback while exposing the corresponding physical-state target to
the archive and Viewer. Pass `target_space="lifted"` only for a target already
defined in the KROM's real split-coordinate representation.

SINDy eigenfunctions may instead be represented in arbitrary ordered model
observables `y=h(x)`. Analytic observable expressions retain exact derivative
observables, SINDy uses `ydot=J_h(x)f(x)`, and KROM/KLQR control uses
`B_psi=nabla_y(Psi)J_hB`. Reddy--Szeri exposes `radial_state`, `c_k2`, and
`c_k5` coordinate presets. See
[`docs/architecture/OBSERVABLE_COORDINATE_SINDY_AND_KLQR.md`](docs/architecture/OBSERVABLE_COORDINATE_SINDY_AND_KLQR.md).

## Reliability defaults

Archive publication and manifest writes are atomic. Archive read-modify-write
operations use per-record locks, while readers remain lock-free. Failed jobs
are retried only when their `execution.retry_limit` permits it (the default is
zero extra retries). Low-level dense DMD, Hankel-DMD, and Koopman-with-SINDy
calls retain a conservative `max_memory_mb=1024` preflight. The explicitly
large `trajectory-hdmd-sindy` preset declares an 8192 MiB workstation budget;
its CLI planning command reports the estimate before execution. Pass `None`
only when an external resource manager has reserved sufficient memory.

## RPE presets

`dynsyslib/lib/rpe.py` contains ten named liquid/gas starting presets. Their source
notes and limitations are in `dynsyslib/lib/PRESET_SOURCES.md`; several deliberately carry
cautions because classical RPE does not validate a non-Newtonian or strongly
viscoelastic liquid.

The default `air_water_20c` RPE uses `R0=1e-5 m`, `kappa=1.3333`,
`p0=1e5 Pa`, `pV=2330 Pa`, `rhoL=998 kg/m^3`, `sigma=0.0725 N/m`, and
`muL=0.001 Pa s`. Its natural-frequency nondimensionalization is

`f0 = sqrt(3*kappa*pG0 - 2*sigma/R0) / (2*pi*R0*sqrt(rhoL))`,

where `pG0 = p0 - pV + 2*sigma/R0`. The characteristic time is `T0=1/f0`,
dimensionless time is `tau=t/T0=f0*t`, radius is `r=R/R0`, and dimensionless
wall velocity is `u=(T0/R0)*dR/dt`. Dimensionless archives record the complete
scaling under `metadata["nondimensionalization"]`; the Viewer reads that
provenance instead of reconstructing a scale from mutable defaults.
