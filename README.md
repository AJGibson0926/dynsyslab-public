# dynsyslab

`dynsyslab` provides a native C++ dynamical-system library, single and composite
numerical solves, durable archives, post-processing, read-only viewing, and a
durable job queue with campaign planning. Everything runs in the compiled core
(`dynsyslab-cli`, `DynSysLab.exe`, `DynSysLabViewer.exe`); Python is optional.

**New here? Read the [DynSysLab overview](docs/OVERVIEW.md)** — what the suite does, its modules and experiments, and how to analyse and control a dynamical system, with architecture diagrams.

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
`hindmarsh_rose`, `hodgkin_huxley`, `fput`, and `kuramoto` are native ODE models selected by `system.name` in a solve job; `dynsyslab-cli systems`
lists them with their family, state dimension, and presets. Every model supplies
its right-hand side, Jacobian, initial state, period, and invalid-state guard; the
campaign machinery never branches on equations.

Mackey--Glass is the native `dde` model because a delay equation requires a history
buffer and a different integrator. The 2-D/3-D Ising and Potts models are native
`lattice` models with Metropolis dynamics; they produce temperature/field phase maps
rather than a deterministic bifurcation diagram or Lyapunov exponent.

Campaign axes are encoded in the common compute-job contract and executed by
the compiled workflow, either directly or through its durable SQLite queue.
`dynsyslab-cli queue submit|list|status|cancel|run-next|serve|recover` operates
that queue (see [the job queue](docs/architecture/JOB_QUEUE.md)). Numerical data remains exclusively in the archive. A
positive Lyapunov screen remains a candidate, not by itself proof of a chaotic
attractor.

## Adding a system

Built-in systems are registered in the native core (`native/src/Systems*.cpp`).
A third-party system is a compiled extension implementing
`native/include/dynsyslab/extension_api.h`: an ODE extension supplies `rhs`,
`jacobian`, `input_matrix`, `invalid`, and observables, with optional delay, field
(PDE), and lattice tables. Forcing belongs to the system's parameters or to a
prescribed controller, not to the integrator. See
[the native solver](docs/architecture/NATIVE_SOLVER.md).

## Standalone solver

The system-agnostic solver is the native `solve` job; see `SOLVER.md`. It dispatches ODE,
PDE, delay-equation, and lattice requests and saves every completed result to a
per-system archive. History checkpoints are opt-in.

Numerical execution is owned by the native C++ engine. Euler, RK4, adaptive
Dormand-Prince RK45, and SUNDIALS/CVODE BDF are available through the same
workflow API. The optional Python client calls this engine; it does not own a
second solver implementation. See [native build and installation](docs/BUILD_NATIVE.md)
and `ADAPTIVE_BACKENDS.md`.

## Independent archives

The archive root (by default `archives/` in a source build, or the installed
archive folder) is the durable data store. Completed solutions go to
`archives/solutions/`; optional in-progress histories go to `archives/history/`.
Archives are read and maintained only through the native workflow and
`dynsyslab-cli archive|storage`; see [archive management](docs/ARCHIVE_MANAGEMENT.md).

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
dependency. The former Python research packages have been retired: every
capability is in the native core (see the
[native port ledgers](docs/architecture/native-port/README.md)).

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
H-DMD-to-SINDy job.

H-DMD's standard default remains the first scalar state component (radius for
RPE and Marmottant, position for Duffing), but the same workflow accepts any
real scalar entry in the archived system's observable catalog:

Set the H-DMD stage's `settings.observable` to a state-component index or a
named observable such as `far_field_acoustic_echo` or
`inertial_kinetic_energy`.

Use separate immutable jobs with identical matrix, window, rank, sampling,
sorting, and linear-algebra settings to compare spectra from several scalar
observables in Viewer.

A history-dependent or externally sampled scalar signal can be joined to the
source as an additional time series; see
[H-DMD scalar observables](docs/architecture/HDMD_SCALAR_OBSERVABLES.md).

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
lifecycle state. The primary human front end is the native Solver and
Post-processor described below.

Sampled solves should set `solver.control_evaluation: "sampled"` with a
`control_sample_time`. This creates an integer-tick zero-order-hold clock
instead of changing a PID/observer state at internal RHS evaluations. Stateless
LQR/SDRE/KLQR laws may instead request `control_evaluation: "continuous_rhs"`
to recompute their state-dependent gain at every internal solver evaluation.
Native archives are written with lossless deflate; choose LZMA (denser) or stored (fastest)
with `--compression` or `DYNSYSLAB_COMPRESSION`, and re-encode existing archives with
`dynsyslab-cli archive recompress`. See [archive management](docs/ARCHIVE_MANAGEMENT.md).

### Native Solver and Post-processor

The C++17 workbench is the primary interactive application. Solver and
Post-processor are separate panels over a shared, versioned native boundary;
post-processing an archived solve does not rerun the dynamics.

Solver job names apply to the entire proposed tree. Each System node preserves
its own model, parameters, ICs, controller, target, campaign, and numerical
settings while the user navigates between nodes. The same `solve_set` contract
is accepted by the compiled planner, runner, and durable queue.

```powershell
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\DynSysLab.exe"
# or: python -m dynsyslab ui   (source checkout)   /   dynsyslab-compute (Python client)
```

Use `dynsyslab-cli template`, `validate`, `plan`, and `run` for automation. See the
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
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\DynSysLabViewer.exe"
# or: python -m dynsyslab viewer   (source checkout)   /   dynsyslab-viewer (Python client)
```

Source builds produce the same executable through CMake (see
[native build and installation](docs/BUILD_NATIVE.md)); the
[Viewer CLI](docs/architecture/VIEWER_CLI.md) exposes its versioned boundary for
tooling and tests. See [the native Viewer architecture](docs/architecture/VIEWER_CPP_WORKBENCH.md)
for the information model, UI modules, supported views, performance policy, and
first-pass limitations.

## Control and estimation

Canonical dimensionless RPE state targets are native expression targets (see
[native observable expressions and
targets](docs/architecture/NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md)). Each
harmonic target supplies the analytic pair `(r, dr/dtau)` from angular-frequency
parameters, avoiding numerical differentiation.

Controlled solves archive three distinct records: the solution-aligned actual
control (`control_values`), the solution-aligned target (`control_target_values`), and the exact
controller-evaluation log. Continuous-RHS logs include repeated and rejected
adaptive-solver evaluations; sampled logs contain only actual clock updates,
not reuse of held values. Viewer keeps
these roles separate: the target is a dashed red overlay on state and phase
plots, while actual controller output is rendered in its own control
time-series panel(s).

## KROM extraction

The post-processing method `krom_extract` (with `selected_modes` or a
`mode_selection` rule) creates a new KROM from selected nonredundant KROM indices.
Existing conjugate partners are included, so each complex pair stays one
real/imaginary group. Eigenvalues, eigenfunctions, sampled values,
coefficients, reconstruction modes, and compatible linear/bilinear control
matrices are copied into an independent model. Extraction provenance records
the source spectrum and all source-spectrum indices; the source KROM is not
mutated.

The subset is published as an immutable derived archive, including coefficient
vectors, sampled eigenfunction values when available, dictionary/state
preprocessing, coordinate-system metadata, and source lineage. Full and
extracted KROM archives are both executable by the native control and rollout
stages; `dynsyslab-cli krom-info <krom.npz>` summarizes one.

The `state_observability` mode-selection rule applies grouped greedy
state-observability selection to polynomial KROMs. Complex conjugate coordinates remain one
indivisible real/imaginary group, control authority is an optional tie-breaker,
and the report records every selected full-spectrum index and cumulative state
`R^2`. The implementation works in the small polynomial-feature space rather
than loading the full eigenfunction-value matrix. Used with `krom_extract`, it
publishes an immutable selected KROM with the complete policy report
(`mode_selection_report`) embedded in its provenance. See
[KROM tools](docs/architecture/KROM_TOOLS.md) and
[mode selection](docs/architecture/MODE_SELECTION.md).

KROM archives expose an axis-aligned training domain when training states are
available. Koopman LQR/MPC defaults to `target_domain_policy="error"`; an
out-of-domain state or target must be explicitly changed to `"warn"` or
`"ignore"` for a labeled extrapolation study.

State reconstruction is a separate validated capability. The post-processing
method `observable_reconstruction` fits an affine real/imaginary-split ridge
decoder with held-out quality gates; failed maps are not published unless the
diagnostic override is explicit. `dynsyslab-cli krom-training-plan <request.json>` plans a
system-agnostic state-domain expansion for any compatible ODE KROM and writes
no campaign: it returns the plan and its screened solve job, which runs only
when you submit it. The screen uses the non-persistent canonical ODE kernel and
checks candidate ICs plus any supplied state probes. RPE's canonical target envelopes
are now a preset adapter over this same generic workflow. See
[`docs/KROM_OBSERVATION_AND_TRAINING_DOMAIN.md`](docs/KROM_OBSERVATION_AND_TRAINING_DOMAIN.md).

A solve job's `controller` section selects a multi-tone `waveform`, `lqr`,
`sdre`, finite-horizon linear `mpc`, Koopman-ROM `klqr`/`kmpc`, or `pid` law.
Targets may be constants, expressions, or archived time series. The `estimator`
section supplies `static_map`, `luenberger`, `lqe`, or Koopman (KLQE)
estimators; the solver gives the controller the estimator's current state at
each evaluation. Without an estimator the controller receives the exact full
state. Every controller accepts `active_time_window: [start, end]` on its own
time coordinate; outside that
interval it has no effect. No controller or estimator changes an ordinary solve. A controller
with `actuation="system_input"` calls a model's explicit input channel (RPE
supports this); `additive_state` controllers return a state-derivative vector.
Native models expose their input channels through an input-aware right-hand
side and an `input_matrix(x)` accessor. For most models this returns a constant matrix
`B`; RPE returns its state-dependent `B(x)`. RPE input is external acoustic
pressure: positive pressure compresses the bubble. New systems record
`pressure_input_convention="external-acoustic-pressure/v2"`. Archives that
predate this field reconstruct with `legacy-positive-acceleration/v1` so old
solutions and analytic KROM control experiments remain reproducible.
When a controller is used, its actual solve-time output is stored as an aligned
`control_values` time series (and `raw_control_values` when interpolation is
requested); Viewer and post-processing read it from the archive.
The native workflow binds that history to the model's declared coordinates.
For physically nondimensionalized bubble systems, prescribed pressure and
clock settings are accepted in Pa, Hz, and seconds and converted to
`p_control/p0` and model time. Archives retain both the model units and the
scale back to physical units; explicit `model_*` controller settings are
available when an experiment is already specified in nondimensional
coordinates.

A KROM is an explicit reduced-order model built from a Koopman spectrum
with callable eigenfunctions. A spectrum always identifies its source dynamical
system; the KROM control mode determines how its lifted input matrix is formed.
KROM capability metadata distinguishes lifted evolution, state reconstruction,
and analytic control. Delay-coordinate H-DMD modes cannot be used as physical
state modes, while SINDy KROMs expose eigenfunction lifting/control without
claiming reconstruction modes. Dimensionless controller adapters accept physical
solver state/time and preserve the system input's physical units.
A `klqr` controller with a state target lifts the canonical model-state
target for feedback while exposing the corresponding physical-state target to
the archive and Viewer. Set `target.space: "lifted"` only for a target already
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
zero extra retries). Dense DMD, Hankel-DMD, and SINDy stages refuse a job
before allocation when their working-set estimate exceeds the execution memory
limit (`execution.max_memory_mib` or `--max-memory-mib`); `dynsyslab-cli plan`
reports the estimate before execution.

## RPE presets

The native bubble models carry named liquid/gas starting presets
(`air_water_0c` ... `co2_honey_25c`); `dynsyslab-cli systems` lists them, and each
system reports its source papers in its
[research metadata](docs/SYSTEM_RESEARCH_METADATA.md). Several deliberately carry
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
