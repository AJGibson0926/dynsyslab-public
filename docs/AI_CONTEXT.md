# dynsyslab AI context and maintainer guide

## Purpose

DynSysLab is a local, modular native C++ laboratory for defining dynamical
systems, solving them numerically, preserving results, processing trajectory
data, and displaying results. It installs as compiled desktop software and does
not require Python. The distributed Python wheel contains only the optional
dependency-free `dynsyslab_client`. The former Python application packages have
been retired; every capability is in the native core (see the
[native port ledgers](architecture/native-port/README.md)). This document is a
compact orientation packet for a new developer or AI assistant.

## Spectral Koopman analysis

The [Lie algebra service](architecture/LIE_ALGEBRA.md) owns exact polynomial geometric analysis in the independent `dynsyslab_lie` library. Preserve finite/infinite proofs versus search exhaustion, and physical-image versus ambient stability. The private study record `KOOPMAN_LIE_REPRESENTATION_RESEARCH.md` formalizes the central research question and Worthmann/Strasser connection.

The [Spectral KROM contract](architecture/SPECTRAL_KROM.md) adds native ResDMD,
mpEDMD, Hermitian/TLS/MultDMD, stochastic residuals, Rigged measures and RKHS
adjoint diagnostics. Scientific algorithms belong to dynsyslab_spectral;
Koopman owns frozen dictionaries/data/archive/reference binding. Foundation
observables carry passive references only. Do not interpret continuous spectral
measures or adjoint RKHS sections as a finite forward eigenmode model.

## Non-negotiable architecture

The user-designated **critically important Koopman reference** is the
exact linear Duffing benchmark (private development record `KOOPMAN_LINEAR_DUFFING_EXACT_REFERENCE.md`)
(reference ID: `koopman-linear-duffing-exact`). Preserve this baseline and its
local asset pointers for future analysis. Its eigenfunctions and Bpsi are
analytic ground truth, not a learned-model success. Keep new variants separate.

Also preserve the user-designated **critically important gradient-closure
benchmark**, `koopman-triangular-gradient-closed` (private record
`KOOPMAN_GRADIENT_CLOSED_TRIANGULAR_REFERENCE.md`), and the **important arctangent
counterexample**, `koopman-arctangent-infinite-gradient` (private record
`KOOPMAN_ARCTANGENT_GRADIENT_REFERENCE.md`). These have exact proofs distinguishing
autonomous eigenfunction closure, closure under original-state partial derivatives,
and finite bilinear input models. Do not treat these properties as interchangeable.

1. **Systems** (`native/src/Systems*.cpp`) define equations, presets,
   Jacobians, observables, state dimensions, and drivers. They do not solve,
   schedule, archive, or visualize.
2. **Solvers** execute one solve request at a time. They support ODE, PDE, DDE,
   and lattice pathways and write completed/checkpoint data only through the
   archive contract.
3. **Archives** own storage, audit, deduplication, and the canonical archive
   index. Applications reach them through the native workflow and
   `dynsyslab-cli archive|storage`, never by editing archive files directly.
   Viewer opens archives read-only.
4. **Campaign design** (in the workflow service) produces only job manifests.
   It never integrates or analyzes.
5. The **queue** executes lightweight job manifests and never owns numerical
   solution arrays. Its coordination state is SQLite; jobs remain JSON files
   and numerical arrays remain archive files.
6. **Viewer** (`viewer/native/`) is display-only relative to the primary
   archive. It may create disposable display artifacts only in
   `viewer_compute_archive`, never in `archives`.
7. **Post-processing** analyzes archived solution data. Its Koopman and
   POD/PCA methods are independent stages; LCS (`ftle`) and SSM (`ssm`) are native
   post-processing methods ([LCS and SSM](architecture/LCS_AND_SSM.md)).
8. The native Solver/Post-processor (`compute_workbench/native/`) and Viewer are
   clients of versioned public C++/C/JSON boundaries. They do not own
   scientific algorithms or data.

## Package map

| Module | Responsibility | Must not do |
|---|---|---|
| systems (`native/`) | System definitions, presets, drivers, Jacobians | I/O, scheduling |
| solvers (`native/`) | One numerical solve plus optional controller/estimator | Campaign design |
| archives (`native/`, `workflow/native/`) | Store arrays, index, audit, storage locations | Integrate equations |
| workflow (`workflow/native/`) | Plan and run jobs, campaigns, experiments, post-processing | Mutate source trajectories |
| queue (`dynsyslab-cli queue`) | Claims, scheduling, retries, pause/recovery | Store arrays |
| Koopman, control, estimation | KROMs, feedback laws, estimators | Alter an uncontrolled solve |
| Viewer (`viewer/native/`) | Read-only rendering | Solve, schedule, change primary archive |
| Solver/Post-processor (`compute_workbench/native/`) | Versioned job editing | Implement scientific algorithms |
| `dynsyslab_client` (Python) | Optional thin client over the C ABI and CLI | Own algorithms or data |
| `operations/` (Python) | `python -m dynsyslab` front door | Own algorithms or data |
| `quality` (Python) | Static architecture/documentation audit | Import or execute numerical engines |
| `tester` (Python) | Test selection, execution, reports | Production archive mutation |

## Data model

The hierarchy is `Batch -> System -> Cluster -> Trajectory -> TimeSeries`.

- A `System` contains one dynamical-system object, integration mesh/time-span,
  and one or more clusters.
- A `Cluster` groups trajectories, usually distinct initial conditions.
- A `Trajectory` has an initial condition and, after solving, a `TimeSeries`.
- Native `TimeSeries` is the authoritative typed view of time values and
  component-major numerical samples. It may own buffers loaded from an archive
  or expose a zero-copy view over solver buffers. It records uniform-grid state
  and `delta_t` without changing the samples.
- Fourier analysis is a lazy, in-memory sidecar on `TimeSeries`. The first
  explicit request performs an O(N log N) transform; the latest requested
  option set is cached and can be cleared. Frequencies, coefficients,
  amplitudes, powers, and phases are derived data and are never written into a
  solution archive.

Archives retain raw solutions only once. Tiered hierarchy archives point to
child solution archives rather than duplicating their numerical arrays.
Viewer's catalog therefore retains archive pointers and metadata rather than
eagerly constructing every `TimeSeries`; it materializes only the selected
record and invokes its Fourier sidecar only for a frequency-domain view.

## Archive format and safety

- Completed solutions are compressed `.npz` records, primarily under
  `archives/solutions/`.
- Hierarchy manifests are compact `.npz` metadata records under
  `archives/hierarchies/`.
- Derived postprocessing results are under `archives/derived/` and retain
  source pointers/provenance rather than copying source solution arrays.
- Queue coordination state is a transactional SQLite store under
  `archives/job_manager/` ([job queue](architecture/JOB_QUEUE.md)).
- Archive publication and manifest writes are atomic. Do not overwrite a
  completed solution to extend it; use the archive extension pathway.
- Never modify the original, protected RPE chaos-campaign archive without the
  user explicitly providing the agreed password in that conversation.

## System-library contract

Built-in ODE systems live in `native/src/Systems*.cpp`; third-party systems are
compiled extensions over `native/include/dynsyslab/extension_api.h`. A typical
system exposes:

- its right-hand side and, where available, its Jacobian;
- initial state, state dimension, and period metadata;
- an invalid-state check for model-domain/collapse guards;
- an input-aware right-hand side and `input_matrix(x)` for
  affine/state-dependent control when meaningful;
- named presets and observables, listed by `dynsyslab-cli systems`.

RPE has both dimensional and nondimensional forms. RPE defaults to its
nondimensional form; archive metadata records which form was used.

The native catalog also carries an explicit system-family field. Mackey--Glass
is the first native `dde` model and lives in `dynsyslab_delay.dll`; it uses a
constant-history method-of-steps workflow and the same immutable hierarchy as
ODE solves. Do not pass it through the ODE factory. Parameter/time/IC campaigns
are native; DDE control remains an explicit error.

Native 2-D/3-D Ising and 2-D Potts systems live in
`dynsyslab_lattice.dll`. They archive scalar thermodynamic histories plus exact
integer initial/final site arrays. Their native random streams are seed-stable
but intentionally not bitwise NumPy-compatible. Lattice controllers and
remain explicit errors; free lattice campaigns are native.

The general structured-field Euler/RK4 API and the built-in one-dimensional
heat equation live in `dynsyslab_pde.dll`. The workflow archives space and full
field histories, enforces a conservative explicit stability limit, and exposes
PDE parameter/time/IC campaigns. Viewer owns the space-time heat-map display;
the PDE module has no UI dependency. An optional System-owned
`system.pde_partitioning` object selects deterministic Cartesian partitions,
boundary-aware ghost regions, and
a bounded local worker pool; MPI and adaptive refinement are not implied. PDE
control is not defined. Legacy `solver.parallel` jobs remain readable.

## Solve lifecycle

A solve is a versioned `dynsyslab.compute-job/v1` document with
`"kind": "solve"` (`dynsyslab-cli template solve` writes one). It chooses a
system, time span, initial conditions, solver options, history policy,
interpolation policy, optional controller/estimator, and an archive
destination. `dynsyslab-cli validate|plan|run` and the queue accept the same
document; a run returns versioned JSON with status and archive paths. See
[SOLVER.md](../SOLVER.md).

History defaults to lightweight completed-solution storage. Checkpoints and
high-frequency history are opt-in through Job/Campaign metadata. A pause is
cooperative and occurs at a safe solve boundary or configured checkpoint.

Controller families are deliberately separate. `lqr`, `sdre`, `mpc`, and
`pid` operate on the dynamical state. `klqr` and `kmpc` require a KROM and
operate on split-real Koopman coordinates. KLQR is state-dependent: at each
continuous adaptive RHS evaluation it recomputes `B_Psi(x)`, solves the CARE,
and updates `K(x)`. PID, local state-dependent LQR, state-space MPC, and KMPC
use exact sampled ticks with zero-order hold. Do not replace KLQR with the
fixed state-space LQR or silently switch its default to sampled evaluation.

## Job and campaign workflow

A job document is a lightweight plan: systems, solver settings,
post-processing stages, dependencies, resources, and tags. It contains no
arrays. Experiments chain stages over one pipeline
([experiments](architecture/EXPERIMENTS.md)).

The native queue (`dynsyslab-cli queue ...`) supports dependency conditions,
duplicate-claim protection, resumable pauses, recovery, retries, resource
admission, and archive reuse/extension lookup. It checks compatibility before
reuse, including time span, dimensional/nondimensional form, solver settings,
and forcing metadata for forced trajectories
([job queue](architecture/JOB_QUEUE.md)).

Campaign axes in a solve job generate coherent collections of trajectories and
campaign metadata for bifurcation sweeps, Lyapunov sweeps, trajectory
ensembles, and parameter sweeps; `dynsyslab-cli campaign ...` inspects and
repairs them. A campaign is not analysis; it is a reproducible data-generation
plan.

## Postprocessing and reduced-order models

- Koopman methods: DMD, EDMD, Hankel DMD/HDMD, HAVOK, Koopman-with-SINDy,
  validation, dictionaries, spectra, and resampling, as post-processing
  stages of a job ([KROM tools](architecture/KROM_TOOLS.md)).
- A Koopman spectrum stores eigenvalues/eigenfunctions and source-system
  metadata. Eigenvalues maintain discrete and continuous representations with
  the `delta_t` conversion provenance.
- A KROM is built from a spectrum. Complex output is the mathematical default;
  controllers and estimators request split real/imaginary state coordinates.
- POD/PCA stages produce POD modes and POD-Galerkin models.
- Uniform-time requirements must be checked before methods that assume them.
  Offer resampling for postprocessing only; never silently alter archive data.

## Viewer rules

The native Viewer browses the archive catalog and renders bounded plot
payloads from the compiled archive/payload services
([Viewer architecture](architecture/VIEWER_CPP_WORKBENCH.md),
[Viewer CLI](architecture/VIEWER_CLI.md)).

Viewer may display trajectories, spectra, Koopman/POD modes, FFTs, phase
portraits, 3-D data, and already-computed campaign data. It must never request
a solve, schedule a job, mutate the primary archive, or infer missing data by
running a new calculation. It may create the minimum required temporary
display artifact in `viewer_compute_archive`.

## Development rules

- Prefer public package interfaces over reaching into internals.
- Preserve source pointers and metadata; do not duplicate large arrays merely
  for convenience.
- Use structured warnings/events with preserved error context instead of broad
  `except: pass` blocks.
- New archive I/O goes through the native archive and workflow APIs.
- Native tests go in the module's `tests/` folder (CTest); Python tests go in
  `tester/tests/`, not a top-level testing package.
- Generated `build/`, `dist/`, egg-info, caches, viewer artifacts, and test
  reports are not source. `project_paths.source_files()` excludes them.
- Keep every source file below 1,000 lines and treat 500 as the normal ceiling.
  Add API/CLI gates before UI bindings and preserve versioned process contracts.
- Run `python -m quality --format markdown` to measure current size,
  documentation, link, and package-boundary debt.

## Verification checklist

Run the installed native gates first:

```powershell
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" doctor
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" capabilities
```

When changing the optional Python client, the front door, or the audit, also
run from the project root:

```powershell
$env:PYTHONPATH = (Resolve-Path '.').Path
python -m dynsyslab_client doctor
python -m quality --format markdown
python -m tester --all
python -m build
```

The audit reports the currently discoverable test-method count so this guide
does not preserve a stale number. If architecture or archive behavior changes,
add a regression test and run the full suite before handoff. The distribution
command requires the standard `build` frontend; CI additionally installs the
wheel outside the checkout and exercises its public commands.
