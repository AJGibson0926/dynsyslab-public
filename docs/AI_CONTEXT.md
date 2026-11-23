# dynsyslab AI context and maintainer guide

## Purpose

DynSysLab is a local, modular native C++ laboratory for defining dynamical
systems, solving them numerically, preserving results, processing trajectory
data, and displaying results. It installs as compiled desktop software and does
not require Python. The distributed Python wheel contains only the optional
dependency-free `dynsyslab_client`; historical Python packages are source-only
reference implementations during migration, not an application runtime or an
alternate owner of canonical algorithms. This document is a
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

1. `dynsyslib` defines equations, presets, Jacobians, state dimensions, and
   drivers. It does not solve, schedule, archive, or visualize.
2. `solver` executes one solve request at a time. It supports ODE, PDE, DDE,
   and lattice pathways; it writes completed/checkpoint data through Archiver.
3. `archives` is the low-level storage implementation. Application packages
   must use `archiver.ArchiveGateway`, not direct `archives` I/O.
4. `archiver` is the archive policy boundary. `VIEWER_ARCHIVE` is read-only;
   `SOLVER_ARCHIVE` and `POSTPROCESSING_ARCHIVE` may write; `MANAGER_ARCHIVE`
   may perform management operations.
5. `campaigner` designs only manifests. It never integrates or analyzes.
6. `job_manager` executes lightweight Job manifests and never owns numerical
   solution arrays. Its coordination state is SQLite; jobs/manifests remain
   JSON files and numerical arrays remain archive files.
7. `viewer` is display-only relative to the primary archive. Its private
   `viewer/_compute` service may create disposable display artifacts only in
   `viewer_compute_archive`, never in `archives`.
8. `postprocessing` analyzes archived or supplied solution data. Its Koopman,
   POD/PCA, LCS, and SSM layers remain independent extensions.
9. `compute_workbench` and the native Viewer are clients of versioned public
   C++/C/JSON boundaries. They do not own scientific algorithms or data.
10. `infrastructure` holds only neutral, lightweight mechanisms needed across
    layers. Producers must not import Viewer to announce source changes.

## Package map

| Package | Responsibility | Must not do |
|---|---|---|
| `dynsyslib` | System definitions, presets, drivers, Jacobians | I/O, scheduling |
| `solver` | One numerical solve plus optional controller/estimator | Campaign design |
| `archives` | Store arrays and low-level indexes | Be imported by new application code |
| `archiver` | Archive capabilities, audit, search, storage locations | Integrate equations |
| `campaigner` | Build reproducible campaign/job manifests | Solve or plot |
| `job_manager` | Claims, scheduling, backends, retries, pause/recovery | Store arrays |
| `postprocessing` | DMD/Koopman/POD etc. | Modify source trajectories |
| `rom` | KROM and POD-Galerkin models | Own archive files |
| `control`, `estimation` | Optional live solver plug-ins | Alter an uncontrolled solve |
| `viewer` | Read-only rendering | Solve, schedule, change primary archive |
| `compute_workbench` | Versioned Solver/Post-processor UI jobs | Implement scientific algorithms |
| `infrastructure` | Neutral events and cross-cutting mechanics | Own scientific state |
| `quality` | Static architecture/documentation audit | Import or execute numerical engines |
| `tester` | Test selection, execution, reports | Production archive mutation |

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
- JobManager coordination state is `archives/job_manager/state.sqlite`, using
  a transactional SQLite store. Legacy `state.json` is only imported once.
- Archive publication and manifest writes are atomic. Do not overwrite a
  completed solution to extend it; use the archive extension pathway.
- Never modify the original, protected RPE chaos-campaign archive without the
  user explicitly providing the agreed password in that conversation.

## System-library contract

ODE system modules live in `dynsyslib/lib/`. A typical system exposes:

- `rhs(t, x)`;
- `jacobian(t, x)` where available;
- `initial_state`, `state_dimension`, and `period` metadata;
- `invalid(t, x)` for model-domain/collapse checks;
- `rhs_with_input(t, x, u)` and `input_matrix(x)` for affine/state-dependent
  control when meaningful;
- named presets registered in the common preset loader.

RPE has both dimensional and nondimensional forms. RPE defaults to its
nondimensional form; archive metadata records which form was used.

The native catalog also carries an explicit system-family field. Mackey--Glass
is the first native `dde` model and lives in `dynsyslab_delay.dll`; it uses a
constant-history method-of-steps workflow and the same immutable hierarchy as
ODE solves. Do not pass it through the ODE factory. Parameter/time/IC campaigns
are native; DDE control remains an explicit error rather than a Python fallback.

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

Use `app.RunConfiguration` for a stable frontend configuration or construct a
`solver.SolveRequest` directly. The request chooses a system, time span,
initial state, solver options, history policy, interpolation policy, optional
controller/estimator, and an archive destination. The solver returns a result
with status and archive path.

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

`JobManifest` is a lightweight plan: systems, solver metadata, postprocessors,
dependencies, schedule, resources, and tags. It contains no arrays.

`JobManager` supports series or dependency-safe parallel backends, duplicate
claims, resumable pauses, recovery manifests, schedule/deadline metadata, and
archive reuse/extension lookup. It checks compatibility before reuse, including
time span, dimensional/nondimensional form, solver settings, and forcing
metadata for forced trajectories.

`Campaigner` generates coherent collections of jobs and campaign metadata for
bifurcation sweeps, Lyapunov sweeps, trajectory ensembles, and parameter
sweeps. A campaign is not analysis; it is a reproducible data-generation plan.

## Postprocessing and reduced-order models

- Koopman layer: DMD, EDMD, Hankel DMD/HDMD, HAVOK, Koopman-with-SINDy,
  validation, dictionaries, spectra, and resampling.
- `KoopmanSpectrum` stores eigenvalues/eigenfunctions and source-system
  metadata. Eigenvalues maintain discrete and continuous representations with
  the `delta_t` conversion provenance.
- `KROM` is built from a spectrum. Complex output is the mathematical default;
  controllers and estimators request split real/imaginary state coordinates.
- POD/PCA layer produces `PODMode` objects and `PROM` collections.
- Uniform-time requirements must be checked before methods that assume them.
  Offer resampling for postprocessing only; never silently alter archive data.

## Viewer rules

Initialize `Viewer`, attach one or more sources, add `PlotRequest` objects,
then call `render()`. Source routing is automatic only when exactly one source
of the required type is attached; otherwise specify `source_name`.

Viewer may display trajectories, spectra, Koopman/POD modes, FFTs, phase
portraits, 3-D data, and already-computed campaign data. It must never request
a solve, schedule a job, mutate the primary archive, or infer missing data by
running a new calculation. It may ask its private compute service to create
the minimum required temporary display artifact.

## Development rules

- Prefer public package interfaces over reaching into internals.
- Preserve source pointers and metadata; do not duplicate large arrays merely
  for convenience.
- Use structured warnings/events with preserved error context instead of broad
  `except: pass` blocks.
- New archive I/O goes through `archiver`.
- New tests go in `tester/tests/`, not a top-level testing package.
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

When changing the optional Python adapter or retained parity references, also
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
