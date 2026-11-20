# DynSysLab suite architecture remediation

**Date:** 2026-08-31

**Source baseline:** `main` at `3d2d5e8`

**Implementation branch:** `refactor/audit-roadmap`

**Scope:** all P1, P2, and P3 work prioritized by the
suite-wide architecture audit (private development record `DYNSYSLAB_SUITE_ARCHITECTURE_AUDIT_2026-08-31.md`),
plus the requested low-overhead Solver progress display and the later
system-agnostic KROM-training generalization.

## Outcome

The architecture audit now passes in absolute strict mode and against an empty
debt baseline. No maintained Python or C++/CMake source exceeds 500 lines, no
hard-size or package-boundary finding remains, the package graph has no cycle,
all narrative documentation is indexed with valid local links, and the public
API docstring floor is enforced at 70%.

This was a behavior-preserving architecture refactor. It did not change an
equation of motion, solver method, KROM algorithm, controller law, estimator,
accepted scientific result, or completed archive. The suite's versioned
archive and UI payload schemas remain compatible.

## Before and after

| Measure | Audit baseline | Remediated suite |
|---|---:|---:|
| Maintained Python | 295 files / 39,266 lines | 328 files / 36,727 lines |
| Maintained native/CMake | 36 files / 3,382 lines | 39 files / 3,558 lines |
| Python files over 1,000 lines | 3 | 0 |
| Python files over 500 lines | 10 | 0 |
| Native files over 500 lines | 0 | 0 |
| Static boundary findings | 4 | 0 |
| Package-cycle components | 1 eight-package component | 0 |
| Module docstring coverage | 83.1% | 84.8% |
| Public-definition documentation | 49.8% baseline measure | 71.0% outward-API measure |
| Discoverable test methods | 405 | 413 |
| Broken/unindexed narrative docs | 0 / 0 after audit indexing | 0 / 0 |

The current public-documentation measure deliberately excludes test classes;
it measures top-level outward definitions. Package-level counts are included in
the audit JSON so weak areas cannot hide behind a suite-wide percentage.

## P1: hard-policy refactors

### Viewer

`viewer/core.py` is now a small façade and coordinator. Source resolution,
record contracts, plot-family rendering, and serialization live in focused
modules under `viewer/` and `viewer/rendering/`. Existing `Viewer`,
`PlotRequest`, source-routing, figure, HTML, and SVG behavior remains behind
the same public entry points.

### JobManager

`job_manager/manager.py` is now a façade over state, planning, execution,
reuse, transition, archive-inspection, and runtime services. Durable state,
claim, pause, repair, scheduling, and campaign behavior retains the established
manager API; the decomposition removes implementation ownership from the
façade without creating a second scheduler.

### Duffing bilinear search

The reusable bounded-search engine moved out of the one-off-sized script and
into `operations.duffing_bilinear`. Candidate selection, regression, exact and
approximate controller evaluation, schedules, paths, and campaign records are
independent modules. The retained script is a thin CLI. Existing ROM imports
have a compatibility gateway rather than a duplicate engine.

### Remaining preferred-size modules

KROM capabilities were separated into coordinate, observation, control, and
operation services while `KROM` remains the discoverable façade. Large
post-processing, estimation, control-candidate, and spectral modules were
similarly divided along existing behavior boundaries. The largest maintained
Python module is now `viewer/ui/catalog.py` at 499 lines.

## P2: package seams

### Neutral scientific contracts

- `analysis_data` owns trajectory/derived-record persistence contracts that
  would otherwise force scientific packages through storage or UI layers.
- `spectral` owns eigenvalue, eigenfunction, value, dictionary, and error
  contracts. Post-processing produces them; ROM consumes them.
- `infrastructure.persistence` provides narrow capability protocols for
  publishing/loading derived data.
- `infrastructure.controlled_data` owns controlled-data interchange instead of
  making system definitions depend on ROM.
- `infrastructure.coordination` and `infrastructure.estimation` provide small
  neutral coordination/estimator protocols.

These seams remove the former eight-package strongly connected component. The
strict audit now rejects any reintroduction of the explicit Viewer-producer,
application-to-low-level-Archives, and DynSysLib-to-ROM boundary violations.

### Ownership corrections

- Archiver no longer depends on JobManager implementation details.
- Solver diagnostics and estimator fallback are neutral operations rather than
  reverse imports into analytical/application packages.
- Observable-spectrum rendering belongs to Viewer; analytical packages return
  typed data.
- Campaigner designs manifests only. Generic numerical screening lives under
  `operations` and calls system/solver contracts.
- Legacy controlled CSV parsing is independent from construction of a KROM
  regression dataset.

## System-agnostic KROM training workflow

The former RPE planner is now a generic campaign workflow:

- `KROMTrainingPlan` accepts arbitrary registered systems, dimensionality,
  presets, constructor parameters, target states, state probes, and component
  envelope requirements;
- numerical screening accepts componentwise lower, upper, and absolute limits;
- initial-condition and probe screens are reusable operations with explicit
  publication gates;
- RPE remains a convenience preset using its canonical targets, not a special
  implementation path;
- the CLI is `dynsyslab krom-training` or the expert alias
  `dynsyslab-krom-training`.

Characterization covers 2-D Duffing, 3-D Lorenz, arbitrary target conversion,
generic screening, and campaign publication.

## Script lifecycle cleanup

The suite now distinguishes three things explicitly:

1. reusable scientific packages and normal CLI/API workflows;
2. retained expert campaign tools that remain repeatable research operations;
3. one-off MATLAB-parity wrappers, case reruns, ad hoc diagnostics, and export
   scripts whose evidence already lives in reports or immutable artifacts.

The third group was removed rather than being presented as permanent product
capability. The retained policy is documented in
[Script lifecycle policy](SCRIPT_LIFECYCLE_POLICY.md). No scientific archives
or case-study reports were deleted.

## P3: outward-facing quality

### One front door, specialist owners

`dynsyslab` remains non-executing without an explicit operation and now routes:

```text
doctor, systems, ui, viewer, legacy-gui, docs, audit
solve, postprocess, archive, jobs, estimate, krom-training
```

The normal routes delegate to package-owned parsers, so their complete help and
safety gates remain authoritative. Existing specialist console scripts remain
compatible for automation and research reproduction.

### Native UI support and diagnostics

Viewer and Solver/Post-processor now link one `dynsyslab_ui_support` static
library instead of compiling shared implementation files independently. The
library owns the native theme, JSON/model parsing, process bridge, navigation,
inspector, and Help/About action. Each application retains its own workbench,
panels, and scientific role.

Help/About reports the UI support version, Viewer and compute schema versions,
and the exact project and canonical archive roots supplied by the Python
launcher. It can open `docs/START_HERE.md` in a source checkout.

### Low-overhead Solver progress

The Solver footer now contains a progress label and native progress bar. The
implementation reuses existing lifecycle/heartbeat signals:

- adaptive Python solves expose in-trajectory progress through the throttled
  heartbeat already emitted by the solver;
- hierarchy solves advance at completed-trajectory archive boundaries;
- bulk operations without a safe fraction remain explicitly indeterminate;
- the UI polls one small atomic JSON record four times per second;
- the Python reporter publishes UI state at most once per second;
- ordinary API/CLI solves keep their established 60-second telemetry heartbeat
  unless a progress output is requested.

No numerical state is copied and no extra solver callback is inserted merely
to animate the UI.

### Documentation and release gates

- `docs/START_HERE.md`, `README.md`, `ARCHITECTURE.md`, and `docs/AI_CONTEXT.md`
  describe the same front door and package boundaries.
- Public API documentation is 71.0%; compute workbench, control, performance,
  infrastructure, scripts, and quality public definitions are fully
  documented, while solver and JobManager are each at 72%.
- The audit enforces a 70% suite floor, reports package-level coverage, and
  excludes test classes from the outward-API denominator.
- CI runs the complete suite on Python 3.10 and 3.12, performs a strict audit,
  builds sdist and wheel artifacts, smoke-tests the installed wheel outside the
  checkout, and compiles both Windows native applications.
- The distribution manifest excludes runtime campaigns, solutions, derived
  data, histories, JobManager state/reports, management reports, bytecode, and
  build output. It includes the compute workbench required by the public front
  door.

## Verification evidence

Final verification on the implementation branch:

- complete Python suite: **413 passed**, 132.071 seconds;
- strict architecture audit: **PASS**;
- files over 500/1,000 lines: **0 / 0** for Python and native source;
- boundary violations: **0**;
- package cycles: **0**;
- documentation: **0 broken links, 0 unindexed narrative documents**;
- Viewer CMake/MSVC Release build: **PASS**;
- Solver/Post-processor CMake/MSVC Release build: **PASS**;
- sdist and wheel build: **PASS**;
- wheel: 343 entries, 685,202 bytes, compute workbench present, no runtime
  archive records or bytecode caches;
- isolated installed-wheel help/import smoke: **PASS** for the front door,
  system registry, solve, post-process, estimation, and core packages.

The native applications were compiled but not left running; this respects the
user's request to interrupt the app while preserving automated UI-contract and
build verification.

## Scientific and compatibility invariants

- Completed archives were not rewritten, moved, deduplicated, or reclassified.
- The protected RPE campaign was not modified or resumed.
- MATLAB remained a read-only evidence source; no MATLAB behavior was copied in
  this architecture refactor.
- Dimensionless dynamics/controller/control-history metadata rules are
  unchanged.
- The hierarchy remains `Batch -> System -> Cluster -> Trajectory ->
  TimeSeries`.
- Viewer remains read-only relative to the primary archive.
- Planning remains non-executing; numerical work still requires explicit run
  gates.

## Remaining non-blocking maintenance

There is no remaining strict architecture debt from the audit roadmap. The
following are sensible incremental work, not refactor blockers:

- raise package-local public documentation in post-processing, DynSysLib, and
  low-level Archives while keeping descriptions substantive;
- add manual accessibility/keyboard/high-DPI testing when the desktop apps are
  next run interactively;
- continue moving reusable expert-script behavior into packages when a second
  workflow needs it, and retire scripts that cease to be repeatable research
  operations;
- increase the supported-Python CI matrix if the project elects to promise
  additional versions beyond the current minimum and primary release runtime.

The suite should now evolve through these measured gates rather than another
large architecture campaign.
