# DynSysLab suite-wide architecture and usability audit

> **Baseline record.** The prioritized remediation described here has now been
> completed. See the [suite architecture remediation report](DYNSYSLAB_SUITE_ARCHITECTURE_REMEDIATION_2026-08-31.md)
> for the implemented boundaries, before/after measurements, and final gates.

**Date:** 2026-08-31

**Scope:** Python suite, native numerical extension, native Viewer, native
Solver/Post-processor, CLIs, tests, documentation, and repository/data layout.

**Scientific-change policy:** No equations, solvers, KROM algorithms, controller
laws, estimators, or archived numerical records were changed during this audit.

## Executive conclusion

DynSysLab has a strong scientific-data foundation: systems are registered,
solves and derived products carry provenance, completed records are intended to
be immutable, post-processing is separable from integration, and both native
applications consume versioned Python/JSON contracts rather than reimplementing
the science. The C++ UI code is notably compact: all 36 maintained native files
are below the preferred 500-line ceiling.

The suite is not yet architecture-policy complete. Three Python files exceed
the hard 1,000-line limit, ten exceed the preferred 500-line limit, four static
package-boundary exceptions remain, and a large import-cycle component joins
eight principal packages. Public API documentation is uneven, general Python
CI is absent, and 23 console scripts make the expert surface more visible than
the normal workflow. These are maintainability and onboarding risks, not
evidence that the accepted scientific results are invalid.

This audit adds a safe `dynsyslab` front door, a human/agent start guide, a
repeatable architecture audit, a complete documentation map, and a neutral
source-event channel. It deliberately leaves the high-risk monolith splits as
explicit follow-up refactors with regression gates.

## Method and measured baseline

The repeatable command is:

```powershell
python -m dynsyslab audit --format markdown
```

It parses maintained Python sources, counts physical lines, measures top-level
public docstrings, extracts package-import edges, checks declared boundary
rules, checks local Markdown links and index coverage, inventories console
scripts, and separately measures maintained C/C++ files. Generated builds,
caches, test dependencies, temporary output, Viewer payloads, and runtime
archives are excluded from source metrics.

Measured after the audit infrastructure was added:

| Measure | Result |
|---|---:|
| Maintained Python | 295 files / 39,266 lines |
| Maintained native C/C++ and CMake | 36 files / 3,382 lines |
| Python files over 1,000 lines | 3 |
| Python files over 500 lines | 10 |
| Native files over 500 lines | 0 |
| Python module docstring coverage | 83.1% |
| Top-level public-definition docstring coverage | 49.8% |
| Statically discoverable test methods | 405 |
| Documentation | 43 Markdown documents including this report |
| Documentation-index omissions | 0 after reorganization |
| Console scripts | 23 |
| Remaining static boundary findings | 4 |

The test-method count is static and can differ slightly from executed unittest
cases produced through parameterization or dynamic loading. Full execution,
not this count, is the release gate.

## What is already working well

### Scientific and data boundaries

- `dynsyslib` provides registered ODE, delay, and statistical-system contracts
  with model-owned presets and observables.
- `SolveRequest` and solver options preserve a single-solve execution boundary;
  campaigns and jobs remain separate planning/scheduling concepts.
- The `Batch -> System -> Cluster -> Trajectory -> TimeSeries` hierarchy matches
  the scientific workflow and can represent free, forced, controlled,
  identification, and campaign data through roles and relationships.
- Archive and derived-result formats retain coordinate, source, parameter,
  method, and control provenance. Dimensionless controller metadata follows
  dimensionless dynamics.
- Post-processing can run on an existing solve. H-DMD observable selection,
  SINDy eigenfunction construction, KROM extraction, control-map regression,
  error maps, KLQE, and Viewer presentation are distinct callable operations.

### UI architecture

- The native Viewer and Solver/Post-processor use bounded, versioned JSON
  payloads and Python execution gates. Scientific interpretation remains in the
  Python services.
- The applications share archive/catalog concepts and native visual elements;
  they do not maintain separate scientific databases.
- Native sources are modular by current size policy. The largest are Viewer
  `Workbench.cpp` (374 lines), `PlotCanvas.cpp` (360), and native-kernel
  `bindings.cpp` (345).
- Catalog and payload generation are lazy enough to avoid eagerly materializing
  every heavy plot in the navigation tree.
- Execution remains explicit. Planning and validation do not silently run a
  solve or decomposition.

### Testing and reproducibility

- The test suite is broad and includes solver, archive, KROM, control,
  post-processing, UI-contract, CLI, and native-backend behavior.
- Case studies retain unusually detailed positive and negative evidence for
  RPE and Duffing KLQR behavior. This is valuable: failure cases remain
  reproducible rather than being overwritten by the latest successful run.
- The optional native numerical wheel has its own cross-platform build/test
  workflow and explicit capability admission.

## Critical size and responsibility findings

### Hard policy exceptions

| File | Lines | Responsibilities currently combined | Recommended split |
|---|---:|---|---|
| `viewer/core.py` | 1,549 | request/source contracts, routing, data shaping, plot-family dispatch, HTML/SVG generation | keep contracts in `viewer/core.py`; move source resolution to `viewer/sources.py`, plot families to `viewer/renderers/`, and serialization to `viewer/serialization.py` |
| `job_manager/manager.py` | 1,154 | claims/state, events, telemetry, transition/repair, scheduling, reuse, execution, interruption | a small `JobManager` façade over state/event, scheduler/reuse, execution, and recovery services |
| `scripts/search_duffing_bilinear_krom.py` | 1,061 | campaign model, candidate selection, regression, control evaluation, artifact/report orchestration, CLI | move reusable search/evaluation policy into `rom.candidates`; retain a thin experiment CLI |

These splits should be behavior-preserving first. Each phase should freeze the
public API and archive fixtures, add characterization tests, move one
responsibility, and compare artifacts before changing algorithms.

### Preferred-ceiling exceptions

The remaining files above 500 lines are `rom/koopman.py` (764),
`rom/candidates/control_identification.py` (596),
`postprocessing/koopman/koopman_sindy.py` (555), `estimation/koopman.py` (543),
`scripts/import_legacy_random_forced_rpe.py` (521),
`rom/control_campaign.py` (520), and
`postprocessing/koopman/workflows.py` (506).

The most important is `rom/koopman.py`: `KROM` presently combines immutable
model state, observation maps, exact and fitted control maps, error records,
eigenfunction evaluation, reconstruction, extraction, binding, and rollout.
Split capabilities into independent services or small owned value objects while
keeping `KROM` as the discoverable façade. Do not create mutually inconsistent
KROM subclasses for every capability combination.

## Dependency and boundary findings

Static imports contain one strongly connected component spanning `archiver`,
`dynsyslib`, `estimation`, `job_manager`, `postprocessing`, `rom`, `solver`, and
`viewer`. Some edges are type/conversion conveniences rather than runtime
coupling, but the component means a local change can propagate surprisingly
widely and makes independent package installation difficult.

High-value seams are:

1. **Shared schemas.** `campaigner/designer.py` and
   `job_manager/manifests.py` import `archives.schema` only for format/key
   validation. Move generic versioned-schema helpers into a neutral `contracts`
   package; storage can consume the same helpers.
2. **Post-processing source metadata.** `postprocessing/workflow.py` imports the
   private `archives.store.serializable_system`. Promote this conversion to an
   Archiver or model-provenance API.
3. **Legacy controlled-data import.** `dynsyslib/legacy_dynbatch_csv.py` imports
   `rom.ControlledDataset`. Parsing legacy system data is legitimate, but ROM
   conversion belongs in a data-adapter/ROM boundary so system definitions do
   not depend on reduced models.
4. **ROM/post-processing cycle.** Spectrum/eigenfunction value contracts should
   live in a small neutral spectral-contract module. Algorithms remain in
   post-processing; executable reduced models remain in ROM.
5. **Reporting adapters.** Estimation/post-processing should return typed plot
   data or reports; Viewer-specific rendering belongs in Viewer adapters.
6. **Campaign screening.** `campaigner/rpe_krom.py` both designs and numerically
   screens candidates. Move screening to a solver/operations diagnostic service
   while leaving manifest construction in Campaigner.
7. **Archiver/job cycle.** Archive inspection may describe job-manager state,
   but it should consume a coordination protocol or read-only DTO, not the
   scheduling implementation.

### Boundary correction completed during the audit

Solver, low-level Archives, Archiver, and JobManager previously imported the
private `viewer._sessions` package to publish refresh events. That reversed the
documented producer-to-Viewer dependency and made headless producers depend on
a UI package. The passive JSON event channel now lives in
`infrastructure.source_events`; Viewer retains a compatibility import. No
primary archive or event format changed.

## API and CLI findings

The reusable scientific workflows generally have Python gates and many have
CLI gates, but discoverability was inverted: 22 specialist commands existed
without a single `dynsyslab` command. The specialist commands remain useful for
automation and reproducible case studies, yet they should be the second layer.

The new front door provides:

```text
dynsyslab doctor       exact root, versions, archive and executable paths
dynsyslab systems      grouped system registry
dynsyslab ui           primary native Solver/Post-processor
dynsyslab viewer       native Viewer
dynsyslab docs         start-document path
dynsyslab audit        architecture-policy report
```

No-argument and `--help` paths do not launch an application. The compatibility
Tk GUI also now requires an explicit `launch` command and is lazily imported;
its former help path could import the desktop graph and block behind a window.

Long term, retain the existing specialist commands as compatible aliases but
organize normal operations beneath the front door (`solve`, `postprocess`,
`archive`, `jobs`, `estimate`). Experimental RPE/Duffing commands should remain
documented case-study tools rather than appear as equivalent top-level product
features.

## Documentation and onboarding findings

Before this audit, all local Markdown links were valid, but 17 of 38 narrative
documents under `docs/` were absent from the documentation index. The root
README had no immediate quick start, omitted Marmottant and Reddy-Szeri from
its registered-system list, emphasized the compatibility GUI before the native
Solver/Post-processor, and the AI guide preserved a stale count of 118 tests.
The agent harness also described an obsolete “future MCP” phase despite an
existing harness.

The new `docs/START_HERE.md` provides one five-minute human workflow and one
agent checklist. The documentation index now categorizes operating guides,
architecture, scientific workflows, RPE cases, Duffing cases, and operations.
Case studies are explicitly evidence rather than defaults. The architecture
contract and agent context now point to the measurable audit rather than
claiming implicit compliance.

Public-definition documentation remains weak at 49.8%. Priority is not to add
boilerplate everywhere; document the outward contracts first:

1. `compute_workbench` (about 7% public-definition coverage), because it is the
   principal UI/API boundary;
2. `campaigner`, `solver`, `job_manager`, and control public requests/results;
3. `performance` and supported experiment scripts;
4. internal helpers only when their behavior is non-obvious.

Docstrings should state coordinate conventions, array shapes, mutability,
failure modes, execution cost, and provenance effects—not merely repeat names.

## Repository and data-layout findings

The checkout also serves as a live scientific workspace. Approximate on-disk
measurements during the audit were:

| Area | Approximate size | Interpretation |
|---|---:|---|
| `archives/` | 3.97 GiB | canonical runtime/scientific data beneath the package root |
| `viewer_outputs/` | 232 MiB | generated/presentational artifacts, some historically tracked |
| `native/` | 28 MiB | source plus local build material |
| `build/` | 14.6 MiB | generated |
| `tmp/` | 11.2 MiB | generated experiments/intermediates |

The canonical archive location is valid for this local workspace, but source
tools and packaging must never walk it as code. `project_paths` now excludes
test dependencies, caches, temporary paths, and Viewer outputs in addition to
build products and archives. Longer term, prefer a configured external data
root for installed deployments and store only small manifests/fixtures in Git.
Do not relocate or delete current scientific data as part of an architecture
cleanup.

## Native application assessment

The native UIs conform better to the size policy than the Python orchestration
layer. Solver and Post-processor panels are separate, and the Viewer has
separate model, navigation, inspector, canvas, process, theme, and workbench
components. The compute application reuses Viewer implementation files by
compiling their `.cpp` files directly. This avoids duplication today but is a
build-time coupling: a shared static UI-support library would provide one set
of compile definitions, tests, and version ownership.

The Python suite is version `0.1.0`, optional native kernels are `0.2.0`, and
both UI CMake projects currently declare `0.1.0`. Independent component
versions are reasonable, but diagnostics and release notes must label them
explicitly; a bare “DynSysLab version” is ambiguous.

Recommended UI improvements after the architectural seams are stable:

- add a small in-app Help/About action opening `docs/START_HERE.md` and showing
  Python/UI/schema versions;
- show whether the catalog is current and which archive root it represents;
- preserve the plan-before-run gate and never hide execution behind navigation;
- test keyboard navigation, high-DPI behavior, empty archives, malformed
  payloads, and very large catalogs;
- build shared native components once through a common CMake target.

## Test and continuous-integration findings

Before this audit, `.github/workflows/native-wheels.yml` was the only GitHub
workflow and focused on the optional native wheel. A change could therefore
pass hosted checks without general Python tests, documentation-link checks,
CLI smoke tests, or an architecture audit.

The new `.github/workflows/quality.yml` now runs the full test suite, rejects
new architecture debt relative to the explicit baseline, publishes the JSON
audit, and builds the Python wheel. Continue extending it with:

1. a Python-version matrix rather than only Python 3.12;
2. source-distribution and installed-wheel smoke execution;
3. native UI compile/smoke jobs on Windows, distinct from native-kernel wheels.

Because three hard violations are grandfathered debt, CI should use a checked
baseline initially and reject regressions. Remove each baseline item when its
behavior-preserving split lands; then enable unconditional strict mode.

## Prioritized remediation sequence

### P0 — completed in this audit

- Establish a safe, discoverable `dynsyslab` front door and `doctor` output.
- Make legacy GUI help non-launching and lazy.
- Add a repeatable JSON/Markdown architecture audit.
- Index all narrative documentation and add a human/agent start guide.
- Remove producer imports of private Viewer session code.
- Expand generated-source exclusions without deleting artifacts.
- Add baseline-aware full-test/audit/wheel continuous integration.

### P1 — hard-policy refactors

1. Characterize and split `viewer/core.py` by source, renderer, and serializer.
2. Characterize and split `job_manager/manager.py` behind the existing façade.
3. Move the Duffing bilinear search engine into reusable ROM candidate services
   and reduce its script to CLI orchestration.
4. Add baseline-aware CI so no new hard violation is possible.

### P2 — package seams

1. Introduce neutral versioned schema and spectral value contracts.
2. Remove direct low-level archive imports from application packages.
3. Separate legacy CSV parsing from ROM controlled-dataset conversion.
4. Move campaign screening out of manifest design.
5. Replace Viewer-specific output imports in analytical packages with typed
   report/plot-data adapters.

### P3 — outward-facing quality

1. Raise public documentation coverage on `compute_workbench`, solver,
   campaigner, JobManager, control, and estimation.
2. Consolidate normal operations beneath the front door while retaining expert
   aliases.
3. Add in-app help/version/archive-root diagnostics.
4. Create a shared native UI-support CMake library.
5. Define release checks for Python suite, native kernels, UI schemas, docs, and
   artifact compatibility.

## Acceptance criteria for future architecture work

A refactor is complete only when:

- no maintained source file exceeds 1,000 lines and new files normally remain
  below 500;
- the baseline architecture audit has no unexplained new boundary violation;
- public APIs and versioned payload schemas remain compatible or have explicit
  migrations;
- full tests and targeted characterization fixtures pass;
- accepted scientific outputs agree within declared numerical tolerances;
- archives are not rewritten, duplicated, or silently reclassified;
- both native applications still load the same catalog and use the same
  execution services;
- human and AI start documents describe the actual commands and versions.

## Verification performed for this audit

- Full Python suite: all 405 tests passed in 86.98 seconds after the final
  changes. The 23 focused entry-point, event, and audit tests also passed in
  isolation.
- All maintained Python files compiled successfully with `compileall`.
- Viewer and Solver/Post-processor rebuilt successfully with MSVC/CMake.
- The native catalog and compute registry generated successfully with 998
  catalog records and 16 registered ODE systems.
- A Python wheel built successfully; its installed metadata exposes the new
  front door and architecture audit and contains the versioned debt baseline.
- The final baseline-aware audit passed with zero syntax errors, zero broken
  documentation links, zero unindexed narrative documents, and no new size or
  boundary regressions. Absolute strict mode correctly remains failing until
  the documented debt is removed.

## Bottom line

DynSysLab should not be rewritten wholesale. Its best qualities—scientific
provenance, explicit coordinate metadata, independent post-processing stages,
typed KROM/control artifacts, and shared archive-centered UIs—are the correct
foundation. The next architecture work should be a sequence of small,
characterized extractions around three Python monoliths and a few shared
contracts. That path improves flexibility and speed of development without
disturbing the control and reduced-order-model results that have already been
carefully validated.
