# DynSysLab package-boundary contract

The suite is a set of independently useful scientific services over one
versioned archive and hierarchy model. The canonical services live in the
independently compilable C++ modules; native applications link them through a
neutral workflow service. Stable C ABIs serve external tools without exposing
C++ ownership across DLL boundaries. Python is an optional thin binding and CLI layer, not
an alternate scientific implementation or an application runtime dependency.

## Intended dependency direction

```text
native UI / compiled CLI / optional language adapters
                         |
                stable C++ and C APIs
                         |
                 neutral workflow
                         /|\
        +----------------+----------------+
        |                |                |
   foundation         Koopman         operations
        |                 \              /
     systems          control / estimation
        |                   |
     solvers             archives
        ^                   ^
        |                   |
 versioned extension DLLs   +--- read-only Viewer
        |
 compiled linear-algebra ABI (when a plug-in needs SVD/QR/eigensystems)
```

`infrastructure` contains only small neutral mechanisms, such as passive
source-change events, that would otherwise force a producer to import a UI.

## Package ownership

The independent **dynsyslab_lie** library owns exact polynomial geometric analysis and invariant-module certificates. It depends only on foundation, not on solvers, control, archives or UI; workflow exposes read-only queries.

The native [pipeline boundary guide](docs/architecture/PIPELINE_BOUNDARY_REFACTOR_2026-10-31.md)
specifies scientific binding, frozen measurement transforms, runtime/fitting
targets and coordinate-driven control/estimation. The implementation preserves
the original feedback laws and exposes the same workflow to UI and CLI clients.

The [measured-output feedback guide](docs/architecture/KOOPMAN_OUTPUT_FEEDBACK.md)
defines the optional KLQE/KLQR loop. Estimation and control remain independent
libraries; workflow coordinates sampled measurements, estimated Koopman values
and actually applied inputs through the solver's existing control policy.

- **dynsyslib** defines systems, presets, observables, Jacobians, and drivers.
  It does not solve, schedule, view, or store numerical outputs.
- **solver** executes one requested numerical solve. It writes completed or
  checkpointed results only through the Archiver-owned archive contract.
- **archives** is private storage/index implementation. Its Archiver-owned
  `management/archive_index.sqlite` is the canonical incremental metadata,
  content-digest, and relative-pointer summary. Producers update it only after
  successful publication; Viewer opens it read-only and never rebuilds it.
  New application code
  uses **archiver** gateways rather than low-level storage APIs.
- **archiver** owns archive capabilities, audit, reuse lookup, campaign
  discovery, storage locations, and durable diagnostics.
- **campaigner** designs manifests. Numerical screening/execution belongs in a
  solver or operations service, not in manifest construction.
- **job_manager** schedules manifests and stores only coordination metadata;
  it never owns numerical solution arrays.
- **postprocessing** transforms supplied or archived solution data into derived
  scientific objects without mutating source trajectories.
- **rom** owns reduced-model contracts and algorithms, not archive files.
- **control** and **estimation** are optional solver plug-ins and must not alter
  an uncontrolled solve.
- **viewer** is read-only relative to the primary archive. Disposable display
  payloads live only under `viewer_compute_archive`.
- **compute_workbench** and native UI code compose the public APIs above through
  versioned JSON contracts; they do not contain scientific algorithms.
- **initial-condition generation** belongs to the neutral workflow service.
  `resolveOdeInitialConditions` is the single C++ API used by execution and UI
  previews; `dynsyslab-cli generate-ics` is its machine-readable CLI gate. The
  editable `cluster_specs` array permits independent deterministic generators
  per Cluster while Batch/System/Cluster/Trajectory remains a UI and archive
  hierarchy concern.
- **quality** measures policy conformance without importing numerical engines.
- **extensions** add systems, controllers, and post-processors as independently
  compiled DLLs through a versioned C ABI. They do not include private core
  headers or force established scientific modules to rebuild.
- **dynsyslab_linalg** owns the template-heavy Eigen implementation behind a
  versioned, row-major C ABI. Extension clients can request SVD, least-squares,
  and eigensystem operations without including Eigen themselves.

## Shared contracts

- Persisted formats require a `format` or `schema` version. A new version must
  provide a migration reader or an explicit compatibility error.
- The hierarchy is `Batch -> System -> Cluster -> Trajectory -> TimeSeries`.
- `System + physical parameters` identifies a dynamical variant. Free/forced,
  campaign membership, training use, and control use are data-role metadata,
  not different systems.
- Dimensionality and coordinate scaling are explicit provenance. A controller
  and control history inherit the dynamics' coordinate convention.
- Completed source records are immutable. Derived results refer to their source
  artifacts instead of copying or rewriting them.

## Module policy

Reusable behavior receives a documented C++ API, a stable C ABI wherever it
crosses a tool or language boundary, and a lightweight compiled CLI
gate before UI coupling. Optional Python bindings expose the same core contract
without reimplementing algorithms. Source files must remain under 1,000 lines;
500 lines is the normal design ceiling. Routine builds compile only public
applications and affected runtime DLLs; test binaries are a separate gate
target. Adding an independent system or method should normally compile one
extension DLL, not rebuild the suite. During migration, run the native CTest
and isolated-package gates; the legacy Python audit is advisory until it is
replaced by a compiled quality command.

See [Start here](docs/START_HERE.md),
[the engineering policy](AGENTS.md), and
[the native-core migration plan](docs/architecture/NATIVE_CORE_MIGRATION.md),
the dated architecture audit (private development record `DYNSYSLAB_SUITE_ARCHITECTURE_AUDIT_2026-08-31.md`).
The audit roadmap's completed implementation and final gates are recorded in
the [architecture remediation report](docs/architecture/DYNSYSLAB_SUITE_ARCHITECTURE_REMEDIATION_2026-08-31.md).
