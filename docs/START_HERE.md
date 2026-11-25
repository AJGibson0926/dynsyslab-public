# Start here: DynSysLab

**Central research question:** finite Koopman representations of drift/control Lie algebras with global state reconstruction and quantified control errors. The [native Lie algebra module](architecture/LIE_ALGEBRA.md) provides exact closure, structure, invariant-module and residual-bound queries.

**Critical Koopman reference:** preserve and revisit the
exact unstable linear Duffing example (private development record `KOOPMAN_LINEAR_DUFFING_EXACT_REFERENCE.md`).
The user designated it critically important on 2026-11-05; it records exact
eigenfunctions, input/output maps and the successful position-only KLQE/KLQR run.

**Critical gradient-closure reference:** preserve the exact nonlinear triangular
system (private development record `KOOPMAN_GRADIENT_CLOSED_TRIANGULAR_REFERENCE.md`,
reference ID `koopman-triangular-gradient-closed`). Its finite eigenfunction
collection closes under every partial derivative and gives exact global state
reconstruction and bilinear control. The complementary user-designated important
arctangent case (`KOOPMAN_ARCTANGENT_GRADIENT_REFERENCE.md`) proves why a finite
autonomous linearization need not have finite gradient closure. Both are indexed
in the private study records.

The [quadratic triangular benchmark](QUADRATIC_TRIANGULAR_BENCHMARK.md) is available
as a native system with exact solutions, analytic eigenfunctions and derivative
operators. Its stable and saddle presets support finite-versus-infinite closure studies.

This is the shortest safe route into DynSysLab for a scientist, developer, or
AI agent. The compiled CLI, native Solver/Post-processor, native Viewer, and
optional Python adapter all operate on the same native system definitions and
canonical archives. There is no separate UI data store or required Python
runtime.

## Five-minute human start

For an ordinary installation, follow [Install DynSysLab](INSTALL.md). The
prebuilt release includes both desktop apps and the CLI; no build tools are
required. Anyone can start from `AJGibson0926/dynsyslab-public` alone. Running
DynSysLab means `dynsyslab-cli quickstart`: the overview, the bundled examples,
their results and figures, as Markdown, without opening a window. All data is generated
locally: the release ships no solver output, and every shipped example creates its own data in
your archive folder.

For source development, install and launch from the repository root:

```powershell
.\tools\Install-DynSysLabNative.ps1 -PinToTaskbar
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\DynSysLab.exe"
```

The installer runs the native release gates and creates a stable Start-menu
entry; `-PinToTaskbar` creates a shortcut to the installed executable. The
installer records the canonical archive root for the compiled apps, and the
app's left tree is the same system/parameter/data hierarchy used by Viewer.

In the workbench:

1. Select a system and physical-parameter preset.
2. Select or generate initial conditions and set the trajectory/cluster count.
3. Select free dynamics or configure a controller/campaign.
4. Use **Plan / validate** before **Run and archive**.
5. In Post-processor, select an existing solution and configure H-DMD, SINDy,
   HAVOK, DMD, POD, PCA, or SVD without rerunning the dynamics.
6. Open Viewer to inspect the full compact `system → physics → major lane → study` data
   tree. Viewer initially selects a fully contextualized training set, preferring
   the newest peer when several are available.
   Numbered IC, KROM, Bψ, error, control, and estimator records retain their
   complete archive labels and provenance in the Inspector.

Use **Import data…** for a tagged legacy DynBatch CSV. The native dialog
previews the hierarchy and preserves the external dynamics identity unless you
explicitly map it to a registered system and parameter set.

Use **KROM control search…** in Post-processor for the retained bounded
control-candidate experiment. It
selects a parent KROM and controlled regression ensemble, screens mode subsets,
and compares exact, constant, and bilinear KLQR control for no more than 25
fitted candidates. The complete job can be planned, saved, queued, run, and
reviewed in Viewer without Python or hand-editing JSON.

Use `DynSysLabViewer.exe` to open Viewer directly. The old Python/Tk GUI is a
frozen source reference and is not part of the native distribution.

## Native CLI and optional Python adapter

The installed, Python-free CLI provides the canonical workflow gate:

```powershell
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" capabilities
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" doctor
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" systems
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" validate job.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" run job.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" fourier trajectory.npz --component 0
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" archive audit --all-ages
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" archive deduplicate --all-ages
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" queue submit examples\native\duffing-free-solve.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" queue list
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" import legacy-csv data.csv --label experiment
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" phase-space examples\native\radial-limit-cycle-poincare-query.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" phase-space examples\native\quadratic-fold-bifurcation-query.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" observable-closure examples\native\rpe-echo-closure-query.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" observable-closure examples\native\duffing-position-global-jet.json
```

`systems` reports each model's family as well as its state dimension and
presets. Native ODE models, the Mackey--Glass `dde` model, the one-dimensional
heat-equation `pde` model, and Ising/Potts `lattice` models use the same solve
job schema and archive hierarchy.
Unsupported family-specific operations are rejected explicitly; no native
command falls back to Python.
`fourier` reads the selected archived `TimeSeries`, performs a lazy native
O(N log N) transform, and returns frequencies and amplitudes as versioned JSON.
It never modifies the archive or stores redundant spectral arrays.
Queue operations return versioned JSON for reliable UI and automation clients;
`queue cancel <job-id>` cancels a queued or running record explicitly.
`phase-space` can compute state-defined first-return maps, projected section
stability, periodic shooting, and pseudo-arclength equilibrium or cycle
branches. See [Poincare sections and bifurcation
sidecar](POINCARE_AND_BIFURCATION_SIDECAR.md) for the numerical contract and
candidate-event limitations.
`observable-closure` v2 computes global certificates and transformed equations from
autonomous equations of motion. Its v1 query performs a bounded numerical search for observable-jet
rank loss, coordinate collisions, and closure obstructions. Its output
describes the scope of that numerical search. See
[Coordinate datasets and dynamics SINDy](architecture/COORDINATE_DYNAMICS_SINDY.md#observable-jet-closure-side-car).

The optional installed Python client is a dependency-free adapter over those
compiled executables and the workflow C ABI:

```powershell
python -m dynsyslab_client --help
python -m dynsyslab_client doctor
python -m dynsyslab_client systems
python -m dynsyslab_client plan job.json
python -m dynsyslab_client run job.json --archive-root C:\data\dynsyslab-archives
python -m dynsyslab_client archive audit --all-ages
python -m dynsyslab_client queue list
```

The client accepts the same saved job documents as the desktop workbench. It
does not validate, schedule, solve, process, or archive independently: those
operations remain inside the compiled CLI/workflow ABI. Historical Python
research packages are source-only migration references and are not included in
the client wheel or either application. The compiled
archive audit is read-only and verifies candidate duplicates byte-for-byte;
deduplication is also a dry run unless `--apply` is supplied. Applied cleanup
moves only unreferenced exact duplicates into recoverable archive-local trash;
the optional `dynsyslab archive audit|deduplicate` front doors delegate to
these compiled operations. See
[native build and installation](BUILD_NATIVE.md) for the unified CMake and
portable package workflow.

## Archive and execution safety

- Planning and validation are read-only. Execution requires an explicit Run
  action or `--execute` gate.
- Completed numerical records are immutable. Extensions and derived products
  create new records with provenance rather than overwriting source data.
- Application packages use Archiver capabilities. The low-level `archives`
  package is storage implementation, not a general application API.
- Viewer is read-only relative to the primary archive. Its temporary payloads
  live under the ignored `viewer_compute_archive` directory.
- Dimensionless systems, targets, controllers, and control histories must share
  the same coordinate convention and record it in metadata.

## AI-agent start

Before changing code, an agent should:

1. Run `dynsyslab-cli doctor` and `dynsyslab-cli capabilities`; use the optional
   Python `doctor` only when working on the adapter itself.
2. Read [the engineering policy](../AGENTS.md),
   [the architecture contract](../ARCHITECTURE.md), and
   [the maintainer context](AI_CONTEXT.md).
3. Run `python -m quality --format markdown` to see current structural
   debt rather than assuming the suite is already policy-complete.
4. Use public API/CLI gates. Do not place algorithms in UI event handlers or
   create a parallel archive/catalog.
5. Preserve unrelated working-tree changes and scientific artifacts.
6. Add focused tests and update the relevant user/architecture document.

MATLAB sources are read-only evidence unless the user explicitly authorizes a
wrapper or adaptation. Do not infer permission to copy legacy behavior.

## Where to go next

- [Solver and Post-processor guide](SOLVER_POSTPROCESSOR.md)
- [Native observable derivatives and echo control](architecture/NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md)
- [Tutorial and Python API examples](TUTORIAL.md)
- [Documentation map](README.md)
- [AI context and package contracts](AI_CONTEXT.md)
- [Native Viewer architecture](architecture/VIEWER_CPP_WORKBENCH.md)
- [Functions on their domains and coordinate maps](architecture/DOMAIN_FIELDS.md)
- [Native Solver/Post-processor architecture](architecture/SOLVER_POSTPROCESSOR_CPP_WORKBENCH.md)
