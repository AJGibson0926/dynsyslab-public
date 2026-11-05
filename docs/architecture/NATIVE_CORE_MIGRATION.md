# DynSysLab native-core migration

**Status:** complete (native 0.4.0 migration gates passed)

## Required outcome

DynSysLab must install and run as compiled desktop software without locating,
embedding, or launching Python. The canonical engine is a modular C++ library.
Viewer and Solver/Post-processor link it directly. A compiled CLI exposes the
same services for automation. Python remains useful only as an optional binding
and convenience CLI over the compiled API.

An installed application may contain compiled DLLs and data resources beside
the executable. “Standalone” means the user does not install Python and the
application never invokes Python source or a Python subprocess.

## Dependency rule

```text
desktop apps      compiled CLI      optional Python binding
       \               |                 /
        +--------------+----------------+
                       |
                 dynsyslab_core
                       |
   systems -> solvers -> archives -> post-processing -> ROM/control/estimation
```

Python may call down through the binding. The core must never call Python, and
no scientific algorithm may be maintained independently in both languages.

## Migration gates

Each slice is accepted only when:

1. its public C++ API is independent of UI and Python;
2. a Python-free CMake configuration builds it;
3. native unit tests pass;
4. fixtures demonstrate parity with the retained Python reference behavior;
5. the optional Python binding delegates to the native implementation;
6. archive schemas and scientific provenance remain compatible;
7. source files remain below 1,000 lines and normally below 500;
8. the native desktop workflow uses the slice directly before its Python
   implementation is retired.

## Work packages

| Package | Native ownership target | State |
|---|---|---|
| Foundation | value types, errors, JSON, paths, versioning, progress | native; stable C ABI and execution context added |
| Systems | drivers, parameters, observables, ODE/DDE/PDE/lattice models | ODE registry, Mackey--Glass DDE, 1-D heat equation, and Ising/Potts lattice systems native; one workflow resolver supplies every family to planning and UI parameter/campaign controls |
| Solvers | fixed/adaptive ODE, controlled solves, batches, PDE/DDE/lattice | Euler/RK4/RK45, CVODE/BDF, serial and locally partitioned structured-field PDE Euler/RK4, DDE method of steps, and Metropolis lattice execution native |
| Archives | NPZ/JSON metadata, hierarchy, immutable publication, indexes | native NPZ/JSON, hierarchy, equation fingerprints, immutable publication, and tagged legacy DynBatch CSV import with SHA-256 provenance |
| Post-processing | SVD/POD/PCA, bagged/robust PCA, Hankel POD, DMD/H-DMD/EDMD, SINDy, HAVOK | native |
| ROM | spectra, eigenfunctions, KROMs, extraction, reconstruction, errors | analytic observable-coordinate KROMs, self-contained independent extraction workflow/UI, C maps, and map errors native |
| Control | targets, LQR/SDRE/KLQR/MPC/PID, control maps | native classical LQR, sampled local SDRE/PID/MPC, continuous state-dependent KLQR, sampled KMPC, and exact/constant/bilinear Koopman maps |
| Estimation | Luenberger/LQE/KLQE and observability | native classical observability/SVD, finite discrete/continuous and stable infinite Gramians, and continuous-time LQE/KLQE |
| Operations | campaigns, jobs, planning, progress, recovery | system-agnostic native cartesian campaign executor, KROM training-domain planner and pre-publication safety screen, compiled registry/resource planner, monotonic progress, durable SQLite queue service shared by UI/CLI, tagged legacy trajectory import, exact-duplicate audit, and recoverable deduplication |
| Viewer/services UI | native catalog, payloads, export, structured setup | metadata-indexed hierarchy/lineage catalog, lazy payload, standard five-panel KROM dashboards, synchronized pairwise comparison, numeric axis controls, native all-family parameter/initial-condition, KROM training-domain, controller, cartesian-campaign, and legacy-import editors, structured KLQR/control-map/evaluation and MPC controls, structured H-DMD → SINDy, KROM extraction, control-regression, decomposition, observability, and KLQE configuration, reusable native job files, durable queue management, bidirectional Post-processor navigation, clipboard and disk export |
| Desktop distribution | no-argument launch, resources, installer, taskbar pin | unified dynamic build and app-local package gate implemented |
| Python | thin optional bindings and CLI, no duplicate algorithms | distributable wheel is dependency-free and contains only `dynsyslab_client`, which delegates workflow/API and launch operations to compiled contracts; historical reference executors are source-only, excluded from the product, and non-authoritative |

Complete native UI parity is a release gate, not a follow-up usability item.
The source-backed, machine-readable ledger in
[`NATIVE_UI_PARITY_CONTRACT.json`](NATIVE_UI_PARITY_CONTRACT.json) maps every
retained human workflow to its compiled service, desktop surface, and regression
evidence. See [`NATIVE_UI_PARITY.md`](NATIVE_UI_PARITY.md) for the policy and
audit command.

## Completion evidence

The canonical product is now the compiled core, workflow services, CLI, Solver/
Post-processor, and Viewer. The optional Python wheel contains only the thin
client. The unified native release gate builds every registered CTest binary,
and the installed-package gate rejects Python artifacts, initializes both
desktop executables from the isolated install tree, executes representative
ODE/DDE/PDE/lattice and control/post-processing workflows, and builds and loads
an independent SDK extension. The native UI parity ledger maps every retained
CLI scientific capability and human workflow to a compiled desktop surface.

Migration completion does not freeze the scientific product. New systems,
solvers, post-processors, displays, and research workflows remain ordinary
incremental feature work behind the same DLL/API/UI boundaries.

The numerical implementation is additionally isolated behind runtime-sized,
Eigen-free contracts so those ordinary changes do not trigger template-heavy
recompiles. Concrete scientific matrices remain runtime objects. See
[`NATIVE_COMPILE_GRAPH.md`](NATIVE_COMPILE_GRAPH.md) for the target graph,
developer rules, and measured clean/incremental build results.

## Retained KLQR case-study parity

The accepted Duffing transfer/release and RPE moderate-multisine experiments
now run as a native CTest gate. It exercises exact, constant, bilinear, and
conjugate-complete bilinear `B_Psi` representations with continuous controller
evaluation inside adaptive BDF integration. The native outputs reproduce the
retained release states, final states, control magnitudes, and RPE tracking
errors. Legacy Fortran-ordered arrays, discrete-to-continuous eigenvalue
conversion, explicit conjugate duals, and independently archived control maps
are handled at the native archive boundary. See
`case-studies/NATIVE_KLQR_PARITY_REPLAY_2026-09-05.md` for the numerical record.

## First completed slice

The former pybind-only kernels were split into a normal `dynsyslab_core`
library and a separate adapter. A Python-free build now provides:

- Euler/RK4 single and batched integration for Duffing, Van der Pol, and
  dimensionless RPE;
- autonomous diagonal KROM rollout;
- 2-D Ising Metropolis execution;
- `dynsyslab_native_cli.exe` capability and self-test commands;
- native CTest coverage.

The optional wheel links exactly this library. Existing Python parity tests
continue to pass through the binding.

## Second completed slice

The canonical system contract and general ODE executor are now native:

- all 16 registered ODE systems have C++ factories, parameters, initial-state
  overrides, Jacobians, and control-input maps;
- the two-state Marmottant implementation retains its discrete hysteresis
  latch and resets it between trajectories;
- RPE, Marmottant, and Reddy--Szeri execute in their established
  nondimensional coordinates, including state-dependent pressure-input maps;
- fixed Euler/RK4 and adaptive Dormand--Prince RK45 execute arbitrary state
  dimensions and trajectory batches;
- controlled integration records a control history beside every solution;
- progress callbacks are throttled to percentage or 250-ms heartbeats;
- native CTest covers system, solver, batch, control-history, and coordinate
  contracts, while the optional Python binding delegates to these APIs.

## Current compiled scientific path

The Solver/Post-processor now executes native solves and native spectral jobs
directly. CVODE supplies variable-order BDF integration through a modular
SUNDIALS runtime and accepts the `BDF` and compatibility `ode15s` method
names. The native application writes immutable trajectory archives containing
the state, target, control history, physical configuration, coordinate system,
and solver provenance. H-DMD and SINDy consume those archives without Python.

Native polynomial KROMs own their eigenvalues, eigenfunctions, physics
configuration, reconstruction data, and any exact, constant, or bilinear
control-map representations. Extracted sub-KROMs are independent objects:
compatible regressed matrices are copied into the reduced object rather than
resolved through a parent at runtime. Continuous-time KLQR preserves the
established law `u = -K(Psi(x) - Psi_target)` and writes the applied control in
the same dimensional or nondimensional convention as the dynamics.

Control-map identification is algorithmically explicit in the compiled core
and desktop editor. Constant maps may be fitted from continuous lifted
residuals, exact pointwise maps, or the Korda--Mezić multistep observable-output
convolution using a selected reconstruction `C`. The latter preserves the raw
discrete `B_d`, its sample interval, and the converted continuous `B_c` used by
KLQR. Bilinear maps independently select controlled-residual or exact-map
projection fits, including an optional conjugate-complete feature basis.

The compiled workflow also owns bounded KROM control-candidate campaigns.
Read-only virtual screening ranks system-agnostic mode subsets before any fit;
up to the hard limit of 25 real attempts then extract independent sub-KROMs,
pool selected controlled datasets, regress constant and bilinear maps, and run
matched exact/constant/bilinear KLQR cases. Progress from nested extraction,
regression, and adaptive solves is mapped monotonically into the campaign
progress bar. Durable per-attempt ledgers support interruption diagnosis, while
the immutable campaign archive drives a native Viewer comparison dashboard.
The Solver exposes every campaign input plus plan, save, queue, and run actions,
so the capability does not depend on raw JSON or a command-line-only path.

Classical state-space control and Koopman control are separate public APIs.
Fixed LQR solves one CARE about a declared linearization point. The local
state-dependent mode refreshes the system Jacobian and input map and solves a
CARE on an exact sampled clock. PID and finite-horizon state-space MPC use the
same per-trajectory sampled-clock boundary. KLQR remains Koopman-space and,
by default, evaluates continuously at every adaptive solver RHS state: it
recomputes `B_Psi(x)`, the CARE, and `K(x)` on that state. KMPC is sampled.
No classical controller replaces or bypasses the KLQR/KMPC KROM pathway.

The native system boundary also owns model/physical time and input scaling.
Prescribed and built-in bubble forcing is supplied in Pa, Hz, and seconds, converted once
to model input/time, and archived with reversible scale metadata. Explicit
model-coordinate settings remain available for nondimensional experiments.
This removes the former native/Python discrepancy without changing KLQR:
KLQR still returns model-coordinate control directly and still evaluates
continuously at every adaptive RHS state by default.

The release uses coarse DLL boundaries for foundation, Koopman, control,
operations, and workflow services. App-local packaging copies SUNDIALS,
libarchive, SQLite, compression, and Microsoft runtime DLLs beside the public
executables. The isolated-package gate rejects Python runtimes and performs a
real solve, proving that the installed application neither locates nor launches
Python.

## ABI, workflow, and scientific identity slice

JSON/path handling and all job execution now live below the UI in a neutral
workflow DLL. A C++ API serves the native applications and compiled CLI. Stable
C ABIs use opaque handles, caller-owned buffers, explicit status codes, and
run-once workflow result handles for optional MATLAB, Python, or third-party
adapters. Cancellation, structured progress, thread limits, and memory limits
cross the same execution boundary.

Instantiated systems retain fully resolved parameters with units,
dimensionless groups, equation-version tags, and deterministic equation
fingerprints. Hierarchy grouping therefore identifies actual equations rather
than just preset names and sparse overrides. A versioned extension ABI can add
ODE systems, numeric observables, and generic controller/postprocessor
capabilities from separately compiled DLLs; a real extension-load test guards
the boundary.

The neutral workflow dispatcher owns job validation, extension loading, and
solve/campaign/post-process routing. UI code, the compiled CLI, and the workflow
C ABI all call this one dispatcher rather than maintaining parallel execution
paths. A job may list extension DLLs explicitly; post-processing stages may
route to an extension capability named `postprocessor:<method>`.

Every ODE, DDE, lattice, and PDE trajectory produced inside a campaign carries
the campaign ID, campaign kind, point index, controller configuration, and
hierarchy-manifest path in its own archive metadata. The campaign and hierarchy
manifests remain authoritative relationship records, while individual archives
are still self-identifying when copied or inspected independently.

Archive management is also available without Python. The operations DLL scans a
bounded selection of canonical NPZ records, groups repeated numerical plans,
and verifies exact solution duplicates by descriptor, shape, and bytes after a
fast fingerprint prefilter. Audit and deduplication are dry-run operations by
default. Explicit apply moves only unreferenced exact duplicates into an
archive-local recovery tree and writes a recovery manifest; it never unlinks a
record or silently rewrites provenance. Reference-relinking remains an explicit
compatibility operation until its native semantics have separate acceptance
tests. Operations is no longer linked through the scientific core aggregate,
preventing archive or job changes from invalidating workflow and Viewer objects.

## Build-isolation and linear-algebra slice

Routine release builds now target only Solver/Post-processor, Viewer, the CLI,
and the compiled linear-algebra service. Test executables are built by the
separate `dynsyslab_test_binaries` gate. Ninja is the default developer
generator, with Visual Studio retained as a supported fallback. On the measured
Windows host, a no-change Ninja build completes in about 0.35 seconds. Cold
builds are still substantial because the current core contains several
Eigen-heavy translation units; an Eigen precompiled-header experiment made the
cold build slower (about 823 seconds versus about 725 seconds) and was removed.

`dynsyslab_linalg.dll` now provides a versioned row-major C ABI for real SVD,
real and complex least-squares (including ridge regularization), and real-matrix
eigensystems. Eigen is compiled inside that DLL. A new system or independent
post-processing method can therefore be an extension DLL that includes only
small DynSysLab ABI headers and calls the compiled numerical service when it
needs dense linear algebra. The helper `dynsyslab_add_extension()` standardizes
such independently compiled targets.

Workflow support code was separated from the algorithm dispatcher so every
maintained native source file remains below the preferred 500-line ceiling.
Execution contexts now reject projected solution buffers that exceed the
configured memory budget, complementing cancellation, progress, and thread
limits. The complete native CTest set and isolated, Python-rejecting packaged
solve pass for this slice.

## Delay-system slice

Mackey--Glass is the first non-ODE family routed through the canonical native
workflow. Its defaults (`beta=0.2`, `gamma=0.1`, delay `17`, exponent `10`, and
constant history `1.2`) preserve the retained reference behavior. A dedicated
`dynsyslab_delay.dll` owns the model and method-of-steps Euler integrator so
subsequent delay-system work does not enlarge the ODE or UI translation units.

The system catalog exposes an explicit `dde` family, while the generic archive
hierarchy groups the resulting trajectories by fully resolved equations just
as it does ODE data. Delay archives retain requested and grid-effective delay,
history convention, dimensionality, native backend provenance, and immutable
hierarchy identity. Default, explicit, uniform, and normally sampled constant
histories are supported. DDE parameter/time/IC campaigns route through the
same cartesian dispatcher as ODEs. DDE control remains an explicit unsupported
operation; it is never routed through an ODE or Python fallback.

The isolated package gate performs a real two-trajectory Mackey--Glass solve
from the installed application, in addition to the ordinary Duffing solve and
independently compiled SDK extension test.

## Lattice-system slice

The 2-D Ising kernel was removed from the foundation DLL and joined by native
3-D Ising and q-state 2-D Potts implementations in `dynsyslab_lattice.dll`.
All use periodic boundaries and single-site Metropolis updates. The catalog
reports them as the separate `lattice` family with magnetization/order
parameter and energy-density observables.

The canonical workflow supports seeded random or explicit initial lattices and
trajectory ensembles. Every solution archive stores the observable time
series, exact integer initial and final lattice snapshots, acceptance
statistics, Hamiltonian parameters, RNG identity, and hierarchy provenance.
Native seeds reproduce a native run exactly; streams are intentionally marked
statistical rather than bitwise parity with NumPy because the canonical engine
uses `std::mt19937_64`.

Lattice parameter/time/IC campaigns route through the same cartesian workflow
as ODE and DDE jobs; the installed-package gate exercises a real Potts
temperature campaign. Lattice control remains an explicit unsupported
operation. No request falls back to the retained Python implementation. The
lattice test is part of the ordinary native CTest gate rather than a
special-purpose migration script.

## PDE-system slice

`dynsyslab_pde.dll` owns a UI-independent structured-mesh/field integration
API with separate internal and output steps, exact output coordinates, Euler
and RK4 stepping, finite-field checks, cancellation, progress, and memory
preflight. Its generic callback contract preserves arbitrary field shapes and
matches the retained decay-field RK4 reference benchmark.

The same DLL owns deterministic Cartesian mesh partition plans and
boundary-aware ghost construction for edge, periodic, and zero-value axes.
Subdomain callbacks receive only their padded local field and write only their
interior derivative. A bounded local worker pool evaluates independent
partitions at each synchronized explicit stage. Tests cover uneven partition
cuts, neighbor identity, overlap rejection, multidimensional fields with
trailing components, and bit-for-bit serial/partitioned heat evolution. MPI
transport and adaptive mesh refinement remain later, optional backends; the
local path does not imply distributed execution.

The first registered PDE is a one-dimensional heat equation with periodic,
zero-flux, or zero-value boundaries and sinusoidal initial-field presets. Its
workflow accepts explicit and seeded random field ensembles, enforces a
conservative explicit stability bound, publishes physical parameters and
spatial coordinates, and participates in the common hierarchy/campaign
dispatcher. Viewer consumes the same archive directly as a downsampled
space-time heat map and final spatial profile. The installed-package gate runs
a real partitioned heat-equation solve without Python.

## Observable and reconstruction slice

Every native ODE exposes a catalog containing its state components and full
state. Duffing, RPE, Marmottant, and Reddy--Szeri additionally expose their
scientific observables; Reddy--Szeri defines the `c_k2` and `c_k5` coordinate
presets. An observable owns an explicit value map and either an analytic
Jacobian or an explicit finite-difference fallback.

H-DMD accepts any named real scalar observable. Its archive retains discrete
and continuous spectra, eigenfunction values along the Hankel snapshots, and
the observable reconstruction row `C`. SINDy may fit polynomial
eigenfunctions in arbitrary real observable coordinates. Exact Koopman control
then evaluates the chain rule
`nabla_x Psi = nabla_y Psi J_h(x)` before applying the system input matrix.
When `h(x)=x`, `J_h` is the identity and the ordinary formulation is recovered.

The generic projection API solves `G ~= Psi C^T` for scalar or vector
observables. Named reconstruction maps are stored inside a KROM archive and
their selected columns are copied into independently extracted sub-KROMs.

## Estimation slice

The native estimator consumes those named reconstruction matrices directly.
It converts complex Koopman modal coefficients to the same real/imaginary
split coordinates used by KLQR, constructs the classical observability matrix,
and retains its singular values and right singular directions. It can instead
construct finite discrete, finite continuous, or stable infinite-horizon
observability Gramians. Gramian eigenvectors and square-root eigenvalues use
the same ordered direction/strength contract as the matrix SVD. Continuous-
time LQE is solved through the dual stabilizing CARE with explicit
detectability, residual, and observer-pole checks.

`runKoopmanLQE` propagates the observer exactly for zero-order-held
measurements, retains lifted estimates, reconstructed physical state,
measurement predictions, innovations, and benchmark errors, and never uses
ground truth in the estimator update. The compiled Post-processor exposes
`observability` and `klqe` stages; its archives render natively in Viewer.
