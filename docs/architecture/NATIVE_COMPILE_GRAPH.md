# Native compilation graph

**Status:** canonical policy and measured implementation (2026-09-05)

## Purpose

DynSysLab's scientific engine is compiled C++, but ordinary scientific and UI
work must not repeatedly instantiate Eigen's template-heavy algorithms. The
native build therefore separates stable runtime data contracts, scientific
orchestration, and the small set of implementation units that use Eigen.

This is a build boundary, not a scientific-data boundary. Concrete matrices
are always runtime objects.

## Runtime matrix invariant

No trajectory, Hankel matrix, KROM matrix, or control/estimation matrix is
baked into a compilation. This includes \(A\), \(B\), \(B_\Psi(x)\), \(C\),
\(Q\), \(R\), covariance, Gramian, SVD, and regression matrices.

Public numerical calls accept contiguous row-major `std::vector` buffers plus
runtime row, column, rank, and output counts. Matrix values and dimensions may
change on every call without rebuilding any code. The compiled numerical DLL
contains algorithms; it does not contain a user's concrete matrices.

Fixed-size Eigen matrix types, Eigen headers, and Eigen expression templates
must not cross a public DynSysLab header. All native DLLs use the same dynamic
MSVC runtime so standard-library values can safely cross their boundaries.

## Target structure

```text
Viewer / Solver / CLI
        |
workflow and payload services
        |
scientific orchestration DLLs
  eigenfunctions | spectral | decomposition | reconstruction
  control maps   | control (incl. estimation)  | linear analysis | PCA
        |
dynsyslab_numerics.dll
  SVD | least squares | eigensystems | nullspace | products
  CARE/LQR | matrix functions/Gramians
        |
Eigen (private implementation dependency)
```

`dynsyslab_numerics` owns the dense numerical interfaces and their Eigen kernels. Its
backends are split by operation and, where beneficial, by real/complex or
finite/infinite algorithm family. High-level modules include only the narrow
operation headers they use:

- `Svd.hpp`
- `LeastSquares.hpp`
- `NormalEquations.hpp`
- `Eigensystem.hpp`
- `Nullspace.hpp`
- `Product.hpp`
- `Control.hpp`
- `MatrixFunctions.hpp`

`LinearAlgebra.hpp` remains an external compatibility umbrella. Native
implementation files may not include it, because changing an aggregate header
would create unnecessary fan-out.

The causal-history extension also confines Eigen to two private targets:
`dynsyslab_observation` uses the FFT correlation kernel for frozen projection,
and `dynsyslab_history_koopman` applies matrix-free Hankel operators and fits
the reduced eigensystem. `HistoryCorrelation.hpp` is a private source header;
neither target exposes Eigen types in its public APIs. These structured kernels
avoid materializing the full delay matrix and are explicitly listed in the
architecture test. Ordinary system, controller, archive, and UI code still uses
runtime buffers and narrow numerical interfaces.

The classical linear-systems kernels (Lyapunov/Stein on the Schur form, QZ
transmission zeros, balanced truncation, DARE doubling, H-infinity) are
`EigenNumerics*` units of `dynsyslab_numerics`, with private helpers in
`EigenNumericsLinearSystems.hpp`. Their public API is the Eigen-free
`core/StateSpace.hpp` family ([LINEAR_SYSTEMS_TOOLKIT.md](LINEAR_SYSTEMS_TOOLKIT.md)).
Numerics has one Lyapunov/Stein solver and one finite-horizon gramian kernel;
the observability gramians in `MatrixFunctions.hpp` call them.

All three Eigen targets use a shared precompiled-header configuration and a three-job Ninja pool.
This avoids redundant parsing and memory contention while allowing ordinary
targets to continue compiling in parallel. The cache variable
`DYNSYSLAB_EIGEN_COMPILE_JOBS` permits a lower value on memory-constrained
builders.

## Measured result

Measurements used Release, MSVC 19.44, Ninja `-j4`, x64-windows, shared native
modules, and the same `dynsyslab_compute`, `dynsyslab_viewer`, and
`dynsyslab_cli` product targets. Times are machine-specific; graph counts are
the more stable architectural metric.

| Measure | Before | After |
|---|---:|---:|
| Clean product build (recorded reference) | 582.988 s | 552.005 s |
| Eigen-exposed compile commands | 93 / 134 (69.40%) | 13 / 147 (8.84%) |
| Targets exposed to Eigen | broad transitive product graph | `dynsyslab_numerics` only |
| No-op product build | — | 0.159 s |

The number of source units rose slightly because the old monoliths were divided
into smaller independent DLL modules and numerical kernels, but Eigen exposure
fell by 86.0% and is confined to one target.

The recorded clean build was 5.3% faster, but clean wall time is not the main
acceptance metric: additional samples on this laptop ranged from 452.612 to
615.932 seconds as thermal state and concurrent machine load changed. A clean
tree must instantiate every Eigen algorithm once. The goal is to stop normal
changes from accidentally becoming clean numerical builds, which the graph and
incremental measurements establish directly.

After correctness was settled, an isolated clean build of the complete
numerical backend took 278.913 seconds with three Eigen compiler jobs and
302.288 seconds with four. Three is therefore the measured default; four was
8.4% slower from memory/cache contention.

Representative settled incremental measurements:

| Edited unit | Wall time | Compiles | Eigen compiles |
|---|---:|---:|---:|
| RPE system implementation | 4.16 s | 1 | 0 |
| SINDy orchestration | 5.10 s | 1 | 0 |
| control-map regression orchestration | 5.21 s | 1 | 0 |
| LQR/KLQR wrapper | 5.39 s | 1 | 0 |
| Solver UI | 7.91 s | 1 | 0 |
| Eigen nullspace backend | 25.11 s | 1 | 1 |
| matrix-product interface | 54.62 s | 9 | 2 |

The product-interface case intentionally recompiles its two implementations
and seven direct consumers. It does not rebuild unrelated SVD, least-squares,
eigensystem, control, estimator, workflow, or UI code.

## Developer rules

1. Add systems and observables using the foundation/system interfaces; never
   include Eigen in a model implementation.
2. Put a new numerical algorithm behind an Eigen-free, runtime-sized public
   interface. Give a template-heavy implementation its own
   `EigenNumerics*.cpp` unit.
3. Include the narrow numerical operation header in internal code. Do not use
   `LinearAlgebra.hpp` inside the source tree.
4. Keep Eigen linkage `PRIVATE`; never expose an Eigen type in installed
   headers, archive schemas, API/CLI contracts, or UI code.
5. Never encode a concrete matrix dimension or value as a template argument or
   generated source merely to improve runtime performance. Optimize the
   generic runtime kernel instead.
6. Avoid editing `EigenPch.hpp` unless the backend genuinely needs a broadly
   shared Eigen facility; that file intentionally invalidates all Eigen units.
7. Measure before changing the Eigen job-pool limit. More compiler processes
   can be slower when memory bandwidth is saturated.

Architecture tests enforce rules 3 and 4. Native parity tests protect DMD,
H-DMD, SINDy, reconstruction, decomposition, control-map regression,
continuous/adaptive KLQR, and estimation behavior.

## Reproducible checks

Inspect the configured graph without rebuilding:

```powershell
./tools/Measure-NativeBuildGraph.ps1 `
  -BuildDirectory build/native-suite-fast
```

Measure one intentional incremental edit; the tool restores the original file
timestamp after the build:

```powershell
./tools/Measure-NativeIncrementalBuild.ps1 `
  -Source native/src/SystemsRPE.cpp `
  -Target dynsyslab_compute
```

Run the canonical native gates:

```powershell
cmake --build --preset windows-fast-release-gates -j4
ctest --preset windows-fast-release --output-on-failure
```

Build directories and generated JSON measurements are local evidence and are
not product artifacts.
