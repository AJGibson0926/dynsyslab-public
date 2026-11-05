# Performance and native execution

DynSysLab keeps its public API, configuration, orchestration, archives, and
Viewer in Python.  Numerical acceleration is attached behind solver-owned or
analysis-owned backend interfaces and must not leak into those boundaries.

`SolverOptions.execution_backend` accepts:

- `auto`: use a compatible optional native kernel when installed, otherwise
  use the Python/SciPy reference implementation;
- `python`: require the Python/SciPy implementation;
- `native`: require native support and fail clearly rather than falling back.

Every numerical result records both the requested and selected backend in its
metadata.  `solver.available_solver_backends()` reports capabilities without
requiring the extension to be installed.

Run `dynsyslab-backends` for a JSON diagnostic of installed solver and ROM
backends plus effective numerical thread environment variables.
The source-tree equivalent is `python -m solver`.

Run the reproducible baseline with:

```text
dynsyslab-benchmark --iterations 5 --warmup 1 --output benchmark.json
```

List or select focused cases when a full numerical sweep is unnecessary:

```text
dynsyslab-benchmark --list-cases
dynsyslab-benchmark --case sindy-lazy-polynomial-values-5000x256 \
  --iterations 5 --warmup 1 --output sindy-values.json
```

`python -m performance` is the equivalent package entry point.

The report includes Python, platform, NumPy, SciPy, and NumPy runtime/SIMD
information, plus structured BLAS/OpenMP thread-pool identities and active
thread counts. Timing comparisons are valid only when their runtime metadata
and workload definitions are comparable.

Compare a candidate run to a versioned baseline with explicit admission or
regression thresholds:

```text
dynsyslab-benchmark --baseline baseline.json \
  --require-speedup native-duffing-rk4-20k-steps=2.0
```

The command exits with status 2 when a gate fails and refuses comparisons
across different platform, machine, Python, NumPy, or SciPy identities.

Admission comparisons within the same run use paired case names:

```text
dynsyslab-benchmark \
  --require-relative-speedup native-duffing-rk4-20k-steps/python-duffing-rk4-20k-steps=2.0
```

This gate also exits with status 2 on failure. It is the appropriate check for
deciding whether an optional implementation may replace its reference pathway.

The dimensionless RPE fixed-step gate uses the same 25,000-step workload as
the standard `tau=25`, `dtau=0.001` free-response trajectory:

```text
dynsyslab-benchmark \
  --require-relative-speedup native-rpe-dimensionless-rk4-25k-steps/python-rpe-dimensionless-rk4-25k-steps=20
```

On the 2026-08-21 Windows reference run, native 0.2.0 completed this kernel in
0.00411 s versus 2.918 s for the Python reference (709.7x). Numerical parity
is tested separately; timing never substitutes for correctness.

Native kernels are optional private implementation details.  Each one must
retain a tested Python fallback, preserve public request/result types, pass
numerical parity tests, avoid filesystem and application state, and produce a
meaningful end-to-end improvement on a representative benchmark before it is
selected by `auto`.

The admitted RPE kernel is explicitly dimensionless. Its native model record
contains the RPE groups, pressure scale, natural-period time scale, and the
physical sinusoidal-driver parameters. Native evaluation converts model time
`tau` back to physical seconds only for driver phase and preserves driver
amplitude in Pa. Dimensional RPE, tabulated/custom drivers, controllers,
in-solve callbacks, pause, and cancellation remain on the Python reference
path.

KROM backend admission must also preserve coordinate contracts. A dimensionless
KROM evolves on model time even when its controller-facing adapter accepts
physical seconds; compiled rollout kernels receive already-normalized model
time and model state. Backend code must not infer units from array magnitudes.

Dense H-DMD already executes in compiled SciPy LAPACK; moving the Python call
itself to C++ does not replace the underlying SVD. `hankel_dmd` and the
standard CLI therefore expose `svd_threads` / `--svd-threads` and
`eigendecomposition_threads` / `--eig-threads` as explicit BLAS boundaries
while inheriting the caller or Job resource policy by default.
On the 2026-08-21 OpenBLAS/12-logical-core machine, controlled 8,000x800 tests
measured `gesvd` at 4.73 s with one BLAS thread versus 6.59 s with 12; `gesdd`
was fastest at two threads (1.61 s). These are machine-specific observations,
not portable defaults. Plan and benchmark on the target runtime before setting
a limit; use `--svd-threads 1` for the measured reference-compatible `gesvd`
configuration on this machine.

H-DMD now selects that compiled work through `linear_algebra_backend`. The
admitted `auto` route is `scipy-lapack`, whose provenance includes the linked
BLAS/LAPACK provider. An optional `torch-cuda` route owns transfer, full SVD,
projection, reduced eigendecomposition, and mode construction on-device; it is
explicit-only, performs a device-memory preflight, and refuses non-`gesvd`
drivers or CPU thread-limit controls. Use the paired full-spectrum benchmark
cases `scipy-hdmd-full-3000x256` and `torch-cuda-hdmd-full-3000x256` with a
relative speed gate before considering admission on a CUDA workstation.

SINDy publication has a separate, non-BLAS performance boundary. The standard
pipeline now stores the 15-term polynomial coefficients and training states,
then reconstructs requested eigenfunction columns lazily. On the full RPE
result this reduced publication from a 419,158,070-byte, 24.30-second archive to
459,218 bytes in 0.036 seconds without dropping an eigenvalue or coefficient.
Moving compression code to C++ would not address the underlying redundant-array
design; coefficient-first storage removes it at the archive contract. Standard
SINDy execution now also defers that matrix in memory, avoiding an additional
800,032,000-byte live allocation while keeping callable eigenfunctions and the
same KROM coefficients.
The versioned `sindy-lazy-polynomial-values-5000x256` case exercises the shared
archive/Viewer evaluator. Its 2026-08-21 OpenBLAS median was 0.00488 seconds for
a 5,000-by-256 output; compare only against reports with matching runtime
fingerprints.

An exact archive-free 20,000x2,000 `gesvd` run with `svd_threads=1` retained
all 2,000 modes and took 126.74 s: 67.30 s SVD, 4.87 s projected operator,
42.29 s eigendecomposition, and 12.28 s modes. Its eigenvalue set was within
`5.6e-10` of the prior full-spectrum record. Reusing the saved 2,000x2,000
operator isolated NumPy eig timings of 14.33, 14.64, 12.79, and 19.76 s for
1, 2, 4, and 12 threads respectively, with cross-thread set distances below
`3.3e-13`. An exact end-to-end candidate run with `--svd-threads 1
--eig-threads 4` completed the numerical kernel in 67.27 s (44.83 s SVD,
2.93 s projection, 10.81 s eigendecomposition, 8.70 s modes), retained all
2,000 modes, and matched the prior eigenvalue set within `5.54e-10`. This is
the measured workstation profile, but remains opt-in because BLAS vendor,
matrix shape, thermal state, and hardware change the optimum.

KROM autonomous rollout follows that rule explicitly.  Its vectorized NumPy
reference uses optimized matrix multiplication.  The currently available C++
comparison kernel is callable with `execution_backend="native"`, but it is not
auto-admitted because representative 64-mode reconstruction benchmarks show
the NumPy pathway is faster.  Both cases remain in `dynsyslab-benchmark` so a
future optimized implementation can earn admission with measured evidence.

## Process and thread resources

Job resource requests are also execution controls.  For local worker
processes, `resources.cores` is translated into OpenMP, OpenBLAS, MKL,
Accelerate, and NumExpr thread limits before Python or a numerical library is
loaded.  Serial in-process Jobs apply the equivalent runtime limit through
`threadpoolctl`.  This prevents a pool of worker processes from each creating
a full-machine BLAS thread pool.  The effective policy is recorded with Job
provenance and backend results.
