# Performance and native execution

All numerical execution is native C++: solvers, post-processing, KROM rollout,
control and estimation run in the compiled core behind the workflow service and
`dynsyslab-cli`. There is no Python execution backend to select; the optional
Python client only forwards jobs to the compiled CLI. Every result records its
solver, method and runtime fingerprint in its archive metadata, and
`dynsyslab-cli capabilities` and `dynsyslab-cli doctor` report what the
installation provides.

## Benchmarks

Benchmarks are native ([job queue](docs/architecture/JOB_QUEUE.md)). Run the reproducible baseline with:

```text
dynsyslab-cli benchmark --iterations 5 --warmup 1 --output benchmark.json
dynsyslab-cli benchmark cases
dynsyslab-cli benchmark --case ode_rk45_duffing --iterations 5 --output rk45.json
```

The report (`dynsyslab.benchmark-report/v1`) keeps every sample, the min/median/max and the
runtime fingerprint (version, commit, compiler, build type, machine). Compare a candidate to a
baseline with explicit thresholds:

```text
dynsyslab-cli benchmark gate baseline.json candidate.json --required-speedup 1.0
dynsyslab-cli benchmark gate baseline.json candidate.json --require ode_bdf_duffing=2.0
```

The gate exits with status 4 when a requirement fails and refuses comparisons across different
machine, platform, compiler or build-type identities unless `--any-runtime` is given. Numerical
parity is tested separately; timing never substitutes for correctness. Compare only reports with
matching runtime fingerprints.

## Coordinate contracts

The RPE model is explicitly dimensionless by default. Its native model record
contains the RPE groups, pressure scale, natural-period time scale, and the
physical sinusoidal-driver parameters. Native evaluation converts model time
`tau` back to physical seconds only for driver phase and preserves driver
amplitude in Pa; the dimensional form is a separate, recorded choice.

A dimensionless KROM evolves on model time even when its controller-facing
adapter accepts physical seconds; rollout kernels receive already-normalized
model time and model state. Kernels must not infer units from array magnitudes.

## Dense linear algebra and memory

Dense H-DMD, DMD and SINDy stages run on the compiled linear-algebra layer
(`dynsyslab_linalg`). The SVD driver is an explicit stage setting
(`svd_driver`). Before allocation, a stage refuses a job whose dense
working-set estimate exceeds the execution memory limit; `dynsyslab-cli plan`
reports that estimate.

SINDy publication is coefficient-first: the standard pipeline stores the
polynomial coefficients, dictionary, state transform and training states, and
reconstructs requested eigenfunction columns lazily. This is a storage policy,
not a truncation: every eigenvalue and coefficient is kept, and the large
sample matrix is neither archived nor held in memory until a column is
requested.

## Process and thread resources

Job resource requests are execution controls. `execution.max_threads` and
`execution.max_memory_mib` in a job, or `--max-threads N` and
`--max-memory-mib M` on `run`, `queue run-next` and `queue serve`, set the
execution limits; the smaller value wins. The thread limit caps native worker
pools and sets `OMP_NUM_THREADS` and the BLAS thread variables, so a pool of
queue workers does not each create a full-machine thread pool. `queue serve
--pool-threads T --pool-memory-mib M` bounds what runs beside already running
jobs. See [the job queue](docs/architecture/JOB_QUEUE.md).

## Historical measurements

These figures were measured on the author's Windows workstation while the
former Python reference implementations still existed. They are
machine-specific observations, not portable defaults; re-measure with
`dynsyslab-cli benchmark` on the target runtime.

- 2026-08-21: native 0.2.0 completed the RPE solve kernel in 0.00411 s versus
  2.918 s for the Python reference (709.7x).
- 2026-08-21, OpenBLAS on 12 logical cores, 8,000x800 SVD: `gesvd` took 4.73 s
  with one BLAS thread versus 6.59 s with 12; `gesdd` was fastest at two
  threads (1.61 s).
- Coefficient-first SINDy storage reduced the full RPE publication from a
  419,158,070-byte, 24.30-second archive to 459,218 bytes in 0.036 seconds, and
  deferring the eigenfunction matrix avoided an 800,032,000-byte live
  allocation.
- An exact archive-free 20,000x2,000 `gesvd` H-DMD run retained all 2,000
  modes in 126.74 s (67.30 s SVD, 4.87 s projected operator, 42.29 s
  eigendecomposition, 12.28 s modes); its eigenvalue set was within `5.6e-10`
  of the prior full-spectrum record.
