# Post-processing package

This package reads completed archives and writes derived results only. It is independent of `solver/` and `dynsyslib/`.

- `koopman/`: native DMD, Hankel DMD, EDMD, and SINDy equation discovery.
- `pca_pod/`: native PCA, snapshot POD, robust PCA, bagged PCA, direct SVD,
  and Hankel Alternative View of Koopman (HAVOK).
- `lcs/`: LCS adapters and an isolated source-snapshot location for LCStool.
- `ssm/`: SSM adapters and isolated source-snapshot locations for SSMtool and SSMLearn.

Fixed-lag Koopman methods (DMD, Hankel DMD, and EDMD) verify archive metadata
and the actual time values before running. Nonuniform input raises
`ResamplingDecisionRequired`, whose message asks whether the user wants an
in-memory resample. A caller may retry explicitly, for example with
`KoopmanResamplingOptions(time_step=0.01, method="pchip")`; available methods
are `pchip` (recommended), `linear`, and `nearest`. This never edits or appends
an interpolated solution trajectory to the raw archive. Only the requested
derived Koopman result is saved, with its in-memory resampling provenance.

The `trajectory-hdmd-sindy` preset uses radius (state index zero) and constructs
a shifted pair of 20,000-by-2,000 Hankel matrices: 20,000 window samples by
2,000 snapshot columns. A full pair consumes 22,000 analysis samples. It
recovers the archived run's integration coordinates in memory. For a
dimensionless RPE run this means `tau`, `r=R/R0`, and `dr/dtau`, even though
the immutable solver archive stores canonical dimensional time and state
values. SINDy evaluates the dimensionless RHS in the same coordinates. The
source trajectory is never rewritten. The preset applies
the established absolute SVD tolerance of `1e-10`, retains every singular
direction above that tolerance (no fixed modal truncation), uses the LAPACK
`gesvd` driver for reference-compatible full SVD behavior, and passes the exact
discrete-time eigenvalue objects into SINDy. SINDy converts each value with
`lambda=log(mu)/delta_t` for its continuous-time generator equation; KROM and
Viewer spectrum displays retain the dimensionless discrete-time `mu` values.

The lower-level `hankel_dmd` API also retains its general all-state
`matrix_layout="delay-state"` mode. The preset deliberately selects the scalar
`matrix_layout="window-snapshots"` contract instead. Derived metadata records
the requested and realized matrix geometry, scalar observable, rank policy, and
spectral-domain transfer.

The dense thread counts are also explicit. By default H-DMD inherits the caller's
or JobManager's BLAS resource policy. Direct calls may set `svd_threads` and
`eigendecomposition_threads`, and the CLI exposes `--svd-threads` and
`--eig-threads`; each requested limit is applied only around its corresponding
factorization and recorded in result metadata. This is preferable to a
machine-specific global default because OpenBLAS/MKL scaling varies strongly
with matrix shape and hardware.

The complete DMD factorization now crosses one typed dense-linear-algebra
boundary. `linear_algebra_backend="auto"` selects the admitted
`scipy-lapack` reference and records SciPy plus the linked BLAS provider (for
example OpenBLAS or oneMKL). `torch-cuda` is an optional explicit-only backend
that keeps host-to-device transfer, full `gesvd`, projection, reduced
eigendecomposition, and mode construction inside one CUDA kernel boundary. It
requires an available CUDA runtime/device, accepts only reference-compatible
`gesvd`, performs a conservative free-device-memory preflight, and never
silently falls back. It is not auto-admitted until the full-spectrum parity and
relative benchmark gates pass on the target machine. Inspect capabilities with
`dynsyslab-postprocess backends` and select explicitly with
`--linalg-backend scipy|torch-cuda`.

`hankel_dmd_ensemble` accepts a pooled `TimeSeries` and builds each member's
shifted pairs independently. It uses the same rank, absolute-tolerance, LAPACK
driver, memory-preflight, and timing kernel as single-trajectory DMD. Its
metadata records the exact paired-snapshot count and zero cross-trajectory
transitions; no synthetic pair is formed between the end of one trajectory and
the beginning of the next.

Native SINDy uses sequentially thresholded least squares with a selectable
observable library (polynomials by default). It derives time derivatives from
the supplied timestamp values, so it can operate on nonuniform trajectories;
its results are stored as independent derived archive entries.

Koopman-with-SINDy fits only one member of each numerically paired conjugate
eigenvalue pair when its dictionary and generator samples are real. The other
coefficient vector is the exact conjugate and therefore does not require a
second sequence of SVDs. Eigenfunction values for the complete spectrum are
then evaluated in one batched matrix multiplication. This optimization is on
by default, can be disabled with `reuse_conjugate_fits=False`, and records its
direct/reused fit counts in `fit_acceleration` metadata.

The trajectory preset also uses `null_space_backend="gram-eigh"`. It
precomputes the three small Gram terms shared by every generator eigenvalue,
then performs sequential sparsification with feature-sized Hermitian
eigenproblems instead of repeatedly factoring the full samples-by-features
matrix. Every subproblem is checked for a separated smallest eigendirection,
small Hermitian residual, and a nonnegative numerical floor. An ambiguous or
ill-conditioned subproblem automatically falls back to the exact SVD. The
public `koopman_with_sindy` function retains `null_space_backend="svd"` as its
reference default; backend, fallback, spectral-gap, and residual diagnostics
are persisted with the result.

Composable `processor` runs expose separate `algorithm_seconds`,
`archive_seconds`, and `total_seconds` values for every attached method. This
keeps dense numerical work distinguishable from compression and publication
costs when profiling large complete-spectrum KROM archives.

The processor itself owns only ordered workflow concerns: source adaptation,
coordinate policy, optional resampling, timing, lineage context, and attachment
metadata. Algorithm execution and derived-publication contracts live in closed,
family-owned method handlers under `method_handlers/`. This keeps Koopman and
decomposition code boundaries explicit and prevents the orchestration layer
from growing another algorithm-specific branch for every method.

`PostprocessingPipeline` is the conventional public workflow class; the
historical lowercase `processor` name remains a compatibility alias.
The standard settings are held by the immutable named
`TRAJECTORY_HDMD_SINDY_POLICY`; every preset receives fresh mutable attachment
copies, and descriptions, run summaries, and result attachment metadata record
the canonical policy name.
`workflow.py` owns the remaining cross-method concerns explicitly.
`SourceResolver` flattens archive paths, completed results, hierarchies, and
campaigns into provenance-bearing trajectories while rejecting cyclic source
graphs. `WorkflowContext` validates and removes workflow controls before method
dispatch, caches the common archive/integration-coordinate view for an ordered
run, and advances the latest spectrum plus its archive lineage atomically.
Resampling accepts only `ResamplingOptions`; an unrelated object is an error,
not a silently ignored request. These types are available through the lazy
top-level `postprocessing` façade without loading numerical method families.
The read-only planner uses this same preparation gate, so planned coordinate
and resampling geometry cannot drift from execution behavior.

Derived publication is immutable by default. Every call receives a unique
publication identifier and UTC timestamp while an analysis signature groups
runs having the same source and settings. `discover_derived` filters records by
analysis and source archive newest-first. Overwriting the deterministic legacy
path requires the explicit `publication="replace"` option. SINDy records the
exact H-DMD derived archive from which its eigenvalue objects came, and Viewer
uses that lineage automatically.

The named `trajectory-hdmd-sindy` policy uses coefficient-first SINDy storage.
It writes the fourth-order polynomial coefficient matrix, training states,
dictionary definition, and complete state-normalization transform. The logical
`eigenfunction_values` array is versioned as an exact virtual array and
`load_derived_array` reconstructs either all columns or only requested columns.
Viewer plots the complete spectrum and all KROM-retained spectrum points, but
reconstructs only its explicitly selected dashboard columns (five by default).
This is display selection, not model truncation. The policy does not truncate the
2,000-value spectrum. Generic/manual processor attachments and
`koopman_with_sindy_archive` retain full materialized storage by default for
compatibility; select coefficient-first explicitly there, or use
`--sindy-storage full` to override the standard CLI policy. Custom dictionaries
that cannot be serialized must use full storage.

Coefficient-first execution also sets `compute_training_values=False`. The
algorithm returns the complete coefficient matrix and callable eigenfunctions
without first allocating their full Cartesian product with the training rows.
For the standard RPE shape this avoids 800,032,000 bytes (about 763 MiB) of
complex values. `koopman_with_sindy` itself retains
`compute_training_values=True` for compatibility; compact archive wrappers
select deferred evaluation automatically. Full publication evaluates the same
serialized polynomial contract before writing. In-memory callers can request
all or selected training columns through
`KoopmanSINDyResult.evaluate_training_eigenfunctions(columns=...)`.

New derived archives use level-1 ZIP deflate by default. This is a publication
policy local to post-processing; it does not change raw solver archives.
`save_derived(..., compression="lzma")` remains available when minimum storage
size matters more than publication latency, and `recompress_derived` can
losslessly switch existing derived records among `deflate`, `lzma`, and
`stored`. The actual codec is written to `storage_compression` metadata.

For the 25,001-by-2,000 RPE eigenfunction matrix, the earlier materialized
archive was 419,158,070 bytes and required 24.30 seconds to publish. Replaying
the same coefficients, states, and metadata through coefficient-first storage
produced a 459,218-byte archive in 0.036 seconds (about 913 times smaller).
Reconstructing the 1,001 retained KROM columns took 0.74 seconds in a local
benchmark. These figures are workload- and machine-specific; `plan` reports
both estimated full and coefficient-first payload sizes before execution.

When an H-DMD result uses the same matrix as both its projected modes and its
trajectory eigenfunction samples, the archive writes that matrix once and
records `eigenfunction_values -> modes` in `storage.array_aliases`.
`load_derived_array` follows this alias transparently. This removes one of the
largest redundant arrays in a complete-spectrum H-DMD record while keeping the
logical archive contract explicit.

Each DMD-family extraction also returns a typed `KoopmanSpectrum`: a source
system reference, paired `KoopmanEigenvalue` objects, associated
`KoopmanEigenfunction` objects, and a separate list of modes. Current methods
derive discrete-time eigenvalues from their fitted one-step operator and only
then compute continuous-time values with the verified archive spacing
`delta_t`; provenance records which representation was originally derived.

`koopman_with_sindy` constructs sparse dictionary representations of
continuous-time Koopman-generator eigenfunctions from supplied eigenvalue
information. It accepts either one or more archived trajectories or a sampled
state-space domain. It never runs HDMD itself. Domain mode is the sole optional
post-processing dependency on `dynsyslib/`; see `DEPENDENCIES.md`.

The LCS and SSM reference suites are MATLAB-based. Their adapters validate an installed upstream snapshot but never import or merge its source into dynsyslab.

Regular SVD operates directly on a samples-by-features trajectory matrix and is
uncentered by default. HAVOK requires an individual, uniformly sampled
trajectory: it delay-embeds a specified scalar observable, performs an SVD of
the Hankel matrix, and fits the retained delay coordinates as a linear system
forced by the final retained coordinate. It never invokes the solver.
