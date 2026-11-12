# H-DMD to SINDy pipeline contract and performance

## Standard trajectory contract

`PostprocessingPipeline.preset("trajectory-hdmd-sindy")` is the standard
pipeline. The historical `processor` class name remains an alias. It is
also exposed by `dynsyslab-postprocess describe|plan|run`.
Its defaults come from the immutable named `TRAJECTORY_HDMD_SINDY_POLICY`;
each constructed pipeline receives independent setting copies and records the
canonical policy name in workflow and attachment provenance.

H-DMD uses one scalar observable (state index zero), a direct one-sample-shifted
pair of 20,000-by-2,000 Hankel matrices, an absolute singular-value cutoff of
`1e-10`, no fixed rank truncation, projected modes `U @ W`, and LAPACK `gesvd`.
The complete spectrum above the cutoff is retained. The source archive's fixed
sampling interval defines discrete eigenvalues `mu`. SINDy reuses those exact
eigenvalue objects and uses `lambda=log(mu)/delta_t` in the continuous generator
equation. Viewer displays `mu` as unitless.

Solver archives always publish canonical dimensional states. If their metadata
says the solve used dimensionless equations, this preset reversibly recovers
the integration coordinates in memory. For RPE those are `tau=t/T0`, `r=R/R0`,
and `u=(T0/R0)dR/dt`. H-DMD and the SINDy dimensionless RHS therefore operate in
the same coordinates as the generating run. No solution archive is modified.

SINDy uses a fourth-order polynomial dictionary (15 features for two states),
state normalization, the originating system RHS, guarded Gram-Hermitian
nullspace solves, conjugate-pair fit reuse, batched eigenfunction evaluation,
and an 8 GiB preflight budget. The exact-SVD nullspace implementation remains
available as a reference fallback.

Derived outputs are atomically and immutably published by default. A UUID,
creation timestamp, source-metadata hash, and analysis signature distinguish
and group records. SINDy records the exact upstream H-DMD archive. Viewer follows
that lineage automatically. Replacement requires an explicit opt-in.

New scattered-eigenfunction archives store sampled coordinates once under the
canonical `eigenfunction_values` array. Readers continue to accept the legacy
`training_eigenfunction_values` name, but writers no longer publish both names
with byte-identical payloads. Existing immutable archives are not rewritten.

The named standard policy publishes SINDy in coefficient-first form. It keeps
all 2,000 eigenvalues and all 2,000 polynomial coefficient vectors, but does not
duplicate their 25,001-sample evaluations in the ZIP payload. The complete
dictionary and state transform define an exact virtual array, and Viewer
reconstructs only the dashboard-selected columns lazily. Manual archive APIs preserve
their historical full-array default; `--sindy-storage full` is the explicit CLI
override and unserializable custom dictionaries require full storage.

Viewer keeps model retention and display selection separate. The spectrum panel
always contains all 2,000 H-DMD eigenvalues and marks all 1,001 nonredundant KROM
modes. The standard dashboard then builds phase, phase-versus-time,
magnitude-versus-time, real-versus-time, and imaginary-versus-time panels for
the first five modes in retained-KROM order. Callers may choose explicit spectrum
indices or request every retained dashboard, but neither choice mutates or
truncates the KROM. The Plotly payload stores repeated time and state coordinate
arrays once; the full t030 five-mode page fell from 26,613,804 to 13,408,468
bytes while preserving all model-spectrum points.

The same policy defers training-value evaluation during SINDy construction.
The full 2,000-by-15 coefficient matrix, state callables, gradients, residual
operators, and KROM are still constructed, but the 25,001-by-2,000 complex
sample matrix is not. This avoids an 800,032,000-byte (about 763 MiB) live
allocation. Direct `koopman_with_sindy` calls retain materialized training
values by default for compatibility.

## Comparison with the established MATLAB computation

The reference source builds, for one radius observable and one trajectory,
`X=hankel(c1,r1).'` and its one-sample-shifted `Y`, both 20,000-by-2,000. It then
runs `[U,S,V]=svd(X,'econ')`, retains `diag(S)>1e-10`, forms
`Atilde=((U'*Y)*V)*diag(1./diag(S))`, eigendecomposes `Atilde`, and uses projected
modes `U*W`. The corrected Python path performs the same matrix-level operations.
Singular/eigenvector signs, phases, and bases inside repeated subspaces are not
expected to be bit-identical across libraries.

That equivalence statement applies to the H-DMD matrix computation, not to the
implementation of every SINDy eigenfunction fit. The MATLAB trajectory path
uses fourth-order polynomials, forms `lambda*Theta-Gamma`, runs one economy SVD
for every supplied continuous generator eigenvalue, takes the final right
singular vector, and zeros coefficients below `1e-6` once. The standard Python
pipeline uses the same fourth-order polynomial function class and generator
equation, but normalizes state, applies sequential thresholding with a configured
threshold of `0.05` for at most 20 iterations, solves guarded feature-sized
Gram-Hermitian eigenproblems, reuses conjugate fits, and evaluates all functions
in a batch. These are deliberate robustness/performance changes, so sparse
coefficient vectors are not expected to be identical to MATLAB even though the
exact H-DMD eigenvalue objects are transferred unchanged into the SINDy KROM.

The earlier Python result was not equivalent. It stacked both RPE states into a
20,000-by-4,000 delay-state matrix and retained rank 2,789. It also consumed
canonical physical archive coordinates despite the dimensionless integration.
That run took about 22 minutes 57 seconds for H-DMD and 13 minutes 6 seconds for
SINDy (about 36 minutes total). It is preserved only as a historical comparison.

An exact, archive-free benchmark on the corrected RPE radius signal and this
machine produced:

| Stage | MATLAB R2024b / oneMKL | Python SciPy / OpenBLAS, `gesvd` |
|---|---:|---:|
| Matrix construction | 2.04 s | <0.01 s as a strided view; copy is charged to SVD |
| Economy SVD | 16.57 s | 52.48 s |
| Projected operator | 3.44 s | 3.29 s |
| Eigendecomposition | 9.01 s | 19.42 s |
| Projected modes | 12.12 s | 9.81 s |
| Total numerical pipeline | 43.18 s | 85.01 s (85.47 s wall) |
| Retained rank | 2,000 | 2,000 |

SciPy `gesdd` retained the same rank and reduced its SVD stage to 39.44 seconds,
but the measured end-to-end result was 81.61 seconds; one run is insufficient
to auto-admit it. It remains an explicit full-spectrum performance option.

The first fully published corrected pipeline run had additional thermal/runtime
variation and measured:

- H-DMD algorithm: 142.06 s. Inside it: SVD 87.20 s, projected operator
  4.29 s, eigendecomposition 36.52 s, projected modes 13.14 s, and about
  0.91 s of spectrum/KROM and orchestration overhead. Deflate publication was
  another 22.14 s.
- SINDy algorithm: 1.68 s. Inside it: coordinate/library construction 0.56 s,
  Gram and fit-plan preparation 0.13 s, 1,001 direct sparse nullspace fits
  0.43 s, batched evaluation 0.12 s, eigenfunction objects 0.01 s, and KROM
  construction 0.24 s. Deflate publication of the complete 25,001-by-2,000
  eigenfunction-value matrix was another 24.30 s.

That 24.30-second publication step was the dominant SINDy wall-time defect, not
the sparse fitting. The materialized record was 419,158,070 bytes. The new
standard coefficient-first payload for the same shapes measured 459,218 bytes
and 0.036 seconds to write, a roughly 913-fold size reduction. Lazy evaluation
of its 1,001 retained columns measured 0.74 seconds. Full storage remains an
opt-in compatibility representation.

This first corrected run predates the final named storage policy. Its H-DMD
settings already match the standard policy, and its SINDy numerical settings
already used Gram-Hermitian fits and conjugate reuse, but its attachment did not
set `archive_storage="coefficient-first"` or
`compute_training_values=False`. That is the significant difference between
the existing 419 MB record and a new run through the current standard preset;
it changes storage and peak memory, not eigenvalues or polynomial coefficients.

The original timed algorithm also spent 0.12 seconds producing and then
retaining that 763 MiB matrix. The new standard run records
`deferred-polynomial-evaluation` instead; this change is primarily a peak-memory
and object-size improvement rather than a large compute-time speedup.

## Phase-specific BLAS thread policy

The pipeline now exposes two independent, opt-in limits:

- `svd_threads` / CLI `--svd-threads` applies only to the dense SVD;
- `eigendecomposition_threads` / CLI `--eig-threads` applies only to the
  reduced operator eigendecomposition.

Both default to `None`, which inherits the caller or JobManager resource
policy. The requested values are recorded in H-DMD metadata. They are separate
because the two LAPACK workloads have different measured scaling on this
machine.

Using the exact immutable t030 source, full 20,000-by-2,000 `gesvd` geometry,
absolute `1e-10` cutoff, and no archive publication, the measured workstation
candidate `svd_threads=1`, `eigendecomposition_threads=4` produced:

| Stage | Seconds |
|---|---:|
| SVD | 44.83 |
| Projected operator | 2.93 |
| Eigendecomposition | 10.81 |
| Projected modes | 8.70 |
| Total numerical kernel | 67.27 |

All 2,000 modes were retained. The nearest-set distance to the prior published
eigenvalues was `5.54e-10`; isolated cross-thread eigensolver comparisons were
below `3.3e-13`. A one-thread-SVD run with inherited eig threads took 126.74 s,
illustrating both the value of phase-specific limits and the substantial
thermal/runtime variability. Therefore `(1, 4)` is a documented opt-in profile
for this workstation, not a portable mathematical default.

The corrected H-DMD and SINDy records resolve to exactly the same 2,000 stored
discrete eigenvalues (`max(abs(mu_hdmd-mu_sindy)) == 0`) and exactly the same
1,001 retained KROM indices. Their representations are intentionally different:
H-DMD stores projected delay-coordinate modes; SINDy stores sparse coefficient
vectors over 15 polynomial features and can evaluate them at arbitrary states.
The standard archive virtualizes, rather than materializes, trajectory values.

## Executable KROM capability boundaries

Spectrum identity does not imply representation identity. `KROM.capabilities()`
therefore reports lifted evolution, state reconstruction, analytic control
coupling, model coordinate system, and state representation independently.

- Scalar H-DMD is `delay-observable`: its rows are delayed radius samples, not
  the physical `(r,u)` state. It may evolve lifted coordinates, but the API now
  rejects attempts to slice its first rows and call them a physical-state
  reconstruction.
- Koopman-with-SINDy is `eigenfunction-only`: it supplies callable sparse
  eigenfunctions and analytic gradients, but no reconstruction modes. It does
  not advertise state reconstruction.
- An ordinary DMD/EDMD model whose modes have exactly the state dimension is
  `state-modes` and may reconstruct state.

For a dimensionless KROM, eigenfunctions and generator eigenvalues are defined
in model state and model time. Controller adapters accept physical solver state
and seconds, transform to model coordinates, scale `A` back to inverse seconds,
and use a dimensionless input Jacobian expressed per physical input unit. RPE's
input remains pressure in pascals. Koopman LQR/MPC return `system_input` values;
they never masquerade as additive state derivatives. Controller coupling
requires explicit analytic eigenfunction gradients by default; numerical
finite-difference gradients require an explicit `gradient_backend` opt-in.

## Acceleration policy

The expensive calls already execute in compiled Fortran/C through SciPy. Moving
the surrounding Python loops to C++ or Fortran would not materially accelerate
this workload. Backend selection is the useful boundary:

1. Prefer a Python environment whose NumPy/SciPy is linked to Intel oneMKL on
   this Intel workstation, then benchmark numerical parity and end-to-end time.
2. Keep SciPy/OpenBLAS as the portable reference fallback and control thread
   counts at the job boundary.
3. Add an optional CUDA backend only after hardware detection and parity tests.
   CuPy calls cuSOLVER; PyTorch exposes `gesvd`, `gesvdj`, and approximate
   `gesvda` CUDA drivers. The approximate driver is not admissible for the full
   reference spectrum.
4. Do not substitute `scipy.sparse.linalg.svds`, PROPACK, ARPACK, randomized
   SVD, or a normal-equations eigensolve while full-spectrum equivalence is
   required. They are appropriate only for an explicitly truncated/approximate
   policy.
5. Do not make a native backend automatic until a reproducible benchmark passes
   both numerical and performance gates.
6. For this workstation's reference-compatible `gesvd` path, start with
   `--svd-threads 1 --eig-threads 4`; re-benchmark after changing NumPy, SciPy,
   BLAS vendor, CPU, or cooling conditions.

The implementation choices above are backed by the current vendor interfaces:

- [SciPy `svd`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html)
  exposes both `gesvd` and the usually faster divide-and-conquer `gesdd`, and
  explicitly documents `gesvd` as the MATLAB/Octave-compatible approach.
- [MATLAB `svd`](https://www.mathworks.com/help/matlab/ref/double.svd.html)
  confirms that `"econ"` computes the reduced factors used here and supports
  threaded, GPU-array, and distributed-array execution.
- [Intel oneMKL `gesvd`](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2026-0/gesvd.html)
  and [`gesdd`](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/gesdd.html)
  provide optimized CPU LAPACK implementations. A separate MKL-linked
  NumPy/SciPy environment is the highest-priority CPU experiment because the
  current runtime is OpenBLAS, whereas the faster MATLAB measurement used
  oneMKL.
- [PyTorch CUDA SVD](https://docs.pytorch.org/docs/stable/generated/torch.linalg.svd.html)
  exposes accurate QR `gesvd`, Jacobi `gesvdj`, and approximate tall-matrix
  `gesvda`; only `gesvd` is currently admissible for reference-compatible full
  spectrum work. [CuPy SVD](https://docs.cupy.dev/en/stable/reference/generated/cupy.linalg.svd.html)
  similarly delegates large matrices to cuSOLVER.
- [MAGMA](https://icl.utk.edu/magma/) supplies LAPACK-compatible hybrid
  CPU/GPU dense factorizations, and
  [ScaLAPACK `PDGESVD`](https://www.netlib.org/scalapack/explore-html/d7/db3/pdsvdtst_8f_a77c8bbe8780697c212eb7c74d633206f.html)
  is the distributed-memory option. Both add deployment complexity and should
  be considered only after the oneMKL and single-GPU gates.

Writing a bespoke C++ or Fortran wrapper around the current call is not expected
to help: SciPy is already entering compiled LAPACK. C++/Fortran becomes useful
only as the host for a different admitted library or for fusing surrounding
operations after profiling proves transfer/allocation overhead material. The
full-spectrum requirement rules out ARPACK/PROPACK `svds`, randomized SVD, and
approximate `gesvda` as drop-in replacements.

This policy is executable through the H-DMD `linear_algebra_backend` boundary.
`auto` remains `scipy-lapack`; the backend provenance records the actual BLAS
library, so a SciPy environment linked to oneMKL is distinguishable from the
portable OpenBLAS environment without changing the mathematical API. The
optional `torch-cuda` candidate owns the entire factorization path rather than
moving only SVD to the GPU and immediately copying its largest factors back to
the CPU. It remains explicit-only and supports `gesvd` only. The versioned
benchmark contains comparable `scipy-hdmd-full-3000x256` and
`torch-cuda-hdmd-full-3000x256` cases when CUDA is available.

The generic `hankel_dmd`, `hankel_dmd_archive`, and GUI H-DMD actions remain
general-purpose entry points; they do not silently imply the complete scalar
trajectory H-DMD-to-SINDy contract. Use
`PostprocessingPipeline.preset("trajectory-hdmd-sindy")` or
`dynsyslab-postprocess` for the
standard combined pipeline. Its defaults are executable policy, covered by
contract tests and visible without computation through `describe` and `plan`.

Primary documentation:

- [MATLAB economy SVD](https://www.mathworks.com/help/matlab/ref/double.svd.html)
- [SciPy dense SVD and LAPACK drivers](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html)
- [Intel oneMKL `gesdd`](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-0/gesdd.html)
- [SciPy partial `svds`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html)
- [PyTorch CUDA SVD drivers](https://docs.pytorch.org/docs/stable/generated/torch.linalg.svd.html)
- [CuPy/cuSOLVER SVD](https://docs.cupy.dev/en/stable/reference/generated/cupy.linalg.svd.html)
