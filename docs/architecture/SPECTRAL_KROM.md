# Spectral KROMs and native Koopman analysis

## Ownership and entry points

Implemented November 2026. The existing observable/Hankel sampler is unchanged in
its mathematical convention. The new estimators consume its X/Y pairs or explicit
snapshot pairs. They never manufacture adjacent pairs across trajectory boundaries.

| Owner | Responsibility |
|---|---|
| `dynsyslab_foundation` | Observable definitions and passive decomposition references; no KROM loading |
| `dynsyslab_koopman` | Frozen dictionaries, spectral result classes, evaluation, archive codecs, reference binding |
| `dynsyslab_spectral` | Fitting, dictionary construction, residuals, invariance diagnostics, measures, wave packets |
| `dynsyslab_numerics` | SVD, eigensystems and least squares; scientific modules call its compiled API |
| Workflow / Post-Processor | Resolve source observables and clocks, prepare/reuse Hankel pairs, call native estimators, publish immutable results |
| Workbench | Configure methods; owns no scientific algorithm |
| Python client | Marshal JSON through the compiled C ABI; no second fitting implementation |

C++ entry points are in
[KoopmanSpectralAnalysis.hpp](../../native/include/dynsyslab/core/KoopmanSpectralAnalysis.hpp).
Typed results are in [SpectralKROM.hpp](../../native/include/dynsyslab/core/SpectralKROM.hpp).
The C entry point is `dsl_workflow_koopman_spectral_json`; it follows the existing
opaque-result lifecycle and error-status protocol. Python exposes
`NativeWorkflow.koopman_spectral(request)`. The compiled CLI accepts
`dynsyslab-cli koopman-spectral query.json`. The optional Python CLI delegates the
same command. The Post-Processor has a **Spectral Koopman** checkbox and configuration
dialog with estimator, dictionary, rank, normalization and explicit structural assumptions.

## Supported calculations

| Method | Actual implemented calculation |
|---|---|
| EDMD / ResDMD | Weighted, SVD-whitened least-squares compression; explicit eigenfunction coefficients and residuals |
| Kernel-section EDMD / ResDMD | Same estimator on frozen Gaussian, Laplacian, Lorentzian or Matern sections; explicit centers and out-of-sample evaluation |
| Weighted EDMD / ResDMD | Explicit nonnegative quadrature weights or normalized smooth Birkhoff weights |
| mpEDMD | Orthogonal Procrustes in the weighted orthonormal trial basis; unitary finite compression |
| Hermitian DMD | Symmetric compression in that basis, for an asserted self-adjoint problem |
| TLS DMD | Joint snapshot-subspace projection before fitting; an explicit joint rank is required |
| MultDMD | Deterministic successor-cell selection minimizing the Gram-normalized objective on a fully occupied indicator partition |
| Stochastic ResDMD | Pooled regression over repeated successors; sample-mean residual, unbiased innovation standard deviation and total pathwise residual |
| Residual pseudospectrum | Minimum residual and minimizing dictionary coefficients at each requested complex point, over the retained trial subspace |
| Principal-angle diagnostics | Angles between source and successor evaluation subspaces, their ranks, and relative projection defect |
| Correlation spectral measure | Biased autocorrelation followed by Fejer smoothing on a requested angular grid |
| Rigged DMD calculation | First- through eighth-order equispaced rational resolvent smoothing of a finite mpEDMD compression; densities and dictionary wave packets |
| SpecRKHS adjoint calculation | Gaussian kernel Gram matrices, rank truncation, adjoint eigenvectors and kernel-norm residual/pseudoresidual calculations |
| Personalized PageRank | Nonnegative graph normalization, dangling-row handling, personalization and converged iteration |

The original exact DMD, Hankel DMD, DMDc and EDMDc remain available. Existing
`edmd` Post-Processor jobs retain their old implementation. Use
`koopman_spectral` with `settings.method="edmd"` for the new common spectral contract.
Radial sections support a kernel dictionary but are **not labeled as the full
dual kernel ResDMD algorithm**. Principal-angle diagnostics are not the entire
TrustKoopman PAD selection/error-bound framework. PageRank is a reusable graph
calculation, not yet a complete automatic EDMD-PPR dictionary-search workflow.

Not ported in this change: the full dual kernel ResDMD construction, kernelized
Rigged DMD extension, neural DeepMDMD training, randomized/fbDMD variants beyond
existing capabilities, certified infinite-dimensional error bounds, complete PAD
selection, and application-specific ENSO workflows. These require distinct
contracts or further validation. The published DeepMDMD repository URL was
unavailable; no substitute code was invented.

## What a spectral KROM contains

A `SpectralKROM` is a second KROM representation, separate from the existing
polynomial/scattered eigenfunction `KoopmanModel`. It retains a finite operator
even when a diagonal eigenmode representation is unreliable:

- Frozen dictionary and normalization; linear, polynomial, Chebyshev, Laguerre,
  Hermite, Fourier-Hermite, radial sections or nearest-center indicator partitions.
- Basis B with columns orthonormal in the training quadrature metric, coefficient
  action A, discrete multipliers, and explicit eigenfunction coefficients.
- Named dictionary readouts for measured coordinates or other supplied outputs.
- Typed spectral measures, residual landscapes and minimizing coefficients,
  invariance diagnostics, wave packets, stochastic diagnostics and RKHS-adjoint
  results. Measures store angular densities, moments, optional atoms and mass.
- Sampling interval, coordinate-provider contract, method, assumptions and
  content identity.

Automatic basis construction and deterministic farthest-point or k-means center selection
use training points only. Training normalization is frozen; evaluation never
relearns it from held-out data. Explicit nodes/exponents/normalization are also accepted.
Dictionary size defaults to a maximum of 4096 features. Dense SpecRKHS analysis
rejects more than 4096 sections; reduce the source/Hankel sample count explicitly.

With row evaluations P, reduced features are Q=P B. EDMD fits Q A approximately
equal to the successor reduced features. Eigenfunctions are P B v when A v=lambda v.
A free readout rollout starts from q0=P(x0)B and applies q(k+1)=q(k)A, using no later
measurements. This also works for defective operators. No logarithm is taken
implicitly; frequencies in spectral measures are radians per sample, and physical
angular frequencies require division by the archived sample interval.

For an ensemble, `action="rollout_batch"` accepts flat sample-major
`initial_points`, `steps`, and the same named `reconstruction`. It returns flat
`trajectory/time/output` values with `trajectories`, `samples`, `outputs`, and
`sample_interval`. The compiled core factors the readout once for the ensemble;
it neither refits nor uses later observations. The query bounds output at ten
million scalars. The existing single-point operation is unchanged. These APIs
advance only by the archived sample interval. For a centered history dictionary,
causal use requires shifting the supplied history to end at acquisition time
and choosing the corresponding readout; the API does not infer this convention.

A spectral measure is not an eigenmode list. Wave packets are finite-width
functions, not certified distributional limits. RKHS adjoint sections are stored
as **adjoint** results and are not accepted as an autonomous forward model.

Existing `KoopmanModel` objects also accept `SpectralData` attachments. Their archives
persist `spectral_data_json`. Mode extraction clears these diagnostics because
they refer to the original trial subspace; it does not silently relabel them as
diagnostics of the extracted subspace. New spectral KROM archives use
`dynsyslab.spectral-krom/v1`, a compact JSON contract, and separate NPZ numeric
arrays. Identity is checked on loading. Spectral KROMs do not automatically become
control-ready polynomial KROMs: no continuous input map or physical-state provider
is inferred.

## Observable references to an existing KROM

`ObservableDecompositionReference` contains artifact locator, model identity,
component name, kind and index. `Observable::withDecomposition` associates one
with an existing observable. Expression serialization preserves nonempty
associations; differentiating an observable clears the old association because a
different function is produced.

Kinds include reconstruction, eigenfunction, spectral-measure and wave-packet.
Measures can be associated with an observable but cannot be evaluated as a
pointwise time signal. The Koopman layer's `existingKromObservable` resolves
references; foundation never imports Koopman or reads model files.
For a spectral KROM, the identity includes the dictionary, operator, diagnostics
and readouts. Legacy KROM references use the existing model identity and a named
reconstruction in the immutable referenced archive.

A Hankel channel can use a resolved reference:

~~~json
{
  "observable": {
    "artifact": "path/to/spectral-krom.npz",
    "model_identity": "identity returned by inspect",
    "component": "coordinates",
    "kind": "reconstruction",
    "output": 0,
    "part": "real",
    "coordinate_contract": "exact contract returned by inspect"
  },
  "taps": [0, 2]
}
~~~

This evaluates the supplied source rows as the KROM's input coordinates.
It requires an explicit matching coordinate contract; equal vector dimensions
alone do not establish compatibility. It does not transform physical state into
an unknown delay/history coordinate system. Complex channels require choosing
real or imaginary part. No time derivative of a referenced model is inferred.

## API examples

A direct data fit (flat arrays are samples by coordinate, row-major):

~~~python
from dynsyslab_client.workflow import NativeWorkflow
api = NativeWorkflow()
fit = api.koopman_spectral({
    "action": "fit",
    "data": {"dimension": 2, "x": X.ravel().tolist(), "y": Y.ravel().tolist()},
    "settings": {
        "method": "resdmd", "rank": 6, "sample_interval": dt,
        "dictionary": {"kind": "polynomial", "degree": 3, "normalize": True},
        "residual_grid": [[0.9, 0.0], [1.0, 0.1]]
    },
    "readouts": {"echo": {"values": echo.tolist(), "outputs": 1}},
    "output_archive": "new-model.npz"
})
prediction = api.koopman_spectral({
    "action": "rollout", "archive": "new-model.npz",
    "reconstruction": "echo", "initial": X[0].tolist(), "steps": 100
})
~~~

Output archives are created atomically and never overwrite an existing path.
Queries accept either `archive` or the returned `model` object. Actions:
`fit`, `inspect`, `evaluate`, `reconstruct`, `rollout`, `residuals`,
`rigged_measure`, `correlation_measure`, `specrkhs`, `pagerank`, `observable`.
Complex arrays use pairs [real, imaginary]. `include_model:false` suppresses the
full fitted model in responses. An `inspect` query returns the exact coordinate contract.

For repeated stochastic data supply `replicate_count` and
`successor_replicates` in replicate by sample by coordinate order in `data`.
At least two independent successors at each identical source point are required.
The empirical mean residual retains finite-replicate sampling error; the reported
innovation variance is unbiased under the asserted independent-replicate model.

For Rigged DMD, choose `method:"mpedmd"`, explicitly set
`assume_measure_preserving:true`, and add:

~~~json
"rigged_measures": [
  {"name": "observable", "coefficients": [1, 0],
   "angles": [-3.14, -1.57, 0, 1.57, 3.14], "bandwidth": 0.08, "order": 4}
]
~~~

This example assumes a two-feature dictionary; supply one coefficient per frozen feature.
Coefficients describe the observable in the frozen dictionary, not in eigenmode
order. The observable must lie in the retained trial subspace; projection is not
silently performed. Higher-order rational kernels can have negative lobes.
`correlation_measure` requires `assume_stationary:true`. `specrkhs` returns typed
associated data without constructing a forward KROM. `pagerank` requires a
nonnegative transition graph; a signed EDMD matrix is not automatically such a graph.

## Pipeline jobs

The [observable and dictionary guide](OSCILLATOR_OBSERVABLES_DICTIONARIES.md)
documents the Pendulum/Duffing scalar observables, reusable Chebyshev, Laguerre,
Hermite, Fourier-Hermite and radial dictionaries, normalization, source formulas,
and the standalone `action:"dictionary"` API. These families are shared across
systems, including Van der Pol. The SpecRKHS adjoint method remains Gaussian-only.

Stages `resdmd`, `mpedmd`, `hermitian_dmd`, `tlsdmd`, `multdmd`,
`stochastic_resdmd`, `specrkhs` or `koopman_spectral` accept the shared
[observable/Hankel contract](OBSERVABLE_HANKEL.md). With no observable selection,
all source coordinates are selected, and a one-block window is used.
`use_previous_hankel:true` consumes the preceding immutable matrix.
Its paired shift determines the sample interval; conflicting intervals are rejected.
The [rotation example](../../examples/native/spectral-koopman-query.json) is a
standalone CLI fit; the [pipeline template](../../examples/native/spectral-koopman-job.json)
shows a shared Hankel comparison and requires a real source path.

Readout fitting is available directly through the C++ API and query `readouts`.
The pipeline attaches a readout for the input coordinates by default.
Each source is fitted independently; use explicit pooled X/Y pairs when fitting
several trajectories jointly.

## Interpretation and validation

Residuals are finite-data evidence, not universal certificates. Training
interpolation can make an empirical residual misleadingly small; use independent
snapshot pairs through `action:"residuals"`, examine conditioning/rank, and refine
sampling/dictionaries. Small residual does not equal spectral distance for a
nonnormal operator. Principal angles test a finite linear trial space, not
existence of some nonlinear closed coordinate law.

mpEDMD and spectral measures require the stated invariant-measure assumptions.
A damped RPE transient should not receive unitary structure merely because it is
an oscillation. Hermitian DMD requires a different, self-adjoint problem.
SpecRKHS changes the function space and norm; its spectrum need not be the L2
spectrum of an empirical trajectory.

Validation commands:

~~~powershell
cmake --build build/native-suite-fast --target dynsyslab_test_binaries --parallel 2
ctest --test-dir build/native-suite-fast --output-on-failure
python -m unittest discover -s tests -p test_native_koopman_spectral.py -v
~~~

The new native tests cover rotation, damped spectra, held-out mismatch, free
rollouts, residual landscapes, unitary and self-adjoint constraints, MultDMD
cycles, stochastic mean/innovation separation, RKHS identity, spectral mass,
archive identities and reference binding. Python tests exercise the C ABI, CLI,
immutable archive readback and shared-Hankel workflow, and independently compare
weighted EDMD to normal equations.

The optional [MATLAB parity script](../../tools/experiments/miscellaneous_computes/references/colbrook/check_native_reference.m)
executes the pinned upstream mpEDMD/Rigged routines on a synthetic rotation.
It uses explicit observable coefficients and requests the spectral-density grid.
Native and upstream operators, order-1/order-4 densities and wave packets agree
within the tests' 2e-10 absolute tolerance. The parity artifact is local build output;
the test skips that comparison if the reference script has not been run.

Validation run on 2026-11-11: all 121 native CTest targets passed; all four
Python/C-ABI/CLI/workflow tests passed against both the build and the isolated
package; the upstream MATLAB parity test passed; the isolated package/SDK smoke
gate passed; the strict architecture audit reported zero boundary violations,
dependency cycles, broken links or unindexed documents. The installed application
was not replaced; the rebuilt application is in `build/native-suite-fast/bin`.

Rigged spectral grids reuse one complex Schur factorization of the frozen
compression for all angles and kernel poles. The full upper-triangular factor
is retained, rather than replacing it by eigenvalues; shifted triangular solves
avoid repeatedly factoring dense matrices. This handles repeated eigenvalues
without assuming numerically returned eigenvectors form a unitary matrix.
The scalar kernel, finite compression, observable projection, and archive
contracts are unchanged. The numerics API exposes `complexSchur` as row-major
unitary and triangular factors. See the
[Eigen factorization contract](https://libeigen.gitlab.io/eigen/docs-3.4/classEigen_1_1ComplexSchur.html)
and the Rigged DMD source cited above.

Source inventory and attribution:
[reference acquisition](../../tools/experiments/miscellaneous_computes/references/colbrook/README.md),
[third-party notices](../provenance/COLBROOK_NATIVE_ADAPTATION.md).
