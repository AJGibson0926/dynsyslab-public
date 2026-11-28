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
Every fit is a krom (`KoopmanModel`, one `krom` data object); its compressed operator is in
[KoopmanCompression.hpp](../../native/include/dynsyslab/core/KoopmanCompression.hpp).
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

## What a spectral fit contains

A spectral fit is a krom like any other method's: one `KoopmanModel`, one `krom` data object.
What is particular to it is its **compressed operator**, `KoopmanModel::compression()`
(`KoopmanCompression`), which it keeps because some analyses need the operator itself rather than
its eigenpairs (it stays valid even when a diagonal eigenmode representation is unreliable):

- **The atomic part, on the spectrum.** The eigenpairs: values with stable mode ids (numbered in the
  operator's eigensystem order), each pair's one-step residual, and the eigenfunctions as an
  expansion on the frozen dictionary (and an interpolation over physical states).
- **The continuous part and the observables, on the spectrum.** Rigged DMD's spectral bands (a named
  angle grid with its kernel) and, per observable, its Koopman modes and its measure on each
  band (density, moments, mass) with its wave packets as an expansion on the same dictionary
  (`SpectralMeasure::packets`); atoms follow from the Koopman modes (`Spectrum::atoms`).
- **Evidence, on the spectrum.** The pseudospectrum (grid, minimum residual) with its minimizers
  (`Pseudospectrum::minimizers`).
- **The operator, as `compression()`.** Frozen dictionary and normalization (linear, polynomial,
  Chebyshev, Laguerre, Hermite, Fourier-Hermite, radial sections or nearest-center indicator
  partitions); basis B with columns orthonormal in the training quadrature metric; the coefficient
  action A; named readouts of dictionary features; sampling interval, coordinate contract, method,
  assumptions. Its identity (`koopmanCompressionIdentity`) is what references to the krom's
  operator-derived observables name. Child kroms and extraction keep it.
- **Diagnostics**: invariance angles, stochastic residuals and RKHS-adjoint results.

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

A `koopman_spectral` step publishes the fit as a krom (see
[EIGENFUNCTION_OBJECTS.md](EIGENFUNCTION_OBJECTS.md#eigenpairs-spectrum-and-pseudospectrum)). Mode
extraction keeps the selected pairs (with their ids and residuals), the operator, its pseudospectrum,
measures and spectral functions, and clears the diagnostics because they refer to the original trial
subspace. The archive stores the operator as `compression_basis`, `compression_action` and
`compression_readout__<name>` with its description in the krom block, and the spectral functions as
`spectral_function__<key>`. There is no second archive kind: the query and periodic routes write
kroms too. A spectral fit's krom has no input maps until a stage fits them.

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

Output archives are created atomically and never overwrite an existing path. A fit is a krom:
`output_archive` writes it, and actions on an existing fit name it as `krom` (or `archive`); inline
models are no longer accepted. Results report the spectrum (`eigenvalues`, `mode_ids`, `residuals`
with the pseudospectrum, `measures` each with its `packets`) and, unless `include_model:false`, the
operator as `compression`. Actions:
`fit`, `inspect`, `evaluate`, `reconstruct`, `rollout`, `residuals`,
`rigged_measure`, `correlation_measure`, `specrkhs`, `pagerank`, `observable`.
Complex arrays use pairs [real, imaginary]. An `inspect` query returns the exact coordinate contract.

For repeated stochastic data supply `replicate_count` and
`successor_replicates` in replicate by sample by coordinate order in `data`.
At least two independent successors at each identical source point are required.
The empirical mean residual retains finite-replicate sampling error; the reported
innovation variance is unbiased under the asserted independent-replicate model.

### Rigged DMD: the `rigged_measure` stage

Rigged DMD (Colbrook, Drysdale & Horning, arXiv:2405.00782, 2024) is a stage of its own that acts
on a krom. The krom's spectral fit must be a unitary compression, which `mpedmd` with
`assume_measure_preserving:true` gives (Colbrook, SIAM J. Numer. Anal. 61, 2023). The stage checks
the operator (A^T A = I), not the method's name. A child krom finds the fit through its parents, so
measures can be added in several steps.

~~~json
[{"name": "mpedmd", "settings": {"assume_measure_preserving": true, "rank": 2,
   "dictionary": {"kind": "linear", "constant": false}}},
 {"name": "rigged_measure", "settings": {"measures": [
   {"name": "observable", "coefficients": [1, 0],
    "angles": [-3.14, -1.57, 0, 1.57, 3.14], "bandwidth": 0.08, "order": 4}]}}]
~~~

`krom` names the krom to analyse (default: the stage's source); it must carry a unitary compression,
inherited from its parents when it is a child. Each item names a spectral band (`band`,
default `name`: the angle grid with the kernel's bandwidth and order; an existing band must match)
and an observable (`observable`, default `name`, given by dictionary `coefficients`; or output
`output` of a fitted feature readout, `reconstruction: {name, output}`, whose name it then takes).
The result is a child krom whose spectrum gains the band, the observable's measure on it (the
smoothed density on the grid, its mass) and, when the observable has none, its Koopman modes, its
exact coefficients on the eigenfunctions. A unitary compression's measure of an observable is purely
atomic, one atom per eigenpair with mass |c_k|^2, the observable's Koopman mode on that unit-norm
eigenfunction, and the masses sum to its squared norm (Colbrook 2023); `Spectrum::atoms` derives
them from the Koopman modes by stable mode id, so they resolve in any child or extraction that keeps
those pairs. Rigged DMD's density is this atomic measure smoothed by its kernel. The measure's
generalized eigenfunctions, one wave packet per angle, are an expansion on the frozen dictionary held
by the measure (`SpectralMeasure::packets`); they survive `krom_extract` and evaluate as an
observable of kind `wave-packet` (`component`: the measure's key `band/observable/output`, or
`index`: its position). The stage's metadata lists the new measures by key. A fitting stage refuses
`rigged_measures` and points to this stage. The `koopman-spectral` query keeps its `rigged_measures` setting and its
`rigged_measure` action.

This example assumes a two-feature dictionary; supply one coefficient per frozen feature.
Coefficients describe the observable in the frozen dictionary, not in eigenmode
order. The observable must lie in the retained trial subspace; projection is not
silently performed. Higher-order rational kernels can have negative lobes.
`correlation_measure` requires `assume_stationary:true`. `specrkhs` returns typed
associated data without constructing a forward KROM. `pagerank` requires a
nonnegative transition graph; a signed EDMD matrix is not automatically such a graph.

### Pooled fits over a family of orbits

`mpedmd` (and the other `koopman_spectral` fits) take `"additional_sources": [archives...]` of the same
system and clock, or `"pool": true`, which makes a one-stage job's `sources` one fit (the first leads,
the rest are its additional sources; `"sources": ["${solve.archives}"]` in a job group pools every
orbit a family solve wrote): every record's delay pairs enter one fit, none crossing a record boundary. Settings
`"source_weights"` weight the records: `"uniform"`, a list, or
`{"canonical_area": {"invariant": "energy", "band": "<band>"}}`, the trapezoid rule in the invariant
divided by each orbit's frequency, dE/Ω, which is the invariant measure's transverse density across a
family of periodic orbits in action-angle coordinates (dI = dE/Ω; Arnold, *Mathematical Methods of
Classical Mechanics*, §50). Each record's weight is split evenly over its pairs, so the time averages
along each orbit discretize the measure along it. The krom's origin lists every record.

An observable given to Rigged DMD (or any measure) must lie in the fit's trial subspace. Its Koopman
modes are solved in two steps, z = basis⁺g on the whitened trial basis's orthogonal columns (the span
test) and then V c = z on the action's eigenvectors, so the whitening's condition number
σ_max/σ_min does not enter the solve. `residual_grid` also takes a lattice
`{"real": [from, to, n], "imag": [from, to, n]}`, which the Viewer draws as the pseudospectrum.

### A system's known spectrum against a fit: `band_projection` and `packet_validation`

When the system catalogue gives a band with generalized eigenfunctions (a shell entry: exact phase
φ and frequency Ω(H)), orbit records measure an observable's spectral measure on that band exactly up
to quadrature. `band_projection` pools its sources, one orbit each: along each record it takes the
weighted Birkhoff averages a_n = ⟨g e^{−inφ}⟩ with the smooth window exp(−1/(s(1−s))), which converge
faster than any power of the record length for a smooth quasi-periodic signal (Das & Yorke,
*Nonlinearity* 31, 2018). These are g's Fourier coefficients in the band's exact phase. Weighted
across orbits by `weights` (as above), they give the measure's atoms at nΩ_j, which the stage smooths
with Rigged DMD's kernel onto `angles` (`reference_density__<obs>`, the density a Rigged DMD fit is
compared with) and turns into wave-packet amplitudes at `packet_centers`. With `"krom"` it also
publishes a child of that krom whose spectrum gains the reference band (default name
`fourier-reference`) with the measures of the named observables (`measures: [{observable, as,
output}]`), beside the fit's own bands: one spectrum holds the exact and the fitted measures of one
observable.

`packet_validation` evaluates a fitted krom's Rigged DMD wave packets along held-out orbit records (on
the fit's own delay windows, its Hankel provenance checked against the settings) and compares their
RMS amplitude on each orbit with the band projection's reference packets
(`packet_rms__<obs>__<k>`, `reference_packet_rms__<obs>__<k>`, `validation_frequency`).
With `"krom"`, the report's metadata `density_errors` lists, for every density the krom already has of
the compared output on the same angle grid (Rigged DMD's), its relative L2 error against the reference
and both masses; `packet_validation`'s metadata `errors` gives the packets' relative L2 amplitude error
over the held-out orbits. The data pack `data-packs/known-spectrum-comparisons-2026-11-28` runs this
comparison on every catalogue system with a shell band.
`workflow/native/tests/band_projection_tests.cpp` runs the whole comparison on the pendulum as one job
group: solve, project, pooled mpEDMD, Rigged DMD, both bands on one spectrum, held-out packets.

## Experiments

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

`dynsyslab_colbrook_reference_tests` (`workflow/native/tests/colbrook_reference_tests.cpp`)
holds upstream parity natively. The pinned upstream `mpEDMDqr` and `riggedDMD` were
run once on the synthetic rotation (128 samples, rotation 0.31, bandwidth 0.08,
257 angles). Their operator and their order-1 and order-4 densities and wave packets,
taken at every 16th angle, are frozen in
[`colbrook_rotation_reference_v1.json`](../../workflow/native/tests/data/colbrook_rotation_reference_v1.json)
together with the upstream commits. The native results match within 2e-10. The order-1
density and packets also equal the closed-form Poisson-smoothed two-atom answer.
The same test checks two invariances of the weighted fit (Colbrook 2023, Sec. 3).
Block-QR compression of the weighted joint rows (Demmel et al. 2012) leaves the
operator, density and packets unchanged. Whitening the dictionary changes only the
coefficients. These checks replace the retired MATLAB reference adapters.

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
reference acquisition ([reference acquisition](../provenance/colbrook/README.md)),
[third-party notices](../provenance/COLBROOK_NATIVE_ADAPTATION.md).
