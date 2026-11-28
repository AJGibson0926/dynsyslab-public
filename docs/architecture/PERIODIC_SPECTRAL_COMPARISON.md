# Periodic observable spectral comparisons

Repeatable mpEDMD/Rigged-DMD comparisons for supplied **physical-time Fourier
orbit banks**. The workflow does not infer integrability, invariant transverse
weights, or infinite Fourier tails. Case-study inputs and numerical results
belong in caller-selected output/research folders, not these modules.

## Ownership and APIs

| Responsibility | Module / API |
|---|---|
| Delay translation, Parseval, phase-quadrature audit | `dynsyslab_math`, `PeriodicSnapshots.hpp` |
| Matched fitting and frozen packet comparisons | `dynsyslab_spectral`, `PeriodicSpectralComparison.hpp` |
| Reused Schur factorization | `dynsyslab_numerics`, `Eigensystem.hpp` |
| Native orchestration and immutable archives | `workflow/native/src/PeriodicSpectralQuery.cpp` |
| Optional thin client | `dynsyslab_client.periodic` |
| Analytic orbit sampling | `dynsyslab-cli analytic` (`periodic`, `histories` operations) |
| Separate observable/family/angle pages | native Viewer (`viewer/native/src/DocumentPayload.cpp`) |

All scientific algorithms execute in C++. The optional Python client only
transports requests. Viewer never fits, renormalizes, clips, interpolates,
or smooths the scientific values.

```console
dynsyslab-cli periodic-analysis job.json
```

Schema: `dynsyslab.periodic-query/v1`. Operations: `spectral_sweep`,
`spectral_compare`, `quadrature_audit`, `delay_fit`, `delay_forecast_audit`.
Results: `dynsyslab.periodic-result/v1`.
C++ workflow entry point: `dynsyslab::workflow::evaluatePeriodicQuery`.

`packet_features` exports a frozen causal packet map for downstream inference.
The separate `feature_identify` operation compares nonlinear orbit/phase fits
using arbitrary linear history features and full transformed noise covariance;
see [periodic feature identification](PERIODIC_FEATURE_IDENTIFICATION.md).

A small public synthetic input is
[`periodic-spectral-comparison.json`](../../resources/examples/periodic-spectral-comparison.json).
It supplies two analytic rotating circles and contains no fitted results.

## Input contract

Every `training` and `validation` row supplies `weight` and `channels`, each
channel being `{angular_frequency, coefficients:[[n,real,imag],...]}`. Channel
indices have one meaning across the bank. All channels on an orbit share its
physical angular frequency and phase convention. Optional `parameter` and
`family` are presentation coordinates. Separate orbit families get separate
Viewer pages: frequency alone need not distinguish branches.

Training masses and validation error weights are explicitly supplied. Each
bank's weights are normalized by their total; observable spectral masses are
**not** normalized. Supply frequency-quadrature weights for a frequency integral,
or appropriate transverse weights for an invariant canonical-area measure.

Common fields are `sample_interval`, `angles` (density grid, rad/sample),
`packet_centers` (independent of fitting), `bandwidth`, `order` (1–8), optional
`channel_names`, `audit_phase_counts`, and caller `source` provenance. Limits
are 16384 density angles and 64 packet centers.

`spectral_sweep` also takes `groups:[{name,channels}]`, `delay_counts`,
`phase_count`, `assume_measure_preserving:true`, optional `rank`, and
`svd_tolerance` (default `1e-10`). Each group/count is fitted from scratch with
matched settings. Signed offsets are
\((j-\lfloor N/2\rfloor)\Delta t\), including current time. Arbitrary offsets
are supported by the public C++ fitter and existing snapshot/fitting APIs.
Two selected channels produce twice as many features as one at the same delay
count; compare both feature count and retained rank.

`require_resolved_quadrature` defaults to true. For degree-\(M\) finite Fourier
polynomials, \(N>2M\) resolves all pairwise products of shifted copies. Explicitly
setting false permits an under-resolution experiment. The audit reports exact
Parseval and sampled masses; it cannot certify omitted harmonics or transverse
orbit quadrature. Snapshot row limits remain enforced, with no silent resolution
changes. Large fits can use the bounded-memory QR path described below.

## Bounded-memory fitting

Operation `delay_fit` fits a single linear history model without calculating
spectral densities or requiring a validation bank. Supply `training`, `channels`,
arbitrary signed `time_offsets`, `phase_count`, `sample_interval`, and `method`
(`edmd`, `resdmd`, or `mpedmd`). The same optional rank, SVD tolerance,
quadrature-resolution flag and storage settings apply. mpEDMD requires explicit
`assume_measure_preserving:true`; ordinary EDMD does not impose unitarity.
The methods reuse the same native fitting implementation and equivalent weighted
feature rows, so controlled comparisons can keep the trial space/data fixed.

The C++ entry is `fitPeriodicDelayModel`; `fitPeriodicDelaySpectrum` retains its
mpEDMD-only contract for packet analysis. The optional Python client exposes
it as `dynsyslab_client.periodic.fit_delay_model`. An optional `output_archive`
creates an immutable standard spectral KROM, with `periodic_fit_provenance_json`
containing the bank digest, source, coordinate contract and fitting settings.
The result always includes model identity, rank, original sample/feature counts
and the phase audit. A fit is a krom, written to `output_archive`; later operations name it as `archive`.
`include_model` (the operator, as `compression`) defaults to true without an archive and
false when one is written. Existing destinations reject before fitting.

This is standard [EDMD](https://arxiv.org/abs/1408.4408) or the existing mpEDMD
algorithm, not a new spectral estimator. A finite least-squares prediction is
not a certificate of nonlinear observability or long-horizon accuracy. Exact
tests use weighted cosine families whose scalar correlation is known in closed
form, separating unconstrained regression from the unitary scalar solution.
Both dense and QR paths are checked, along with archived analytic rotations.

`storage` accepts `dense`, `streaming_qr`, or `auto` (default). Auto selects QR
when a source matrix would exceed two million entries. It samples one orbit at
a time and reduces weighted blocks by Householder QR, retaining only the R
factor. If \(A=[X\;Y]\) denotes the **horizontal** concatenation of feature
matrices, the retained factor satisfies

\[
 R^\mathsf{T}R=A^\mathsf{T}WA.
\]

Only exactly identical time offsets share columns. Splitting the factor back
into source and successor columns preserves the weighted least-squares problem,
including the successor norm needed for residuals. Multiplication by
\(\sqrt{m_c/\sum_i w_i}\), with \(m_c\) compressed rows, compensates for the
existing fitter's normalized uniform weights. These are equivalent **feature
rows**, not physical states; this route applies only to the already evaluated
linear delay dictionary. Applying a nonlinear dictionary after compression
would change the problem and is not supported.

QR introduces floating-point roundoff but no Gram-matrix squaring, row
subsampling, rank truncation, or new fit regularization. The subsequent SVD
tolerance and mpEDMD algorithm are unchanged. Degenerate subspaces may choose
different bases; compare predictions and spectra, not raw basis coefficients.
`coordinate_contract` records storage, original samples, fit rows, and joint
feature count. The model's sample count remains the original count. Residuals
are training-quadrature diagnostics, not independent validation evidence.

The generic numerics API is `numerics::WeightedQrAccumulator`; the periodic
sampler is `math::periodicDelaySamples`. The accumulator accepts nonnegative
weights without normalizing them. Its block-plus-factor workspace and the
sampler each retain a fifty-million-entry limit; fitting retains one million
original rows and 4096 features as limits. No number of phases is silently
reduced to satisfy memory limits. Memory grows with one orbit and the dictionary,
not the number of orbits. This implementation is a straightforward sequential
row-block reduction, not an optimized parallel TSQR implementation.

With `output_directory`, immutable canonical NPZ model archives and
`periodic-spectral-request.npz` preserve the fitted models and complete request
(`metadata_json`). Existing destinations are rejected. A later failure can leave
earlier completed archives; compare those with `spectral_compare` rather than
overwriting. Model identities and SHA-256 bank digests are returned.
`include_models:true` also returns small model definitions inline.

`spectral_compare` takes `models:[{label,archive,channels,time_offsets}]`.
Inline `model` can replace `archive`. It does not refit. Newly fitted models
store versioned channel/offset contracts and reject mismatches. Legacy models
require per-model `allow_unverified_coordinate_contract:true`, retained in the
output metadata. Matching dimensions alone cannot verify coordinate identity.

## Causal free-forecast phase integrals

`delay_forecast_audit` takes an archived or inline periodic delay model,
`validation` Fourier orbits, and one `steps_per_orbit` integer per orbit.
Alternatively, `evaluation_steps` supplies one strictly increasing integer
list starting at zero per orbit. Use exactly one of these two forms. Every
returned row includes `step_indices`; its largest entry is `steps`. Sparse
horizons use cached binary matrix powers across long gaps, avoiding unused
Fourier/time output columns. No eigenbasis, matrix logarithm or interpolation
is used. Powers and ordinary repeated steps can differ by roundoff; this is
tested independently. The shared power cache has a fifty-million-entry limit.
The optional Python client entry point is
`dynsyslab_client.periodic.audit_delay_forecast`. The model must retain its
verified coordinate/clock contract and an unnormalized linear dictionary.
No fitting occurs. Offsets are shifted so their maximum is zero; the newest
coordinate is reconstructed at each native time step, including step zero.

For the supplied finite Fourier series, native Fourier translation and
Parseval give the exact uniform-starting-phase prediction and error mean
squares, up to roundoff. If the initial history coefficients at harmonic n are
the row vector b_n and the propagated linear readout is d_k, the error is

\[
 E_{k,c}=\sum_n\left|b_n d_{k,c}-a_{n,c}
                   e^{in\Omega k\Delta t}\right|^2.
\]

`reference_mean_square` has one entry per selected channel;
`prediction_mean_square` and `error_mean_square` use time/channel order.
`persistence_error_mean_square` is the independent baseline that holds the
true current value fixed. Its Fourier error coefficients are
`a_n * (1 - exp(i*n*Omega*k*dt))`; a small short-step error should be compared
with this baseline before attributing it to learned dynamics.
The result records model identity and a validation-bank digest. This removes
starting-phase quadrature error, not Fourier truncation, coefficient error,
time discretization, or transverse quadrature error. Real channels must be
conjugate symmetric. Coefficient/Fourier workspaces are limited to fifty million
entries per array and returned time/channel entries to ten million.

## Reference and errors

\[
 \rho_g(\vartheta)=\sum_jw_j\sum_n|a_{j,n}|^2
 K_\epsilon(\vartheta-n\Omega_j\Delta t),\qquad
 (P_\vartheta g)_{j,n}=a_{j,n}K_\epsilon(\vartheta-n\Omega_j\Delta t).
\]

Fitted packets use the original frozen dictionary and projected current-time
coordinate readout. They are compared against the original supplied observable,
so projection error is retained. Fourier translation and Parseval evaluate
packet RMS and complex-field error without another phase grid. Packet arrays
are center-major, then orbit. Summary packet norms use normalized validation
weights; density error is discrete grid L2. A zero reference norm yields a null
relative error, while absolute errors remain available.

Higher-order kernel densities may have negative lobes. A finite-width packet is
not an unsmoothed delta distribution. Amplitude agreement does not imply phase
agreement, off-orbit prediction, or closed-loop control.

## CLI, client and Viewer

```console
dynsyslab-cli periodic-analysis sweep.json > comparison.json
```

The optional client (`dynsyslab_client.periodic.spectral_sweep`,
`spectral_compare`, `quadrature_audit`) sends the same request to the compiled
CLI and returns its JSON result.

The native Viewer opens the saved result directly: `DynSysLabViewerCli render --source
comparison.json --page N --output page.png` (or `.svg`, or `.html` for every page,
interactive); `--view parameter` uses orbit parameters as the abscissa.

Fourier orbit banks for registered analytic systems come from the native
analytic sampler (`dynsyslab-cli analytic`, operation `periodic`; client
`dynsyslab_client.periodic.sample`), which preserves requests, resolved
parameters, equation fingerprints and coefficient-source descriptions. Native qualification rejects
unavailable periodic analytic solutions. FFT quadrature of an analytic orbit
is distinct from an exact closed-form Fourier formula. Independently derived
Fourier banks can be supplied directly.

Viewer uses raw signed-symlog densities and linear RMS amplitudes, with one page
per observable, family and center. Native Viewer owns wheel zoom, drag, rendering
and export; display payloads remain disposable and separate from primary data.

## Sources and tests

- Colbrook, *The mpEDMD Algorithm for Data-Driven Computations of Measure-Preserving
  Dynamical Systems*, [SIAM (2023)](https://doi.org/10.1137/22M1521407).
- Colbrook, Drysdale and Horning, *Rigged Dynamic Mode Decomposition*,
  [primary preprint](https://arxiv.org/abs/2405.00782), including Appendix B.
- Fourier translation, Parseval and periodic quadrature: Katznelson,
  *An Introduction to Harmonic Analysis*, chapter I.
- Binary integer matrix powers: Golub and Van Loan, *Matrix Computations*.
- Demmel, Grigori, Hoemmen and Langou, *Communication-optimal Parallel and
  Sequential QR and LU Factorizations*, [SIAM (2012)](https://doi.org/10.1137/080731992).

Native tests cover exact rotations, changed validation amplitude, redundant
dictionaries, repeated invariant eigenvalues, complex delay readouts and an
aliased harmonic. No private dataset is required by these tests.
Streaming tests additionally cover weighted blocks, a nearly dependent column,
irregular and regular delay grids, changed orbit weights, residuals, packet
fields and free rollouts against the dense fitter. Forecast phase integrals
are checked against independently evaluated phase ensembles for irregular
and regular delays, multiple orbits and channels. An exact rotation checks
the compiled client/archive contract.
Sparse evaluation is checked against ordinary repeated propagation. The
hold-current baseline is checked independently on an exact rotation. A
two-frequency scalar Galerkin example verifies that shrinking native one-step
errors need not imply accurate prediction over a fixed physical interval.
