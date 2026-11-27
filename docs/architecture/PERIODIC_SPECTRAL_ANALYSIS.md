# Periodic spectral analysis and estimation

These tools compose existing analytic system properties, Fourier analysis and
spectral fitting. They do not introduce a new dynamical-system registry or
archive store. Case-study data and scientific conclusions are not bundled.

Solutions expressed in an auxiliary time can use the native
[periodic-clock Fourier transform](PERIODIC_CLOCK_FOURIER.md). It includes the
clock Jacobian and phase correction before temporal coefficients are compared.

## Fourier coefficient covariance checks

`core/CoefficientShift.hpp` exposes `checkPeriodicCoefficientShift`.
For an autonomous periodic flow, a Fourier coefficient obeys
\[
 c_n(\Phi^s x)=e^{in\Omega s}c_n(x).
\]
The checker samples the registered analytic family, advances to each supplied
shift, independently restarts that state as an initial condition, and compares
its Fourier coefficients with this relation. It reports coefficient errors,
reference norms, initial-state consistency, period changes and endpoint errors.
This numerical check does not prove a global expansion or bound omitted tails.

Use the analytic CLI with `operation:"coefficient_shift"`, or
`dynsyslab_client.periodic.check_coefficient_shift`. Supply `system`,
`initial_state`, `observables`, `samples`, `maximum_harmonic` and `shifts`
in native model time; `initial_time` and `period` are optional.
The coefficient convention is the ordinary \(c_ne^{in\phi}\) convention.
At a repeated temporal frequency, the temporal coefficient is the aggregate,
not an identified basis of the eigenspace.

The covariance is a standard Fourier/Koopman relation; see
[Mezić (2005)](https://doi.org/10.1007/s11071-005-2824-x).

## Ownership and APIs

For arbitrary delay windows on registered analytic flows, use
`core::sampleAnalyticHistories` (`core/AnalyticHistories.hpp`), or the compiled
`analytic` query with `operation:"histories"`. Its optional thin Python adapter
is `dynsyslab_client.periodic.analytic_histories`.

Supply `system`, `initial_states` (rows of state vectors), `time_offsets`, and
named `observables`. Optional `initial_time` anchors every row, `family` selects
an analytic solution family, and `output_archive` writes an immutable NPZ.
The returned values have point/time/channel order and complex128 storage;
vector observables concatenate their channels. Each value is
\(g(\Phi^{t_0+\tau,t_0}(x))\), using the evaluator's physical/native clock.
Offsets may be negative, irregular, or span multiple periods. This operation
does not require a periodic orbit or infer an invariant sampling measure.
Unavailable analytic families and invalid values are rejected, never replaced
with numerical integration. Caller-supplied IC blocks limit memory; each block
is bounded to fifty million complex entries. Requests above 100000 point/time
pairs require an archive instead of inline JSON. Archives preserve ICs, offsets,
equation fingerprints, full requests and the selected solution family per row.

Parity with the pinned upstream mpEDMD and Rigged-DMD routines is held by
the native `dynsyslab_colbrook_reference_tests`, against frozen upstream output (see
[Spectral KROM](SPECTRAL_KROM.md)). That test also checks that compressing weighted
joint history rows by block QR, which keeps every weighted product (Demmel et al.,
doi:10.1137/080731992), leaves the fitted operator, density and wave packets
unchanged. It checks too that whitening the dictionary changes only coefficients.
The former MATLAB reference adapters were retired with the native port.

A small Gram defect after whitening tests the supplied quadrature only. It does not
certify dictionary completeness, the sampling measure, or a spectral approximation.

For refinement, hold the evaluated observable, target time, kernel, validation
states and error norm fixed. Vary spatial quadrature independently of dictionary
size. Longer delays can make products of dictionary functions harder to
integrate, especially across orbits with rapidly varying frequencies. More
dictionary columns on the same grid need not improve the result. Compare full
complex fields as well as amplitudes and scalar densities. This distinction
between data count, dictionary size and smoothing is explicit in
[Rigged DMD, sections 3 and 4.3](https://arxiv.org/html/2405.00782v2#S3), and
[mpEDMD, section 3](https://doi.org/10.1137/22M1521407).

Delay construction follows [Arbabi–Mezić (2017)](https://doi.org/10.1137/16M1105230);
equation-specific analytic evaluators retain their own references. The QR
reduction follows [Demmel et al. (2012)](https://doi.org/10.1137/080731992).

- Foundation: `AnalyticTrajectory.period` is a known period in native model time,
  zero when unavailable. An equilibrium does not receive an arbitrary period.
  `samplePeriodicTrajectory` samples a supplied evaluator and observables on
  one uniform physical-time mesh, omitting the repeated endpoint. It reports
  an endpoint closure residual. This does not certify a fundamental period.
- Mathematics: `RationalCircleKernel`, `spectralPacket`, `smoothOrbitSpectra`
  and `predictFrequencyBands` own spectral filtering and finite Fourier-bank
  diagnostics. The native Rigged DMD implementation shares these kernel poles
  and residues. `uniformFourierSamples` uses the existing inverse FFT.
- Mathematics: `PeriodicAtlas` interpolates frequency and Fourier coefficients
  across a scalar orbit parameter. Phase origins and channel order must agree.
  Extrapolation is rejected. This is a nonlinear frequency-conditioned model,
  not a finite constant-matrix Koopman invariant subspace.
- Estimation: `EnsembleBayesFilter` owns deterministic sequential importance
  weighting with supplied transition/readout callbacks. `filterPeriodicAtlas`
  composes it with a frozen unforced atlas and periodic lookup tables. It has
  no process-noise model, resampling or hidden access to the true state.
- Workflow and compiled CLI own JSON parsing and named-system binding. Optional
  `dynsyslab_client.periodic` and `dynsyslab_client.mathematics` are thin clients.

## Uniform-time orbit samples

`dynsyslab-cli analytic query.json` accepts the existing analytic schema with
`operation: "periodic"`, `system`, `initial_state`, `samples`, and a list of
scalar registered `observables`. Optional fields are `family`, `initial_time`,
`period` and `maximum_harmonic`. If the system evaluator has no period metadata,
an explicit positive period is necessary. The analytic family must be applicable.

The result contains sample-major `values`, complex `[real,imag]` measurements,
channel names, period, angular frequency, closure error, resolved parameters,
dimensionless groups and the equations fingerprint. Optional `fourier_series`
contains one object per channel with coefficients `[harmonic,real,imag]` for
the convention \(c_n e^{in\phi}\), \(\phi=\Omega(t-t_0)\).
Require \(2n_{\max}<N\); increase the sample count to check aliasing.

```python
from dynsyslab_client.periodic import sample
result = sample({
    "system": {"name": "rpe", "parameters": {"mu_pa_s": 0}},
    "initial_state": [1.1, 0], "samples": 1024,
    "observables": ["radius", "far_field_acoustic_echo"],
    "maximum_harmonic": 64,
})
```

An elliptic auxiliary clock is not interchangeable with native model time.
The current analytic RPE family exposes its inverted-clock physical period.

## Spectral mathematics queries

Use `dynsyslab-cli mathematics query.json`, schema `dynsyslab.math-query/v1`.

| Operation | Fields and result |
|---|---|
| `spectral_filter`, action `kernel` | `bandwidth`, optional `order` (1–8), `angles`, optional integer `lags`; returns kernel values and time multipliers. |
| `spectral_filter`, action `packet` | Same kernel settings, `series`, `sample_interval`, `angle`; returns filtered Fourier coefficients. |
| `spectral_filter`, action `density` | Same kernel settings, `orbits` with `series` and explicit quadrature `weight`, `sample_interval`, `angles`; returns density, total mass and zero-harmonic mass. |
| `frequency_band_prediction` | `series`, `initial_phase`, `elapsed_times`, `band_width`; returns predictions and an error bound for the supplied coefficients only. Zero width retains exact finite-series frequencies. |
| `periodic_atlas` | `knots` with `parameter` and Fourier `channels`, requested `parameter`, `phase`, optional `elapsed`; returns interpolated series and predicted values. |
| `fourier`, action `to_samples` | `series`, `samples`; returns uniform-phase values using the native inverse FFT. |

The circle kernel has unit integral. Higher orders may have negative sidelobes.
The `kernel` action also accepts `tail_after_lag`. It returns an absolute
coefficient-sum bound for all kernel Fourier lags beyond that symmetric history
window, including the \(1/(2\pi)\) normalization. Multiplying this bound by
the observable norm bounds truncation of the exact unitary time-domain filter.
It does not include learned-operator, quadrature, conditioning or measurement
errors, and is not a floating-point enclosure.
Weights are not automatically normalized. Neither a finite orbit quadrature
nor a finite matrix spectrum certifies that the original spectrum is pure point
or continuous. Scalar spectral densities do not contain the spatial phase
information of spectral packets. Sampled spectra contain aliases.

## Causal periodic-family filter

`dynsyslab-cli periodic-analysis query.json` takes schema
`dynsyslab.periodic-query/v1`, operation `filter`, and:

- `knots`: the frozen atlas, as above;
- increasing `parameters`: particle parameter grid within the atlas;
- increasing `time`, flat sample-major `measurements`;
- `measurement_channels`, `output_channels`: zero-based channel indices;
- positive `measurement_variances` per measured channel;
- optional `parameter_prior`, `phase_count` (default 256), `table_samples`
  (default 4096), and `frozen_frequency` for a diagnostic constant-frequency
  baseline (default zero means each orbit uses its own frequency).

The prior is uniform over initial phases, with the specified parameter masses.
Nonuniform parameter grids require explicit prior masses for a desired density.
Every correction uses only observations at or before its output time. Outputs
include decoded means, parameter mean/variance, effective sample sizes,
conditional log evidences and final weights. The table uses linear periodic
interpolation; refine both the Fourier bank and table for sharp waveforms.

This is an **unforced periodic-family estimator**. It does not support changes
in orbit parameter caused by control or process noise. Finite particles can
lose effective sample size or miss the correct phase/frequency. A small
posterior variance is conditional on this discrete model and is not a calibrated
guarantee under mismatch. Callbacks in the generic ensemble API permit separate
transition models, but no controlled atlas estimator is inferred automatically.

## Efficient snapshot input

`periodic-analysis` operation `snapshot_pairs` constructs those arrays natively
from `orbits`, each with `channels` (Fourier series) and an explicit `weight`.
Supply `phase_count`, signed `time_offsets`, `sample_interval`, and a new
`output_archive` path. Negative offsets are past delays; positive offsets are
future values and cannot be used for a causal online initialization. Optional
`phase_origin` shifts all phase grids. `split_complex:true` makes separate real
and imaginary columns; otherwise genuinely complex samples are rejected.
Rows are orbit then phase, columns offset then channel then real/imaginary.
Weights are orbit mass divided by phase count, without automatic normalization.
Optional `source` and `measure` fields preserve provenance in the canonical NPZ.
The optional Python adapter is `periodic.snapshot_pairs`.

This computes exact point values of the **supplied finite series** using folded
FFT indices. High harmonics are not discarded when phase sampling is sparse.
This differs from inferring coefficients from an unresolved sampled function.
The mathematics operation `fourier`, action `sample_values`, exposes the same
evaluation with `series`, `samples`, and optional `phase_origin`.

`analytic_fourier` provides Jacobi `sn`, `cn`, `dn` coefficients (DLMF 22.11)
and `log_cosine`, meaning \(\log(1+k\cos\phi)\). Inputs are `parameter`,
`maximum_harmonic` and optional `angular_frequency`. Outputs include `series`,
`parameter_period` and `uniform_tail_bound` for the entire omitted series in
exact arithmetic. The bound is not a floating-point enclosure. `sn`/`cn` use
\(\phi=2\pi u/(4K)\); `dn` uses \(\phi=2\pi u/(2K)\).

The existing `koopman-spectral` query also accepts `data: {"archive": "pairs.npz"}`.
Use the standard NPZ reader/writer with equal nonempty sample-by-coordinate
`x` and `y` arrays, optional length-sample `weights`, and `metadata_json` whose
schema is `dynsyslab.koopman-snapshot-pairs/v1`. Inline and archived pairs are
mutually exclusive. This avoids decimal JSON copies of large snapshot arrays.
The declared coordinate contract, method and sample interval remain in the fit
settings. No trajectory boundaries are inferred from raw arrays.

When SVD removes redundant dictionary directions, raw observable coefficient
vectors need not lie in the retained coefficient subspace even if their
functions are representable. A `rigged_measures` item can explicitly choose
`reconstruction: {"name":"coordinates","output":0}` instead of `coefficients`.
The selected fitted readout already lies in the retained subspace. This choice
is recorded in measure/packet assumptions. It measures that projected
observable; inspect readout error and held-out projection error separately.

## Window identification and free prediction

Operation `predict` evaluates an existing atlas without fitting: supply `knots`,
`parameter`, `initial_phase`, and `time`. The phase is anchored at `initial_time`
(default zero); `table_samples` and `frozen_frequency` are optional. This allows
dense display clocks without changing the data or coefficients used in a fit.
Python exposes `periodic.predict_orbit`.

`periodic-analysis` operation `identify` uses the same `knots`, measured `time`,
`measurements`, channels, variances and table options as the filter. Provide
`initial_candidates`, each `[parameter,phase]`, and optional `maximum_iterations`.
The native bounded least-squares solver refines parameter and initial phase
against that supplied measurement window. Optional `prediction_times` then
produces a free forecast; no measurements at those additional times enter the
objective. The return reports convergence, iterations, evaluations and squared
residual, not a claim of a unique/global minimum. Table/atlas interpolation,
noise and narrow collapse peaks can affect long forecasts. This remains an
unforced constant-parameter model.

## Algorithm sources and verification

For matched delay sweeps, frozen-model reference comparisons, phase-quadrature
qualification and separate Viewer pages, use
[PERIODIC_SPECTRAL_COMPARISON.md](PERIODIC_SPECTRAL_COMPARISON.md).

### Exact finite-Fourier delay readouts

`math::periodicDelayReadout(channels, timeOffsets, coefficients)` evaluates a
linear readout of an arbitrary signed-delay dictionary in Fourier coefficient
space. Columns are delay-major, then channel. Complex coefficients and channels
are supported. All channels must describe the same orbit frequency. The Fourier
translation theorem gives

\[
 b_n=\sum_{j,c} C_{jc}a_{c,n}\exp(i n\Omega t_j).
\]

Consequently `meanSquare` computes its squared phase RMS by Parseval without
phase-grid aliasing or a large snapshot matrix. The compiled mathematics query
`operation="periodic_delay_readout"` accepts `channels`, `time_offsets`, and
complex `coefficients`. Alternatively, `readouts` batches coefficient vectors;
`summary_only=true` suppresses returned coefficient series. A readout's optional
`reference_channel` compares it to that channel's spectral packet using the
top-level `angle`, `sample_interval`, `bandwidth`, and `order`, returning both
`reference_mean_square` and `error_mean_square`. These are finite-series
identities; omitted input Fourier tails remain the caller's responsibility.

Uniform offsets use a phase recurrence restarted every 32 terms; arbitrary
offsets use direct phases. Native tests compare complex multichannel readouts
against direct time shifts and independent phase quadrature (Katznelson,
*An Introduction to Harmonic Analysis*, chapter I).

Controlled planar orbit-family transport is described separately in
[PERIODIC_CONTROL_TRANSPORT.md](PERIODIC_CONTROL_TRANSPORT.md). It retains both
amplitude and phase input effects and uses the existing native integrator.

- Colbrook, Drysdale and Horning, *Rigged DMD*, especially Appendix B:
  [primary paper](https://arxiv.org/abs/2405.00782).
- Mezić, spectral expansions on periodic families:
  [Nonlinear Dynamics 41 (2005)](https://doi.org/10.1007/s11071-005-2824-x).
- Doucet, de Freitas and Gordon, *Sequential Monte Carlo Methods in Practice*:
  [Springer (2001)](https://doi.org/10.1007/978-1-4757-3437-9).
- Korda, Putinar and Mezić, finite spectral projections:
  [primary paper](https://arxiv.org/abs/1710.06532).

Native tests check kernel moments/mass/time multipliers, an independent
two-sided Poisson projector, Fourier quantization error bounds, physical-time
period closure, exact Gaussian Bayes weights, held-out parameter interpolation
and estimator prefix causality. These tests are numerical diagnostics and
contract checks, not global certification of a learned model.
