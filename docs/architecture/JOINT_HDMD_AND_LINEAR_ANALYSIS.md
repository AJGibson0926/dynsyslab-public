# Joint H-DMD reconstruction and linear-system diagnostics

The native postprocessor supports one shared Hankel decomposition of any number
of real scalar observable channels. It retains the existing KROM, observable,
reconstruction, archive, and KLQE classes. No second estimation pipeline is used.

See the [observable Hankel contract](OBSERVABLE_HANKEL.md) for shared matrices, channel taps, derivative links, and signed POD coordinates.

## Joint H-DMD

An `hdmd` job can use these settings:

```json
{
  "observables": ["radius", "far_field_acoustic_echo"],
  "matrix_size": 3200,
  "sampling_window": 12751,
  "sample_stride": 1,
  "snapshot_stride": 245,
  "window_alignment": "end",
  "channel_scales": [1, 1],
  "rank": 3200,
  "svd_tolerance": 1e-10,
  "reconstruction_observables": ["state"]
}
```

`matrix_size` is the number of fitted window pairs; `sampling_window` is the
number of delays per channel. The delay-major row index is `delay * channels +
channel`. `sample_stride` decimates the measurement grid before forming delays.
`snapshot_stride` spaces the starts of training windows on that grid; the paired
X-to-Y shift remains ONE grid sample. Thus widely spaced fitting windows do not
change the inferred eigenvalue dt or discard samples inside each history.

For start index j and lag l, X contains y[(j*snapshot_stride+l)*sample_stride],
and Y is shifted by sample_stride original samples. The recorded eigenfunction
values contain the X columns plus the last Y column. Their training times are
explicit; when snapshot_stride differs from one the final interval is shorter.
They must not be assumed to be uniformly spaced at the DMD pair dt.

The default is `window_alignment: "end"` for every window setting, the legacy
keys included (before 0.9.1 the legacy keys defaulted to `"start"`). End
alignment labels each window by its newest sample, so the value at time t uses
only the history up to t. This is causal once a complete history is available,
and the first window length of a record has no value; no zero-padding or future
measurements are supplied to an observer. `"start"` labels a window by its
oldest sample, so the value at t uses samples up to t + T and the last window
length of a record has no value. Recorded training window
coordinates are reference projections, not guaranteed exact eigenfunctions.

Positive channel scales divide inputs before fitting; the corresponding C rows
are multiplied back into physical observable units. There is no implicit mean
subtraction. Every channel receives a named reconstruction `reconstruction__NAME`
and the names are retained in the existing KROM map. Catalog aliases resolve to
canonical names. Additional named observables, including vector `state`, are
fitted against the same Psi using the existing least-squares reconstruction.

The scalar `observable` setting and C++ overloads remain compatible. The public
C++ sampled-channel overload accepts a samples-by-channels table and unique
channel names without requiring equations of motion. Another overload evaluates
a vector of supplied scalar Observable objects. All channels share the same
modal ordering and normalization. `hankel_singular_values` records the complete
snapshot SVD spectrum; requested rank remains capped by the numerical threshold.

## Fitting a new output to an existing KROM

The C++ `fitAndSetObservableReconstruction` overloads accept either (a) paired
Psi and arbitrary measured output arrays, or (b) a KROM evaluator, states, and an
Observable. They fit G approximately equal to Psi times C-transpose and retain
C under its name. A nonnegative ridge penalty is optional. No output is silently
refitted during KROM subset extraction.

The compiled workflow method `observable_reconstruction` creates a new immutable
KROM archive. Settings include `krom`, `observable`, and optional `ridge`. For a
catalog observable, the source is a matching state trajectory. For measured data,
provide `observable_array` and optionally `eigenfunction_array` (default
`eigenfunction_values`) in one aligned source archive. The caller is responsible
for supplying Psi in the specified KROM's exact ordering and normalization; array
shape agreement alone does not prove identity. Existing model training pointers
are retained and the reconstruction source is recorded separately.

## Observability and controllability

Public C++ surfaces are `Observability.hpp`, `Controllability.hpp`, `Pbh.hpp`,
`ModalObservability.hpp`, and `EstimationSelection.hpp`, exposed alongside the
existing Estimator and Controller APIs. Their shared implementation is a small
native linear-analysis module; neither optional package depends on the other.

`analyzeObservability(A,C,...)` and `analyzeControllability(A,B,...)` provide:

- Classical stacked observability / horizontal controllability matrices and
  their SVDs, with state-space directions ordered strongest first.
- Finite discrete sums and finite continuous Gramians.
- Infinite discrete and continuous Gramians, rejecting non-strictly-stable A.
- Rank, threshold, effective rank, Gramian eigenvalues, and conditioning. Overall
  condition number is infinite when rank deficient; resolved conditioning is
  reported separately for the retained numerical subspace.

For matrix inputs, a discrete Gramian with sampleTime=0 treats A as the supplied
transition. Positive sampleTime interprets A as a continuous generator and uses
exp(A*dt). KROM observability overloads take a reconstruction name and require an
explicit sampleTime for discrete analysis. KROM controllability accepts constant
continuous maps or a stored discrete input map at its recorded sample interval.
A state-dependent B is not silently treated as an LTI matrix; pass a stated
operating-point matrix for local analysis.

PBH tests provide mode-by-mode singular-value margins, observability or
controllability, and detectability or stabilizability, in either time domain.
These dense per-pole tests are intended for moderate models. Large diagonal
Koopman banks can use `analyzeModalObservability`: it assembles finite/infinite
complex Gramians analytically in O(modes squared) before the spectral
decomposition, rather than forming thousands of powers of a dense generator.
It accepts explicit modal scales. Its complex-coordinate directions must not be
confused with independent real degrees of freedom when conjugate pairs exist.

The `linear_analysis` compiled workflow method accepts a KROM as source and
`measurement` and/or `control_map` names. Alternatively the source can contain
raw A plus C and/or B arrays. Settings are `analysis_method`, `blocks`, `horizon`,
`sample_time`, and `relative_tolerance`; output includes matrices, Gramians,
spectra, directions, and metadata. This diagnostic does not require an LQE solve.

`analysis_method: "modal-matrix"` (K-ROM source) forms the observability matrix in
discrete modal coordinates, `O = [C; CM; …; CM^(n−1)]` with `M = diag(exp(λΔt))`, one
complex coordinate per eigenfunction (a conjugate pair is two modes, not four split real
states), and takes its direct complex SVD (`modalObservabilityMatrix`,
`ModalObservability.hpp`). Settings: `measurement`, `sample_time` (default: the K-ROM's
sample interval; required when it records none), `blocks` (default: the mode count),
`output_rows`. Arrays: `modal_observability_matrix`, `modal_singular_values`,
`modal_directions` (complex, column k = direction k), `modal_direction_magnitudes`,
`modal_multipliers`, `continuous_eigenvalues`. The singular values depend on Δt and the
block count: n blocks span only nΔt of output, so with a data-rate Δt neighbouring
directions differ by tiny phase increments and the spectrum falls off steeply.

`analysis_method: "eigenfunction-observability"` (K-ROM source) answers how well each
eigenfunction can be determined from the measurement when every other eigenfunction is
also unknown (`eigenfunctionObservability`, `EigenfunctionObservability.hpp`). With the
finite-horizon modal gramian `G` over `[0, horizon]` and output noise intensity `V`, it
reports per eigenfunction: `signal = sqrt(G_nn)` (how strongly ψₙ alone moves the output),
`distinguishable = 1/sqrt((G⁻¹)_nn)` (the part of ψₙ's output no combination of the other
eigenfunctions can imitate), `isolation` (their ratio, 0–1), the Cramér–Rao `bound_std =
sqrt(V (G⁻¹)_nn)`, `pair_bound_std` (the worst direction of a conjugate pair, which a real
output determines jointly), and with `process_covariance` W > 0 the KLQE's steady-state
`kalman_std` of the conjugate-consistent estimate ½(ψₙ + conj ψ_partner). Settings:
`measurement`, `horizon` (default 1), `noise_intensity` (1e-4), `regularization` (relative
Tikhonov floor on G, 1e-12), `process_covariance` (0 = skip). Rows of V from the
observability SVD carry the same information: `(G⁻¹)_nn = Σₖ |V_nk|²/σₖ²`. The case-study
driver `observability/eigenfunction_observability.py` plots both panels as stems for any
list of K-ROMs (`--krom name=path`).

The fuller classical toolkit (balanced truncation, Kalman decomposition, norms,
zeros, design) and its `linear_systems` method are described in the
[linear-systems toolkit](LINEAR_SYSTEMS_TOOLKIT.md).

## Subset selection and real realization

`makeRealModalRealization` verifies sampled real/conjugate relationships,
normalizes retained coordinates by RMS, and transforms every C map consistently.
It produces an ordinary sampled KROM. It does not refit C, and records source and
partner indices and normalization divisors. This avoids treating both complex
conjugates as four independent real observer unknowns. The equivalence asserted
is the reconstruction of real outputs on the verified coordinate relation.
An explicit real-eigenvalue tolerance can remove numerical imaginary residue;
the sampled real-up-to-phase check must also pass. Record that tolerance when
publishing a realization. This is a numerical correction, not a stability repair.

`assessEstimationSubset` reports state truncation MSE in caller-supplied accuracy
units, relative echo reconstruction MSE, autonomous propagation MSE, and a
state-weighted posterior-noise benchmark. It uses the normalized finite-window
observability Gramian with a declared unit Gaussian prior and independent
measurement noise. Weak directions retain prior uncertainty, rather than being
assigned zero uncertainty through a bare pseudoinverse. The benchmark is
propagated to the observation-window endpoint before decoding x. The component
scores remain available and weighting is explicit. This is not a certificate of
observer convergence under model mismatch.

The case-study protocol documents the required ranking and validation comparisons
in `examples/native/case-studies/rpe-joint-hdmd-klqe-20261105/PROTOCOL.md`.

The RPE experiment report (private development record `JOINT_HDMD_KLQE_SELECTION_2026-11-05.html`)
records the 3,188-mode bank, failed subset observers, and conditioning comparisons.

Definitions follow the standard [observability matrix](https://www.mathworks.com/help/control/ref/statespacemodel.obsv.html),
[controllability matrix](https://www.mathworks.com/help/control/ref/statespacemodel.ctrb.html),
and [Gramians](https://www.mathworks.com/help/control/ref/statespacemodel.gram.html).
