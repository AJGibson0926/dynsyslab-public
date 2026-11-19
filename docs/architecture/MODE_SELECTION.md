# Eigenfunction subset selection for echo-only KLQE → KLQR

## The design problem

Online, the estimator measures the echo and nothing else. Its measurement equation is the truncated
Koopman expansion of the echo,

    y = e ≈ e₀ + Σₙ νₙ ψₙ = e₀ + C_e ψ,

over the retained eigenfunctions ψₙ. ψ is never measured or evaluated online; KLQE infers it from
the echo through ż = Λz + B_ψ u. The state is read out the same way, x ≈ x₀ + C_x ψ.

Truncation error is unavoidable, so the quality of the controller is decided by *which*
eigenfunctions are retained. One subset must satisfy four demands at once:

1. the echo is reconstructed well (C_e);
2. ψ is well observable from the echo;
3. the state x = (r, ṙ) is reconstructed well (C_x);
4. KLQR controls well with a constant B_ψ regressed from forced data.

Changing the subset changes C_e, C_x and which B_ψ entries are present. The `mode_selection`
post-processing method scores subsets against all four demands and returns the Pareto front.

## Method

Sources: unforced records (readouts, invariance) and forced records with one scalar input (B_ψ).
`krom`: the candidate pool.

- **Pools.** Polynomial coordinate pools ψ(e, e′) are evaluated on each coordinate record's
  *measured* columns (a `coordinate_transform` archive with `source_states`), never through the
  analytic coordinate expressions. State pools, including interpolated H-DMD eigenfunctions, are
  evaluated on recorded states.
- **Units.** A real mode or a conjugate pair. Growing and marginal modes are skipped.
- **Per unit (derivative-free, horizon h = `invariance_horizon`, default 0.25):**
  - invariance: ‖ψ(t+h) − e^{λh}ψ(t)‖ / ‖ψ(t+h) − ψ(t)‖ on unforced records;
  - process-noise intensity for the predicted KLQE: var(ψ(t+h) − e^{λh}ψ(t)) / h;
  - constant B_ψ: ψ(t+h) − e^{λh}ψ(t) = B_ψ ∫₀ʰ e^{λ(h−s)} u(t+s) ds, least squares over forced
    records, with its unexplained fraction.
- **Per subset**, in scaled real coordinates (Re ψ, Im ψ)/scale:
  - `echo_error`, state readout errors: least-squares readouts with intercepts on unforced records;
  - `predicted_state_error`: steady-state Kalman covariance P of (Λ, C_e) with the process noise
    above and measurement noise (readout residual² + `sensor_noise`²)·var e, carried into the
    state, √(ε_x² + (C_x P C_xᵀ)ᵢᵢ / var xᵢ);
  - `input_map_residual`: worst unit;
  - `state_reach`: √λ_min of the normalized output controllability gramian D⁻¹ C_x W_c C_xᵀ D⁻¹;
  - diagnostics: observability-gramian balance, PBH margin, Hankel singular values of
    (Λ, B_ψ, C_e) and (Λ, B_ψ, C_x).
- **Search.** The `pool_size` units that explain the most echo alone; every subset of up to
  `max_units` of them. Pareto front over echo error, predicted state error, B_ψ residual and state
  reach; a rank-sum composite picks one subset.
- **Output.** A `dynsyslab.mode-selection/v1` archive (units, subset metrics, Pareto flags) and,
  with `write_krom`, a ready K-ROM of the chosen subset (`derived/koopman-mode-selection`) with
  reconstructions `<measurement>` and `state` and the constant map `constant`, usable directly by
  an echo-only estimator (`measurement`, `output: state`) and KLQR (`control_representation:
  constant`).

| Setting | Default | Meaning |
|---|---|---|
| `krom` | — | candidate pool |
| `measurement` | `far_field_acoustic_echo` | coordinate column or catalog observable |
| `time_window`, `forced_time_window` | all | rows of unforced / forced records |
| `max_samples_per_record`, `sample_stride` | 20000, auto | thinning |
| `invariance_horizon` | 0.25 | h of the integral form |
| `pool_size`, `max_units` | 10, 4 | search size |
| `sensor_noise` | 0.01 | sensor noise relative to the echo's standard deviation |
| `stability_margin`, `ridge` | 1e-6, 1e-10 | |
| `write_krom` | false | write the chosen subset as a K-ROM |
| `krom_units` | composite pick | unit indices to write instead |
| `eigenfunction_polynomial_degree` | none | replace a scattered (H-DMD) pool by a polynomial fit of the same eigenfunction values first |
| `readout_records` | `"unforced"` | records that fit C_e and C_x: `"unforced"`, `"forced"`, `"all"`, or a list of source indices |

**Readouts where the controller operates.** C_e and C_x are a static map ψ → (e, x), so any
offline record with known e and x can fit them. Fitted over a wide range (e.g. a collapse
decay), a one-pair readout is right on average but biased on a particular orbit, and an
echo-only KLQE then shifts ẑ to cancel the echo-readout error, biasing the lifted state the
controller acts on. `readout_records` fits them on the controller's operating orbit instead.

**Analytic-coordinate pools on state records.** A pool whose coordinates are analytic
observables (recentred (r, ṙ), or (e, e′) with the analytic e′) is evaluated through that map
on the recorded true state when a record is a plain state record (offline only; recorded in
`records[].e_dot`). Measured or causal-history coordinates still need coordinate records.

**Smooth pools.** H-DMD eigenfunctions exist only at training states; a scattered K-ROM
interpolates between them by inverse distance over the 8 nearest. That surface is continuous
but rough: its slope jumps whenever the neighbour set changes, which is worst between the turns
of a single decaying training trajectory. With `eigenfunction_polynomial_degree: d` the pool is
first refitted (least squares on the training values, all monomials of total degree ≤ d in the
normalized state; same eigenvalues) and saved in `derived/koopman-polynomial-fit`; scores,
readouts, B_ψ and the written K-ROM then all use the smooth surface. The archive records the
per-mode fit errors in `eigenfunction_polynomial_fit`. On the 2026-11-16 RPE pool, degree 5 fits
the fundamental pair to 3e-4 and removes the closed-loop drive chatter (case study).

The written K-ROM rewrites each pair's lower eigenfunction as exactly conj(upper) before attaching
the readouts and B_ψ. H-DMD normalizes every mode separately, so its stored partner is
c·conj(ψ) with an arbitrary complex c; without the rewrite the conjugate readouts and map are
inconsistent and the controller sees the wrong reference (found in the first closed-loop runs).

## Where a numerically differentiated e′ enters

Recorded for every run in the archive's `e_dot_entry_points` and per record (`records[].e_dot`):

| Place | Offline / online |
|---|---|
| `coordinate_transform`: e′ of each (e, e′) record, a numerical derivative of the recorded echo | offline |
| SINDy ψ(e, e′) with `input_space: coordinates`: ψ̇ = ∇ψ · (e′, e″) from the record's measured derivatives | offline |
| ψ(e, e′) evaluated on records (readouts, invariance, B_ψ) | offline |
| ψ̇ | never differentiated numerically: the integral form replaces it |
| echo-only KLQE / KLQR | no e′ online |

H-DMD pools on (r, ṙ, e) involve no e′ at all.

## Tests

`dynsyslab_mode_selection_tests`: exact exponential eigenfunctions with a conjugate pair and a real
mode, a constant-input forced record (B_ψ recovered to 2e-4), exact readouts, intercepts, and a
finite predicted error with positive reach and observability.

## Results

See `examples/native/case-studies/rpe-echo-mode-selection-20261116/README.md`.
