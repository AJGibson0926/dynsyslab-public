# K-ROM control and estimation tools

Native replacements for the Python `rom/`, `control/`, `estimation/` and K-ROM training code
(port phase E). Every item is a job key, a post-processing method or a CLI command over the C++
core; sources and papers are cited in the code comments.

## Estimators in the solve loop

`estimator.kind` is `klqe` (Koopman, `EstimatorFactory.cpp`) or a classical estimator on a
linear(ized) plant (`core/LinearObserver.hpp`, `ClassicalEstimator.cpp`):

| kind | Keys |
|---|---|
| `lqe` | `A`, `B` (default: Jacobian and input matrix at `linearization_state`), `C` (default identity, or the sensors' Jacobian), `measurement_observables`, `process_covariance` W, `measurement_covariance` V. Gain from the filter Riccati equation (Kalman & Bucy 1961). |
| `luenberger` | as `lqe` with a given gain `L` (Luenberger 1964). |
| `static_map` | `x_hat = M y + offset`. |

The observer is advanced exactly over each sample with the input and measurement held. Traces are
archived as `estimator_time`, `state_estimate`, `innovation`, `estimator_gain`, `observer_poles`;
`controller.settings.feedback = "estimated_state"` runs any controller on the estimate.
`estimator.design_archive` reloads a published `klqe` design (gain, covariance, poles) after
checking its `krom_identity` and measurement matrix against the bound K-ROM.

## Targets

`target.kind`: `constant`, `waveform`, `trajectory` (`interpolation` `linear` or `pchip`,
Fritsch & Carlson 1980), `expression` (`components`: one expression of `t` per state, with
`parameters`), `rpe_canonical`. `target.space` is `state`, `observable` or `lifted` (z_ref in the
K-ROM's real split coordinates; KLQR/KMPC). `controlTargetEnvelope` gives bounds containing a
target (analytic for harmonic targets: offset +- sum |a_k|, velocity +- sum |a_k w_k|).

## Control

- Multi-input `lqr`, `mpc`, `pid`: an explicit n x m `B` (`settings.B` or `input_matrix`) with
  `actuation: "additive_state"`; `R` m x m; `control_limit` per input (`ControlFactoryClassical.cpp`).
- Multi-input `klqr`/`kmpc`: `actuation: "additive_state"` with an n x m `B`; the exact map is
  grad psi(x) B (`KoopmanModel::setAdditiveInputMatrix`), `R` is m x m and a stored map needs m
  columns. History and online-measurement coordinates stay single-input.
- `sdre` with `sdc_A` (n x n expressions of the state names and parameters, f(x) = A(x) x) and
  `sdc_B` (Cloutier 1997).
- KLQR/KMPC `target_domain_policy`: `warn` (default; archives `domain_report_counts`,
  `domain_report_maximum_excursion`, `domain_report_bounds`), `error` (refuses), `ignore`. The
  domain is `settings.training_domain` or the K-ROM's `training_minimum`/`training_maximum` or
  `training_states` (`KromDomainPolicy.cpp`).

## Post-processing methods

| Method | What it does |
|---|---|
| `krom_extract` `mode_selection.rule` | adds `state_observability` (greedy grouped R^2 to the state, optional control-authority weight, per-step report), `nearest_frequency` (`target_rate` or `target_frequency`, `band`) and `validation_score` (anchored error on `validation_sources`). The report is `mode_selection_report`. |
| `koopman_validation` | `windows` (re-anchored relative L2 per window), `binned_defect` (one-step defect by state bin; with `axes` [i, j], `bins` or `x_edges`/`y_edges` and `horizon`, the log10 \|psi\|-weighted finite-horizon defect over occupied 2-D cells with the occupied mask), `segments` {`starts`, `horizon`} (relative L2 of psi(t+h) against exp(lambda h) psi(t) per segment), `history` (psi and exp(lambda t) psi_0). |
| `observable_reconstruction` `quality` | held-out gate: `validation_fraction` with `split` interleaved or tail, or `validation_source`; `minimum_r2`, `maximum_normalized_rmse`, `maximum_condition_number`; refuses a failing readout unless `allow_quality_failure` (then `validated: false`). |
| `input_map` `fit: "symmetric_pulse"` | B_psi from paired +-u one-sample responses (`symmetric_pulse` input, `pair_id`/`polarity`), fitted constant or bilinear. Also reads `pulse_branches` solves: unpaired branches are skipped, pairs must last one sample, and the solve must use `control_evaluation: "sampled"`. |
| `pulse_response` | paired causal sensitivities G(t) = (psi(x_+) - psi(x_-)) / (u_+ - u_-) of `pulse_branches` solves, per anchor and mode (WF/PulseControlPostprocess.cpp): immediate and peak sensitivity, time to peak, the semigroup error of exp(lambda (t - t_1)) G(t_1) overall and per `bands` [from, to], the zero branch's drift over the pulse effect and per band, the even (nonlinear) over odd response, the change of G with amplitude (`amplitude_checks`), the post-pulse semigroup error of wider pulses (`width_checks`) and, with `reference_sources` (the same branches solved with a finer step), `reference_step_sensitivity_change`. `evaluate`: `state` (default) or `coordinates`. Arrays `time`, `causal_sensitivity` (anchors x samples x modes). Korda & Mezic 2018. |
| `control_map_audit` | for each stored input map (`maps`, default all), the representation error against the exact map grad psi(x) B along each source path: `rmse`, `relative_frobenius`, and the median / p95 / maximum pointwise relative error (`relative_floor`). For each record in `training_sources`, the persistent-excitation design audit: input and state ranges (with p01/p99), the rank of the constant design [u], and the rank, shape and 2-norm condition number of the bilinear design [u_j psi(x), u_j] every `audit_stride` samples (core/ControlMapAudit.hpp; Willems et al. 2005). The exact map enters only this audit. |
| `control_regression` | error grids now report median and p95 (`relative_floor` in the core). |
| `control_error_compare` | ranks control-map candidates (`candidates: [{krom, map, label}]`, `grid`, `metric`: relative mean/median/p95/maximum, absolute mean/p95/maximum) without selecting a winner. |
| `krom_rollout` | psi(t) = exp(Lambda (t - t0)) psi(x0) at any increasing times, decoded by `reconstruction` (a `time_series` object). |
| `pod_galerkin` | POD modes, reference, reduced A and B from trajectory snapshots (Holmes, Lumley & Berkooz 1996). |
| `krom_refit` | refit each eigenfunction at its fixed multiplier, gauge-align to the original (or, with `gauge_reference`, to another K-ROM's mode of nearest eigenvalue, e.g. a hand-derived eigenfunction), prune small terms, report anchored errors before and after; `orient` fixes each phase so the largest-magnitude value on the training runs is real positive, and `coefficient_summary` (a relative floor) lists each mode's terms with \|c\| >= floor max\|c\|, largest first. |
| `reconstruction_transfer` | moves a readout C to another K-ROM after per-mode scale alignment, with the alignment residuals. |
| `control_score` | adds `tone` (amplitude, phase, residual, harmonic content), `settling_band`, `drive_lowpass_hz`, `estimate_rmse`, `innovation_lag1_autocorrelation`, `compare_to` run differences, and `estimation` scores (`components`, `innovation_lags`, `shadow_force`): per-component estimate error rms, bias and relative error, innovation relative size and whiteness at time lags, and the force error of a shadow controller (archived `commanded_input` against `applied_input`); with `estimation` the `reference` is optional (workflow/native/src/ControlScoreEstimation.hpp). |
| `linear_analysis` | `analysis_method` `modal-gramian` and `klqe-assessment` (finite-window KLQE subset benchmark on a `trajectory`). |

## Queries and commands

- `krom-audit` action `control_map_path`: along `trajectory` (or `states`), the relative error,
  direction cosine and norm ratio of a stored map against the exact one, the truth-coupled frozen
  KLQR poles, and a pass/fail gate.
- `dynsyslab-cli krom-training-plan <request.json>`: the training-domain expansion plan (current
  domain from the K-ROM, union of domain and target requirements, probes, margin, r > 0 guard,
  fingerprint) and a solve job whose preflight screens every planned state
  (`workflow/KromTrainingPlan.hpp`). The RPE preset is `"preset": "rpe_targets"`.
- `krom_control_search`: `screening.schedule` (explicit candidates, `schedule_mode` prepend or
  only), `screening.closure_weight` (reorder by bilinear-closure quality),
  `screening.map_variants` (extra exact-map cases on controlled, prior-exact, arc-length or tube
  samples; default `bilinear_map`), and `resume_from` a ledger.

## Tracked presets

`examples/native/presets/rpe-clean-isolation-klqr.json` (the clean-isolation KLQR settings) and
`examples/native/presets/rpe-krom-training-plan.json` (the RPE training plan with the screening
preset r in [0.05, 15], |dr/dtau| <= 250, tau 0.25, DOP853).
