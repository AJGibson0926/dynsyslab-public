# Estimator module (KLQE) beside the Controller module

The native workflow has two plant-facing modules with matching shapes:

| | Controller module | Estimator module |
|---|---|---|
| Factory | `ControlFactory.hpp` → `buildController` | `EstimatorFactory.hpp` → `buildEstimator` |
| Job block | `"controller": {kind, settings, target}` | `"estimator": {kind, krom, measurement, …}` |
| Runtime | `ScalarControlPolicy` (input law) | `KoopmanLQEObserver` (lifted estimate) |
| Uses the KROM for | `A`, `B_psi` (`control_representation`) | `A`, `C_y`, `C_x`, optional `B_psi` |
| Archive traces | `control_values`, targets | `estimator_time`, `lifted_estimate`, `state_estimate`, `measurement_values`, `measurement_estimate`, `innovation`, `estimator_applied_input`, `estimator_gain`, `observer_poles` |

Core numerics live in `native/include/dynsyslab/core/Estimation*.hpp`:
- `designKoopmanLQE*`: continuous LQE (CARE dual), detectability and pole checks.
- `runKoopmanLQEFromLifted`: offline exact-ZOH filtering of a measurement record.
- `KoopmanLQEObserver`: the online observer.
- Observability and Gramians: `Observability.hpp`, `ModalObservability.hpp`, `Pbh.hpp`.

The workflow factory only binds those numerics to a solve job and its simulated
sensors.

## Solve-job API

```json
"estimator": {
  "kind": "klqe",
  "krom": "<KROM archive>",
  "measurement": "far_field_acoustic_echo",
  "measurement_observables": ["far_field_acoustic_echo"],
  "output": "state",
  "process_covariance": 1.0,
  "measurement_covariance": 1e-4,
  "initial_estimate": "zero",
  "control_representation": "none",
  "sample_time": 0.001
}
```

- `measurement` names the KROM reconstruction `C_y`.
- `measurement_observables` are the simulated sensors. The default is the
  measurement name.
- `output` names the decoded reconstruction `C_x`.
- `process_covariance` and `measurement_covariance` each accept a scalar, a
  diagonal or a full square matrix.
- `process_covariance` may instead be shaped from the K-ROM's own prediction error:
  `{"from_records": [solution archives], "residual_horizon": h, "window": [t0, t1],
  "scale": 1, "floor": 1e-6, "add": 0}` gives `W = scale·E[r rᵀ]/h + add·I` with the
  one-step residual `r = ψ(x(t+h)) − e^{Ah}ψ(x(t)) − Γ(h)u(t)` on the records' states
  and inputs (zero-order hold); eigenvalues below `floor·max` are raised to keep W
  positive definite. The estimator then distrusts the K-ROM dynamics only where they
  mispredict. `process_covariance_input_weight: q²` adds `q² B Bᵀ` (constant map):
  model error entering through the input. Both are recorded as
  `estimator_process_covariance` in the archive evidence (with the relative residual).
  The records are offline data; use the operating orbit, never the run being scored.
- `initial_estimate`: use `"zero"`, or give `initial_lifted`.
- `control_representation`: `"none"`, or the name of a KROM input map.
- `balanced_truncation`: optional reduced-order gain lifted to the full KROM
  (see [LINEAR_SYSTEMS_TOOLKIT.md](LINEAR_SYSTEMS_TOOLKIT.md)). The observer still
  estimates every lifted coordinate.

Composition with the controller (`Solve.cpp`):

- **No controller.** The estimator observes the free plant. The archive stays
  unforced: no `control_values` are written.
- **Prescribed input or state-feedback controller.** The input is held on the
  estimator clock, and the estimator is told the applied input.
- **`controller.kind = "klqr"` with `settings.feedback = "estimate"`.** This is
  output feedback, KLQE → KLQR. The controller acts on the estimator's lifted
  estimate. Both must use one KROM (`estimator.krom` equals the controller
  `krom`, or is omitted). The legacy form `controller.settings.estimator` is
  still accepted; give the estimator in one place only.
- **Any controller with `settings.feedback = "estimated_state"`.** Separation of
  estimation and control: the estimator (its own KROM, its own sensors) decodes
  x̂ = x₀ + C_x ẑ through `output` (which must give the full state), and the
  controller's ordinary state-feedback law is evaluated at x̂ instead of the
  plant state. A KLQR then evaluates its own eigenfunctions ψ_ctrl(x̂), so the
  estimator and the controller may use different KROMs and coordinates. The
  estimator observes first each sample; the controller clock must equal the
  estimator clock. Errors in x̂ propagate through ψ_ctrl: a steep ψ_ctrl
  amplifies them (measured in the 2026-11-17 RPE study, where a biased ṙ̂ was
  amplified about 4× into ψ). Use `active_time_window` to let the estimator
  converge before control starts.
- **Shadow KLQE → KLQR: `settings.applied_controller`** (with `feedback =
  "estimate"`). Any non-estimating controller spec (waveform, full-state KLQR,
  ...) drives the plant on the true state; the KLQE → KLQR runs beside it,
  predicts with that applied drive, and its own force is archived as
  `commanded_input` but never applied (`applied_input` is the drive the plant
  got). This separates estimation quality from the feedback loop. Use the same
  KLQR on the true state as `applied_controller` to make the force error
  exactly −K(ẑ − ψ(x)): KLQR has no feedforward, so an exact inverse drive is
  not what it would command even with a perfect state.

The estimator reads the plant **only** through its sensors. It never evaluates
an eigenfunction at the true state.

## Offline estimation (post-processing `klqe`)

`klqe` filters a stored measurement record with a KROM: `sources` are solution
archives, and the measurement is evaluated from the archived states (or read
from an array). Settings:

- `krom`, `measurement`, `output_reconstruction`;
- `process_variance` / `process_diagonal`;
- `measurement_variance` / `measurement_diagonal`;
- `initial_estimate`: `zero`, `truth` (ψ at the true first state; benchmark
  only), or the legacy default;
- `time_window` and `sample_stride` to select rows.
- `balanced_truncation`: as in the solve-job API.

The output archive holds the estimate, the truth, the innovation, the
observability analysis and the gain.

## Building a KROM for KLQE: x = C_x Ψ and y = C_y Ψ from one H-DMD

This mirrors the MATLAB suite (`dynsyslab (5-15-25)`, `KLQE.m`):

- One Hankel DMD over the observables [r, ṙ, e].
- `extract(method="observable", threshold, caps)` selects the modes.
- `Kf = lqr(A', C_y', Vd, Vn)'` gives the gain.

1. **Joint H-DMD on unforced data.** Use
   `hdmd {"observables": ["radius", "radial_velocity", "far_field_acoustic_echo"],
   "window_alignment": "end", "reconstruction_observables": ["state"], "time_window": [...]}`.
   Every channel gets its own `C` row (its Koopman modes). `time_window`
   restricts the fit to part of a free record.
2. **Mode selection.**
   `krom_extract {"mode_selection": {"rule": "observable_energy",
   "reconstructions": [...], "thresholds": [...], "caps": [...],
   "max_real_part": 0, "exclude_constant": true}}` ranks modes by
   `|C_o,j| · rms(ψ_j)` per observable. It keeps only non-growing modes and
   completes conjugate pairs.
3. **Readouts with an intercept.**
   `"refit_reconstructions": {"observables": ["state", "far_field_acoustic_echo"], "intercept": true}`
   refits `C` on the selected modes' training samples as `h = h0 + C Ψ`. The
   intercept is stored as `reconstruction_offset__NAME` and is used by every
   estimator path (offline, online and output feedback).

The intercept is essential. The equilibrium radius appears in H-DMD as a
nearly-constant mode with `λ ≈ 0`. The echo has zero mean, so that mode is
unobservable from it. An observer started at zero therefore never recovers it,
and the radius estimate goes to 0. That was the failure in the
5 November joint-H-DMD KLQE attempt. Dropping the constant mode and carrying
`x0` as an intercept removes it.

## Results (RPE, viscous, 2026-11-14)

The KROM is the late unforced baseline (fine record, τ 50–80, dt 1e-3, 2,000
delays). Five modes were selected: the λ ≈ −0.060 ± 6.283i pair, the
−0.121 ± 12.566i harmonic pair, and one real decaying mode. The intercepts were
x0 = (0.99995, 0.0002).

The observer runs on the echo alone with a zero initial estimate. Q = 1 and
R = 1e-4. The ICs below are unforced trajectories; 1.25, 1.65, 2.15 and 2.5 were
never used for fitting.

| IC (window) | r RMSE / r amplitude | ṙ RMSE / ṙ amplitude | echo relative error |
|---|---|---|---|
| 2.33 (τ 50–80, training) | 0.0015 / 0.026 | 0.0013 / 0.16 | 0.65% |
| 1.25 (τ 20–80) | 0.0006 / 0.015 | 0.0006 / 0.093 | 0.55% |
| 1.65 (τ 30–80) | 0.0020 / 0.025 | 0.0020 / 0.15 | 1.0% |
| 2.15 (τ 40–80) | 0.0025 / 0.032 | 0.0036 / 0.20 | 1.4% |
| 2.5 (τ 45–80) | 0.0046 / 0.040 | 0.0063 / 0.25 | 2.0% |

These errors are scored after the first 5 τ of each window.

**Limits.** This is a near-equilibrium estimator. Including the large-amplitude
collapses (for example IC 1.65 over τ 0–80, or IC 2.5 from τ 0) gives errors of
30–100% of the oscillation amplitude. A linear modal model of those collapses,
fitted on the whole decay, has growing modes and fails selection.
