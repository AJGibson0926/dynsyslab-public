# KLQE → KLQR with a regressed B_ψ on the three RPE targets (2026-11-15)

**Goal.** Track the mild (1 + 0.05 sin, 100 kHz), strong (1.5 + 0.4 sin, 100 kHz) and hard four-tone (50/130/210/290 kHz) radius targets of the viscous RPE. The controller is KLQR acting on a KLQE estimate built only from echo sensors, and B_ψ is a constant regressed from forced data.

**Method.** Training follows the data-role rules: eigenfunctions are fitted on unforced data only, and B_ψ is regressed on three forced records (rpe-forced-bpsi a005p05, a02p05 and a02p2) with the `input_map` step (`representation: constant`, `fit: exact_map`). Control is sampled at 0.001 τ, as KLQE feedback requires, and runs for 300 τ. Scores are on τ 100–300 against the exact inverse drive.

## Results

| Target | Estimator → controller | Eigenfunctions | Drive error | Correlation | Radius RMSE | Chatter |
|---|---|---|---|---|---|---|
| mild | **KLQE (ψ measured) → KLQR**, R 1e-2, process 100, measurement 1e-6 | (e, e'), r₀ = 0.3747 early fit | **0.7%** | 1.000 | 0.0002 | none |
| | full-state KLQR (reference) | same | 0.6% | 1.000 | 0.0002 | none |
| | KLQE (echo-pair readout) → KLQR, R 1e-2, process 1, measurement 1e-4 | (r, ṙ), r₀ = 0.3747 early fit | 4.7% | 0.999 | 0.0022 (r̂ error 0.0024) | none |
| strong, 10 τ ramp | **KLQE (ψ measured) → KLQR**, R 1e-2, process 100, measurement 1e-6 | (e, e'), r₀ = 0.3747 | **0.3%** | 1.000 | 0.0022 | none |
| | full-state KLQR (reference) | same | 0.3% | 1.000 | 0.0022 | none |
| hard | **KLQE (ψ measured) → KLQR**, R 1e-2, process 10, measurement 1e-5 | (e, e'), r₀ = 0.3747 | **3.1%** | 1.000 | 0.0018 | none |
| | full-state KLQR (reference) | same | 1.9% | 1.000 | 0.0013 | none |

The strong KLQE run also works without the ramp (0.3%). Its estimator start-up then reaches ±1.7 MPa for about 0.05 τ.

## What made it work

1. **Measure the eigenfunctions, not a linear readout.** The KLQE sensors are the echo and its time derivative (e, ė). The measurement is ψ evaluated on those sensor values, in split real form: `measurement_transform: "eigenfunctions"`. The measurement map is then the identity (`lifted_split` reconstruction), so it is exact at every amplitude. The filter fuses that measurement with the model ż = Λz + B_ψ u.

2. **(e, e') coordinates.** A constant regressed B_ψ in (e, e') coordinates suffices for all three targets: full-state 0.3–1.9%, KLQE 0.3–3.1%. In (r, ṙ) coordinates the same regression reproduces the analytic ∇ψ·g at equilibrium (1.096i against 1.027i). Full-state control with it reaches 0.9% (mild), 1.7% (strong) and 14.5% (hard). That map varies about 4× between r = 1 and r = 1.5.

3. **Trust the measurement on the hard target.** The best settings are process covariance 10 and measurement covariance 1e-5. At process 100 with measurement 1e-6 and R = 1e-2 the hard loop chatters.

## What did not work, and why

- **Joint H-DMD [r, ṙ, e] K-ROMs with a constant regressed B_ψ**, both the late window (5 modes) and an early r₀ = 0.3747 window. Even full-state, the drive carries a ±1–3 kPa sample-to-sample jitter. The best mild KLQE run gave 56% drive error. Strong and hard gave 70–100%, or diverged.
- **(r, ṙ) eigenfunctions with the K-ROM's fitted linear echo readout.** Three problems:
  - **ψ is not zero at equilibrium.** ψ(1, 0) = 0.108, which the linear model ż = Λz cannot hold, so the estimator is biased. Re-centring ψ at equilibrium removes this.
  - **Global readouts.** Readouts fitted over the whole excursion are about 4× off near equilibrium, so the loop under-drives at about 37% amplitude. Exact local readouts at equilibrium fix this: C_x = J⁻¹ and C_y = (∂e/∂x)·J⁻¹, with S = ρR₀³/τ_c² = 0.10998 Pa·m.
  - **Velocity is nearly unobservable from e alone near equilibrium.** e ≈ S(−39.5 δr − 0.12 ṙ). Adding ė, with readout (∂e/∂x)·Df·J⁻¹, makes ṙ observable.

  With all three fixes the mild target reaches 4.7%. Strong and hard still diverge, because a local linear readout cannot represent large excursions. That is why the eigenfunction measurement was needed.

- **Input feedthrough.** A real hydrophone echo contains the drive directly (∂e/∂u = −p₀R₀ r). The simulated `far_field_acoustic_echo` sensor, however, is a function of state only. The new optional `measurement_feedthrough` is therefore correct physics, but must stay off for this sensor.

## Code added (uncommitted)

- **Estimator measurement transform.** `measurement_transform: "eigenfunctions"` in `EstimatorFactory` and `EstimatedControl`. The sensors must be the K-ROM's coordinate observables, the measurement reconstruction must be p = n, and y = `evaluateCoordinatesSplit(sensor values)`.
- **Measurement feedthrough.** `measurement_feedthrough` (D, p × m) and `measurement_feedthrough_scale_output` are wired through `KoopmanLQEObserver::setMeasurementFeedthrough`. The observer uses y = y₀ + C_y z + D u.
- **Tests.** The estimation, feedback, estimation-runtime, control-factory and pipeline-job tests all pass.

## Experiment archives (hand-built, prefixed `experiment-`)

| Archive | Purpose |
|---|---|
| `koopman-krom-extract/experiment-ee-ic0p37-split-krom.npz`, `experiment-ee-late-split-krom.npz` | (e, e') K-ROMs plus the `lifted_split` identity readout |
| `koopman-krom-extract/experiment-rr-ic0p37-recentred*.npz` | (r, ṙ) K-ROMs: re-centred at x*, anchored or local readouts, echo-pair readout |

Their B_ψ maps were regressed by the normal `input_map` step. The readout and re-centring edits were made with scratch scripts; making them native K-ROM options is the natural next step.

## Plots

`viewer_outputs/klqe-klqr-regressed-bpsi-20261115/{mild,strong,hard}/` contains:

- drive and radius start-up and steady-state time series;
- log10 drive and radius spectra over τ 100.8–300;
- a phase portrait of each KLQE run;
- KLQE estimate against measured ψ (mild (r, ṙ) run: r̂ against r).

Colours: KLQE → KLQR (e, e') in blue, full-state KLQR in red, KLQE (r, ṙ) in gold (dashed), exact or target in black.
