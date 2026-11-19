# Overnight 2026-11-13/14: ξ eigenfunctions, switchable B_ψ, data roles, KLQE, KLQE→KLQR

All runs use the viscous RPE (air-water, R0 = 10 µm, p0 = 100 kPa). Control targets the
20 kPa, 100 kHz periodic orbit; its driving sine is withheld from the controller and used only
for scoring.

## 1. ξ (causal echo-history) eigenfunctions for the 20 kPa target

The approach is coherent. Along the unforced flow, ξ (the Hankel-SVD coordinates of the echo
history) is a function of the state, so ψ(ξ) can be an eigenfunction. The KROM keeps the frozen
basis, the target ξ_ref comes from the target orbit's own echo history, and B_ψ is either
`history_exact` or regressed. It runs in the existing 20 kPa pipeline with only parameters
changed.

What it needed:

- **Late-regime spectrum.** The 25 τ baseline is dominated by collapses. H-DMD on all of it
  gives λ = −0.21 ± 4.82i, but the small-amplitude regime (where the target lives) rings at
  −0.06 ± 6.28i. H-DMD gained a `time_window` setting (SINDy already had one), so both can fit
  τ 40–80 of the fine baseline. Free closure of ψ(ξ) dropped from 52% to 1–2%.
- **Linear (degree-1) ψ(ξ).** Higher degrees blow up on the target's ξ, which lies off the free
  manifold.
- **A regressed B_ψ.** The analytic `history_exact` map fails in every configuration
  (|B| drops to ~1e-4 of its median once per cycle).

## 2. B_ψ is one switchable pipeline step (`input_map`)

- **Representations:** `exact` | `history_exact` | `constant` | `bilinear`.
- **Fits:** `held_input` (ξ KROMs, measurement only, optional affine drift nuisance),
  `controlled_residual`, `multistep_output` and `exact_map` (state KROMs).
- **Pipeline:** `20kpa.json` now has `B_psi` and `forced_training` parameters. The controller
  reads `${params.B_psi.representation}`; nested parameter fields are new. KLQR gained
  `control_limit`.
- **Results (ξ: 50-delay, 2-latent, degree 1; R = 1000; ±50 kPa limit):** every regressed form
  trained on forced data of 20 kPa or less recovers **18.4 kPa, 1.2–1.4° phase, correlation
  0.9997** with radius RMSE 0.003. The 100 kPa forced data gives the wrong sign. Full table:
  `docs/architecture/INPUT_MAP_B_PSI.md`.

## 3. Training-data roles are enforced

Eigenfunction/spectrum methods (H-DMD, SINDy, EDMD, ResDMD, mpEDMD, HAVOK, …) refuse forced
sources. DMDc and EDMDc accept them. B_ψ regressions require them. Coordinate changes now carry
`forcing_type: forced` forward. See `docs/architecture/TRAINING_DATA_ROLES.md`, which includes
the reasons and references for each method.

## 4. KLQE: estimating Ψ and x from the echo

- **Recipe:** the MATLAB suite (`dynsyslab (5-15-25)`, `KLQE.m`) through native code:
  1. joint H-DMD on [r, ṙ, e] of the unforced baseline;
  2. `observable_energy` mode selection (non-growing modes only);
  3. C_x and C_y refit **with an intercept**;
  4. a continuous Kalman gain.
- **Failure fixed:** Codex's 5 November failure was the equilibrium radius sitting in a λ≈0 mode
  that the echo cannot observe. Reconstructions now carry an offset through archive, extraction,
  design, and offline and online observers.
- **Results:** with a zero initial estimate on unforced held-out ICs, r error is 2–12% of the
  oscillation amplitude and echo error 0.6–2% (`docs/architecture/ESTIMATOR_MODULE.md`). The
  estimator is near-equilibrium only.
- **Estimator module:** `EstimatorFactory` sits beside `ControlFactory`. A solve job takes a
  top-level `"estimator"` block. It runs estimator-only, beside any input policy, or as KLQE →
  KLQR output feedback (`controller.settings.feedback = "estimate"`).

## 5. KLQE → KLQR with regressed B_ψ

Both runs use the same 5-mode KROM and B_ψ regressed from forced data; the target is the 20 kPa
orbit over 100 τ.

| Feedback | Best setting | Fundamental (true 20 kPa) | Phase | Correlation | Radius RMSE | Estimator r error |
|---|---|---|---|---|---|---|
| full state (Ψ(x)) | constant, controlled-residual; R = 5000; ±25 kPa | 19.8 kPa | 0.9° | 0.89 | 0.0008 | — |
| echo only (KLQE) | constant, controlled-residual; ±25 kPa | 28 kPa | 8° | 0.88 | 0.030 | 0.036 |
| echo only (KLQE) | constant, exact-map fit; ±25 kPa | 15–16 kPa | 24° | 0.83–0.85 | 0.016 | 0.008 |

It does not fail catastrophically once the actuator is limited to about the target amplitude.
At ±50 kPa the loop drives the bubble out of the estimator's near-linear domain (estimate error
0.09 in r), and it overdrives at 53 kPa. See section 6 for the lower-gain grid.

## 6. Second grid (lower controller gain, tighter limits)

Settings: process covariance 100, measurement covariance 1e-6, R = 5000.

| B_psi fit | Limit | Q | Fundamental | Phase | Correlation | Radius RMSE | Estimator r error |
|---|---|---|---|---|---|---|---|
| controlled-residual | ±20 kPa | 5000 | **22.6 kPa** | **−15°** | **0.87** | **0.018** | 0.020 |
| controlled-residual | ±20 kPa | 50 | 22.2 kPa | −29° | 0.85 | 0.023 | 0.018 |
| exact-map | ±20 kPa | 5000 | 15.1 kPa | +25° | 0.85 | 0.016 | 0.008 |
| exact-map | ±20 kPa | 50 | 13.0 kPa | +39° | 0.75 | 0.023 | 0.007 |

- **Gain does not help.** Lowering the controller gain (Q) does not improve tracking.
- **The two B_psi fits trade off.** The exact-map fit estimates better but under-drives. The
  controlled-residual fit drives harder with larger estimate error.
- **Overall:** the best echo-only loop recovers the withheld 20 kPa driver within about 13% in
  amplitude, 15° in phase and correlation 0.87, with the bubble tracking the orbit to
  RMSE 0.018. This is far from the full-state result but not catastrophic.
- **Bilinear output feedback:** most runs hit the 240 s timeout; the bilinear maps are kept for
  full-state use.

## 7. Example pipeline

`examples/native/pipelines/rpe-klqe/echo-to-state.json` runs unforced training → joint H-DMD →
observable_energy K-ROM with intercepts → offline KLQE on held-out IC 1.65 (r RMSE 0.0014
over τ 35–80) → online estimator-only solve from r = 1.4. The online r RMSE is 0.12 during the
collapses (τ 0–10), 0.022 over τ 10–30 and 0.001 after τ 30. The run takes about 5 minutes; the
H-DMD dominates.

## Tests

The native suite passes: 122/122. The new checks are in `pipeline_job_tests` (nested
parameters, data roles, intercept round trip, estimator-only solve). Three existing tests that
trained eigenfunction methods on forced data were changed:

- `native_solve_tests` now uses an unforced solve.
- `observable_pipeline_tests` now asserts refusal.
- `hankel_observables_tests` exposed the bug below and needed no change.

Coordinate changes of sources without provenance wrote `forcing_type: "unknown"`, which was then
classified as forced. "unknown" now means no information.

## Code touched

- **Workflow:** `InputMap.cpp` (new), `EstimatorFactory.{hpp,cpp}` (new), `Solve.cpp`,
  `Plan.cpp`, `Pipeline.cpp`, `Postprocess.cpp`, `KromPostprocess.cpp`,
  `KromDatasetSupport.{hpp,cpp}`, `HankelAnalysis.{hpp,cpp}`, `HdmdPostprocess.cpp`,
  `CoordinateTransformPostprocess.cpp`, `ControlFactory.cpp`, `EstimationPostprocess.cpp`.
- **Core:** `Reconstruction.hpp`, `KoopmanArchive.cpp`, `KoopmanModel.cpp`,
  `Estimation.{hpp,cpp}`, `EstimationRuntime.{hpp,cpp}`.
- **Tests:** `pipeline_job_tests.cpp` covers nested parameters, data-role refusal/acceptance,
  the intercept round trip and estimator-only solves.
