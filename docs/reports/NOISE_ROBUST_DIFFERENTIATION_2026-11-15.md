# Noise-robust derivatives for echo-only control (2026-11-15)

**Question.** The echo-only KLQE→KLQR needs time derivatives of a measured, noisy echo:
- ė as a coordinate;
- ψ̇ for the eigenfunction fit;
- ψ̇ for the B_ψ regression;
- ė online for the estimator.

Finite differences amplify noise. Which methods work, under which noise, and where in the pipeline?

**Tools added.** They are documented in `docs/architecture/NUMERICAL_DIFFERENTIATION.md`:
- 12 offline and 6 causal differentiators, plus Hampel despiking, in one core library;
- 8 noise models;
- `coordinate_transform` kind `measured_derivatives`;
- `derivative` options for coordinate SINDy and `input_map`;
- online `measurement_noise` and `measurement_derivative` for the KLQE;
- a `derivative_benchmark` post-processing step;
- an H-DMD `rank_rule: optimal_hard_threshold`.

## 1. Benchmark on the training echo

**Setup:**
- **Data:** the unforced r₀ = 0.3747 record, τ 0–10 (10,001 samples at Δt = 0.001 τ).
- **Reference:** the exact ė.
- **Noise:** levels are relative to the echo RMS, with 2 noise realizations each.
- **Tuning:** each method at its best cutoff from a grid of 2–256 cycles/τ.
- **Error metric:** relative RMS error of ė.
- **Archive:** `archives/derived/derivative-benchmark/…T060626Z.npz`.
- **Plots:** `viewer_outputs/derivative-benchmark-20261115/`.

The first six method columns are offline; the last five are causal (usable online).

| Noise | central | Savitzky–Golay | Kalman RTS | Tikhonov integral | TV-regularized | RTS + despike | causal backward | causal SG | causal Butterworth | causal Kalman | Levant |
|---|---|---|---|---|---|---|---|---|---|---|---|
| none | 0.18 | 0.04 | 0.0006 | 0.08 | 0.18 | 0.0006 | 3.8 | 6.2 | 8.9 | 0.52 | 20 |
| white 0.1% | 3.4 | 0.79 | 1.0 | 1.1 | 1.4 | 1.0 | 7.8 | 6.6 | 9.4 | 4.3 | 20 |
| white 1% | 34 | 4.0 | 3.7 | 4.0 | 3.9 | 3.7 | 68 | 20 | 22 | 18 | 49 |
| white 5% | 170 | 11 | 11 | 12 | 13 | 12 | 341 | 48 | 43 | 48 | 73 |
| white 20% | 681 | 30 | 26 | 27 | **16** | 28 | 1366 | 94 | 72 | 95 | 120 |
| red (AR1, 0.05 τ) 1% | 6.7 | 3.3 | 3.2 | 3.4 | **2.9** | 3.2 | 10 | 9.5 | 11 | 8.6 | 20 |
| pink 1% | 17 | 4.3 | 4.2 | 4.3 | **3.6** | 4.2 | 29 | 16 | 19 | 15 | 32 |
| drift 1%/τ | 1.1 | 0.70 | 0.72 | 0.72 | 0.96 | 0.72 | 4.1 | 6.3 | 9.0 | 2.0 | 20 |
| quantized 8 bit | 27 | 6.0 | 5.9 | 6.0 | 5.5 | 5.9 | 46 | 21 | 22 | 19 | 36 |
| quantized 12 bit | 2.1 | 0.56 | 0.62 | 0.65 | 1.4 | 0.62 | 5.6 | 6.3 | 9.1 | 3.3 | 19 |
| impulses (0.2%, 100% RMS) | 102 | 9.9 | 10 | 11 | 13 | **1.6** | 203 | 35 | 35 | 43 | 31 |
| hydrophone mix (white 2% + pink 1% + 12 bit) | 69 | 8.8 | 8.4 | 8.6 | 8.3 | 8.3 | 137 | 32 | 33 | 32 | 62 |

All values are ė error in %. The full table, which also includes Gaussian, Holoborodko, weak-form, zero-phase Butterworth, Whittaker, spectral and TV-denoise, is in `best_table_wide.md`. The per-noise ranking is in `best_table_ranked.md`.

**What works:**
- **Any tuned linear smoother.** Savitzky–Golay, the Kalman RTS smoother, the Whittaker smoothing spline, zero-phase Butterworth, spectral and the Tikhonov integral formulation all land within about 10% of each other once the cutoff is tuned. They beat central differences by 4× at 0.1% noise and by 10–25× at 1–20% noise. They are near-equivalent linear estimators, so the cutoff matters far more than the choice of method.
- **TV-regularized differentiation (Chartrand's integral form)** is the best choice for heavy white noise (16% against 26–30% at 20%) and for red and pink noise. It costs 10–100× more compute.
- **Hampel despiking before smoothing** is essential for impulsive noise: 10% falls to 1.6%. It costs nothing on the other noise types.
- **Rule of thumb for the cutoff:** about 128 cycles/τ at ≤0.1% noise, 64 at 1%, 32 at 5–20%. The echo's collapse spikes need the high bandwidth.

**What does not work:**
- **Plain finite differences:** 34% error at 1% noise.
- **TV denoising followed by finite differences:** it staircases the smooth echo.
- **Levant's sliding-mode differentiator:** it chatters at 20% even without noise. It needs a third-derivative bound, and the echo collapse makes that bound huge.
- **Online (causal) differentiation is inherently worse.** The best causal methods have 4–5× the error of the offline ones (18–22% against 3.7% at 1%), because any causal smoother lags. Lag-compensated, causal Butterworth reaches 13–14%. With little noise, the causal Kalman filter is best (0.5% clean, 4.3% at 0.1%).

## 2. Noisy training data through the whole pipeline (offline noise, clean online sensors)

**Pipeline.** Each training step was run on noisy measured coordinates:
1. **Coordinates:** noisy echo → robust (e, ė).
2. **Eigenvalues:** H-DMD.
3. **Eigenfunctions:** coordinate SINDy with measured derivatives.
4. **B_ψ:** regressed from the forced records' noisy measured coordinates.
5. **Target:** lifted offline.

The KLQE→KLQR then ran with the clean formula sensors, which isolates the effect of training noise. On clean data this measured-coordinate pipeline reproduces the earlier result: 0.76% / 0.38% / 2.82%.

**Finding 1: the spectrum broke first, not the derivatives.**
- With any noise, even 0.1%, H-DMD at full numerical rank picked a noise pair at about −5 ± 35–40i instead of −0.06 ± 5.94i.
- Fix: `rank_rule: "optimal_hard_threshold"`, the Gavish–Donoho singular-value threshold for unknown white noise.
- Result: λ = −0.081 ± 5.78i with the same value at every noise level tested, both 0.1–20% white and the mixtures.

**Finding 2: the eigenfunction fit then broke, through noise in its inputs.**
- The degree-4 polynomial ψ fitted to noisy (e, ė) drifts away from the clean ψ, with a normalized mismatch of 0.3–0.76.
- With degree-4 ψ, 1% noise lost the bubble on all three targets, whether ė came from RTS or central differences.
- Lowering the degree to 3 fixes most cases:

| Training noise (degree-3 ψ, ė by RTS unless noted) | mild | strong | hard |
|---|---|---|---|
| white 1% | 0.53% | 0.23% | 2.26% |
| white 1%, central differences | 0.72% | fails | fails |
| white 5% | 0.42% | 0.14% | 22% |
| white 20% | 0.34% | fails | fails |
| hydrophone mix | 0.42% | fails | fails |
| quantized 8 bit | 0.45% | 0.19% | 2.16% |
| impulses + despike | drive spikes (2100% raw, 0.66% low-passed; radius error 0.0002) | 0.23% | 30% low-passed |
| impulses without despike | B_ψ regression returns non-finite values | | |

- Degree 2 was best on mild and strong (0.36% and 0.14% at 1% noise) but lost the bubble on hard.
- Some degree-3 combinations still failed, including clean training on mild.
- Two remedies were tried: an estimator warm-up (control from 1 τ) and a 1000 kPa drive cap. Neither fixed the failures:
  - the warm-up alone changed nothing;
  - with the cap, the radius stayed on target for mild and some hard runs (radius error 0.0002–0.002), but the capped drive chatters (99–2600% raw drive error);
  - strong and hard trained on 20% noise or the mix still lost the bubble.

**Status: not robust yet.** Whether a noise-trained model controls well depends on the combination, not on the noise level.

## 3. Noisy sensors online

**Formula ė, noisy e and ė (the clean-trained measured K-ROM):**
- Raw drive error is large at any noise level. The high-gain loop passes sensor noise into the drive.
- The drive's low-frequency content is still right. After both drives are low-passed at 2 cycles/τ (above all target tones), the errors are:

| Online noise | measurement covariance | mild (low-passed drive / radius RMS error) | hard |
|---|---|---|---|
| white 0.1% | 1e-4 | 1.1% / 0.0003 | 4.8% / 0.003 |
| white 1% | 1e-4 | 3.0% / 0.0009 | 4.8% / 0.003 |
| white 5% | 1e-4 | 17% / 0.005 | 5.3% / 0.007 |

**Numerical ė online (sensor = noisy echo only, R = 1, the gain where numerical ė is stable):**

| Online noise | backward2 | causal Kalman (64) | causal SG (64) | causal Butterworth (64) |
|---|---|---|---|---|
| none | 7.9% | 7.9% | 7.9% | 8.1% |
| white 0.1% | 16% | 18% | 18% | 18% |
| white 1% | **542%, radius lost** | 16% (radius 0.0056) | 17% | 17% |

Values are low-passed drive error. The robust causal differentiators are what keep the loop alive at 1% noise. A plain backward difference loses the radius.

## Recommendations

1. **Offline derivatives:** use the Kalman RTS smoother (`kalman`, order 3) or Savitzky–Golay with a cutoff of about 64 cycles/τ, despiking first. Switch to `tvr` when the noise is heavy or colored and compute time allows. Never use raw differences.
2. **Spectrum from noisy data:** H-DMD with `rank_rule: optimal_hard_threshold`.
3. **Eigenfunctions from noisy data:** use a low polynomial degree (3). The next step is a fit that is robust to noise in its inputs, not just its targets. Candidates:
   - weak-form SINDy, which integrates the dynamics against test functions;
   - total least squares;
   - selecting the degree by validation on held-out data.
4. **Online:** causal Kalman for low noise; causal Butterworth or SG at a cutoff of 64 for ≥1% noise. Pair them with a KLQE measurement covariance matched to the noise (about 1e-4 in ψ units at 1%). A drive low-pass or a lower-gain R would remove the high-frequency drive noise.
5. **Open problem:** closed-loop robustness of noise-trained models and the start-up transient. The estimator starts from zero and the first drive spikes reach thousands of kPa.

## Sources

- F. van Breugel, J. N. Kutz, B. Brunton, "Numerical differentiation of noisy data: a unifying multi-objective optimization framework," and PyNumDiff ([PMC10530632](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10530632/), [GitHub](https://github.com/florisvb/PyNumDiff)).
- R. Chartrand, "Numerical differentiation of noisy, nonsmooth data" (TV-regularized differentiation, integral form).
- D. Messenger, D. Bortz, "Weak SINDy" ([PMC10795802](https://pmc.ncbi.nlm.nih.gov/articles/PMC10795802)).
- M. Mboup, C. Join, M. Fliess, algebraic differentiators ([arXiv 1209.1192](https://ar5iv.arxiv.org/html/1209.1192)).
- A. Levant, robust exact differentiators ([review](https://www.tau.ac.il/~levant/Levant,Livne,Yu-SM_Dfr_Review_Ifac2017.pdf)).
- P. Holoborodko, smooth noise-robust differentiators.
- Modulating-function method ([arXiv 1810.08929](https://arxiv.org/pdf/1810.08929)).
- L. Condat, "A direct algorithm for 1-D total variation denoising."
- M. Gavish, D. Donoho, "The optimal hard threshold for singular values is 4/√3."
