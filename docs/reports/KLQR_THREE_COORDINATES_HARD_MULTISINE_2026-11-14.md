# KLQR on a hard four-tone radius target in (r, ṙ), (e, e') and ξ eigenfunction coordinates (2026-11-14)

This follows `KLQR_THREE_COORDINATES_SINE_TARGETS_2026-11-14.md`. It uses the same pipeline (`examples/native/pipelines/rpe-echo-jet-klqr-rediscovery/20kpa.json`) and the same settings: early-time training, Q = 5000, R = 1e-4, continuous evaluation, no cap and no ramp. Only the target changes.

## Units

τ = 1/f0 = 3.0124 µs, with f0 = 331.97 kHz for R0 = 10 µm, p0 = 100 kPa, pV = 2.33 kPa, ρ = 998, σ = 0.0725 and μ = 0.001. The drive coefficient is Eu = p0 τ²/(ρ R0²) = 9.0926, and one drive unit is p0 = 100 kPa.

## Physical limits on large targets

Exact inverse drive u = (r̈ − f(r, ṙ)) / (−9.0926 / r):

| Candidate | Tones (kHz) | r range | Exact drive (kPa) | Far-field minimum p0 + u |
|---|---|---|---|---|
| offset 4, amplitudes [3/2, 1], ω [4, 3] rad/τ | 211, 158 | 1.56–6.44 | −1681 … +2197 | **−1581 kPa** |
| offset 3.5, amplitudes [3/2, 1/3, 1], ω [½, 1, 2]·2π rad/τ | 166, 332, 664 | 1.32–5.68 | −7753 … +9954 | **−7653 kPa** |

Large-offset targets like these need megapascal liquid tension, so they are not physical. They also run r to 5.7–6.4, outside any training decay used here. The static constraints are tight: the liquid at the wall falls below pV once r ≳ 1.95, and the far field needs tension for a static r ≳ 1.9.

## Target used (designed here)

The target was chosen by a random search that maximised excursion × nonlinearity × drive RMS. It is subject to three constraints:

- far-field absolute pressure ≥ pV + 3 kPa;
- wall liquid pressure p_G + p_V − 2σ/R − 4μṘ/R ≥ pV + 1 kPa;
- r ≥ 0.6.

The target is then

r_T(τ) = 1.12 + 0.16 sin(0.946368 τ + 1.6) + 0.35 sin(2.460557 τ + 5.2) + 0.07 sin(3.974746 τ + 0.8) + 0.05 sin(5.488935 τ + 1.3)

The four tones are **50, 130, 210 and 290 kHz**, with ṙ the exact derivative. The common period is 100 µs (33.2 τ).

- **Radius:** r ∈ [0.613, 1.627] and ṙ ∈ [−1.28, 1.28]. This is well inside the training range.
- **Exact inverse drive:** −93.0 … +656.4 kPa, mean +28.1 kPa, AC RMS 151.5 kPa.
- **Physicality:** the far-field absolute minimum is 7.0 kPa and the wall-liquid minimum is 9.4 kPa. Both are above pV = 2.33 kPa.
- **Nonlinearity:** 40% of the drive variance lies outside the four target tones. The drive's line at 290 kHz is only 2 kPa, but it has 62 kPa at 80 kHz (130 − 50) and 76 kPa at 260 kHz (2 × 130). The full set of lines is in the spectra.

The previous "strong" target spanned −95 … −52 kPa with one tone, so this one is about 17× harder in drive span.

## Attempts and where B_ψ comes from

All three attempts use the same eigenfunction training. The fine baseline record (`rpe-viscous-echo-long-training-five-ics-trajectory-0__20261002T001350Z.npz`, **unforced**) is fitted with H-DMD and SINDy on τ 0–25, both with stride 10, and eigenfunction SINDy also uses τ 0–25. Q = 5000.

| Attempt | Coordinates | ψ fit | B_ψ | Evaluation |
|---|---|---|---|---|
| 1 | (r, ṙ) identity | SINDy degree 4, γ `system_rhs` | `exact`: ∇ψ · g(x), g = (0, −9.0926/r). Analytic, **no data** | continuous, R = 1e-4 |
| 2 | (e, e') `echo_jet2` | SINDy degree 4, γ `system_rhs` | `exact`: ∇_(e,e')ψ · J_(e,e')(x) · g(x). Analytic chain rule, **no data** | continuous, R = 1e-4 |
| 3 | ξ: echo_jet2 → history (decimation 10, 50 delays, 2 latents) | SINDy degree 1, γ `finite_difference` | `constant`, held-input least squares on the pipeline's **forced** record: 60 τ from (1, 0), random pulses ±20 kPa, 2 τ, seed 15, `fit_start` 5 | sampled, 0.001 τ |

## Results (scored on τ 100–300)

| Variant | Attempt | Drive error vs exact | Correlation | Radius RMSE (τ ≥ 100) | Radius RMSE (τ ≤ 10) | Peak command | Steady peak (exact 656) |
|---|---|---|---|---|---|---|---|
| no ramp, no cap | (r, ṙ) | **1.05%** | 1.0000 | 0.0013 | 0.014 | **12.3 MPa** (−12 344 kPa at τ≈0) | 652 kPa |
| | (e, e') | **0.62%** | 1.0000 | 0.0007 | 0.007 | **16.8 MPa** (−16 844 kPa at τ≈0) | 652 kPa |
| | ξ, R = 1000 | 77% | 0.946 | 0.217 | 0.218 | 133 kPa | 132 kPa |
| 10 τ ramp, no cap | (r, ṙ) | 1.05% | 1.0000 | 0.0013 | 0.078 | 652 kPa (dips to −124 kPa at τ≈3.5) | 652 kPa |
| | (e, e') | 0.62% | 1.0000 | 0.0007 | 0.0004 | 652 kPa | 652 kPa |
| | ξ, R = 1000 | 77% | 0.946 | 0.217 | 0.131 | 137 kPa | 132 kPa |
| 10 τ ramp + cap [−95 kPa, +1 MPa] | (r, ṙ), (e, e'), ξ | identical to ramp, no cap | | | | (r, ṙ) touches −95 kPa briefly at τ≈3.5 | |
| cap [−95 kPa, +1 MPa], no ramp | (r, ṙ) | 1.05% | 1.0000 | 0.0013 | 0.022 | 652 kPa (cap at start-up, 0.04% of samples) | 652 kPa |
| | (e, e') | 0.62% | 1.0000 | 0.0007 | 0.013 | 652 kPa (cap at start-up, 0.05% of samples) | 652 kPa |

The cap's −95 kPa floor keeps the far-field absolute pressure at or above 5 kPa > pV.

ξ R sweep (no ramp, no cap, sampled):

| R | Outcome |
|---|---|
| 100, 200, 300 | diverges ("left its validity domain") |
| 500 | 70% error, corr 0.84, r 0.62–1.35 |
| 1000 | 77% error, corr 0.95 |
| 10 000 | 95% error, corr 0.93 |

Drive line amplitudes on τ 100.8–300 (6 common periods, 1.67 kHz bins), in kPa:

| Line (kHz) | DC | 50 | 80 | 130 | 160 | 210 | 260 | 290 |
|---|---|---|---|---|---|---|---|---|
| exact | 28.4 | 91.4 | 61.9 | 135.3 | 5.0 | 31.9 | 76.5 | 2.0 |
| (r, ṙ) | 28.1 | 90.8 | 61.4 | 134.7 | 5.0 | 31.5 | 76.0 | 2.1 |
| (e, e') | 28.1 | 90.9 | 61.5 | 134.7 | 4.9 | 31.7 | 76.0 | 2.0 |
| ξ | 1.1 | 23.7 | 11.6 | 43.1 | 1.4 | 8.6 | 11.5 | 3.6 |

## Findings

- **(r, ṙ) and (e, e') with exact B_ψ solve the hard target.** In steady state the controller pressure matches the exact inverse drive to within 1% across all lines, including the intermodulation lines at 80 and 260 kHz that are not in the target. The residual is a ~0.3 kPa mean offset (27.7 vs 28.0 kPa).
- **Start-up without a ramp is extreme.** The target jumps at τ = 0 from (1, 0) to (1.069, 0.666), a velocity step, and the first commands are −12 to −17 MPa. Either the 10 τ ramp or the physical cap [−95 kPa, +1 MPa] alone removes the spike without changing the steady state.
- **(e, e') is the cleaner of the two.** During the ramp, ψ(r, ṙ) (early-trained) commands a −124 kPa dip and a ~250 kPa bump near τ ≈ 3.5 while r is still near equilibrium. This is the same weakness of early-trained ψ(r, ṙ) near r = 1 seen on the mild sine target. (e, e') follows the exact ramp drive throughout.
- **ξ with a regressed constant B_ψ fails on this target.** It produces the right waveform shape (correlation 0.95) but only about ¼ of the amplitude, misses the +28 kPa mean, and tracks r over 0.80–1.17 instead of 0.61–1.63. A lower R diverges. B_ψ is one constant matrix regressed from ±20 kPa pulses about R0, while this target needs −93 … +656 kPa, so the constant map cannot represent the state-dependent g(x) = −9.09/r.

## Plots

`viewer_outputs/klqr-coords-hard-multisine-20261114/{no-ramp-no-cap, ramp10-cap}/` contains:

- drive and radius start-up (τ 0–15, plus drive from τ 1);
- drive and radius steady state (τ 266.8–300, one common period);
- log10 drive and radius spectra over 0–1000 kHz and 0–400 kHz.

Colours: exact or target in black, (r, ṙ) in red, (e, e') in blue, ξ in gold (dashed). Data files (`time_series.npz`, `spectra*.npz`) are next to the PNGs.


## Training-IC dependence: chatter on the mild target and all nine reruns

### Chatter with the r₀ = 2.333 baseline

With early training on the baseline decay from r₀ = 2.333, the (r, ṙ) controller chatters on the mild target (1 + 0.05 sin, 100 kHz). It fires ±2.8 MPa bursts once per period, its slope reverses 15× more often than the exact drive's, and r settles near 0.8.

The cause is the eigenfunction fit. On the violent collapse from 2.333 (to r = 0.23), the degree-4 ψ(r, ṙ) folds at r ≈ 0.89: det ∂(Re ψ, Im ψ)/∂(r, ṙ) changes sign, so ψ is not one-to-one near equilibrium. For example, ψ(0.85, 0) ≈ ψ(1.05, 0) ≈ 0.03. KLQR can then null ψ(x) − ψ(x_ref) on the wrong side of the fold, and the gain blows up at each crossing. The strong and hard targets stay clear of the fold and track exactly.

### The r₀ = 0.3747 training IC

Training on the decay from r₀ = 0.3747 at rest instead (τ 0–25, dt 0.001, stride 1; pipeline override `training.initial_conditions.values=[[0.374736182055959, 0]]`) changes the fit:

- **Eigenvalue:** λ = −0.0161 ± 5.9426i, against −0.211 ± 4.818i from the 2.333 baseline.
- **No fold:** det stays between 0.046 and 0.14 over r 0.7–1.15 and |ṙ| < 0.6.
- **Linearisation reference:** at r* = 1 the eigenvalues are −0.0604 ± 6.2829i. Neither training recovers this; H-DMD on a finite-amplitude transient weights the slower nonlinear frequency.

All nine runs on this training IC use Q = 5000, R = 1e-4 and continuous evaluation, with no cap. For ξ: R = 1000, sampled, history decimation 1 on the 0.001 τ record (the same 0.001 τ history clock as before) and a constant B_ψ from the forced record. Scores are on τ 100–300.

| Target | Controller | Drive error | Correlation | Radius RMSE | Peak command (start-up) |
|---|---|---|---|---|---|
| mild | (r, ṙ) | **0.09%** (was: chatter) | 1.0000 | 0.0000 | 10.3 MPa |
| | (e, e') | 0.05% | 1.0000 | 0.0000 | 9.1 MPa |
| | ξ | 45% (was 68%) | 0.9995 | 0.016 | 13 kPa |
| strong, no ramp | (r, ṙ) | start-up fails: Riccati solution not real | | | |
| | (e, e') | 0.03% | 1.0000 | 0.0002 | 301 MPa |
| strong, 10 τ ramp | (r, ṙ) | **0.01%** | 1.0000 | 0.0001 | 95 kPa |
| | (e, e') | 0.03% | 1.0000 | 0.0002 | 95 kPa |
| | ξ | 38% (was 58%) | 0.52 | 0.30 | 84 kPa |
| hard | (r, ṙ) | **0.20%** (was 1.05%) | 1.0000 | 0.0001 | 104 MPa |
| | (e, e') | 0.19% (was 0.62%) | 1.0000 | 0.0001 | 106 MPa |
| | ξ | 60% (was 77%) | 0.93 | 0.14 | 282 kPa |

**No chatter.** In every run, (r, ṙ) and (e, e') match the exact drive's slope-reversal rate and total variation.

**Start-up.** Without a ramp, start-up peaks are larger with this training (10–300 MPa). The strong target with no ramp fails in (r, ṙ) at start-up. The 10 τ ramp removes the spikes and gives a 95 kPa peak, equal to the steady state.

**Plots.** `viewer_outputs/klqr-coords-all9-ic0p3747-20261114/{mild,strong,hard}/` holds radius and drive time series (start-up and steady), log10 spectra, and one phase portrait per controller. The strong set uses the 10 τ ramp for all three controllers.

## Best result per coordinate system and target

### Settings explained

**Training trajectories.** Both are unforced viscous RPE decays from rest, used only for the eigenfunctions.

- **A: r₀ = 2.333 baseline.** The fine record (dt 1e-4, τ 0–80), subsampled with stride 10. It is a violent collapse to r = 0.23 that then rings down toward r = 1.
- **B: r₀ = 0.3747.** It starts compressed and expands to r = 1.86, over τ 0–25 at dt 1e-3, stride 1. It is a milder excursion that spends more of its time near equilibrium.

**Training windows.**

- **Early fit:** H-DMD, state SINDy and eigenfunction SINDy all use τ 0–25 of the decay. This covers large excursions.
- **Late fit:** H-DMD uses τ 50–80 and both SINDy fits use τ 40–80. Only the 80 τ trajectory A has this window, and it contains only small oscillations about r = 1.

**Eigenfunctions.**

- **(r, ṙ):** degree-4 SINDy on the state, with γ from the system right-hand side.
- **(e, e'):** the same fit on the far-field echo and its derivative (`echo_jet2`).
- **ξ:** causal echo history with 50 delays and 2 latents, a degree-1 fit and finite-difference γ.

**B_ψ.**

- **(r, ṙ) and (e, e'):** exact and analytic, ∇ψ·g(x) with g = (0, −9.09/r). No data is used.
- **ξ:** one constant matrix, regressed on the forced record: 60 τ from rest, random ±20 kPa pulses lasting 2 τ, `fit_start` 5.

**Controller.** KLQR with Q = 5000·I.

- **(r, ṙ) and (e, e'):** R = 1e-4. The Riccati equation is re-solved continuously inside the integrator.
- **ξ:** sampled on its 0.001 τ history clock, with R = 1e3 or 1e4. Smaller R diverges or chatters.

**Toggles.** No cap. The 10 τ ramp (quintic from (1, 0)) is used only where noted.

### Results (scored on τ 100–300)

| Target | Coordinates | Training | Ramp | R | Drive error | Correlation | Radius RMSE | Start-up peak |
|---|---|---|---|---|---|---|---|---|
| mild 1 + 0.05 sin | (r, ṙ) | A, late | no | 1e-4 | **0.01%** | 1.0000 | 0.0000 | 158 MPa |
| | (e, e') | A, late | no | 1e-4 | **0.00%** | 1.0000 | 0.0000 | 193 MPa |
| | ξ | A, late | no | 1e3 | **8.7%** | 0.9997 | 0.0030 | 77 kPa |
| strong 1.5 + 0.4 sin | (r, ṙ) | B, early | 10 τ | 1e-4 | **0.01%** | 1.0000 | 0.0001 | 95 kPa |
| | (e, e') | B, early | 10 τ | 1e-4 | **0.03%** | 1.0000 | 0.0002 | 95 kPa |
| | ξ | A, late | 10 τ | 1e3 | **6.1%** | 0.976 | 0.051 | 94 kPa |
| hard four-tone | (r, ṙ) | B, early | no | 1e-4 | **0.20%** | 1.0000 | 0.0001 | 104 MPa |
| | (e, e') | B, early | no | 1e-4 | **0.19%** | 1.0000 | 0.0001 | 106 MPa |
| | ξ | A, late | no | 1e4 | **37%** | 0.982 | 0.077 | 419 kPa (exact 656) |

None of these nine chatter: each matches the exact drive's slope-reversal rate.

### Why these are the best

- **Fit where the target lives.** Exact-B_ψ coordinates do best when the eigenfunctions are fitted where the target lives. Near equilibrium (mild), the late fit wins. For large excursions (strong, hard), an early fit is needed. Trajectory B beats A there because A's collapse makes ψ(r, ṙ) fold near r ≈ 0.89.
- **ξ always needs the late fit.** Its eigenfunctions are linear (degree 1) and its B_ψ is one constant regressed about r = 1, so it is only accurate near equilibrium. With an early fit it falls to 38–77% error on all three targets.
- **ξ needs a tamed start-up.** Without the ramp or a cap, late-fit ξ diverges at start-up on the strong and hard targets. On the hard target it diverges even with the ramp unless R ≥ 1e4.
- **R trades accuracy against chatter for ξ.** On the hard target, R ≤ 3e3 chatters against any cap and R = 1e4 is best. A constant B_ψ cannot supply the state-dependent gain the 656 kPa compression peaks need.

### Plots

`viewer_outputs/klqr-coords-best-20261114/{mild,strong,hard}/` holds radius and drive time series (start-up and steady), log10 spectra and phase portraits. Legend entries name each run's fit. These runs store the control output every 0.001 τ (`--set control.time.output_step=0.001`) so that start-up transients are resolved. Their scores are identical to the 0.02 τ runs.

Start-up plots use the Viewer inset policy. A short off-scale episode (the spike plus any settling, such as ξ ringing as it locks on) moves into an inset with its own axes. That inset uses a symmetric-log ordinate when the episode spans decades beyond the settled level, and the main panel keeps the meaningful range. The legend stays top right and the inset bottom right unless they would hide the inset's own episode or a prominent extreme.
