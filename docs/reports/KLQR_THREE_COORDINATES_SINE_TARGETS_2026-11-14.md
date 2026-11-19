# KLQR on prescribed sinusoidal radius targets in (r, ṙ), (e, e') and ξ eigenfunction coordinates (2026-11-14)

## Setup

- **Targets.** Two prescribed state waveforms at 100 kHz (1.8927 rad/τ), with ṙ the exact derivative. They are kinematic targets, not orbits produced by sinusoidal forcing.
  - **Mild:** r = 1 + 0.05 sin(1.8927 τ).
  - **Strong:** r = 1.5 + 0.4 sin(1.8927 τ). It is strongly nonlinear: its exact inverse drive has a −77.5 kPa mean and lines of 4.8, 17.7 and 2.7 kPa at 100, 200 and 300 kHz.
- **Exact inverse drive.** The reference is u = (r̈_T − f(r_T, ṙ_T)) / (−9.0926 / r_T), the pressure that makes the RPE follow the target exactly.
- **Pipeline.** `examples/native/pipelines/rpe-echo-jet-klqr-rediscovery/20kpa.json`; only the coordinates change between attempts.
  - Training: the baseline RPE decay from r = 2.333 (fine record).
  - H-DMD and state SINDy produce the spectrum.
  - Eigenfunction SINDy fits ψ in the chosen coordinates.
  - The input map step produces B_ψ, and KLQR closes the loop.
- **B_ψ per attempt:**
  - **(r, ṙ):** `exact`, ∇ψ·g(x), analytic, no forced data.
  - **(e, e'):** `exact`, ∇_(e,e')ψ · J_(e,e')(x) · g(x), analytic.
  - **ξ:** a 50-delay echo history with 2 latents. B_ψ is `constant`, held-input regressed from a 60 τ ±20 kPa random-pulse forced record. ξ always runs on its sampled history clock (0.001 τ).

## What went wrong in the first attempt, and the MATLAB-matched settings

The first runs used settings carried over from the earlier ξ work:

- R = 1000;
- a ±50/±150 kPa cap;
- a 10 τ target ramp;
- **sampled** control (0.001 τ) for every coordinate system;
- **late-time** training: H-DMD on τ 50–80, SINDy on τ 40–80.

With those, (r, ṙ) chattered bang-bang between the limits on the strong target.

The MATLAB suite (`dynsyslab (5-15-25)`, `CTRL.m`, `KLQR.m`, `extract_ROM.m`, `KLQR_ctrl_tuning.mat`) does it differently:

- **Control law.** It solves `lqr(A, B(x), Q, R)` at every evaluation inside ode15s, so control is continuous and the Riccati solve is state-dependent. The native KLQR re-solves the same way.
- **Weights and cap.** The tuning table gives Bubble Q = 5000·I and R = 10⁻⁴; `hardcap = 0`, so there is no cap.
- **No ramp.** The reference is used as given.
- **Early-time training.** H-DMD on r uses 20,000 rows × 2000 columns over the whole trajectory, and SINDy fits the eigenfunctions on the whole trajectory with polynomial order 4 and tolerance 10⁻⁶.

Matching those settings (continuous evaluation, R = 10⁻⁴, no cap, no ramp) removes the chatter: (r, ṙ) eigenfunctions give near-exact drivers. The chatter was a consequence of sampled evaluation combined with the R and cap choices. It was not a flaw in the eigenfunctions.

## Results with MATLAB settings (Q = 5000, R = 10⁻⁴, no ramp, no cap; ξ at R = 1000)

All scores are taken over τ ≥ 100.

| Training | Target | Attempt | Drive error vs exact | Drive correlation | Radius RMSE | Mean drive (exact) |
|---|---|---|---|---|---|---|
| early (τ 0–25) | mild | (r, ṙ) | fails (+126 kPa offset) | −0.09 | 0.18 | +126 (+1.4) |
| | | (e, e') | 0.3% | 1.0000 | 0.0001 | +1.37 (+1.37) |
| | | ξ | 68% (6.5 of 19.9 kPa) | 0.9993 | 0.024 | +0.28 (+1.37) |
| early (τ 0–25) | strong | (r, ṙ) | 0.6% | 0.9998 | 0.0030 | −77.2 (−77.5) |
| | | (e, e') | 0.2% | 0.9999 | 0.0016 | −77.4 (−77.5) |
| | | ξ | 58% | 0.37 | 0.45 | −34.6 (−77.5) |
| late (τ 40/50–80) | mild | (r, ṙ) | 0.01% | 1.0000 | 0.0000 | +1.37 (+1.37) |
| | | (e, e') | 0.00% | 1.0000 | 0.0000 | +1.37 (+1.37) |
| | | ξ | 8.7% | 0.9997 | 0.0030 | +1.17 (+1.37) |
| late (τ 40/50–80) | strong | (r, ṙ), (e, e') | start-up fails: Riccati solution not real | | | |
| | | ξ | diverges | | | |

### Findings

- **Exact coordinates track exactly within their training range.** (r, ṙ) and (e, e') with exact B_ψ reproduce the exact inverse drive to within 0.01–0.6% wherever the training covers the target: late training for the mild target, early training for the strong one.
- **(e, e') is the most robust.** It also works on the mild target with early training, where early-trained ψ(r, ṙ) fails near equilibrium.
- **Start-up without ramp or cap is extreme.** The target jumps at τ = 0, and the first commands reach 50–200 MPa before settling. With late training and the strong target, the start-up leaves the region where B_ψ is well-conditioned and the Riccati solve fails. With the ramp and R = 10⁻³ the same case was exact: correlation 1.0000, mean −77.49 kPa.
- **ξ with a regressed constant B_ψ needs a much larger R.** With R ≤ 30 it diverges. Even at R = 1000 it is only good near equilibrium with late training. Its map was regressed from ±20 kPa pulses about R₀.

## Toggles added

- **Target ramp.** Any state target may carry `ramp` = `{duration, from, derivative_pairs}`; `null` or `false` uses the target as given.
- **Control cap.** Every feedback controller accepts `control_limit` as a number, as `[low, high]`, or as `null`/`false` for no cap.
- **Pipeline parameters.** Both are exposed as `target_ramp` and `control_limit` in `base-r.json` and `20kpa.json`.
- **Null parameters.** A pipeline parameter set to `null` is now a defined value. Parameters merge without merge-patch deletion.

## Plots

These Viewer plots show drive and radius, as time series and as spectra over 60 target periods:

- `viewer_outputs/klqr-coords-matlab-settings-{early_mild,early_strong,late_mild}-20261114/`
- `viewer_outputs/klqr-coords-sine-20261114/` and `viewer_outputs/klqr-coords-sine-strong-20261114/`: the first, sampled, late-trained, capped and ramped runs.
