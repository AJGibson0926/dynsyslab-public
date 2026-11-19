# Eigenfunction subset selection for echo-only KLQE → KLQR (2026-11-16)

First runs of the [`mode_selection`](../../../../docs/architecture/MODE_SELECTION.md) method on
two candidate pools. The estimator measures the echo only, y = e ≈ e₀ + C_e ψ; ψ is inferred, never
measured. Results are in `archives/derived/mode-selection/`; selected K-ROMs in
`archives/derived/koopman-mode-selection/`.

```powershell
dynsyslab-cli run examples\native\pipelines\rpe-echo-mode-selection\echo-pool.json --archive-root archives
dynsyslab-cli run examples\native\pipelines\rpe-echo-mode-selection\hdmd-pool.json --archive-root archives
```

## Pools

| Pool | Eigenvalues | Eigenfunctions | Data |
|---|---|---|---|
| ψ(e, e′) (`echo-pool.json`) | first 12 stable pairs of the echo-only H-DMD of the noisy-echo record (r₀ = 0.3747, 8-bit quantized, RTS e′) | degree-3 SINDy in the measured (e, e′) coordinates, ψ̇ from the measured derivatives | unforced: that record; forced: a005p05, a02p05, a02p2 coordinate records (same noise/RTS) |
| H-DMD (`hdmd-pool.json`) | 2026-11-14 joint H-DMD of [r, ṙ, e] on the late baseline, t 50–80 (the working echo→state KLQE) | the H-DMD eigenfunctions themselves, interpolated in state space | unforced: t 50–80 of that record; forced: a005p05, a02p2 (both inside r 0.88–1.12) |

## Results

Per unit (integral form, h = 0.25; invariance = ‖ψ(t+h) − e^{λh}ψ(t)‖ / ‖ψ(t+h) − ψ(t)‖):

| | ψ(e, e′) pool | H-DMD pool |
|---|---|---|
| invariance, leading pairs | **1.0–1.3** for every mode | **0.13–0.16** |
| echo explained by the best single pair | 66% | 95% (λ = −0.06 ± 6.283i) |
| constant B_ψ, unexplained fraction | 0.92–0.99 | 0.58 (fundamental), 0.97–1.0 (harmonics) |

Best subsets:

| Pool | Subset | Echo error | Predicted state error (r, ṙ) |
|---|---|---|---|
| ψ(e, e′) | any single pair | 64–72% | 78% (52–60%, 92–98%) |
| ψ(e, e′) | 4 pairs (fundamental + harmonics) | **0.4–0.6%** | 100–500 (diverges) |
| H-DMD | fundamental pair | 21% | 20% (22%, 18%) |
| H-DMD | fundamental + 2nd harmonic (`…192709Z`, written K-ROM) | 19% | **18.5%** (21%, 16%) |

## Findings

1. **The ψ(e, e′) polynomials are not Koopman eigenfunctions along these records.** ψ(t+h) is no
   better predicted by e^{λh}ψ(t) than by ψ(t). This holds for every fit of this family in the
   repository, clean or noisy, degree 2–4, including the pair behind the 2026-11-15 KLQE → KLQR
   results (SINDy's own residual Dψ·q̇ − λψ is 99% of Dψ·q̇). Those runs worked because they
   measured ψ(e, e′) directly every step and never relied on ż = Λz to predict it. An echo-only
   KLQE does rely on it, so this pool cannot serve, whatever subset is chosen: harmonics make the
   echo readout nearly exact but the predicted state error diverges.
2. **The H-DMD eigenfunctions are genuine** (13–16% horizon miss, about what the interpolation
   costs), and the fundamental pair alone carries 95% of the echo. Predicted echo-only state
   error is about 18–20%. The working echo→state KLQE measured 5–12% on held-out records, so the
   prediction is conservative.
3. **Constant B_ψ is the binding constraint.** Even inside the training region it explains 42% of
   the fundamental's forced response, with values that differ between the two forced records
   (amplitude dependence), and almost none of the harmonics'. A forced record that leaves the
   region (a02p05, r 0.73–1.37) wipes out the regression.
4. **Operating region.** The H-DMD eigenfunctions cover r ≈ 0.92–1.09. The strong and hard targets
   reach r ≈ 1.9. Eigenfunctions, readouts and B_ψ must be trained over the region the controller
   will visit.
5. **Harmonics** add echo and state information but are poorly actuated by a constant B_ψ. That
   is the trade-off the selection was built to expose.

## Retraining over the target range, and the best echo-only closed loop

Target ranges: mild r 0.95–1.05, strong 1.1–1.9, hard about 0.62–1.63 (ṙ up to ±1.3).

**Eigenfunctions.** H-DMD of [r, ṙ, e] on the r₀ = 2.0 free decay (`range-pool.json`, windows tried):

| Window (τ) | Region | Leading stable pair | Invariance |
|---|---|---|---|
| 0–40 | r 0.34–2.0 | none near ω ≈ 6; four heavily damped pairs | ≈ 1.0 |
| 10–50 | r 0.42–1.77 | −11.2, −21.7 ± 5.6i | 0.54–0.57 |
| 20–60 | r ≈ 0.7–1.35 | −0.76 ± 15.6i (no ω ≈ 6 pair) | 0.54–0.97 |
| 30–70 | r ≈ 0.8–1.2 | −0.185 ± 6.283i | **0.07** |

H-DMD of free decays yields genuine eigenfunctions only up to amplitudes of roughly 0.1–0.2 around
equilibrium; the collapse regime breaks it. **The strong (r to 1.9) and most of the hard target lie
outside every genuine pool**, so neither can be attempted honestly with this family.

**Forced data.** Open-loop drives fitted to the 2026-11-15 full-state drives (`forced/*.json`) stay on
the mild orbit but run away on strong (r to 3.2) and hard (r to 17): holding the bubble expanded is
unstable without feedback. B_ψ is therefore regressed only on in-range forced records.

**Pool used** (`selection-v0-fundamental.json`, pipeline `range-pool.json` with its defaults): H-DMD of
the small-amplitude decay r₀ = 1.1, ṙ₀ = 0.3 over τ 0–42 (amplitude 0.115 → 0.01). Fundamental pair
λ = −0.060 ± 6.283i: invariance 0.05, echo explained 97.8%, constant B_ψ explains 81% of the forced
response.

**Bug found on the way.** H-DMD stores each conjugate partner as c·conj(ψ) with an arbitrary c, while
the written K-ROM gave it conj(ν) and conj(B_ψ). The controller then saw a wrong reference (it barely
drove; the bubble sat at r = 1). The K-ROM writer now rewrites the partner as conj(ψ).

**Best echo-only KLQE → KLQR run** (`closed-loop/echo-only-klqe-klqr-mild-best.json`, archive
`echo-only-klqe-klqr-mild-best__20261116T200921Z`; figure
`viewer_outputs/echo-only-klqe-klqr-20261116/best-mild.png`):

- Online: the estimator measures `far_field_acoustic_echo` alone (y = e₀ + C_e ψ̂); KLQR acts on ψ̂ with
  the constant B_ψ. Q = 1, R = 3, W = 1e5, V = 1e-4, 10 τ target ramp, drive cap 0.5, 0.001 τ clock.
- Scores on τ 40–100: radius RMSE **0.0039** (7.9% of the 0.05 amplitude), radius estimate error
  **0.0015**, drive error 62% raw / **13.5% low-passed** at 2 cycles/τ, against the 2026-11-15 full-state
  drive.
- Same K-ROM with full-state feedback (`fullstate-same-krom-mild.json`): radius 0.0022, drive 12.8%
  low-passed. Echo-only estimation costs little; the remaining error is the one-pair model with a
  constant B_ψ (amplitude-dependent input effect, ripples at the natural frequency near the troughs).
- Sweep (mild, fundamental K-ROM): R 1–10 and W ≥ 1e3 all give radius 0.0038–0.0043; larger R loses
  tracking (R = 1000: 0.0175); adding the second 6.283 pair changes nothing.

### Drive chatter came from interpolating the eigenfunctions

Both loops above chattered (raw drive error 62% echo-only, 69% full-state). The H-DMD
eigenfunctions are interpolated between training states (inverse distance, 8 neighbours), and
the training states are one decaying spiral. Along the mild target, the interpolated reference
ψ_ref has 7.2% of its energy above five times the target frequency; a degree-5 polynomial fit of
the same values has 1.3% and differs from it by 8%.

`selection-v1-fundamental-poly5.json` reruns the selection with
`eigenfunction_polynomial_degree: 5` (fit error 3e-4 on the fundamental pair). The fundamental
pair's invariance miss drops from 0.05 to 0.0003, echo explained rises to 98.5%, and the B_ψ
residual falls to 0.16. Same closed-loop settings (Q 1, R 3, W 1e5, V 1e-4, ramp 10, cap 0.5):

| Mild, τ 40–100 | radius RMSE | raw drive error | radius estimate error |
|---|---|---|---|
| echo-only, interpolated (`echo-only-klqe-klqr-mild-best`) | 0.0039 | 62% | 0.0015 |
| full-state, interpolated (`fullstate-same-krom-mild`) | 0.0022 | 69% | — |
| **echo-only, smooth** (`echo-only-klqe-klqr-mild-poly5`) | **0.0036** | **11%** | 0.0015 |
| **full-state, smooth** (`fullstate-poly5-krom-mild`) | **0.00098** | **2.5%** | — |

The chatter is gone (raw and low-passed drive errors are equal). Echo-only is now limited by the
one-pair state readout: W from 1e3 to 1e7 leaves the estimate error at about 0.0015 and the
radius RMSE at 0.0036; the orbit is a slightly offset, smooth ellipse. Figures:
`viewer_outputs/echo-only-klqe-klqr-20261116/interpolation-test/`.

## Strong and hard: eigenfunctions over the range (2026-11-17)

The smooth mild K-ROM, unchanged, fails on strong and hard once the target leaves its training
region (r 0.89–1.12): at τ ≈ 20–25, with full-state feedback too (`closed-loop/sweep-*`).

Pools scored over the range (`range/selection-*.json`; unforced: the IC 0.37 decay, r 0.37–1.86,
and validation trajectory 0; forced: the 2026-11-15 mild/strong/hard full-state runs):

| Pool | invariance | echo explained | full-state KLQR (Q 5000, R 0.01, ramp 10) | echo-only |
|---|---|---|---|---|
| SINDy ψ(e, e′), IC 0.37 (the 2026-11-15 K-ROM) | 0.66 | 84% | strong 0.0027, hard 0.0016, mild 0.0002 | diverges on strong/hard; mild 0.022 |
| SINDy ψ(r, ṙ), IC 0.37 (3 variants) | 0.17 | 45% | diverges | diverges |
| H-DMD of the IC 0.37 decay (τ 0–25, 5–25), interpolated or degree 5/7 | 0.23–4 | ≤ 32% | — | — |

**Echo-only bias on mild (ψ(e, e′)).** A passive KLQE (full-state control, estimator observing)
shows a +0.079 radius-estimate bias. +0.059 of it is the state readout itself: x₀ + C_x ψ of the
*true* ψ is biased on the mild orbit because C_x was fitted over the whole collapse range. The
echo readout is biased the same way (0.108 of an echo rms of 0.15), and the KLQE cancels that
by shifting ẑ (−0.012 in Re ψ, larger than ψ's oscillation on the orbit).

- An input-disturbance (bias) state, z′ = Λz + B_ψ(u + d), with and without cancelling d̂ in the
  control, did not help (radius 0.022 → 0.025–0.041) and was not kept.
- Readouts fitted on the target orbit (`readout_records`, the 2026-11-15 run of each target)
  remove the estimator bias: radius estimate error 0.08 → **0.007** (`readout-ee-mild-echo-W1e4`).
  The closed loop still settles 0.04 low in mean radius (radius RMSE 0.041): the echo barely
  depends on the mean radius, and with this pair's poor invariance ẑ settles off ψ along a
  direction C_x hardly sees. Strong and hard still diverge echo-only.
- Integral action on the estimated radius error (general LQI, KLQR `integral_action`,
  `docs/architecture/INTEGRAL_AND_PID_CONTROL.md`) removes the remaining offset:
  `closed-loop/lqi-ee-mild-Qi3e6.json`, Q_i = 3e6: radius RMSE **0.0022**, drive 5.3% of the exact
  inverse (Q_i 1e5: 0.0055; 1e7 unstable). A PID correction on r̂ (negative gains, since the drive
  compresses the bubble) reaches 0.015 (`pid-ee-mild-Ki-30.json`). Strong and hard still diverge
  echo-only with LQI, slower ramps and caps (`lqi-ee-strong-Qi1e5-ramp40.json`,
  `lqi-ee-hard-Qi1e4-ramp40.json`).

Summary of the whole study: `docs/reports/ECHO_ONLY_KLQE_KLQR_LIMITS_2026-11-17.md`.

## Shadow runs: echo-only KLQE → KLQR beside the ideal drive (2026-11-18)

Question: is the failure a feedback spiral (a poor estimate drives the plant off target, which
makes the estimate worse), or is the estimator already wrong on the correct trajectory?

New option `controller.settings.applied_controller` (`docs/architecture/ESTIMATOR_MODULE.md`): the
plant receives a known correct drive; the echo-only KLQE → KLQR runs beside it, predicts with the
applied drive, and its force is archived (`commanded_input`) but never applied. The ideal drive is
the **same KLQR on the true state** (KLQR has no feedforward, so the exact inverse drive is not what
it would command even with a perfect state). With the same controller in the shadow, the force error
is exactly −K(ẑ − ψ(x)): estimation alone. Exception (*): the smooth H-DMD K-ROM cannot track strong
or hard on the true state, so those two use the 2026-11-15 ψ(e, e′) full-state drive.

`closed-loop/shadow_runs.py` (jobs `shadow-*.json`, scores `shadow-scores.json`),
`closed-loop/shadow_figures.py` (Viewer CLI, `viewer_outputs/shadow-klqe-klqr-20261118/`). Scores on
τ 40–100; estimate errors as RMS (mean).

| Target | K-ROM (estimator = controller) | plant tracking r (ideal drive) | r̂ − r | ṙ̂ − ṙ (ṙ RMS) | KLQR force error / ideal drive RMS | correlation |
|---|---|---|---|---|---|---|
| mild | ψ(e, e′), orbit readouts | 0.0002 | 0.0014 (−0.0005) | 0.18 (+0.17) (0.067) | 153× | 0.98 |
| mild | smooth H-DMD | 0.0010 | 0.0014 (+0.0007) | 0.016 (−0.006) (0.065) | 3.6× | −0.11 |
| strong | ψ(e, e′), orbit readouts | 0.0027 | 0.028 (+0.010) | 1.78 (−1.38) (0.53) | 218× | 0.17 |
| strong | smooth H-DMD * | 0.0022 | 0.28 (−0.23) | 0.88 (−0.66) (0.53) | 5600× | 0.19 |
| hard | ψ(e, e′), orbit readouts | 0.0016 | 0.11 (0.00) | 0.76 (−0.47) (0.68) | 72× | −0.42 |
| hard | smooth H-DMD * | 0.0013 | 0.16 (−0.11) | 0.72 (−0.36) (0.68) | 195× | −0.25 |

Findings:

1. **The estimator is wrong on the correct trajectory.** With the plant held on the target by the
   ideal drive, so no feedback can push it off, the strong and hard estimates still miss: ṙ̂ errors
   are as large as ṙ itself or larger, biased by −0.4 to −1.4, and r̂ errors reach 0.03–0.28. No
   feedback spiral is needed to explain the failure; feedback can only add to it.
2. **The KLQR force needs an estimate as accurate as the tracking error.** The true-state KLQR
   produces the ideal drive from tracking errors of 0.0002–0.003. On the estimate it commands
   −K(ẑ − ψ_ref), and ẑ − ψ(x) is 1–3 orders larger than ψ(x) − ψ_ref, so the shadow force is
   3.6× to over 100× the ideal drive in every case, including mild. The ψ(e, e′) mild force has the
   right shape (correlation 0.98) at 150× the size: a periodic estimate error locked to the orbit
   (the ṙ̂ bias of +0.17).
3. **Why mild still closes the loop.** Mild with the smooth H-DMD K-ROM tracked at 0.0036 in closed
   loop although its shadow force is 3.6× off. In closed loop the plant settles slightly off target,
   where the estimated tracking error yields a workable drive: the loop absorbs a small estimation
   error at the cost of tracking (0.0036 vs 0.0010 full-state). On strong and hard the error is too
   large to absorb.
4. **It is the velocity.** Every case has its largest relative error in ṙ̂. A linear echo readout of
   a few eigenfunctions does not recover ṙ over a wide orbit.

## Classical loop margins: is this LQG fragility? (2026-11-18)

Doyle (1978): LQR on the true state guarantees margins at the plant input (gain
[1/2, ∞), phase ≥ 60°); an LQG loop, with a Kalman filter in place of the state, has none.
`loop-margins/*.json` runs `linear_systems` on each echo-only K-ROM, with y = echo,
u = drive through the constant B_ψ, and the closed-loop runs' Q, R, W, V, plus an LTR
sweep q² = 1e2–1e8 (`loop-margins/run_margins.py`, table in `loop-margins/margins.json`,
Viewer figures in `viewer_outputs/loop-margins-20261118/`: `bode`, `nyquist`, `loop_bode`,
`sensitivity`, `pole_zero`). α = min|1 + L(jω)|; the disk margins follow from α.

Each K-ROM is one oscillatory pair in split form (four states, two pole/zero
cancellations), so the echo channel is second order with one real zero and relative degree 1.

| K-ROM | echo-channel zero | LQR α | LQG input α | LQG disk gain margin | LQG disk phase margin | LTR q² = 1e8, α |
|---|---|---|---|---|---|---|
| ψ(e, e′), orbit readouts (mild, strong, hard) | **+20.5 to +22.0 (right half plane)** | 1 | 0.65–0.66 | [0.60, 2.9] | 38° | 0.43 (worse) |
| ψ(r, ṙ), R 0.01 | −193 | 1 | 0.998 | [0.50, 413] | 60° | 0.99 |
| ψ(r, ṙ), R 1 | −193 | 1 | 0.96 | [0.51, 25] | 57° | 1.0 |
| smooth H-DMD | −98 | 1 | 0.98 | [0.50, 58] | 59° | 1.0 |

The LQR and Kalman loops sit exactly at their guarantee (α = 1) in every case.

1. **Doyle's effect is real, but only for ψ(e, e′).** Its echo channel is non-minimum
   phase (a zero at about +21), so the LQG input loop loses margin (disk gain [0.60, 2.9],
   38°), and loop transfer recovery makes it worse (α 0.66 → 0.43), as theory predicts for
   a right-half-plane zero. The ψ(r, ṙ) and smooth H-DMD models are minimum phase; their
   LQG loops keep almost LQR margins and LTR recovers them fully.
2. **Margins do not explain the failures.** ψ(r, ṙ) has near-LQR margins yet its echo-only
   estimate was just as poor in the shadow runs (mild r̂ off by 0.21, strong ṙ̂ error 1.0),
   and the shadow runs had no feedback loop at all. The margins describe tolerance to
   perturbations of the K-ROM; the K-ROM itself is far from the plant over the range
   (constant B_ψ explains about 37% of the forced response, invariance misses of 0.17–0.66).
3. **Where fragility can matter.** For ψ(e, e′) the tolerated input-gain range [0.60, 2.9]
   is comparable to the input-map error (B_ψ explaining 37% of the response is roughly a
   factor 2.7). So LQG fragility can add to the ψ(e, e′) closed-loop failures on top of the
   estimator bias; it is a contributor there, not the root cause.

## Observability: what the echo tells the estimator (2026-11-18)

Three views of the same question (`observability/run_observability.py`; numbers in
`observability/observability.json`; Viewer figures in
`viewer_outputs/echo-observability-20261118/`).

**1. The physical echo, along the true orbits** (`empirical_observability`: the real RPE
perturbed in r and ṙ at 56 points over τ 40–95 of each shadow run's plant trajectory, echo
integrated over one period τ = 1, states scaled by the orbit RMS; noise intensity 1e-4 as in
the KLQE):

| Orbit | gramian eigenvalues (scaled, median) | condition | best relative std, r / ṙ | instantaneous \|∂e/∂r\| : \|∂e/∂ṙ\| (scaled) |
|---|---|---|---|---|
| mild | 0.011, 0.0010 | 10.6 | 9.5% / 32% | 163 : 1 |
| strong | 0.30, 0.035 | 6.3 | 3.4% / 4.1% | 9 : 1 |
| hard | 1.0, 0.13 | 4.8 | 1.3% / 2.3% | 19 : 1 |

- One echo sample is almost blind to ṙ (the instantaneous echo is one number of two states,
  and its ṙ sensitivity is 9–163 times smaller); ṙ is learned only through the dynamics.
- Over one period the echo does determine both r and ṙ on strong and hard (condition 5–6,
  best possible errors of a few percent). The weak direction is mostly ṙ.
- Near equilibrium (mild) the signal is 30–100 times weaker and ṙ is the poorly observed
  direction. So "the echo tells us little" holds instantaneously and near equilibrium; on the
  large orbits the information is present in the echo history.

**2. The K-ROM the KLQE uses** (`linear_systems` with `readouts: ["state"]`, echo output):

| K-ROM | gramian singular values | effective condition | Kalman covariance anisotropy (observable plane) | Kalman-predicted std ṙ̂ / r̂ |
|---|---|---|---|---|
| ψ(e, e′) orbit readouts (mild, strong, hard) | ≈ 300–520, 300–520, ~0, 0 | 1.005 | 7,500–9,700 | 128–287 |
| ψ(r, ṙ) | 47.9, 47.6, 0, 0 | 1.005 | 2,900 | 266 |
| smooth H-DMD | 1.5e-3, 1.4e-3, ~0, ~0 | 1.02 | **94** | 61 |

- Each K-ROM is one oscillation in split form: two observable directions (amplitude and
  phase), perfectly conditioned; the other two are the duplicate conjugate copies,
  unobservable by construction and carrying none of r or ṙ. Absolute gramian sizes only
  reflect eigenfunction normalization.
- The tuned KLQE (W = 1e4–1e5 ≫ V = 1e-4) trusts each echo sample over the model's dynamics,
  so its steady-state covariance pins one observable direction and leaves the quadrature
  direction loose by a factor of 3,000–10,000. That loose direction is the one ṙ̂ reads, so
  the filter itself predicts ṙ̂ to be 60–290 times less certain than r̂.
- The smooth H-DMD K-ROM, the only echo-only loop that worked (mild), has an anisotropy of 94:
  its eigenfunctions are invariant, so it could be tuned to trust its dynamics.

**3. Conclusion.** The echo is a weak instantaneous sensor of ṙ and a weak sensor overall near
equilibrium, but its history over one period carries enough to fix r and ṙ on the strong and
hard orbits. Extracting it requires dynamics that are accurate over the orbit. The K-ROMs'
dynamics are not, so the KLQE is tuned to distrust them (large W), which reduces it to
inverting the instantaneous echo, exactly the measurement that cannot see ṙ. That matches
the shadow runs: the largest errors are in ṙ̂, and they appear with the plant on the correct
trajectory. Caveats: the empirical gramian is local (it cannot detect distinct states with
similar echo histories), the Cramér–Rao bound assumes an exact model, and the runs had no
sensor noise; V sets an information scale rather than an actual noise level.

## KLQE tuning sweep (2026-11-19)

Can a different W fix the echo-only estimate? 198 shadow runs (`tuning/run_tuning.py`, results
`tuning/tuning-results.jsonl`, summary `tuning/tuning-summary.json`, Viewer figures
`viewer_outputs/klqe-tuning-20261119/`): plant on the ideal drive, echo-only KLQE → KLQR in
shadow, V = 1e-4, for each target and K-ROM:

- scalar W = w·I, w = 1e-9 … 1e6 (W/V from 1e-5, trust the model, to 1e10, trust each sample);
- W shaped from the K-ROM's own one-step prediction residual on the 2026-11-15 full-state run of
  the same target (τ 20–100), ×1, ×10, ×100 (new `process_covariance.from_records`);
- W = I + q² B_ψB_ψᵀ, model error entering through the input, q² = 1e2, 1e4, 1e6 (new
  `process_covariance_input_weight`).

Errors on τ 40–100 as RMS / orbit RMS (ṙ̂ / r̂), and the shadow-force error / ideal drive:

| Target, K-ROM | study W | best low W (trust model) | best input-direction W | residual-shaped W (×100) |
|---|---|---|---|---|
| strong, ψ(e, e′) | 3.3 / 0.10, force 218× | 1.1 / 1.5 (W 1e-6) | **1.5 / 0.09, force 97×** | 1.6 / 0.09, 107× |
| strong, ψ(r, ṙ) | 1.9 / 0.66, 113× | 1.0 / 1.3 | **0.70 / 0.66, 82×** | 0.90 / 0.66, 90× |
| hard, ψ(e, e′) | 1.1 / 0.39, 72× | worse (14 / 5.7 at 1e-9) | **0.75 / 0.39, 48×** | 0.79 / 0.39, 53× |
| hard, ψ(r, ṙ) | 0.84 / 0.51, 25× | worse | **0.55 / 0.51, 27×** | 0.65 / 0.51, 26× |
| mild, smooth H-DMD | 0.25 / 0.04, 3.6× | 0.26 / 0.11 | **0.16 / 0.04, 3.1×** | 0.26 / 0.06, 4.4× |
| mild, ψ(r, ṙ) | 1.2 / 5.9 (r̂ 0.21 high) | 0.66 / 5.8 | 0.28 / 5.9 | 0.75 / 5.8 |

Findings:

1. **Scalar W only slides along a trade-off.** Above W/V ≈ 10 nothing changes (the filter inverts
   each echo sample). Trusting the model (W ≈ 1e-6) roughly halves or thirds the ṙ̂ error on strong
   but multiplies the r̂ error by 10–15, and on hard it makes both worse because the K-ROM
   dynamics drift off the orbit. No scalar W makes both small.
2. **Structured W helps, modestly.** Process noise along the input direction (or shaped from the
   K-ROM's residual) is the best tuning everywhere: ṙ̂ error falls 1.5–2.7× on strong and hard,
   the shadow-force error 1.4–2.3× (not at all for ψ(r, ṙ) on hard), with r̂ unchanged. The force
   is still 25–100 times the ideal drive and ṙ̂ still errs by 55–150% of its RMS, so retuning does
   not make strong or hard workable.
3. **Biases survive every W.** The ψ(r, ṙ) mild r̂ offset (0.21, 5.9× the orbit RMS) is identical at
   all 22 settings: it comes from the readout and the missing mean-radius direction, not the gain.
4. **The innovations are never white.** The normalized innovation autocorrelation at 0.1 τ is
   0.84–0.98 in every run: the K-ROM's error is deterministic and correlated, so no linear-Gaussian
   (W, V) is consistent with it. Covariance tuning cannot substitute for a better model.
5. One run (hard, ψ(e, e′), W = 1e5) stopped at design time: "continuous LQR Riccati solution is
   not real" in the shadow KLQR; its neighbors (1e4, 1e6) ran normally.

Data roles: eigenvalues and eigenfunctions OFF-r + OFF-e (H-DMD on r, ṙ, e of an unforced record);
readouts OFF-r + OFF-e; B_ψ from forced records (one is an offline full-state KLQR run); target
offline; online ON-e only.

## e′ in these runs

The ψ(e, e′) pool uses the RTS-smoothed numerical e′ of the quantized echo in its coordinates and
SINDy targets. No ψ̇ is differentiated numerically in the selection (integral form), and the
H-DMD pool uses no e′. Each archive lists this in `e_dot_entry_points` and `records[].e_dot`.

## Superseded runs

`…191729Z`, `…191832Z` (pool SINDy without `input_space: coordinates`, so ψ̇ came from the RPE
right-hand side evaluated on (e, e′), wrong), `…191956Z`, `…192247Z` (derivative-based invariance,
inflated by differentiating ψ), `…192415Z`, `…192450Z`, `…192451Z` (three forced records, and
single-record probes).
