# Echo-only KLQE → KLQR on the RPE: what works and where it stops (2026-11-16/17)

Setting: viscous RPE bubble, single scalar drive. Online, the estimator measures only the
far-field echo e (y = e₀ + C_e ψ, linear readout of a few Koopman eigenfunctions, no ė, no
other sensor). KLQR acts on the estimate ẑ with a constant B_ψ. Targets in state space:
mild r = 1 + 0.05 sin(1.893 τ) (100 kHz); strong r = 1.5 + 0.4 sin (held expanded, r 1.1–1.9);
hard, a four-tone signal with r 0.62–1.63. Scores on τ 40–100.

## Results

| Target | Eigenfunctions | Full-state KLQR (radius RMSE / drive) | Echo-only KLQE → KLQR |
|---|---|---|---|
| mild | H-DMD small-amplitude pair, interpolated | 0.0022 / 69% (chatter) | 0.0039 / 62% (chatter) |
| mild | same pair, degree-5 polynomial fit | 0.00098 / 2.5% | 0.0036 / 11% |
| mild | SINDy ψ(e, e′), IC 0.37 | 0.0002 / 0.12% | 0.022 → 0.041 (biased, see below) |
| mild | same + orbit readouts + **integral action (LQI)** | — | **0.0022 / 5.3%** |
| mild | same + orbit readouts + PID correction | — | 0.015 |
| strong | SINDy ψ(e, e′), IC 0.37 | 0.0027 / 0.4% | diverges (all variants) |
| hard | SINDy ψ(e, e′), IC 0.37 | 0.0016 / 0.5% | diverges (all variants) |
| strong, hard | any small-amplitude K-ROM | fails beyond r ≈ 1.25 | fails beyond r ≈ 1.25 |

Drive errors are relative to the exact inverse drive (mild) or the 2026-11-15 full-state drive.

## What we learned

1. **Control and estimation need different properties of the eigenfunctions.** Full-state KLQR
   tolerates poor invariance: ψ(e, e′) misses invariance by 0.66 over the range, yet controls
   all three targets almost perfectly when ψ is measured. KLQE does not: it predicts ẑ with
   ż = Λz + B_ψu, so it needs eigenfunctions that are invariant along trajectories (small-amplitude
   H-DMD: 0.0003 after smoothing). No pool found is both invariant and controllable over the
   strong/hard range.
2. **The linear echo readout limits the range.** One pair explains 98.5% of the echo near
   equilibrium, 84% (ψ(e, e′)) or 45% (ψ(r, ṙ)) over the collapse range. Fitted over a wide range
   the readouts are right on average and biased on any given orbit; the KLQE then shifts ẑ to
   cancel the echo-readout error. Fitting readouts on the operating orbit fixes the estimate
   (radius estimate error 0.08 → 0.007).
3. **The echo barely sees the mean radius.** With a correct estimate the loop still settled 0.04
   low. An input-disturbance (bias) state in the KLQE did not help. Integral action on the
   estimated radius error (LQI) did: 0.041 → 0.0022 on mild. A PID integrator on the same signal
   helps less (0.015) because it only removes the mean; LQI redesigns the whole gain.
4. **How eigenfunctions are evaluated matters.** Interpolating H-DMD values (inverse distance,
   8 neighbours) made ψ_ref rough and the drive chatter (62–69% raw drive error). A polynomial fit
   of the same values removed it (2.5% full-state, 11% echo-only) and improved invariance.
   H-DMD and polynomial (SINDy-style) forms of one eigenfunction are compatible once readouts
   and B_ψ are refit on the same form.
5. **Training data bound the reach.** Free decays reach the strong/hard region only during
   violent collapses (ṙ to ±10); H-DMD finds no genuine fundamental there and polynomial fits
   cannot represent what it finds. Strong (held expanded) is open-loop unstable.
6. **Constant B_ψ** explains 81% of the forced response at small amplitude and ~37% over the range.

## Limits of echo-only KLQE → KLQR as built

- Works: small-amplitude tracking around equilibrium (mild, ±5%), with smooth invariant
  eigenfunctions, readouts fitted on the operating orbit and integral action.
- Fails: offset or large-amplitude targets (beyond about r ≈ 1.25 here), even though full-state
  KLQR with suitable eigenfunctions handles them. The bottleneck is estimation: invariance,
  a linear echo readout from few modes, and an echo that is blind to the mean radius.
- Likely needed for more: training data covering the operating region (forced or multi-IC
  eigenfunction fits), state-dependent input maps, more modes with a richer (still linear)
  readout, or separate estimator/controller coordinates.

## Separate eigenfunctions for estimation and control

New general mode `controller.settings.feedback = "estimated_state"`: the KLQE (own KROM, echo
only) decodes x̂ = (r̂, ṙ̂), and the controller evaluates its own eigenfunctions at x̂. Controller:
SINDy ψ(r, ṙ) from IC 0.37 with its constant map, Q 5000 (full-state: mild 0.0003 at R 0.01;
strong 0.007, hard 0.009 at R 1).

Estimators ranked passively (full-state ψ(r, ṙ) control, KLQE observing; radius / ṙ errors):

| Estimator | mild | strong | hard |
|---|---|---|---|
| smooth H-DMD (small amplitude) | r 0.0015, ṙ 0.017 | r 0.28 | r 0.14 |
| ψ(e, e′), orbit readouts | r 0.0014, ṙ 0.18 (bias +0.17) | r 0.030, ṙ 1.78 (bias −1.37) | r 0.099, ṙ 0.70 |
| ψ(e, e′), range readouts | r 0.079 | r 0.29 | r 0.12 |

Separated closed loops (`closed-loop/separated-*.json`):

| Target | Estimator → ψ(r, ṙ) KLQR | Radius RMSE | Drive |
|---|---|---|---|
| mild | smooth H-DMD | **0.0064** | 15% |
| mild | ψ(e, e′) orbit readouts | 0.041 | 141% |
| strong, hard | either (control from τ = 5, ramp 20, cap) | diverge at τ ≈ 7–11, when the target passes r ≈ 1.2–1.3 | — |

- **Error propagation through ψ_ctrl.** Relative ψ(r, ṙ) error is about 2× the relative state error
  with the H-DMD estimator (ψ 0.51 for ṙ̂ 0.25) and 4–13× with the ψ(e, e′) estimator, whose ṙ̂ is
  biased (its ẑ offset barely moves r̂ but moves ṙ̂). The velocity estimate, not the radius, decides.
- The ψ(e, e′) estimator does not even converge at rest (echo ≡ 0 at equilibrium): x̂ = (0.84, 0.44)
  while x = (1, 0). The H-DMD estimator converges there and follows the ramp until r ≈ 1.2, then lags.
- Separation works where an accurate echo-only estimator exists (mild), and the controller can
  then use eigenfunctions chosen for control. It does not extend the reach: the estimator's own
  range (small amplitude) becomes the limit.

## Shadow runs: estimation alone, on the ideal trajectory (2026-11-18)

The plant receives the ideal drive (the same KLQR on the true state); the echo-only KLQE → KLQR runs
beside it and its force is recorded, not applied (`controller.settings.applied_controller`). The
estimator is wrong on the correct trajectory: on strong and hard, ṙ̂ errors are 0.7–1.8 against a
ṙ RMS of 0.5–0.7, and r̂ errors 0.03–0.28. The force is wrong in every case, by 3.6× (mild, smooth
H-DMD) to over 100×, because the true-state KLQR makes the ideal drive from tracking errors of
0.0002–0.003, far below the estimate errors. The failure needs no feedback spiral: the estimate is
already too poor on the trajectory the controller is trying to hold. Details and table:
`examples/native/case-studies/rpe-echo-mode-selection-20261116/README.md` (shadow runs).

## Classical loop margins (2026-11-18)

`linear_systems` now reports transfer functions, Bode data and the margins of the LQR,
Kalman, LQG and LTR loops. LQR and Kalman loops keep their guarantee (min|1 + L| = 1) on every
K-ROM. The ψ(e, e′) echo channel has a right-half-plane zero (about +21): its LQG input loop
keeps only a disk gain margin of [0.60, 2.9] and 38°, and loop transfer recovery makes it worse.
ψ(r, ṙ) and smooth H-DMD are minimum phase and their LQG loops keep near-LQR margins, yet their
echo-only estimates fail just the same. LQG fragility is a contributor for ψ(e, e′), not the
cause: the estimator is already wrong with the loop open, because the K-ROM is far from the plant
over the range. Table: case-study README (classical loop margins).

## Observability (2026-11-18)

The real RPE's echo, perturbed along the true orbits (`empirical_observability`), is almost
blind to ṙ at any instant (∂e/∂ṙ 9–163 times weaker than ∂e/∂r, orbit-scaled), but over one
period it determines r and ṙ to a few percent on strong and hard (gramian condition 5–6);
near equilibrium the signal is 30–100 times weaker. The K-ROMs see their single oscillation
perfectly (condition 1.005), but the tuned KLQE (W ≫ V, needed because the K-ROM dynamics are
wrong over the range) pins one observable direction and leaves the quadrature direction loose
by 3,000–10,000, so it predicts ṙ̂ 60–290 times less certain than r̂. The smooth H-DMD K-ROM
that worked on mild has an anisotropy of only 94. The echo's information is there on the
large orbits; the estimator cannot use it without accurate dynamics. Details: case-study
README (observability).

## KLQE tuning (2026-11-19)

198 shadow runs over scalar W (W/V from 1e-5 to 1e10), W shaped from the K-ROM's own prediction
residual, and input-direction W + q²B_ψB_ψᵀ. Scalar W only trades ṙ̂ error for r̂ error (trusting
the model helps ṙ̂ on strong but ruins r̂, and ruins both on hard). Structured W is best everywhere
and cuts ṙ̂ error 1.5–2.7× and the shadow force error up to 2.3× on strong/hard, but the force stays
25–100 times the ideal drive; biases (ψ(r, ṙ) mild r̂ +0.21) do not move, and the innovations stay
strongly correlated (autocorrelation 0.84–0.98 at 0.1 τ). Tuning cannot replace a model that is
accurate over the orbit. Table: case-study README (KLQE tuning sweep).

## Framework additions made along the way

`mode_selection` (eigenfunction subset selection, `eigenfunction_polynomial_degree`,
`readout_records`, analytic-coordinate pools), general PID compensator and LQI with KLQR
adapters (`integral_action`, `pid_correction`), `feedback: "estimated_state"` (any controller on
the KLQE state estimate, `docs/architecture/ESTIMATOR_MODULE.md`), `applied_controller` (shadow KLQE → KLQR), classical transfer functions and loop margins (`LoopMargins.hpp`, `TransferFunction.hpp`, Viewer `bode`/`nyquist`/`loop_bode`/`sensitivity`), empirical observability (`EmpiricalObservability.hpp`, `empirical_observability`, `linear_systems` `readouts`), KLQE process covariance from K-ROM residuals and input direction (`KlqeProcessCovariance.hpp`), Viewer CLI for archive
comparisons.
Details: `docs/architecture/MODE_SELECTION.md`, `docs/architecture/INTEGRAL_AND_PID_CONTROL.md`,
`examples/native/case-studies/rpe-echo-mode-selection-20261116/README.md`.
