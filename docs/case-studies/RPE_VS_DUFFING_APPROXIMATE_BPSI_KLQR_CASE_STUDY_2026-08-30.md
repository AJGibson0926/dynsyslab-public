# RPE versus Duffing: approximate `B_psi` KLQR case study

**Date:** 2026-08-30  
**Classification:** important retained KLQR control case study  
**Status:** read-only comparative diagnosis completed  
**Scope:** explain why constant and bilinear Koopman input-map approximations succeeded for RPE but failed or remained inconclusive in the bounded Duffing search  

## Executive finding

The different outcomes are not primarily explained by a Python regression bug,
by a generally inferior Duffing SINDy eigenfunction, or by failure of KLQR
itself. The RPE experiment is a smooth, local, same-region tracking problem in
which the fitted input maps preserve the effective direction of analytic
`B_psi` well enough for the Riccati gain to compensate for most magnitude
error. The Duffing experiment is a discontinuously engaged, cross-separatrix
transfer-and-release problem whose success depends strongly on one accurately
directed engagement impulse. Most Duffing approximations rotate the lifted
input vector sufficiently to weaken, reverse, or destabilize the feedback when
the gain designed with the approximation acts on the true plant.

The most important distinction is therefore not raw matrix-norm error by
itself. It is the combination of:

1. input-vector direction and polarity;
2. target displacement in lifted coordinates;
3. truth-coupled closed-loop stability margin;
4. whether the controlled data cover the actual engagement and transfer tube;
5. whether the retained eigenfunctions provide an approximately closed
   representation of `B_psi` over that tube; and
6. whether integration is robust to the engagement discontinuity and nearby
   off-path states.

## Compared retained experiments

### RPE positive case

- Minimal one-nonredundant-complex-eigenfunction SINDy KROM.
- Continuous eigenvalue:

  ```text
  lambda = -0.21061644578341268 + 4.818634561556464 i
  ```

- Dimensionless state and pressure input.
- Initial state `(1,0)`.
- Target `rpe_multisine_equilibrium_moderate_identification`:

  ```text
  r_ref(tau) = 1 + 0.09 sin(4 tau) - 0.12 sin(3 tau)
  v_ref(tau) =     0.36 cos(4 tau) - 0.36 cos(3 tau)
  ```

- The reference begins exactly at equilibrium and moves smoothly.
- Continuous-RHS KLQR, SciPy BDF, `Q=diag(5000,5000)`, `R=1e-4`.
- No input cap, feedforward term, or target shaping.
- Exact, constant, complex-bilinear, conjugate-complete bilinear, and
  persistent-data approximations all completed successfully.

### Duffing negative and inconclusive cases

- Parent degree-4 SINDy KROM derived from 30 free BDF trajectories.
- Candidate sub-KROMs containing one, two, or three nonredundant complex
  eigenfunctions.
- Dimensionless Duffing system:

  ```text
  x_dot = v
  v_dot = -0.5 v + x - 4 x^3 + u
  ```

- Initial state `(0.5,0)` at the right equilibrium.
- Constant target `(-0.25,0)` on the opposite side of the separatrix.
- Control inactive until `t=1`, active through `t=10`, then released through
  `t=20` so the state should fall into the left equilibrium.
- Continuous-RHS KLQR, SciPy BDF, `Q=12500 I`, `R=0.01`.
- No input cap or feedforward term.
- The bounded campaign retained 24 durable extraction-regression-control
  attempts and 100 approximate-map controller cases.

## Quantitative comparison

| Property | Successful RPE | Representative exact-controllable Duffing KROMs |
|---|---:|---:|
| Maximum initial/target lifted displacement | `0.0151063` | `0.399421` to `0.585054` |
| Directional excursion of exact `B_psi` from engagement direction | at most `25.57 deg` | approximately `154` to `180 deg` |
| Best standard target/path approximation error | `0.08623` | ordinary bilinear medians approximately `0.98` to `1.05` |
| Constant approximation error | target-path `0.30384` | global medians approximately `12.6` to `19.4` |
| Stored exact-control peak | `1.40209` | commonly `165` to `761` for successful candidates |
| Recorded internal controller-call peak | modest | as large as approximately `1.5e6` to `1.8e6` |
| Typical approximate/true frozen-state pole behavior | strongly stable | often unstable over most of the successful exact path |
| Control objective | continuing same-region tracking | threshold basin crossing followed by release |

The RPE lifted target excursion is roughly 26 to 39 times smaller than the
Duffing engagement error. RPE therefore starts with zero tracking error and
continues through small corrections. Duffing presents KLQR with the full
cross-well error as soon as the active window opens.

## RPE outcome

| Input representation | Radius RMSE, second half | Velocity RMSE, second half | Target-path `B_psi` error | Direction error, p95 | Relative control-history difference from exact |
|---|---:|---:|---:|---:|---:|
| Exact analytic | `0.00121273` | `0.00412142` | `0` | approximately `0 deg` | `0` |
| Pulse-data constant | `0.00118033` | `0.00397416` | `0.303837` | `24.7409 deg` | `0.003986` |
| Pulse-data complex bilinear | `0.00117205` | `0.00393333` | `0.210358` | `21.1021 deg` | `0.003627` |
| Pulse-data conjugate-complete bilinear | `0.00118979` | `0.00400489` | `0.0862298` | `4.07411 deg` | `0.001291` |
| Persistent-data constant | `0.00114386` | `0.00381339` | `0.398742` | `29.4127 deg` | `0.004276` |
| Persistent-data bilinear | `0.00111899` | `0.00372411` | `0.391077` | `46.7299 deg` | `0.006960` |

All five approximate maps met the declared tracking gate. Their realized
control histories differed from exact control by less than one percent, and
their state histories differed by less than two tenths of one percent in
relative L2 norm.

This demonstrates that a 30 to 40 percent Frobenius error need not cause a
large closed-loop error when it is dominated by an input-vector rescaling that
CARE can compensate. Direction and stability margin matter more than the raw
norm alone.

## Duffing campaign outcome

Across the 100 approximate-map controller cases:

| Status | Count | Meaning |
|---|---:|---|
| Complete but unsuccessful | `25` | Solver completed, but target/release gates failed |
| Numerical failure | `28` | Usually BDF reported that the required step was below floating-point spacing |
| Wall-time timeout | `47` | No completed scientific outcome was obtained before the campaign limit |

The best completed constant-map release error was `0.447450`; the best
completed bilinear-map release error was `0.534932`. The required release
error was at most `0.02`.

Eighteen exact-map cases completed, and eleven retained KROM selections passed
the full exact-control gate. The campaign therefore separates two issues:

1. some eigenfunction subsets are unsuitable even with exact `B_psi`; and
2. for subsets that work with exact `B_psi`, the tested approximations usually
   destroy the successful feedback geometry.

## Engagement-command evidence

The transfer is controlled mainly by a highly localized engagement impulse.
Most of the remaining active interval uses a small input; for successful exact
cases the 99th percentile is generally below approximately `1.7` despite the
hundreds-scale stored peak.

Representative static engagement calculations were:

| Retained Duffing modes | Exact command | Constant approximation | Bilinear dynamics | Direct bilinear map |
|---|---:|---:|---:|---:|
| `6` | `-380.3` | `-139.8` | `+281.4` | `+286.5` |
| `6,8` | `-663.4` | `+236.4` | `+874.8` | `+881.2` |
| `8,9` | `-638.0` | `-185.1` | `-35.7` | `-65.3` |

For modes `6` and `6,8`, several approximations reverse the sign of the
initial action. For modes `8,9`, the sign is retained but the action is much
too weak. In either case the approximation fails at the event that determines
whether the state crosses the separatrix.

## Input-vector geometry and truth-coupled poles

For a candidate input map, KLQR computes a gain using the candidate matrix:

```text
K_candidate = R^-1 B_candidate^T P_candidate.
```

The physical plant nevertheless responds through analytic chain-rule
`B_psi`. A relevant local diagnostic is therefore not only the nominal matrix

```text
A - B_candidate K_candidate,
```

but also the truth-coupled matrix

```text
A - B_exact K_candidate.
```

For representative Duffing candidates, the ordinary global constant and
bilinear maps had median direction errors from approximately `66` through
`134` degrees, with maxima approaching `180` degrees. Their truth-coupled
frozen systems were unstable over roughly 69 to 100 percent of the exact
successful path, and maximum real pole parts reached order `1e2` to `1e3`.

An input-direction error above 90 degrees means the fitted lifted input vector
is more opposed than aligned with the true vector. CARE may stabilize the
model it was given while producing positive or badly misdirected feedback on
the true lifted dynamics.

For RPE, the pulse-data conjugate-complete map had a p95 direction error near
`4.07` degrees and a worst target-path truth-coupled pole near `-4.26`. Even
the coarser RPE approximations generally retained sufficient direction and
margin for the nonlinear integration to correct their errors.

## Eigenfunction geometry

The successful one-mode RPE and Duffing KROMs both use a single complex
eigenfunction and its implicit conjugate in split-real coordinates. Their
state-to-lift geometry differs sharply over the tested control paths.

For the RPE target:

- exact `B_psi` stays within approximately `25.57` degrees of its engagement
  direction while oscillating repeatedly;
- the one-mode lift Jacobian does not change orientation; and
- the controller remains inside a local region surrounding equilibrium.

For Duffing mode `6`:

- exact `B_psi` reverses by almost `180` degrees;
- the two-by-two lift Jacobian determinant changes sign four times along the
  successful transfer;
- the Jacobian condition number reaches approximately `10356`; and
- the retained coordinate therefore folds or becomes nearly singular over
  the cross-well path.

A global affine map `B_psi(Psi)=M Psi+b` cannot faithfully describe that
folded input geometry.

The two-mode Duffing embeddings are substantially better conditioned, so
folding is not the entire explanation. Their exact input vectors still rotate
by more than 150 degrees, and the derivative of the retained eigenfunctions
does not lie in the tested low-order affine span of those retained coordinates.
Adding modes improved the state embedding but did not automatically supply
input-map closure.

## Controlled-data comparison

### RPE identification data

The successful RPE candidate set included symmetric opposite-pulse
experiments at 600 anchor states:

- 400 anchors in a narrow tube around the intended target;
- 200 anchors in a supporting rectangle; and
- paired `+1` and `-1` dimensionless pressure commands at each anchor.

The estimator

```text
Bpsi_hat(x) = [Psi(x_plus(dt))-Psi(x_minus(dt))]
              / [(u_plus-u_minus) dt]
```

cancels autonomous evolution to leading order. It recovered held-out analytic
`B_psi` samples with relative Frobenius error `0.0107513`. The remaining error
was therefore mostly the capacity floor of the constant or bilinear family.

### Duffing identification data

The retained Duffing random-forced dataset contains:

- 100 trajectories;
- 20,000 transitions at `dt=0.01`;
- states spanning approximately `[-1,1]^2`; and
- canonical input spanning approximately `[-0.5,0.5]`.

It is globally broad but not aligned with the transfer endpoints. In the
retained realization, no controlled-state samples fall within Euclidean radius
`0.05` of either `(0.5,0)` or `(-0.25,0)`. Only approximately 23 percent fall
inside a broad rectangle surrounding the successful transfer corridor.

The continuous lifted-residual regression also mixes control response with
autonomous SINDy eigenfunction residual. The Korda--Mezic constant-output
regression for the minimal Duffing KROM depends on an observable
reconstruction `C` that failed validation and had a condition number near
`5.3e16`.

These issues damage the standard controlled-data fits, especially the
constant matrix. They are not the complete cause: direct oracle fits made
against analytic exact `B_psi` also retained order-one error and failed. That
isolates insufficient approximation-family capacity as a primary limitation.

## Hypotheses ruled out or weakened

### The Duffing SINDy eigenfunctions are simply much worse

This is not supported by the free-trajectory validation:

| Quantity | RPE retained mode | Duffing retained mode `6` |
|---|---:|---:|
| One-step relative error | `0.005324` | `0.008771` |
| Horizon `0.1` relative error | `0.41185` for RPE at its retained `0.1` audit | `0.08511` for Duffing |
| Horizon `1.0` relative error | `0.70950` | `0.60123` |

The grids and systems differ, so these values are not a universal ranking, but
they are sufficient to reject a dramatic eigenfunction-fidelity gap as the
main explanation. Duffing is slightly better at the one-unit comparison.

### The Python matrix regression is algebraically wrong

Prior MATLAB parity work reproduced supplied MATLAB `bb`, `B1`, and `B2`
arrays to approximately `1e-13` relative error or better when both frameworks
received the same aligned lifted values, derivatives, inputs, and generator.
Furthermore, direct oracle map fits bypass controlled-data derivative
regression and still fail. The retained evidence therefore does not support a
general Python regression-algebra defect.

### Raw Frobenius error alone predicts controller success

RPE disproves this. Maps with 30 to 40 percent target-path error still produce
nearly exact trajectories because the error is sufficiently gain-compensable.
Duffing demonstrates the missing variables: direction, engagement action,
truth-coupled poles, and off-path robustness.

## Stability-tube comparison

For the RPE pulse-data conjugate-complete map, perturbing the target path by a
Euclidean state radius of `0.1` retained a worst truth-coupled pole near
`-3.78`. Instability appeared only sparsely at substantially larger probes.

For the most promising Duffing arc-length conjugate-complete candidates:

| Tube radius | Arc-7 p95 direction error | Arc-7 unstable fraction | Arc-9 p95 direction error | Arc-9 unstable fraction |
|---:|---:|---:|---:|---:|
| `0.025` | `14.37 deg` | `0` | `16.07 deg` | `0` |
| `0.05` | `24.39 deg` | `0` | `26.28 deg` | `0` |
| `0.10` | `46.10 deg` | `0` but nearly marginal | `48.05 deg` | approximately `0.001` |
| `0.20` | `118.44 deg` | `0.353` | `120.00 deg` | `0.358` |
| `0.30` | `149.33 deg` | `0.470` | `149.38 deg` | `0.474` |

The RPE controller has considerably more stability reserve. The Duffing arc
fits are accurate on the successful path but become fragile under moderate
departure from it.

## Important timeout qualification

The final seven-node and nine-node conjugate-complete arc fits should not be
classified as demonstrated physical-control failures. Each was terminated by
the campaign's 60-second wall-time limit before a completed solution existed.

Their static path diagnostics were substantially better than those of the
ordinary failed fits:

- path p95 direction errors approximately `5` to `6` degrees;
- path p95 relative errors approximately `16` to `17` percent;
- stable truth-coupled frozen poles on the exact path; and
- stability throughout a narrow tube of radius `0.025`.

They nevertheless remain numerically and off-path fragile. The exact modes
`8,9` solution itself required approximately 30 seconds and recorded internal
controller-call magnitudes near `1.76e6`. A 60-second timeout can therefore be
reached through stiffness and repeated CARE evaluation without proving that
the resulting physical trajectory would fail.

These candidates are **inconclusive** and merit a longer, instrumented rerun.

## Recommended clean isolation experiment

The next experiment should preserve the KLQR law and use the exact-successful
two-mode Duffing KROM with retained indices `(8,9)`. This model has a much
better-conditioned state embedding than the one-mode model.

1. Generate Duffing symmetric opposite-pulse data analogous to the successful
   RPE experiment.
2. Place anchors at the engagement state, target, successful exact transfer
   path, and a declared surrounding tube.
3. Retain a held-out subset and exact command histories.
4. Fit the conjugate-complete real-linear input map first.
5. Before full control, require:

   - correct engagement-command sign;
   - comparable engagement-command magnitude;
   - small input-vector direction error;
   - adequate truth-coupled pole margin over the entire path tube;
   - smooth CARE gain variation under bounded perturbations; and
   - acceptable oracle capacity floor.

6. If the oracle fit fails those gates, do not generate more data of the same
   kind. Expand the `B_psi` representation beyond a global affine function of
   the retained eigenfunctions.
7. Rerun the arc-7 and arc-9 candidates with a longer wall-time allowance and
   progress/partial-trajectory telemetry.

This sequence does not require changing the established control law
`u=-K(Psi-Psi_ref)`. It separates data-estimation error, representation error,
gain sensitivity, and integration overhead before another expensive campaign.

## Reproducibility and retained evidence

- RPE control-map comparison report:
  `docs/case-studies/RPE_APPROXIMATE_BPSI_KLQR_COMPARISON_2026-08-29.md`
- RPE designed controlled-data report:
  `docs/case-studies/RPE_DESIGNED_CONTROL_IDENTIFICATION_DATA_2026-08-29.md`
- RPE machine-readable comparison:
  `viewer_outputs/rpe-minimal-krom-control-map-comparison.json`
- RPE Viewer comparison:
  `viewer_outputs/rpe-minimal-krom-control-map-comparison.html`
- RPE error-map Viewer:
  `viewer_outputs/rpe-minimal-krom-control-map-error-comparison.html`
- Duffing campaign report:
  `docs/case-studies/DUFFING_BILINEAR_KROM_SEARCH_2026-08-29.md`
- Duffing durable campaign ledger:
  `archives/experiments/duffing_bilinear_krom_search_20260829/campaign-report.json`
- Duffing campaign Viewer:
  `viewer_outputs/duffing-bilinear-krom-search-final.html`
- Accepted exact Duffing reference:
  `docs/case-studies/DUFFING_MINIMAL_EXACT_KLQR_TRANSFER_RELEASE_2026-08-27.md`
- Duffing fresh one-mode refit:
  `docs/case-studies/DUFFING_MINIMAL_FRESH_CONTROL_REFIT_2026-08-29.md`
- Duffing two-mode refit:
  `docs/case-studies/DUFFING_TWO_MODE_FRESH_REFIT_KLQR_2026-08-29.md`

## Final conclusion

RPE succeeded because the selected minimal lift, target, controlled-data
design, and KLQR tuning together created a forgiving local problem. Its input
maps are not exact, but their errors are mostly compatible with Riccati gain
compensation and the true closed loop retains substantial stability margin.

Duffing failed because the chosen control target converts input-map error into
a basin-selection error at one abrupt engagement event. The exact input vector
rotates or folds strongly across the transfer, ordinary affine approximations
do not close on the retained eigenfunctions, the random controlled data do not
sample the decisive endpoints, and the resulting gains are frequently
destabilizing when coupled to the true plant. The remaining promising
path-fitted candidates are too fragile and computationally unresolved to call
successful, but they identify a precise next experiment rather than a need to
replace KLQR or copy a different control law.
