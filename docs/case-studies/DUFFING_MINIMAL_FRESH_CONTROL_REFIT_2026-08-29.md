# Duffing Minimal KROM: Fresh Control-Matrix Refit

Date: 2026-08-29  
Status: completed reproducibility rerun

## Question

Would fitting the constant and bilinear input maps directly on the accepted
one-eigenfunction Duffing KROM, rather than restricting the corresponding
20-mode parent matrices, recover the successful exact-`B_psi` KLQR transfer?

## Inputs held fixed

- Accepted saved one-eigenfunction KROM; no spectrum or eigenfunction refit.
- Continuous eigenvalue:
  `-0.4911019907065521 + 1.0449812139069314 i`.
- Same retained random-forced Duffing control source:
  `korda_mezic_duffing_random_forced_full_v3.mat`.
- Same 100 trajectories, 200 transitions per trajectory, and canonical
  DynSysLab input scaling.
- Ridge parameter `1e-12`.
- Constant method: Korda--Mezic multistep observable-output regression with the
  saved one-mode observable reconstruction `C`.
- Bilinear dynamics method: continuous lifted-dynamics residual regression.
- Bilinear map method: direct least-squares approximation of analytic
  `B_psi(x)` on every fifth controlled-state sample.

This is a genuine reduced-coordinate refit. None of the matrices were selected
from the 20-mode parent fit.

## Freshly fitted matrices

### Constant

```text
B = [[2.7981481638409096 + 13.418950080738032 i]]
Bd = [[0.02721354788587258 + 0.1340038241701857 i]]
```

### Bilinear lifted-dynamics fit

```text
B1 = [[[-0.2355918763463184 - 0.05272876764245868 i]]]
B2 = [[ 0.12495768227120568 + 0.08906604516391588 i]]
```

### Direct bilinear control-map fit

```text
B1 = [[[-0.4890439517530329 - 0.08866487745010618 i]]]
B2 = [[ 0.24327894986701135 + 0.17344924103261639 i]]
```

Every fresh matrix was bit-for-bit identical to the corresponding matrix from
the earlier 2026-08-28 one-mode regression. The maximum absolute difference was
zero for all six stored arrays. This is expected because the algorithms,
controlled data, KROM, sampling, and regularization are deterministic and
unchanged.

## Relation to MATLAB parity work

The MATLAB comparison verified the array-level regression algebra: supplied
with the same aligned lifted values, derivatives, input samples, and generator,
Python reproduced MATLAB's `bb`, `B1`, and `B2` to approximately `1e-13`
relative error or better.

That work added explicit compatibility boundaries and corrected or documented
data ordering, input scaling, conjugation, and continuous/discrete conversion.
It did not replace the canonical regression used here. Two legacy MATLAB
inconsistencies remain opt-in compatibility cases rather than production
defaults: a different Duffing cubic coefficient in one controlled derivative
path, and conjugated lifted values paired with nonconjugated derivatives.

## Approximation errors over the free-training rectangle

| Refit | Median relative Frobenius error | Mean | 95th percentile | Maximum |
|---|---:|---:|---:|---:|
| Constant | 12.5652 | 18.1258 | 44.2157 | 4159.41 |
| Bilinear dynamics | 0.999780 | 1.02806 | 1.21042 | 30.6581 |
| Direct bilinear map | 1.02350 | 1.09836 | 1.47949 | 58.7244 |

The reduced refit improves these global errors relative to principal
restriction of the parent bilinear matrices, but errors near order one remain
too large to consider either bilinear representation an accurate surrogate for
exact `B_psi`.

## Unchanged KLQR retry

- initial state `(0.5,0)`;
- target `(-0.25,0)`;
- active interval `[1,10]`, release at `10`, final time `20`;
- `Q = 12500 I_2`, `R = [[0.01]]`;
- no input cap;
- `u = -K(Psi(x)-Psi(target))`;
- dimensionless dynamics and input;
- continuous-RHS control evaluation;
- SciPy BDF with `rtol=1e-8`, `atol=1e-10`, `max_step=0.01`;
- Viewer output step `0.0001`.

| Representation | Status | Release state | Target error at release | Peak `|u|` |
|---|---|---:|---:|---:|
| Fresh constant | Complete | `(0.3448926,-0.0003143)` | `0.594893` | `119.172` |
| Fresh bilinear dynamics | Failed | unavailable | unavailable | unavailable |
| Fresh direct bilinear map | Failed | unavailable | unavailable | unavailable |

Both bilinear solves terminated with SciPy's message:
`Required step size is less than spacing between numbers.` The constant case
remained in the right basin and relaxed back toward the right equilibrium after
release. These results exactly reproduce the earlier one-mode-refit control
outcome; intervening extraction and controller-boundary work did not change it.

## Interpretation

Fresh regression is the correct way to obtain a self-contained one-mode
bilinear approximation, but it cannot create representational capacity that the
selected coordinate does not possess. In this case, both direct reduced fits
still approximate the exact input map poorly enough to generate an ineffective
or singular closed loop. The exact analytic `B_psi = nabla(Psi) B` remains the
successful and authoritative one-mode control representation.

## Artifacts

- Regression report:
  `archives/experiments/duffing_minimal_fresh_refit_20260829/report.json`.
- KLQR report:
  `archives/experiments/duffing_minimal_fresh_refit_klqr_20260829/report.json`.
- Viewer error maps:
  `viewer_outputs/duffing-minimal-krom-fresh-refit-control-matrix-errors.html`.
- Viewer KLQR result:
  `viewer_outputs/duffing-minimal-krom-fresh-refit-klqr.html`.

