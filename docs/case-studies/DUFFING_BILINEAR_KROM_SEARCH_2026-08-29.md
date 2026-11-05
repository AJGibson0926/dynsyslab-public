# Duffing bilinear KROM control search

## Bottom line

The bounded campaign completed **24 durable attempts** and found **no successful bilinear control approximation**. One earlier pilot was interrupted before it entered the durable ledger; counting it conservatively gives 25 physical extraction-regression-control sequences, exactly the user-specified ceiling.

This is a negative result about the retained approximation families, not about KLQR or the selected KROMs in general. Exact analytic `B_psi` passed the unchanged transfer/release gate in attempts 1, 2, 5, 8, 10, 19, 20, 21, 22, 23, 24. The same KROMs failed when their input maps were replaced by a constant matrix, a controlled-dynamics bilinear fit, a global pointwise bilinear map, a corridor-local map, a conjugate-complete real-linear map, or the path-conditioned candidates.

The final native Viewer dashboard is [duffing-bilinear-krom-search-final.html](viewer_outputs/duffing-bilinear-krom-search-final.html). All 24 per-attempt error-map dashboards were generated under `viewer_outputs/duffing-bilinear-krom-search`.

## Frozen experiment contract

- Existing Duffing degree-4 SINDy KROM and existing 100-member random-forced controlled dataset; neither was regenerated or expanded.
- Initial state `(0.5, 0)`, constant target `(-0.25, 0)`, control active on `[1, 10]`, release at 10, final time 20.
- Dimensionless Duffing dynamics `x_dot=v`, `v_dot=-0.5 v + x - 4 x^3 + u`.
- KLQR law `u=-K(Psi(x)-Psi(target))`, `Q=12500 I`, `R=0.01`, no control cap, continuous RHS control evaluation.
- SciPy BDF, `rtol=1e-8`, `atol=1e-10`, `max_step=0.01`; screening output interval `0.001`.
- Success gate: release target error at most 0.02 and final left-equilibrium error at most 0.10.

## Tooling and efficiency

The campaign is available through `dynsyslab-duffing-bilinear-search` and the public `rom.control_campaign` API. It enforces a hard maximum of 25 attempts, writes an atomic resume ledger after every attempt, isolates every controller solve in a timed process, archives each extracted independent KROM and fitted control representation, and saves every full phase-space error field.

The Korda-Mezic constant-matrix multistep design was algebraically reduced from a nested convolution to an O(horizon) recurrence. Representative fresh constant fits now take roughly 0.4 seconds; pointwise bilinear fits take about 1-4 seconds. The long wall time comes from stiff or unstable closed-loop integration, not matrix regression.

## Attempt ledger

| # | retained KROM indices | exact | constant | bilinear dynamics | direct bilinear map | successful cases |
|---:|---|---|---|---|---|---|
| 1 | `6` | complete; release=0.000111593 | complete; release=0.594893 | failed | failed | exact |
| 2 | `6,8` | complete; release=7.17388e-05 | failed | complete; release=1.12765 | complete; release=1.12178 | exact |
| 3 | `6,7,8` | timeout | complete; release=0.523057 | timeout | timeout | none |
| 4 | `5,8` | timeout | complete; release=1.16897 | timeout | timeout | none |
| 5 | `6,9` | complete; release=8.96067e-05 | complete; release=0.918706 | timeout | timeout | exact |
| 6 | `5,6` | timeout | complete; release=1.38613 | timeout | timeout | none |
| 7 | `6,7` | timeout | complete; release=0.687271 | timeout | complete; release=1.10453 | none |
| 8 | `9` | complete; release=0.000788364 | complete; release=0.897467 | complete; release=0.941171 | timeout | exact |
| 9 | `8` | complete; release=1.0384 | complete; release=1.26338 | failed | failed | none |
| 10 | `8,9` | complete; release=7.22201e-05 | complete; release=0.44745 | timeout | timeout | exact |
| 11 | `7,8` | timeout | failed | timeout | timeout | none |
| 12 | `0` | complete; release=0.545948 | complete; release=1.14737 | failed | failed | none |
| 13 | `1` | complete; release=0.551921 | timeout | failed | failed | none |
| 14 | `5` | complete; release=0.594724 | timeout | failed | complete; release=0.597871 | none |
| 15 | `2` | complete; release=0.549841 | failed | failed | failed | none |
| 16 | `7` | timeout | complete; release=0.792147 | failed | failed | none |
| 17 | `4` | complete; release=0.549138 | timeout | failed | failed | none |
| 18 | `3` | complete; release=0.541703 | complete; release=0.541502 | failed | timeout | none |
| 19 | `6` | complete; release=0.000111593 | complete; release=0.594893 | timeout | failed | exact |
| 20 | `9` | complete; release=0.000788364 | complete; release=0.897467 | timeout | timeout | exact |
| 21 | `6,8` | complete; release=7.17388e-05 | failed | timeout | timeout | exact |
| 22 | `8,9` | complete; release=7.22201e-05 | complete; release=0.44745 | timeout | timeout | exact |
| 23 | `6` | complete; release=0.000111593 | complete; release=0.594893 | timeout | failed | exact |
| 24 | `8,9` | complete; release=7.22201e-05 | complete; release=0.44745 | timeout | timeout | exact |

## Evidence and diagnosis

The smallest release error among completed bilinear solves was `0.534932` for attempt 21 (`bilinear-conjugate-map`), still far outside the 0.02 success gate.

1. **Exact controllability is a necessary independent gate.** Several subsets fail even with exact `B_psi`; those KROMs cannot be rescued by a better regression. Modes 6 and 9 and the pairs `(6,8)`, `(6,9)`, and `(8,9)` are exact-controllable for this target.
2. **Global affine closure is absent.** Across exact-controllable KROMs, global bilinear-map median relative errors remain approximately order one. This means `B_psi(x)` is not in the affine span of the retained eigenfunctions over the training rectangle.
3. **Adding a conjugate term is insufficient.** The candidate `M Psi + N conjugate(Psi) + b` correctly spans a general real-linear map on split coordinates, but its global median errors also remain near one and its feedback solves fail.
4. **Naive locality can make the global map worse.** The rectangular target-corridor dataset is nonuniform and highly correlated; local fits often have worse conditioning and larger out-of-corridor errors.
5. **Time-weighted path fitting hid the engagement defect.** For mode 6 the reference-path median relative error fell to about 0.059, but the error at engagement was about 2.69 because most samples represented the long target dwell. Both the path and path-tube controllers timed out.
6. **Arc-length collocation improved pathwise error but not closed-loop robustness.** For `(8,9)`, seven- and nine-node conjugate-complete arc fits reduced reference-path maximum relative error to roughly 0.18, yet all closed-loop solves timed out. Once the approximate controller leaves the exact reference path, its errors and CARE conditioning amplify each other.
7. **Regression error alone is not a sufficient selection score.** Useful pre-screens must combine exact-control success, target displacement in lift space, controllability/CARE conditioning, engagement-point error, path-tube error, and sensitivity of the gain to perturbations in `B_psi`.

## Selection criteria learned

A future candidate should be admitted only if it passes these gates in order: (a) exact `B_psi` control succeeds; (b) local controllability rank is complete with a usable smallest singular value; (c) the fitted map is accurate at the initial state and target, not merely in aggregate; (d) path-tube error stays small; (e) CARE solutions and gains vary smoothly under bounded map perturbations; and only then (f) full nonlinear integration is run.

The polynomial closure score, eigenvalue frequency spread, and global median map error were individually weak predictors. Mode 8 demonstrates why: it fails alone under exact control but succeeds when paired with mode 9, so complementary target geometry matters. Conversely, adding modes can reduce controllability conditioning enough to make exact control impractical.

## Considerations for later work (not performed)

- Expand the controlled dataset with better state-space coverage and persistent excitation, especially around the transfer tube and engagement state, while preserving a held-out validation set.
- Build a larger or control-aware parent KROM whose eigenfunction span is closer to derivative closure under the Duffing input vector field. The derivative of a degree-4 polynomial eigenfunction is degree 3, but arbitrary small subsets of learned eigenfunctions need not span those derivatives.
- Seek eigenfunctions with different eigenvalues/frequencies only after exact-control and conditioning gates; frequency diversity by itself did not predict success here.
- Consider regularized local bases, orthonormalized lift coordinates, and gain-sensitivity-aware regression objectives. Do not infer that simply adding modes will help: larger exact KLQR systems became much more expensive and sometimes ill-conditioned.
- Preserve exact `B_psi` as the objectively correct default. Constant and bilinear maps remain optional approximations with explicit validation metadata.

## Reproduction

```text
dynsyslab-duffing-bilinear-search --dry-run --max-attempts 24
dynsyslab-duffing-bilinear-search --resume --max-attempts 24 --case-timeout-seconds 60
dynsyslab-duffing-bilinear-report
```

Durable campaign record: `archives/experiments/duffing_bilinear_krom_search_20260829/campaign-report.json`

Every attempt directory contains its independent KROM archive references, fitted matrix arrays, controller records, and `control-map-errors.npz`. Failed solves are retained as explicit statuses; they are never converted into missing or successful data.
