# RPE KLQR with exact and approximated control maps

**Date:** 2026-08-29  
**Status:** matched comparison completed  
**Target:** `rpe_multisine_equilibrium_moderate_identification`  
**Invariant contract:** same minimal KROM, equilibrium start `(1,0)`, dimensionless dynamics/control, continuous-RHS KLQR, BDF, unchanged Q/R, no cap, no feedforward, no target shaping

## Results

| Control map | Status | r RMSE, second half | v RMSE, second half | RMS u* | Target-path Bpsi error | B direction angle p95 | Control difference from exact |
|---|---|---:|---:|---:|---:|---:|---:|
| exact analytic B_psi | complete | 0.00121273 | 0.00412142 | 0.455482 | 0 | 1.20742e-06 deg | 0 |
| pulse-data constant | complete | 0.00118033 | 0.00397416 | 0.456261 | 0.303837 | 24.7409 deg | 0.00398624 |
| pulse-data complex bilinear | complete | 0.00117205 | 0.00393333 | 0.456299 | 0.210358 | 21.1021 deg | 0.00362712 |
| pulse-data conjugate-complete bilinear | complete | 0.00118979 | 0.00400489 | 0.455854 | 0.0862298 | 4.07411 deg | 0.00129146 |
| persistent-data constant | complete | 0.00114386 | 0.00381339 | 0.456695 | 0.398742 | 29.4127 deg | 0.00427602 |
| persistent-data bilinear | complete | 0.00111899 | 0.00372411 | 0.457825 | 0.391077 | 46.7299 deg | 0.00695951 |

## Interpretation

6 of 6 matched solves completed. Approximate maps meeting the predeclared two-times-exact tracking gate: pulse-data constant, pulse-data complex bilinear, pulse-data conjugate-complete bilinear, persistent-data constant, persistent-data bilinear.

All five approximations produce essentially the same closed-loop response as exact `B_psi`: their realized control histories differ by less than one percent and their state histories differ by less than two tenths of one percent in relative L2 norm. The phase-space matrix error therefore does not translate directly into tracking error in this case.

The diagnostic reason is structural. This minimal KROM has one complex eigenfunction, represented by a two-state real block, and one scalar input. CARE is recomputed from the selected local input vector. Much of the approximation error is a rescaling of that vector, which the Riccati gain largely cancels. Direction error is more important than raw matrix-norm error. The conjugate-complete bilinear fit preserves direction most closely and produces the closest trajectory to exact control.

The persistent-data bilinear map is a warning against treating this result as universal: it reaches much larger direction errors and its gain coupled to exact `B_psi` has locally unstable frozen-state poles over a small part of the target path, although the integrated nonlinear trajectory remains stable and accurate. Closed-loop validation remains mandatory.

The aggregate phase-space heat-map error, pathwise error, input-vector geometry, nominal frozen-state poles, truth-coupled frozen-state poles, and realized trajectory differences are retained in the machine-readable summary.

## Artifacts

- Viewer control comparison: `viewer_outputs/rpe-minimal-krom-control-map-comparison.html`
- Viewer error maps: `viewer_outputs/rpe-minimal-krom-control-map-error-comparison.html`
- Machine-readable summary: `viewer_outputs/rpe-minimal-krom-control-map-comparison.json`
- KROM error bundle: `archives/experiments/rpe_control_identification_designed/rpe-minimal-krom-error-evaluations.npz`
