# Duffing optimized-eigenvalue / SINDy control-regression case study

Date: 2026-08-28  
Status: reproducible candidate experiment; no change to the default exact-control workflow.

## Purpose

This experiment exercises the newly separated pipeline by deriving discrete eigenvalues with trajectory-boundary optimization, passing those values unchanged into SINDy to construct analytic polynomial eigenfunctions, reconstructing the full state observable, and comparing constant and bilinear approximations with exact state-dependent `B_psi(x)`.

It is a diagnostic of representation quality. It is not evidence that approximate control matrices should replace exact `B_psi` by default.

## Reproduction

Run from the DynSysLab root:

```powershell
python scripts/run_duffing_optimized_sindy_control_regression.py
```

The script is a reusable CLI and writes a JSON report, compressed numerical arrays, a typed control-error archive, a static heat-map figure, and a native Viewer page.

## Free-data and spectral pipeline

- System: retained DynSysLab Duffing system, `delta=0.5`, `alpha=1`, `beta=4`, `input_gain=1`.
- Free ensemble: 30 trajectories, 801 samples per trajectory, analysis interval `delta_t=0.01`.
- Observable passed to optimization: full state `(x, xdot)`.
- Initial eigenvalue groups: retained Korda/Mezic candidate artifact.
- Optimized full spectrum: 20 modes, divided into ten modes per observable output.
- Conjugate structure: preserved.
- KROM after conjugate reduction: 10 retained eigenfunctions.
- Eigenfunction method: normalized degree-four SINDy polynomial dictionary.
- Eigenfunction representation: analytic for all retained modes.
- Gradient backend: analytic for all retained modes.

Optimization converged for both observable components:

| Observable | Iterations | Relative trajectory residual |
|---|---:|---:|
| `x` | 29 | 0.0052676 |
| `xdot` | 24 | 0.0212608 |

The SINDy stage consumed the optimized spectrum. It did not re-estimate the eigenvalues.

## Observable reconstruction `C`

The full state was reconstructed from the nonredundant complex KROM by least squares with a deterministic interleaved 20% holdout. Validation produced `R² = 1.0` for both state components, with RMSEs of approximately `1.78e-14` for `x` and `1.03e-13` for `xdot`.

This nearly exact reconstruction is expected because the degree-four full-state polynomial dictionary contains the coordinate observables. It does **not** imply that the lifted regression is well conditioned: the reconstruction fit condition number was approximately `1.72e17`.

## Controlled dataset

- 100 accepted trajectories.
- 200 transitions per trajectory; 20,000 total transitions.
- `delta_t=0.01`.
- One held input over each fixed classical RK4 interval.
- Source commands span `[-1, 1]` and are converted by `u_dynsyslab = 0.5 u_source` to reconcile the source Duffing convention with DynSysLab's canonical equation.
- The exact source control array is retained with SHA-256 `bfcf2e6d87b00ffe98e51176f1cbaf370bc38b2ae565db52c1e0edc25940fb24`.

## Regressed control representations

### Constant matrix

The multistep output regression produced a full-rank `40000 x 20` design matrix, but its condition number was approximately `2.91e15`. Output RMSE was `0.7198`. Across the free-training rectangle, relative error versus exact `B_psi` had median `463.6`, mean `522.6`, 95th percentile `1104.0`, and maximum `3539.4`.

This approximation is unusable for this KROM. The failure is a substantive diagnostic: nominal rank did not protect the regression from severe numerical ill-conditioning.

### Bilinear matrix from controlled dynamics

The `B1 Psi(x) + B2` lifted-dynamics regression used a `20000 x 11` design matrix with condition number `9.28e6`. Lifted RMSE was `0.8023`. Relative error had median `0.5684`, mean `0.6002`, 95th percentile `0.8171`, and maximum `1.0042`.

This representation is much closer to exact `B_psi` than the constant matrix and is comparatively uniform over the rectangle, but its typical error remains large.

### Bilinear matrix fitted directly to exact `B_psi`

The direct pointwise map regression used 4,000 training states, a `4000 x 11` design matrix, and condition number `3062`. Its fit residual relative to exact map samples was `0.2524`. On the diagnostic grid, relative error had median `0.3391`, mean `0.4058`, 95th percentile `0.9042`, and maximum `1.1181`.

This was the best average and median approximation. Its higher worst-case value indicates localized poor regions, particularly where the exact map is weak or the polynomial representation is stressed. It still does not justify replacing exact `B_psi` for control.

## Error-map definition and domain

The rectangular diagnostic domain is the axis-aligned envelope of the free eigenfunction-training states:

- `x` in `[-0.9926372474, 0.9848858579]`
- `xdot` in `[-0.9964836952, 0.9967542211]`
- grid: `61 x 61`

Every point compares the approximation to the analytic exact map using the Frobenius matrix norm. Both absolute and relative fields are retained. The plots use logarithmic relative-error color scales because the constant representation spans several orders of magnitude.

## Interpretation

The experiment validates the new plumbing:

- eigenvalue extraction and eigenfunction construction are independently selectable;
- optimized eigenvalues can feed SINDy without being silently replaced;
- analytic function metadata drives analytic `nabla psi` evaluation;
- `C` reconstruction, constant regression, bilinear regression, exact-map regression, and phase-space error fields all operate on the same KROM contract;
- Viewer can load and render the typed error archive directly.

The numerical outcome is deliberately mixed. The composable pipeline worked, while the constant-matrix model did not. The direct bilinear map is promising as an approximation study, but exact state-dependent `B_psi` remains the standard controller coupling.

## Retained artifacts

- `archives/experiments/duffing_optimized_sindy_control_regression/report.json`
- `archives/experiments/duffing_optimized_sindy_control_regression/duffing_optimized_sindy_control_regression.npz`
- `archives/experiments/duffing_optimized_sindy_control_regression/control_approximation_error_fields.npz`
- `archives/experiments/duffing_optimized_sindy_control_regression/duffing_optimized_sindy_control_matrix_error_heatmaps.png`
- `viewer_outputs/duffing-optimized-sindy-control-matrix-errors.html`

## Timing on the recorded run

| Stage | Seconds |
|---|---:|
| Free training | 3.33 |
| Optimized eigenvalues + SINDy | 2.12 |
| Observable reconstruction | 25.63 |
| Constant regression | 63.43 |
| Bilinear dynamics regression | 20.78 |
| Bilinear exact-map regression | 10.65 |
| Error fields | 21.62 |
| Total | 148.77 |

The constant regression is the dominant cost and should be a future target for numerical reformulation or lower-level acceleration if this candidate method is retained.
