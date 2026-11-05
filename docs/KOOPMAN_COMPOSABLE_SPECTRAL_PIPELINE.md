# Composable Koopman spectral and control-identification pipeline

Status: candidate architecture implemented behind explicit native C++ services,
desktop workbenches, archive contracts, and an optional Python client. The
existing exact-control path remains the default. These additions are intended
to be removable or replaceable without altering solver or controller contracts.

## Design boundary

DynSysLab now treats the following as distinct stages:

1. **Eigenvalue extraction** estimates a spectrum from data and does not construct eigenfunctions.
2. **Eigenfunction construction** accepts an existing spectrum and constructs functions without changing the supplied eigenvalues.
3. **Observable reconstruction** identifies the Koopman-mode matrix `C` for a named observable.
4. **Control-map identification** optionally approximates the objectively state-dependent map `B_psi(x)` with constant or bilinear representations.
5. **Approximation diagnosis** compares those approximations with exact `B_psi(x)` on a stated training rectangle.

This separation permits, for example, optimized eigenvalues followed by SINDy eigenfunctions, H-DMD eigenvalues followed by SINDy, or externally supplied eigenvalues followed by fitted-boundary functions.

## Public tools

### Pipeline selection

`postprocessing.koopman.KoopmanPipelineSpecification` selects an eigenvalue method independently from an eigenfunction method. `run_composable_koopman_pipeline(...)` currently supports:

| Stage | Available methods |
|---|---|
| Eigenvalues | `hdmd`, `optimization`, `provided` |
| Eigenfunctions | `hdmd-interpolation`, `sindy`, `boundary-functions` |

The stage outputs retain the source spectrum and method artifact so provenance is not lost when stages are chained. `eigenvalues_from_hdmd(...)`, `provided_eigenvalues(...)`, and `construct_eigenfunctions(...)` also expose the stages directly.

### Optimization eigenvalues

`optimize_trajectory_eigenvalues(...)` minimizes the trajectory reconstruction residual after analytically eliminating the trajectory-specific boundary coefficients. It accepts one initialized discrete eigenvalue group per observable output and can preserve conjugate-pair structure. Its output is deliberately an eigenvalue-only `KoopmanSpectrum`; no implicit eigenfunction fit occurs.

The canonical implementation is now
`dynsyslab::core::optimizeTrajectoryEigenvalues(...)`. The native
Post-processor stage `eigenvalue_optimization` accepts a seed spectrum plus a
compatible free-trajectory ensemble and records both discrete and continuous
spectra. In `DynSysLab.exe`, **Configure eigenvalue optimization…** exposes all
required inputs. In `DynSysLabViewer.exe`, **Eigenvalue optimization** displays
the seed/optimized spectra and convergence evidence. No Python process is
started by any of these operations.

### Fitted-boundary functions

`fit_boundary_eigenfunctions(...)` accepts any compatible eigenvalue spectrum. It fits per-trajectory boundary values, propagates them by `mu**k`, and constructs scattered state-space interpolants. It also constructs the associated observable reconstruction: the modal coefficients `nu` become the columns of `C`.

### SINDy functions

`koopman_with_sindy(...)` already accepts an external spectrum. The composable pipeline now uses that contract directly. Discrete source eigenvalues are converted to continuous generator eigenvalues by `lambda = log(mu) / delta_t`, and the exact source eigenvalue objects remain traceable in metadata.

## Eigenfunction representation and gradients

Every `KoopmanEigenfunction` now identifies its representation:

- `analytic`: a callable analytic expression, including SINDy polynomials and user-supplied functions;
- `scattered-interpolation`: sampled state/value pairs represented by a nearest-neighbour inverse-distance interpolator;
- `sampled-values`: values exist, but no state-space evaluator exists;
- `unknown`.

`KoopmanEigenfunction.nabla_psi(...)` centralizes derivative dispatch:

- analytic representations use their analytic gradient;
- scattered interpolants use fourth-order finite differences of the interpolant;
- an explicitly requested analytic backend remains strict and fails when no analytic gradient exists.

`KROM.B_psi(...)` uses this representation-aware dispatch and applies the coordinate-observable Jacobian. Identity state observables therefore recover the ordinary state-gradient control map. Nontrivial observables apply the chain rule; their analytic Jacobian is used when available and a numerical Jacobian is the fallback.

## What H-DMD now means

An H-DMD result contains three linked products:

1. discrete/continuous eigenvalues;
2. modal/eigenfunction values on the training trajectories, now exposed as scattered-interpolation eigenfunctions;
3. the observable reconstruction matrix `C` obtained from the H-DMD mode values `nu`.

Because those interpolated eigenfunctions can evaluate `nabla psi` numerically, an H-DMD KROM can now be used with exact `B_psi(x)` when the system supplies a control matrix and the state lies within a meaningful interpolation domain. This is a capability statement, not a claim that interpolation-derived gradients are as accurate or robust as analytic SINDy gradients.

The standard operating procedure remains:

- H-DMD uses the first state-component observable (`x` for Duffing, radius for RPE);
- SINDy uses the complete dynamical state and may consume the H-DMD eigenvalues;
- exact state-dependent `B_psi(x)` is the default controller representation.

## Control identification

The regression routines operate on any KROM that has the data required by the selected method:

- `fit_constant_control_matrix(...)`: multistep output regression from a controlled dataset and a compatible `C` reconstruction;
- `fit_bilinear_control_matrix(...)`: continuous lifted-dynamics residual regression from controlled state/input trajectories;
- `fit_bilinear_control_map(...)`: direct regression of exact pointwise `B_psi(x)` onto `B1 Psi(x) + B2`, isolating representation error from trajectory-derivative error.

The exact representation is still authoritative. Constant and bilinear matrices are named, optional KROM representations and must be selected explicitly or attached by the regression result.

When a KROM is spectrally extracted, every retained control representation is
made self-contained. Constant maps restrict to `B[S, :]`. Bilinear maps use the
closed principal restriction `M[:, S, S]` and `b[S, :]`. The resulting KROM has
no live reference to its source KROM and needs only its own retained
eigenfunctions to evaluate the input map. This restriction is deterministic;
if its accuracy is inadequate, fitting a new input map directly on the reduced
coordinates is a separate, explicit operation and never occurs silently.

## Error fields and native Viewer display

`control_approximation_error_fields(...)` evaluates a stored constant or bilinear representation against exact `B_psi(x)` over a two-dimensional rectangle. At each grid point it computes

`absolute_error(x) = ||B_psi_exact(x) - B_psi_approx(x)||_F`

and

`relative_error(x) = absolute_error(x) / max(||B_psi_exact(x)||_F, relative_floor)`.

`save_control_approximation_error_fields(...)` and `load_control_approximation_error_fields(...)` provide the typed archive format `dynsyslab.control-approximation-error-fields/v1`.

Viewer has a native `control_error` plot kind and a `ControlErrorSource`. A request may display `relative`, `absolute`, or `exact_norm` values, with optional logarithmic coloring. Plotly and SVG renderers both produce equal-axis phase-space heat maps; Plotly adds hover values and a color bar.

Optional Python-client example:

```python
from rom import control_approximation_error_fields
from viewer import ControlErrorSource, PlotRequest, Viewer

fields = control_approximation_error_fields(
    krom,
    ([-1.0, -1.0], [1.0, 1.0]),
    grid_shape=(81, 81),
    representations=("constant", "bilinear"),
)

viewer = Viewer().attach(ControlErrorSource(fields))
page = viewer.build(PlotRequest("control_error", control_error_measure="relative"))
viewer.save_html(page, "viewer_outputs/control-errors.html")
```

## Important limitations

- Error maps are meaningful only over a declared region supported by training data. They do not validate extrapolation.
- Relative error becomes sensitive where exact `B_psi` is near zero; the archive records the denominator floor.
- Scattered-interpolation gradients may be noisy near sparse regions and boundaries. Their use is explicit in eigenfunction and error-field metadata.
- A low reconstruction error for `C` does not imply a well-conditioned lifted representation. Condition numbers and held-out validation remain mandatory diagnostics.
- A fitted constant or bilinear control map is an approximation, never a replacement for the exact mathematical definition when exact gradients are available.

## Verification

The contracts are covered by `tester/tests/test_composable_koopman_pipeline.py`, `tester/tests/test_koopman_control_representations.py`, `tester/tests/test_krom_rollout.py`, and `tester/tests/test_viewer_contracts.py`.
