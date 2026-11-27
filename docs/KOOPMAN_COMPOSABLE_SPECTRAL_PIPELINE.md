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

## Native post-processing methods

Each stage is a native post-processing method in a versioned job document or
experiment stage, run with `dynsyslab-cli run` or the Post-processor:

| Stage | Methods |
|---|---|
| Eigenvalues | `hdmd`, `eigenvalue_optimization`, or a supplied spectrum archive |
| Eigenfunctions | `sindy`, `boundary_eigenfunctions` (H-DMD also supplies sampled eigenfunction values) |
| Observable reconstruction | `observable_reconstruction` |
| Control-map identification and error fields | `control_regression` |
| Spectral extraction | `krom_extract` |

Stage outputs retain the source spectrum and method artifact so provenance is
not lost when stages are chained.

### Optimization eigenvalues

`dynsyslab::core::optimizeTrajectoryEigenvalues(...)` minimizes the trajectory
reconstruction residual after analytically eliminating the
trajectory-specific boundary coefficients. It accepts one initialized discrete
eigenvalue group per observable output and can preserve conjugate-pair
structure. Its output is deliberately an eigenvalue-only spectrum; no implicit
eigenfunction fit occurs. The `eigenvalue_optimization` stage accepts a seed
spectrum plus a compatible free-trajectory ensemble and records both discrete
and continuous spectra. In `DynSysLab.exe`, **Configure eigenvalue
optimization…** exposes all required inputs. In `DynSysLabViewer.exe`,
**Eigenvalue optimization** displays the seed/optimized spectra and convergence
evidence.

### Fitted-boundary functions

`boundary_eigenfunctions` accepts any compatible eigenvalue spectrum
(`spectrum`). It fits per-trajectory boundary values, propagates them by
`mu**k`, and constructs state-space interpolants. It also constructs the
associated observable reconstruction: the modal coefficients `nu` become the
columns of `C`.

### SINDy functions

The `sindy` stage accepts an external spectrum. Discrete source eigenvalues are
converted to continuous generator eigenvalues by `lambda = log(mu) / delta_t`,
and the exact source eigenvalue objects remain traceable in metadata.

## Eigenfunction representation and gradients

Every native eigenfunction carries an explicit function representation
(analytic expression, interpolant, or sampled values); see
[eigenfunction objects](architecture/EIGENFUNCTION_OBJECTS.md). Analytic
representations use their analytic gradient; an interpolant's gradient is
numerical, and a request that requires an analytic gradient fails when none
exists. `B_psi(x)` applies the coordinate-observable Jacobian, so identity
state observables recover the ordinary state-gradient control map and
nontrivial observables apply the chain rule.

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

The `control_regression` stage operates on any KROM (`krom`, `selected_modes`)
that has the data required by the selected method:

- `constant_method: "multistep-output"`: multistep output regression from a controlled dataset and a compatible `C` reconstruction;
- `bilinear_method: "controlled-residual"`: continuous lifted-dynamics residual regression from controlled state/input trajectories;
- `bilinear_method: "exact-map"`: direct regression of exact pointwise `B_psi(x)` onto `B1 Psi(x) + B2`, isolating representation error from trajectory-derivative error.

The exact representation is still authoritative. Constant and bilinear matrices are named, optional KROM representations and must be selected explicitly or attached by the regression result.

When a KROM is spectrally extracted, every retained control representation is
made self-contained. Constant maps restrict to `B[S, :]`. Bilinear maps use the
closed principal restriction `M[:, S, S]` and `b[S, :]`. The resulting KROM has
no live reference to its source KROM and needs only its own retained
eigenfunctions to evaluate the input map. This restriction is deterministic;
if its accuracy is inadequate, fitting a new input map directly on the reduced
coordinates is a separate, explicit operation and never occurs silently.

## Error fields and native Viewer display

`control_regression` also evaluates the stored constant and bilinear
representations against exact `B_psi(x)` over a two-dimensional rectangle
(`bounds`, `grid_columns`, `grid_rows`). At each grid point it computes

`absolute_error(x) = ||B_psi_exact(x) - B_psi_approx(x)||_F`

and the relative error against `||B_psi_exact(x)||_F`, and archives both
alongside the regression. The native Viewer displays these fields as
equal-axis phase-space heat maps.

## Important limitations

- Error maps are meaningful only over a declared region supported by training data. They do not validate extrapolation.
- Relative error becomes sensitive where exact `B_psi` is near zero; the archive records the denominator floor.
- Scattered-interpolation gradients may be noisy near sparse regions and boundaries. Their use is explicit in eigenfunction and error-field metadata.
- A low reconstruction error for `C` does not imply a well-conditioned lifted representation. Condition numbers and held-out validation remain mandatory diagnostics.
- A fitted constant or bilinear control map is an approximation, never a replacement for the exact mathematical definition when exact gradients are available.

## Verification

The contracts are covered by the native workflow CTest gates, including
`dynsyslab_eigenvalue_optimization_postprocess_tests`,
`dynsyslab_boundary_postprocess_tests`, and `dynsyslab_krom_postprocess_tests`.
