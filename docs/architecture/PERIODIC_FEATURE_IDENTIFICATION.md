# Periodic history features and nonlinear identification

This interface tests information retained by a linear feature map of measured
past histories. It does not assign autonomous linear dynamics to those features.
Numerical results and system-specific derivations belong in caller-owned data.

## Ownership

- `dynsyslab_spectral`: `periodicDelayPacketFeatures` in
  `PeriodicSpectralComparison.hpp` exports frozen Rigged-DMD packet maps.
- `dynsyslab_control`: `identifyPeriodicHistoryFeatures` in
  `PeriodicFeatureIdentification.hpp` performs batched nonlinear estimation.
- `dynsyslab_math`: the existing Fourier atlas/table owns interpolation.
- `dynsyslab_numerics`: existing SVD and bounded least-squares implementations.
- Workflow exposes both operations through `dynsyslab-cli periodic-analysis`;
  `dynsyslab_client.periodic` is an optional thin Python transport.

## Frozen packet map

Request schema `dynsyslab.periodic-query/v1`, operation `packet_features`:

```json
{
  "schema": "dynsyslab.periodic-query/v1",
  "operation": "packet_features",
  "archive": "caller-selected-model.npz",
  "channel": 0,
  "packet_centers": [0, 0.4, 0.8],
  "bandwidth": 0.05,
  "order": 6
}
```

An inline `model` can replace `archive`. A verified periodic-delay contract,
unnormalized linear dictionary, mpEDMD operator and coordinate readout are
required. The reference observable is the selected channel at the **newest**
delay position. All offsets are shifted to make that position zero, so input
histories are causal. No model is refitted or written.

The result includes parent model identity, original coordinate contract,
sample interval, centers, bandwidth/order and `features`. That object contains
`channels`, `time_offsets`, `outputs` and row-major `coefficients`. Input order
is offset/channel; output order is center/(real, imaginary). For row history
\(h\), features are \(hC\). Centers are radians per model sample. A finite
regularized packet is neither a singular distribution nor necessarily an
individual Koopman eigenfunction.

## Estimation from histories

Operation `feature_identify` takes:

- `knots`: the standard ordered periodic atlas `{parameter, channels}` rows;
- `features`: the same map object, or any supplied real linear map;
- `histories`: flattened sample/input real past measurements;
- `noise_standard_deviations`: one positive scale per independent input;
- optional `whiten` (default true), `relative_tolerance` (default 1e-10);
- optional `forecast_times`, nonnegative times relative to acquisition end;
- optional `table_samples` (8192), `phase_seeds` (256), `starts` (8);
- optional `parameter_subdivisions` (1): uniformly subdivide each training-knot
  interval for coarse initialization without changing the atlas or training data;
- optional `optimization`: maximum iterations, finite-difference step, step
  tolerance and gradient tolerance, using the native bounded solver;
- caller `source` metadata.

Offsets must be nonpositive with at least one zero. The atlas must have at
least two ordered parameter knots. Output channels must be real. Coarse starts
come from the configured atlas-interval grid and uniform phases, not true test parameters. Distinct
seeds are separated by 0.02 of the parameter range or 1/32 turn before local
refinement. This finite search is not a global identifiability certificate.

The returned estimates include inferred parameter/phase, frequency, current
channels, forecast values in time/channel order, residual, convergence reason,
best-start iteration/evaluation counts, original-history relative error, and
local Jacobian singular values. The Jacobian uses parameter divided by its
atlas range and phase divided by 2*pi; these local diagnostics are not global
uniqueness or stability bounds. Feature/atlas digests and effective coefficients
make the applied transform reviewable.

Forecasts advance the inferred orbit using its interpolated frequency and
Fourier readout. They do **not** propagate the fitted mpEDMD matrix. Extrapolation
outside the atlas is rejected. Process noise, forced motion, online recursion,
parameter priors and runtime compensation are not supplied by this operation.

## Correlated feature noise

For independent measurement noise with standard deviations in diagonal \(D\),
let \(DC=USV^T\) be its retained thin SVD. The full transformed covariance is
\(C^TD^2C\). Whitening gives the equivalent effective input map
\[
 B=D^{-1}U,
 \qquad \| (h-\widehat h) B\|_2^2.
\]
Only singular values above the explicit relative threshold are retained. This
handles redundant output features without pretending that they have independent
noise. At full input rank the cost equals the original standardized-history
cost. Rank reduction projects the residual onto the selected input subspace.
With `whiten:false`, the supplied map's Euclidean metric is used directly;
it is not a Gaussian likelihood for independent transformed errors.

Whitening does not create extra information. Small regression residuals or
local conditioning alone do not establish forecast quality. Keep history
window, model, observations and scoring times matched when comparing maps.

## Sources and validation

- Colbrook, Drysdale and Horning, *Rigged Dynamic Mode Decomposition*,
  [primary paper](https://arxiv.org/abs/2405.00782).
- Golub and Van Loan, *Matrix Computations*: SVD and generalized least squares.
- Nocedal and Wright, *Numerical Optimization*, chapter 10: nonlinear least squares.

Native tests independently compare exported packets to an analytic Fourier
kernel on a rotating circle, verify newest-time causality, recover off-knot
orbit/phase parameters, and prove numerically that a full-rank redundant
transform preserves the noise-weighted metric. Zero-rank maps and future
offsets reject. Python contracts exercise the compiled packet and estimation
gates. No private dataset is required by these tests.
