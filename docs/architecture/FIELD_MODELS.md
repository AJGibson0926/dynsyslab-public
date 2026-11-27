# Field models

A field model is a data-driven vector field $x' = f(x)$, or a map $y = g(x)$, written as a JSON
definition with a `kind`. `core::fieldFromDefinition` (native/include/dynsyslab/core/FieldModels.hpp)
turns every kind into an `Observable`, so fitting, rollouts, scoring and archives treat all kinds
alike. These kinds replace the side-car model families of the former miscellaneous-computes
toolkit.

All regression uses one routine, `thresholdedLeastSquares`: sequentially thresholded ridge least
squares on normalized moments (SINDy: Brunton, Proctor and Kutz 2016, PNAS 113:3932-3937; ridge:
Hoerl and Kennard 1970, Technometrics 12:55-67).

| kind | field | fit |
|---|---|---|
| `polynomial_model` | a dictionary and its coefficients (any `core/ModelLibraries.hpp` family) | `fitMap` |
| `sum` | the sum of a baseline and a fitted residual | `fitMap` with a baseline |
| `compact_radial` | Wendland weights $(1-r)^4(1+4r)$ over the $k$ nearest training centers, normalized to a partition of unity, with an optional affine tail and an asinh coordinate warp (Wendland 1995, Adv. Comput. Math. 4:389-396; Shepard 1968) | `fitCompactRadialField` |
| `local_tangent` | sparse quadratic fields in two-dimensional local PCA charts around training centers, blended by the same weights, with transverse damping toward each chart (Zhang and Zha 2004, SIAM J. Sci. Comput. 26:313-338) | `fitLocalTangentField` |
| `paired_phase` | each real pair is a complex amplitude $u_j$ with $u_j' = (a_j + i\,w_j)\,u_j$; $a_j$ and $w_j$ are linear in saturating amplitude features and, with `phase`, in the phase of $u_j$ relative to $u_0^{h_j}$ (harmonic locking; Pikovsky, Rosenblum and Kurths 2001). Weak-form window rows are optional (Messenger and Bortz 2021, J. Comput. Phys. 443:110525). Winding numbers $h_j$ default to the first record's winding ratios. | `fitPairedPhaseField` |
| `second_order` | imposed kinematic rows $x_{2i}' = x_{2i+1}$ and a fitted acceleration map | `fitSecondOrderField` |
| `pair_network` | second-order taps; tap pairs with equal lag separation share one acceleration map, mixed by directed weights | `fitPairNetwork` |
| `interior` | outer taps by a map; interior taps as critically damped responses $z'' + 2\gamma z' + \gamma^2 z = f$ to an endpoint drive, optionally with all-tap endpoint feedback | `fitInteriorField` |

## Waveform fit and objective

`fitDampedDrive` fits the drive $f = \gamma^2 \varphi(u)$ of $x'' + 2\gamma x' + \gamma^2 x = \gamma^2\varphi(u)$
from waveforms. Each library column is passed through the critically damped filter, which is
discretized exactly for inputs held at interval midpoints. Each record's free transient is
subtracted, and the filtered columns are regressed on the observed positions and velocities.
This is output-error identification on filtered regressors (Ljung 1999, section 7.6). The
`interior` fit uses it with `waveform_drive`.

The waveform objective has two parts:

- `waveformWeights` returns per-sample, per-channel weights. They combine a global-RMS and a
  fixed-interval local-RMS normalized squared error, with trapezoid quadrature, so irregular
  clocks are handled.
- `resolvedDirections` rotates proposed coefficient directions by the SVD of a training response
  Jacobian. It keeps only the resolved combinations, and zero coefficient rows stay zero.

## Calibration

`calibrateField` refits the nonzero coefficients of a polynomial field against short free
trajectories of the training records:

- shots start at evenly spaced samples of each record;
- the solver is bounded damped Gauss-Newton (Nocedal and Wright 2006, chapter 10);
- a failed shot carries a finite penalty;
- each coefficient change has a small penalty, expressed in units of that coefficient.

An incomplete or worse candidate is rejected, and the uncalibrated field is kept. Residuals are
RMS-scaled by default. When global weights are given, they are weighted by the waveform objective
instead.

## Pipeline method `field_model`

`field_model` is a pooled post-processing stage. Every source of the job is a training record.
Records listed in `validation_sources` are only scored.

```json
{"name": "field_model", "settings": {
  "kind": "polynomial",
  "rates": "finite_difference",
  "library": {"family": "polynomial", "degree": 3, "threshold": 0.01},
  "options": {},
  "validation_sources": ["<held-out archive>"],
  "rollout": {"horizon": 2.0, "starts": 4},
  "calibrate": {"duration": 1.0, "starts_per_record": 2}}}
```

- `kind` is `polynomial`, `compact_radial`, `local_tangent`, `paired_phase`, `second_order`,
  `pair_network` or `interior`.
- `options` holds the kind's own settings: `stride`, `neighbors`, `asinh_scale`, `affine_tail`,
  `speed_scale`, `samples`, `decay`, `harmonics`, `phase`, `weak_windows`, `instant_weight`,
  `lags`, `weights`, `tolerance`, `gamma`, `outer`, `drive`, `feedback` and `waveform_drive`.
  `ridge` and `threshold` apply to every kind.
- `rates`: by default, rates are fourth-order central differences of each record, with two
  samples dropped at each end. With `"rates": "stored"`, the record's `coordinate_derivatives`
  are used instead.

The stage publishes a `dynamical_model` object containing:

- `field_model_json`, the definition;
- `training` and `validation` sections, each with the rate relative $L^2$ error, the free-rollout
  relative RMS error and the rollout completion, per record and in total;
- with `calibrate`, the calibration record. Calibration applies to `polynomial` fields only.
- with `controls` (`maximum_step`, `restart_time`, `compare_from`, `visits` [[t1, t2], ...]), the replay controls on the held-out records, or on the training records when there are none:
  - `step_change_over_record_rms`: the rollout with the step halved;
  - `original_on_common` and `restart_on_common`: a restart from the state at `restart_time`, compared with the original rollout on common times;
  - `visits`: the separation of the record's states at pairs of visit times, by cubic Hermite interpolation with the rates (Fritsch and Carlson 1980).

Experiment `select` stages can rank candidates by the validation section.

Tests: native/tests/field_models_tests.cpp, workflow/native/tests/field_model_workflow_tests.cpp.
