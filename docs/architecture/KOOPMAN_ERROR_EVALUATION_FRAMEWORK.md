# Koopman Error-Evaluation Framework

**Status:** canonical API and persistence contract

## Purpose

DynSysLab now treats error analysis as part of a Koopman model rather than as
temporary plotting code. The framework separates two scientifically different
questions:

1. Does a fitted eigenfunction obey autonomous Koopman evolution on free data?
2. Does a controlled-data approximation reproduce the chain-rule input map
   of that fitted eigenfunction over a specified state-space region?

This separation is essential. The analytic map

```text
B_psi(x) = nabla(Psi_fit)(x) B(x)
```

is exact with respect to the fitted function `Psi_fit`. It does **not** prove
that `Psi_fit` is an exact eigenfunction of the nonlinear dynamical system.
Autonomous error and input-map approximation error are therefore recorded
independently and can be compared without one concealing the other.

## Free-trajectory validation

`FreeTrajectoryDataset` labels every dataset as either `training` or `test`,
records its coordinate system and provenance, and computes a SHA-256 content
fingerprint from all time and state arrays. A pooled ensemble is expanded into
individual trajectories before evaluation.

Each `KoopmanEigenfunction` can call `evaluate_errors(...)`. Each resulting
`EigenfunctionErrorEvaluation` is attached to that function by a unique ID.
The record contains:

- one-step semigroup error for every adjacent pair;
- finite-horizon semigroup errors at caller-specified horizons;
- an anchored long-term error, propagated from each trajectory's first sample;
- generator error using the trajectory derivative of `Psi`;
- training-domain coverage, including the fraction of test samples outside the
  fit region and the largest normalized excursion beyond it;
- dataset role, fingerprint, coordinate system, source provenance, eigenvalue,
  and numerical conventions.

The finite-horizon comparison is

```text
Psi(x(t + h))  versus  exp(lambda h) Psi(x(t)).
```

For nonuniform trajectories, the observed value is interpolated at the exact
requested horizon. No eigenvalue or eigenfunction coefficient is refit during
validation. The anchored comparison uses

```text
Psi(x(t))  versus  exp(lambda (t - t0)) Psi(x(t0)).
```

The generator comparison is

```text
d Psi(x(t))/dt  versus  lambda Psi(x(t)).
```

The standard `relative_l2` denominator is the observed left-hand side. The
evaluation metadata additionally records generator residual relative to
`||lambda Psi||`, which is the convention used in earlier RPE diagnosis.

`KROM.evaluate_eigenfunction_errors(...)` applies the same dataset to every
retained eigenfunction. It converts physical time/state values to the KROM's
model coordinates when required and then applies the KROM's declared
coordinate observable before calling each eigenfunction. This prevents mixed
dimensional/dimensionless validation.

## Control-map validation

`ControlApproximationEvaluation` retains one candidate approximation together
with:

- the complete `ConstantControlMatrix` or `BilinearControlMatrix` object;
- controlled-dataset summary and SHA-256 content fingerprint;
- regression diagnostics and method provenance;
- the evaluation rectangle and grid;
- exact-map norm, absolute matrix error, and relative matrix error at every
  phase-space grid point;
- aggregate minimum, median, mean, p95, and maximum values.

The pointwise error uses the Frobenius norm:

```text
e_abs(x) = ||B_psi,fit(x) - B_psi,exact(x)||_F
e_rel(x) = e_abs(x) / max(||B_psi,exact(x)||_F, relative_floor).
```

`ControlRegressionResult.evaluate(...)` is the direct bridge from a regression
to this record. `KROM.evaluate_control_approximation(...)` performs the same
operation for any independently supplied representation. Evaluating a
candidate never replaces the KROM's active control representation.

Multiple evaluations of the same family can coexist on one KROM. Thus two
constant matrices or two bilinear maps fitted from different controlled
datasets remain separately executable, auditable, and comparable. Ranking is
allowed only when all candidates share the same grid and the same exact
`B_psi` truth field; comparisons over different regions are rejected rather
than presented as though they were commensurate.

## Persistence and Viewer

`save_krom_error_evaluations(...)` writes one compressed bundle containing:

- all attached free-trajectory evaluation summaries by eigenfunction;
- every attached candidate control representation;
- every full-resolution phase-space error field;
- a KROM signature containing the method, coordinate system, retained
  dimension, and continuous eigenvalues.

`load_krom_error_evaluations(...)` validates that signature before attaching
anything to a KROM. A bundle cannot silently attach to a different spectrum.

Viewer's native `control_error` display accepts both old standalone error-field
archives and the new multi-evaluation bundles. Each retained candidate becomes
its own labeled heat map; relative, absolute, and exact-map-norm displays and
linear/log color scaling remain available.

The installed CLI entry point is `dynsyslab-krom-errors`:

```text
dynsyslab-krom-errors free \
  --krom MODEL.npz --trajectory FREE_SOLUTION.npz \
  --role test --label outside-training-test \
  --horizon 0.001 --horizon 0.1 --horizon 1 \
  --output validation.npz

dynsyslab-krom-errors control \
  --krom CONTROL_EQUIPPED_MODEL.npz \
  --controlled-data CONTROLLED_DATA.npz \
  --representation bilinear --label designed-data-bilinear \
  --bounds 0.76 1.24 -0.9 0.9 --grid-shape 81 81 \
  --output validation.npz

dynsyslab-krom-errors compare \
  --krom MODEL.npz --bundle first.npz --bundle second.npz \
  --metric relative.mean
```

The equivalent Python API is preferred inside reproducible pipelines.

## Retained RPE application

The first production bundle is:

```text
archives/experiments/rpe_control_identification_designed/
    rpe-minimal-krom-error-evaluations.npz
```

It is attached to the established one-complex-mode RPE SINDy KROM with

```text
lambda = -0.21061644578341268 + 4.818634561556464 i.
```

The free training archive contains 25,001 samples. Its finite-time relative
semigroup errors are:

| Horizon `h` | Relative L2 error |
|---:|---:|
| 0.001 | 0.005324 |
| 0.01 | 0.050670 |
| 0.1 | 0.411850 |
| 0.5 | 0.751224 |
| 1.0 | 0.709503 |
| 5.0 | 0.966914 |
| 10.0 | 1.015348 |

The full anchored long-term relative error is `1.12737`. This confirms the
earlier diagnosis: the polynomial SINDy expression has a small one-step defect
but does not preserve ideal exponential eigenfunction evolution over long
horizons. The controlled-map approximation study must therefore not be read as
an autonomous eigenfunction validation.

The trajectory-gradient generator residual is `0.84237` relative to observed
`dPsi/dtau` and `1.09486` relative to `lambda Psi`. The latter is the direct
analogue of the earlier generator-PDE residual normalization; the small
difference from the earlier spot calculation comes from using the archived
trajectory derivative uniformly in the formal evaluator.

Five controlled-data candidates were evaluated over the common rectangle
`0.76 <= r <= 1.24`, `-0.9 <= dr/dtau <= 0.9` on a `61 x 61` grid:

| Candidate | Mean relative map error |
|---|---:|
| pulse-data conjugate-complete bilinear | 0.101221 |
| pulse-data complex bilinear | 0.271089 |
| pulse-data constant | 0.377757 |
| persistent-data bilinear | 0.483733 |
| persistent-data constant | 0.486811 |

The pulse and persistent datasets have different retained content
fingerprints. The ranking therefore records both the approximation-family
effect and the controlled-data-design effect while keeping them identifiable.

## Standard operating policy

- Every new eigenfunction construction should evaluate its source free data as
  `role="training"` before promotion to a canonical KROM.
- New free trajectories must be labeled `role="test"`; tests outside the
  training rectangle are explicitly reported as extrapolation.
- Test data never refit coefficients.
- Every controlled-data regression promoted for control should create a
  `ControlApproximationEvaluation` on a declared region of interest.
- Competing controlled datasets must be retained as separate evaluations, not
  overwritten in the KROM's active `constant` or `bilinear` slot.
- Numerical comparison is valid only over an identical grid and exact-map
  truth field.
- Viewer displays are projections of retained evaluation records; they are not
  the storage boundary.
