# RPE designed controlled data for Koopman input-map identification

**Status:** first designed dataset completed and retained

**Scientific objective:** determine whether constant or bilinear control maps
can replace the analytic, state-dependent `B_psi(x)` in RPE Koopman control.

**KROM used for the audit:** the established clean-isolation one-complex-mode
SINDy KROM with continuous eigenvalue
`-0.21061644578341268 + 4.818634561556464i` and its conjugate dual in split-real
controller coordinates.

## Result

The new controlled data are sufficiently accurate to identify the input map.
The symmetric-pulse estimator recovers held-out analytic `B_psi` samples with
relative Frobenius error `0.0107513`. The remaining approximation error is
almost entirely the capacity limit of the selected constant or bilinear map,
not inadequate excitation or an inaccurate command history.

| Approximation fitted only from controlled data | Held-out relative error | Moderate target-path relative error | Exact-map oracle floor on target |
|---|---:|---:|---:|
| Constant | 0.3369 | 0.303858 | 0.302510 |
| Complex-linear `M Psi + b` | 0.2419 | 0.210411 | 0.208864 |
| Conjugate-complete `M Psi + N conjugate(Psi) + b` | 0.0915 | 0.086198 | 0.085483 |

The data-to-oracle gaps on the target are only `0.00135`, `0.00155`, and
`0.00071`, respectively. More data of the same type will not make the standard
constant or complex-linear forms exact. A more expressive representation, a
more locally weighted fit, or a different eigenfunction is required to lower
those floors.

The conjugate-complete form is not an arbitrary nonlinear extension. For one
nonredundant complex eigenfunction it is the general real-linear bilinear map
on the corresponding two split-real coordinates. Its substantially lower
error is therefore directly relevant to the controller's actual real
representation.

## Dataset 1: persistent-excitation rollouts

The long dataset contains 24 dimensionless RPE trajectories, each with 6,000
fixed-RK4 transitions at `delta_tau=0.001`. Every stored input is the exact
interval-held command in `p_control/p0`; no pressure history is inferred from
state derivatives.

Initial conditions are independently stratified over

```text
0.78 <= r(0) <= 1.22
-0.8 <= dr/dtau(0) <= 0.8.
```

Commands combine a balanced zero-mean multilevel random telegraph signal with
a random-phase multisine at angular frequencies `1.5`, `3.0`, `4.8`, and `7.0`
per `tau`. Peak schedules are `0.45`, `0.60`, `0.75`, and `0.90 p0`, and hold
durations are 10, 20, 40, and 80 integration steps.

All 24 trajectories passed the state-quality gate. The realized dataset has:

| Quantity | Value |
|---|---:|
| radius range | 0.611488 to 1.569710 |
| velocity range | -3.289489 to 3.307845 |
| input range | -0.877177 to 0.898777 |
| input RMS | 0.403964 |
| input mean | approximately zero |
| sampled bilinear design rank | 2 of 2 |
| sampled bilinear design condition number | 140.04 |

This dataset is dynamically rich and is retained for Korda--Mezic-style
trajectory regression, derivative/noise studies, and validation. Nevertheless,
the direct dynamics-residual fits have target-path errors of approximately
`0.40` for a constant map and `0.39` for a bilinear map. The cause is not weak
forcing: autonomous SINDy eigenfunction invariance error contaminates
`dPsi-A Psi` and is comparable to or larger than the lifted control response.

## Dataset 2: symmetric opposite-pulse experiments

The pointwise identification dataset contains 600 dimensionless anchor states:

- 400 samples distributed through a narrow tube around the proposed moderate
  equilibrium-centered multisine target;
- 200 samples throughout the supporting rectangle
  `0.76 <= r <= 1.24`, `-0.9 <= dr/dtau <= 0.9`.

Each anchor launches two otherwise identical trajectories with commands
`+1 p0` and `-1 p0`. Each member contains two transitions at
`delta_tau=0.001`, and the complete archive therefore contains 1,200 short
controlled trajectories. The first-step estimator is

```text
Bpsi_hat(x0) = [Psi(x_plus(delta_tau)) - Psi(x_minus(delta_tau))]
               / [(u_plus - u_minus) delta_tau].
```

This cancels autonomous evolution to leading order and does not require the
analytic gradient or exact `B_psi` during data generation or regression. The
median measured lifted response norm is `4.2743e-5`, with a range from
`2.8569e-5` to `7.8876e-5`, which is cleanly resolved in double precision.

Exact `B_psi` is evaluated only after the fits are complete, as a validation
standard and to calculate the representation-capacity floors in the table.

## Target covered by the design

The data are centered on the retained canonical target
`rpe_multisine_equilibrium_moderate_identification`:

```text
r_ref(tau) = 1 + 0.09 sin(4 tau) - 0.12 sin(3 tau)
v_ref(tau) =     0.36 cos(4 tau) - 0.36 cos(3 tau).
```

It begins at equilibrium `(1,0)` and lies inside the state support of both new
datasets.

## Reproduction artifacts

- `archives/experiments/rpe_control_identification_designed/rpe-persistent-excitation.npz`
- `archives/experiments/rpe_control_identification_designed/rpe-symmetric-pulses.npz`
- `archives/experiments/rpe_control_identification_designed/rpe-control-map-candidates.npz`
- `archives/experiments/rpe_control_identification_designed/summary.json`

The command `dynsyslab-rpe-control-data` regenerates the standard experiment.
The API is exposed through `postprocessing.koopman` and `rom`; all generated
command histories are canonical `ControlledDataset` objects and round-trip
through `save_controlled_dataset` / `load_controlled_dataset`.

The formal eigenfunction/control-map error bundle is retained as
`rpe-minimal-krom-error-evaluations.npz` in the same directory. Its metric,
comparison, persistence, and Viewer contracts are documented in
`docs/architecture/KOOPMAN_ERROR_EVALUATION_FRAMEWORK.md`.

## Next decision

The next experiment should hold the target, KROM, KLQR tuning, solver, and
dimensionless conventions fixed and compare exact, constant, complex-linear
bilinear, and conjugate-complete bilinear control maps. If the 8.6% target-path
map error of the conjugate-complete form still prevents control, the failure
will no longer be attributable to weak or unaligned controlled data.
