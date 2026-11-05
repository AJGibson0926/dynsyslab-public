# Duffing Two-Mode Fresh-Refit KLQR Experiment

Date: 2026-08-29
Status: completed and retained

## Purpose

Test whether adding one control-relevant Koopman eigenfunction to the previously
successful one-eigenfunction Duffing KROM improves controlled-data constant or
bilinear approximations of exact `B_psi`, and then retry the unchanged
transfer-and-release KLQR target.

## Independent two-mode extraction

The anchor was the previously successful nonredundant eigenfunction at parent
retained index 6:

```text
lambda_1 = -0.4911019907065521 + 1.0449812139069314 i
full-spectrum index = 4
```

Every other parent eigenfunction was paired with this anchor and scored using

```text
score = ||Psi(x_initial)-Psi(x_target)||_2
        * min_probe(sigma_min(Ctrb)/sigma_max(Ctrb))
```

where the probes were the initial state, target state, and their midpoint. The
selected complement was retained index 8:

```text
lambda_2 = -0.09147271023857485 + 1.4125559217978911 i
full-spectrum index = 3
lifted target displacement = 0.5850544757
worst normalized controllability singular value = 0.0927038928
score = 0.0542368274
```

This was more than twice the next-best score. The four-dimensional split-real
controllability matrix had full rank at all three probes.

The extracted KROM is self-contained: two complex eigenfunctions, their
diagonal autonomous generator, and no executable parent dependency.

## Fresh controlled-data regression

The same retained random-forced Duffing dataset was used:

- 100 controlled trajectories;
- 200 transition intervals per trajectory;
- `dt=0.01`;
- canonical DynSysLab Duffing input scaling;
- ridge `1e-12`.

No control matrices were inherited from the 20-mode KROM. Constant and both
bilinear representations were refitted directly on the two-mode coordinates.

### Constant matrix

```text
B = [[11.214902007912611 + 18.274802472555040 i],
     [15.524324419153897 +  9.274103208177033 i]]
```

### Bilinear lifted-dynamics fit

```text
B1 = [[-0.958811141145932 + 0.171397329156388 i,
        0.805075627410817 + 0.020923271197538 i],
      [-0.575038972833576 + 0.565989845943990 i,
        0.558710946921742 - 0.248323250428407 i]]

B2 = [[0.233794718212895 + 0.103553769166983 i],
      [0.152322904049353 - 0.064884724782971 i]]
```

### Direct bilinear exact-map fit

```text
B1 = [[-1.969158745868648 + 0.425507021303188 i,
        1.678822023608356 - 0.014623423173693 i],
      [-1.190186056632744 + 1.170164645585496 i,
        1.183533620449188 - 0.499280544161519 i]]

B2 = [[0.474778904952747 + 0.194576199038943 i],
      [0.298012917888076 - 0.118712712430362 i]]
```

The off-diagonal `B1` entries are internal input-map couplings between the two
retained coordinates. They are not off-diagonal autonomous dynamics and do not
refer to removed parent coordinates.

## Shared-actuator interpretation

Free Koopman evolution remains diagonal:

```text
dot(psi_1) = lambda_1 psi_1
dot(psi_2) = lambda_2 psi_2
```

Under control, one physical scalar command enters both equations:

```text
dot(Psi) = A Psi + B_psi(Psi) u
```

Therefore control designed from one eigenfunction can move another whenever
the corresponding row of `B_psi` is nonzero, but it does not automatically make
that second coordinate track its target. The two-mode KLQR explicitly balances
both lifted errors through one actuator. Bilinear off-diagonal terms further
mean that input effectiveness in one coordinate depends on the current value
of the other coordinate.

## Approximation errors

Relative Frobenius error against analytic exact `B_psi` over the free-training
rectangle:

| Representation | Minimum | Median | Mean | 95th percentile | Maximum |
|---|---:|---:|---:|---:|---:|
| Constant refit | 1.0473 | 19.3807 | 22.3264 | 49.3245 | 148.4942 |
| Bilinear dynamics refit | 0.3726 | 0.9799 | 0.9951 | 1.2723 | 2.0391 |
| Direct bilinear map refit | 0.09595 | 1.0065 | 1.0551 | 1.6657 | 3.6636 |

Both bilinear errors remain approximately order one across a typical point in
the training rectangle. The direct-map fit is locally better at its best point
but is not globally accurate enough to replace exact `B_psi`.

## Unchanged KLQR target

- initial state `(0.5,0)`;
- target `(-0.25,0)`;
- controller active on `[1,10]`, then released;
- final time `20`;
- `Q=12500 I_4`;
- `R=[[0.01]]`;
- no control cap;
- dimensionless dynamics and input;
- `u=-K(Psi(x)-Psi(target))`;
- continuous-RHS controller evaluation;
- SciPy BDF, `rtol=1e-8`, `atol=1e-10`, `max_step=0.01`;
- output `dt=0.0001`.

## Results

| Input representation | Status | State at release | Target error | Peak `|u|` | Controller calls |
|---|---|---:|---:|---:|---:|
| Exact analytic `B_psi` | Complete | `(-0.2500717, 1.27e-7)` | `7.17e-5` | `896.02` | 9,297 |
| Fresh constant `B` | Failed | unavailable | unavailable | unavailable | unavailable |
| Fresh bilinear dynamics | Complete | `(0.8776528, 4.32e-10)` | `1.12765` | `1497.66` | 197,606 |
| Fresh direct bilinear map | Complete | `(0.8717753,-3.72e-10)` | `1.12178` | `1422.24` | 51,215 |

The constant case failed with `Required step size is less than spacing between
numbers.` Both bilinear cases held the physical state near `x=0.87` rather than
transferring it to the target. Their nearly identical failure states show that
trajectory-derivative residuals are not the primary cause: even the direct fit
to exact `B_psi` fails similarly.

## Conclusion

The two-mode spectral model is valid for the control task because exact
`B_psi` achieves essentially perfect target acquisition and correct post-release
settling. Adding the most control-relevant second eigenfunction does not rescue
the currently defined constant or affine-bilinear statistical approximations.
Their order-one input-map errors are amplified by KLQR and produce a stiff,
expensive, and inaccurate closed loop.

This supports three conclusions:

1. Exact analytic `B_psi=nabla(Psi)B` remains the authoritative representation.
2. More eigenfunctions do not automatically improve a fixed low-order affine
   approximation of the input map.
3. Any future approximation improvement should target input-map basis capacity,
   locality, weighting, or nonlinear representation—not hidden parent features.

## Artifacts

- Machine-readable report:
  `archives/experiments/duffing_two_mode_fresh_refit_klqr_20260829/report.json`.
- Matrix archive:
  `archives/experiments/duffing_two_mode_fresh_refit_klqr_20260829/duffing_two_mode_fresh_control_regression.npz`.
- Viewer dashboard:
  `viewer_outputs/duffing-two-mode-krom-fresh-refit-klqr.html`.
- Reproduction workflow:
  `scripts/run_duffing_two_mode_refit_klqr.py`.
