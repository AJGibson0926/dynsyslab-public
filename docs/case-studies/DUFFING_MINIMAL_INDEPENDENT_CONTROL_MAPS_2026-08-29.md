# Duffing Minimal KROM: Independent Control-Map Restriction

Date: 2026-08-29  
Status: retained reproducible comparison

## Purpose

This case establishes the extraction contract for fitted Koopman control maps.
A KROM extracted from another KROM is a fully independent model. It does not
retain a parent object, parent eigenfunctions, parent coefficients, or an
external parent-feature evaluator.

For selected nonredundant Koopman coordinates `S`, extraction performs:

- constant map: `B_child = B_parent[S, :]`;
- bilinear map: `M_child = M_parent[:, S, S]` and
  `b_child = b_parent[S, :]`.

This is a closed principal restriction. It is deterministic and introduces no
new fit. A fresh regression on the reduced coordinates is a distinct optional
operation and must be requested explicitly.

## Retained eigenfunction

The previously successful Duffing control mode was reused without re-deriving
the KROM:

- parent retained index: 6 (zero based);
- full-spectrum index: 4 (zero based);
- continuous eigenvalue:
  `lambda = -0.4911019907065521 + 1.0449812139069314 i`;
- discrete eigenvalue at `dt=0.01`:
  `mu = 0.9950466881230384 + 0.010398429461133073 i`;
- reduced complex dimension: 1;
- reduced split-real dimension: 2.

## Meaning of the bilinear principal block

For one input, the parent ten-coordinate bilinear representation is

`B_Psi(Psi) = M_parent Psi + b_parent`,

with `M_parent` shape `(10,10)` and `b_parent` shape `(10,1)`. The controlled
equation for selected mode `j=6` is

`B_Psi,j = b_j + sum(k=0..9) M[j,k] Psi_k`.

An independent one-mode KROM has only `psi = Psi_6`. It therefore retains

`B_psi(psi) = b_6 + M[6,6] psi`.

All `M[6,k] Psi_k` terms for `k != 6` are absent. They are not evaluated through
the parent and are not folded into the offset. Consequently this is an
independent closed approximation, not an exact marginalization of the parent
row.

## Extracted matrices

### Constant regression

Complex form:

```text
B = [[-1191.3966800807323 - 4082.4242013252033 i]]
```

Split-real form:

```text
B = [[-1191.3966800807323],
     [-4082.4242013252033]]
```

### Bilinear lifted-dynamics regression

Complex form:

```text
B_psi(psi) = B1 psi + B2
B1 = [[-2.877493944495333 + 0.09444285560424044 i]]
B2 = [[ 0.2651558084825485 + 0.68845675772793400 i]]
```

Split-real form, with `z = [Re(psi), Im(psi)]^T`:

```text
B_psi,R(z) = B1_R z + B2_R

B1_R = [[-2.877493944495333, -0.09444285560424044],
        [ 0.09444285560424044, -2.877493944495333]]

B2_R = [[0.2651558084825485],
        [0.6884567577279340]]
```

### Direct bilinear control-map regression

Complex form:

```text
B_psi(psi) = B1 psi + B2
B1 = [[-5.809384704526216 + 0.5152942885346369 i]]
B2 = [[ 0.5119804070198101 + 1.2783331234596593 i]]
```

Split-real form:

```text
B1_R = [[-5.809384704526216, -0.5152942885346369],
        [ 0.5152942885346369, -5.809384704526216]]

B2_R = [[0.5119804070198101],
        [1.2783331234596593]]
```

Every array-alignment check against the selected parent row/column was exactly
zero. The written reduced archives contain only the constant matrix, the
`(1,1,1)` bilinear operator, and the `(1,1)` bilinear offset; no parent control
feature arrays are present.

## Unchanged KLQR experiment

- initial state: `(0.5, 0)`;
- target: `(-0.25, 0)`;
- controller active on `[1,10]` and released at `t=10`;
- final time: `20`;
- `Q = 12500 I_2`;
- `R = [[0.01]]`;
- no control cap;
- law: `u = -K(Psi(x)-Psi(target))`;
- dimensionless dynamics and control;
- continuous-RHS control evaluation;
- SciPy BDF, `rtol=1e-8`, `atol=1e-10`, `max_step=0.01`;
- Viewer/output step: `0.0001`.

## Results

| Representation | Release state | Target distance at release | Final distance to left equilibrium | Peak `|u|` | Solve time |
|---|---:|---:|---:|---:|---:|
| Exact `B_psi` reference | `(-0.2501116,-0.0000019)` | `0.0001116` | `0.02018` | `431.82` | retained reference |
| Restricted constant `B` | `(-0.2452139,-0.0026970)` | `0.005494` | `0.02037` | `1057.17` | `14.87 s` |
| Restricted bilinear dynamics | `(0.8938695,0.0000000)` | `1.143870` | `0.04384` | `929.84` | `356.89 s` |
| Restricted bilinear map | `(0.8956284,-0.0000000)` | `1.145628` | `0.04506` | `984.93` | `387.69 s` |

The constant restriction transfers successfully, although it has a larger
initial control spike and worse transient error than exact `B_psi`. Both
one-mode bilinear principal restrictions fail to transfer the state during the
controlled interval. The close agreement of their failed release states is
strong evidence that the missing cross-feature terms are important to the
parent ten-coordinate bilinear formulas.

The failure does not justify retaining the parent in a reduced KROM. The
correct response, if a self-contained bilinear one-mode approximation is
required, is to explicitly regress a new one-mode map and assess its error.
Alternatively, select a larger closed sub-ROM that includes the eigenfunctions
responsible for the important cross terms.

## Approximation-error fields

Over the parent training rectangle, the median relative Frobenius errors were:

- constant restriction: `3943.00`;
- bilinear-dynamics restriction: `1.7080`;
- direct bilinear-map restriction: `2.8001`.

These global field statistics and the control results answer different
questions. The constant fit is globally poor yet happens to support this
specific transfer. The bilinear restrictions are closer globally but place the
closed-loop system on an ineffective branch for this target.

## Reproduction artifacts

- workflow: `scripts/run_duffing_independent_control_klqr_transfer_release.py`;
- machine-readable report:
  `archives/experiments/duffing_independent_control_klqr_transfer_release/report.json`;
- native Viewer dashboard:
  `viewer_outputs/duffing-minimal-krom-independent-control-approximations.html`.

