# Successful Duffing minimal-KROM exact-KLQR transfer and release

Date: 2026-08-27  
Classification: successful retained reference example

## Acceptance statement

The Viewer result was visually accepted by the user on 2026-08-27. The case is retained as the positive reference for the Duffing scalar-H-DMD/full-state-SINDy/exact-`B_psi` workflow.

It is reproduced by:

```powershell
python scripts/run_duffing_exact_klqr_transfer_release.py
```

The authoritative machine record is `examples/duffing-minimal-exact-klqr-transfer-release/manifest.json`. It contains absolute artifact locations, SHA-256 hashes, solver settings, timings, eigenvalues, and metrics. The Viewer page is `viewer_outputs/duffing-minimal-exact-klqr-transfer-release.html`.

## System and target

```text
x_dot = v
v_dot = -0.5*v + x - 4*x^3 + u
```

All state, time, and input coordinates are dimensionless. The controller starts at the right equilibrium `(0.5,0)`, is inactive until `t=1`, targets the non-equilibrium point `(-0.25,0)` until `t=10`, then releases the system. The uncontrolled system subsequently falls into the left well toward `(-0.5,0)`.

## Training and KROM

- 30 exact retained unit-circle initial conditions, seed lineage `2141442`.
- Free SciPy BDF integrations from `0` through `8`, sampled every `0.01`.
- Scalar `x` observable for H-DMD.
- H-DMD window-snapshot geometry: matrix size 400, sampling window 400, rank 20.
- Terminal-amplitude-energy spectrum sorting.
- Full `(x,v)` state for degree-4 polynomial SINDy.
- SINDy threshold 0.05 with normalized coordinates and Gram-eigh null-space solve.
- H-DMD `C` transfer was rejected after scale-alignment validation; a SINDy-native scalar-`x` `C` was fitted and the failed diagnostic retained.

The minimal model contains one nonredundant complex eigenfunction and its implicit conjugate:

```text
lambda = -0.49110199219222445 + 1.0449812136970078 i
mu     =  0.9950466881082771  + 0.010398429458889748 i
```

This is a two-state split-real realization. A second nonredundant eigenfunction is unnecessary for this accepted case: the single pair passes the target-at-release and post-release-basin gates.

## Controller

```text
law           u = -K(x) (Psi(x)-Psi(target))
coupling      exact B_psi = grad(Psi) B_x
Q             12500 I_2
R             0.01
active        [1,10]
integration   SciPy BDF, continuous RHS control evaluation
rtol / atol   1e-8 / 1e-10
max_step      0.01
output_dt     0.0001
```

There is no feedforward pressure term and no arbitrary control cap. The initial high-gain engagement impulse is intentionally retained as part of the established baseline.

## Acceptance evidence

The accepted run reached approximately

```text
x(10) = (-0.2501115763, -0.0000019367)
x(20) = (-0.4837885510, -0.0120155367)
```

with distance `0.0201788058` from the left equilibrium at the final output time. The 99th percentile of `|u|` during the active interval was `0.3123375429`; the isolated engagement spike produced a peak near `431.81611`.

## Reproduction gate

A future run qualifies as this example only if it:

1. uses the manifest equation, coordinates, IC, target, active window, tuning, solver tolerances, and exact input representation;
2. uses one retained nonredundant complex pair with the recorded eigenvalue identity or explicitly records a newly derived successor;
3. reaches the non-equilibrium target at release within a declared numerical tolerance;
4. enters the left basin after release and approaches `(-0.5,0)`;
5. archives training, H-DMD, SINDy, minimal KROM, controlled solution, control history, and Viewer output;
6. records hashes rather than overwriting an earlier archive lineage.

The broader architecture and MATLAB/Korda reconciliation are documented in `docs/DUFFING_KOOPMAN_ARCHITECTURE_AND_CONTROL_REPORT_2026-08-27.md`.
