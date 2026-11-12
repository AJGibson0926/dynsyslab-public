# Native KLQR parity replay

**Date:** 2026-09-05
**Status:** automated native regression gate
**Scope:** accepted Duffing transfer/release and RPE moderate multisine cases

## Purpose

This gate proves that the compiled engine can reuse the retained KROMs and
control-map regressions directly. It does not refit eigenfunctions, regress new
input maps, or invoke Python. Each controller is evaluated at every adaptive
BDF right-hand-side call and recomputes the state-dependent CARE gain.

The jobs are retained under `examples/native/case-studies`. The CTest target
`dynsyslab_klqr_case_study_parity_tests` executes them through the same
`dynsyslab_workflow` service used by the desktop app and compiled CLI, reads
the resulting NPZ archives through the native archive library, and checks the
scientific metrics below.

## Compatibility corrections required by the replay

The accepted research archives predate the native archive schema. Their exact
semantics are now supported at the archive boundary:

1. NumPy Fortran-ordered arrays are converted to row-major storage on import.
2. Legacy KROM `eigenvalues` are discrete-time values `mu`; the loader computes
   `lambda=log(mu)/delta_t` before constructing the continuous generator.
3. Explicit conjugate rows are reduced to the metadata-declared independent
   modes. The native KROM forms the corresponding split-real realization
   internally, so one nonredundant complex mode remains two real states.
4. Polynomial center, scale, multiplier, degree, and constant-feature policy
   are recovered from archive metadata without changing the coefficients.
5. Self-contained legacy control maps and independently archived evaluation
   bundles can supply constant, bilinear, and conjugate-complete bilinear
   representations. The loaded KROM owns the map; it retains no runtime parent
   reference.

The high-gain Duffing constant map also exposed an ill-conditioned Hamiltonian
subspace. Native CARE now balances the Riccati variable, refines the result
with Newton/Lyapunov steps when useful, requires stable closed-loop poles, and
uses an invariant-subspace-condition-aware residual tolerance. This preserves
the residual check without rejecting a stabilizing result solely because
double-precision cancellation makes a fixed `1e-6` cutoff unattainable.

## Duffing parity

All cases use the retained one-complex-mode KROM, dimensionless Duffing
dynamics, initial state `(0.5,0)`, target `(-0.25,0)`, active window `[1,10]`,
release at `10`, final time `20`, `Q=12500 I2`, `R=0.01`, continuous-RHS KLQR,
and BDF with the retained tolerances.

| Control map | Native release state | Retained release state | Native final state | Retained final state |
|---|---|---|---|---|
| Exact | `(-0.2501115763,-0.0000019364)` | `(-0.2501115763,-0.0000019367)` | `(-0.4837885526,-0.0120155376)` | `(-0.4837885510,-0.0120155367)` |
| Constant | `(-0.2452140101,-0.0026969756)` | `(-0.2452138856,-0.0026970207)` | `(-0.4835196299,-0.0119669236)` | `(-0.4835196200,-0.0119669181)` |
| Bilinear dynamics | `(0.8938695581,-3.48e-8)` | `(0.8938695394,1.96e-9)` | `(-0.5404897202,0.0168131240)` | `(-0.5404897285,0.0168130656)` |
| Bilinear map | `(0.8956283586,-6.73e-9)` | `(0.8956283868,-7.11e-11)` | `(-0.5397237133,0.0212783436)` | `(-0.5397237012,0.0212784038)` |

The native exact case retains the successful transfer. The constant
approximation retains its less accurate but successful transfer. Both
bilinear approximations faithfully reproduce the known failure to cross the
barrier; parity is not reclassified as scientific success.

The native peak control magnitudes were `431.814`, `1057.148`, `929.855`, and
`984.960`, respectively, consistent with the retained values `431.816`,
`1057.171`, `929.844`, and `984.928`.

## RPE parity

The RPE cases use exact Air-in-Water physical parameters, dimensionless
`(r,dr/dtau)` dynamics and `p_control/p0`, equilibrium initial state `(1,0)`,
the moderate canonical multisine target, the accepted minimal KROM,
`Q=5000 I2`, `R=1e-4`, no cap or feedforward, continuous-RHS KLQR, and BDF.

| Control map | Native second-half radius RMSE | Retained | Native second-half velocity RMSE | Retained |
|---|---:|---:|---:|---:|
| Exact analytic | `0.0012126532` | `0.00121273` | `0.0041212349` | `0.00412142` |
| Pulse-data constant | `0.0011802973` | `0.00118033` | `0.0039740032` | `0.00397416` |
| Pulse-data conjugate-complete bilinear | `0.0011897613` | `0.00118979` | `0.0040046969` | `0.00400489` |

All three remain successful and agree with the retained metrics well inside
the automated `5e-6` absolute RMSE tolerance.

## Reproduction

From a configured native build:

```powershell
ctest --test-dir build/native-suite-fast `
  -R klqr_case_study_parity --output-on-failure
```

No Python interpreter, Python module, or Python subprocess participates in
this test. Every output archive also records the complete solver and time
configuration, native execution identity, controller configuration,
coordinate system, and input scaling.
