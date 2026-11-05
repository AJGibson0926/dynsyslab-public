# Duffing Koopman architecture, controlled-data regression, and exact-KLQR report

Date: 2026-08-27  
Status: retained implementation and reproducibility report

## Executive result

DynSysLab now has explicit contracts for physical/model observables, observable reconstruction through a Koopman mode matrix `C`, and three independently selectable Koopman input representations:

1. the exact state-dependent map `B_psi(x) = grad_y(Psi(h(x))) J_h(x) B_x(x)`;
2. a fitted bilinear map `B_psi(Psi) = M_j Psi + b_j` for each input channel;
3. a fitted constant lifted input matrix.

The representations coexist in one KROM. Adding an approximation does not overwrite the exact map, and KLQR/KMPC must record which representation they use. KROM extraction subsets every available representation and its observable reconstruction consistently.

The retained Duffing standard operating procedure uses the scalar first state `x` for H-DMD and the full state `(x, xdot)` for SINDy. Using free BDF training data and the exact `B_psi`, a single nonredundant complex eigenfunction plus its implicit conjugate is sufficient for the requested transfer-and-release experiment. It starts at `(0.5, 0)`, reaches `(-0.25, 0)` at control release, and subsequently falls into the left equilibrium near `(-0.5, 0)`. The user visually accepted the Viewer result; it is retained as a successful reference example.

This work does **not** copy the MATLAB architecture wholesale. The MATLAB suites were treated as read-only experimental evidence. Their equations and algorithms were reconciled with DynSysLab contracts, and only the requested behavior was implemented with explicit validation, metadata, persistence, and tests.

## Sources examined and authority boundaries

Two separate local MATLAB code bases were inspected:

- Cleanup suite: `~/Documents/MATLAB/dynsyslab (5-15-25) (cleanup begins)`
- Korda--Mezic companion suite: `~/Documents/MATLAB/KordaMPC`

The second directory is the actual companion implementation, not an analogous file embedded in the cleanup suite. The public upstream code is [Milan Korda's KoopmanMPC repository](https://github.com/MilanKorda/KoopmanMPC/). The controlled-data construction and Koopman control context correspond to Korda and Mezic's work on [optimal Koopman eigenfunctions](https://arxiv.org/abs/1810.08733) and their earlier [linear predictors for nonlinear dynamical systems](https://arxiv.org/abs/1611.03537).

Local MATLAB code was not modified during this implementation. Exact parity inputs should come from an experimental wrapper/export, rather than making the Python package depend on or silently parse a MATLAB working directory.

## Equation reconciliation

### Canonical DynSysLab form

DynSysLab retains the transparent canonical equation

```text
x_dot = v
v_dot = -delta*v + alpha*x - beta*x^3 + g_u*u
```

For the successful control case:

```text
delta = 0.5
alpha = 1
beta  = 4
g_u   = 1
```

so

```text
x_dot = v
v_dot = -0.5*v + x - 4*x^3 + u.
```

The unforced equilibria are `(0,0)` and `(+-0.5,0)`; the two nonzero equilibria are stable for the damped system.

### Korda--Mezic source notation

`KordaMPC/Duffing.m` writes the controlled vector field as

```matlab
0.5 * [2*x(2,:);
       -delta*2*x(2,:) - 2*x(1,:).*(beta + alpha*(2*x(1,:)).^2) + u]
```

with `delta=0.5`, source `beta=-1`, and source `alpha=1`. Expanding it gives

```text
x_dot = v
v_dot = -0.5*v + x - 4*x^3 + 0.5*u_source.
```

Thus the autonomous dynamics are exactly the same. Only the command convention differs:

```text
u_DynSysLab = 0.5*u_source.
```

The source's random command range `[-1,1]` is therefore the canonical DynSysLab range `[-0.5,0.5]` when `g_u=1`. DynSysLab does not adopt the source's rescaled equation notation. The `korda_mezic_optimal_2020` preset instead makes the source convention explicit with `input_gain=0.5`, while the retained exact-KLQR example uses canonical `input_gain=1`.

### Cleanup-suite parameter issue

The cleanup code contains more than one Duffing convention. Its general equation has the canonical form, but a generic default parameter set equivalent to `-0.5*v + x - x^3` is inconsistent with code that hard-codes the nonzero equilibria as `+-0.5`. The requested stabilization experiment actually uses the cubic coefficient `4`, for which `+-0.5` is correct. DynSysLab therefore records the actual experiment as a named preset and does not change the generic Duffing defaults implicitly.

## Free-training data: random ICs, seed, and integration

The Korda source uses:

- top-level seed `2141442`;
- 100 trajectories;
- `randn(2,1)` followed by division by its Euclidean norm;
- initial conditions on the unit circle;
- fixed classical RK4 with `deltaT=0.01`;
- 800 stored states, covering source sample times `0` through `7.99`.

The exact 100-state table is retained in `data/duffing_korda_mezic_initial_conditions.csv`. Its numeric little-endian float64 SHA-256 is:

```text
decabdb25e40936dec8ab5281c4ea299d474ef0483b9b4fcfe66c0f974713cc7
```

The checksum is computed over the numeric 2-by-100 MATLAB orientation, not over CSV bytes, so harmless line-ending changes cannot masquerade as a data change.

The current Duffing H-DMD/SINDy SOP intentionally follows the prior parity decision to use SciPy BDF, the closest available solver family to `ode15s`, with `rtol=1e-10`, `atol=1e-12`, output spacing `0.01`, and 30 exact retained ICs. It includes the endpoint at `t=8`, giving 801 states. This is not silently claimed to be byte-identical to the Korda fixed-RK4 data; both the solver and endpoint convention are recorded in metadata.

## Controlled training data and its random signal

The actual Korda controlled-data section uses:

- the same 100 initial conditions as the free trajectories;
- 200 state transitions per trajectory;
- one scalar command per transition;
- `Ubig = 2*rand([200,100])-1` in source command coordinates;
- a fixed RK4 step of `0.01`, holding each command for the entire step;
- 201 states and 200 commands per accepted trajectory;
- rejection if any state component has absolute value above `3` or is NaN.

An important reproducibility detail is that the top-level seed alone does not uniquely specify `Ubig`. Earlier sections consume random numbers, and the eigenvalue-optimization branch calls `rng(31)`. Exact parity therefore requires exporting the actual `Ubig` from the reference run. `sample_duffing_controlled_inputs` creates new, explicitly NumPy-generated experiments; it is deliberately not labeled MATLAB-stream compatible. `generate_duffing_controlled_training` accepts either the normalized DynSysLab control cube or MATLAB's `Ubig` layout, applies the `0.5` command conversion, performs the fixed-RK4 solve, applies the state-quality gate, and returns a typed `ControlledDataset`.

The cleanup suite's analogous random-control calls are not equivalent when they draw a new random value at every adaptive ODE right-hand-side callback. That signal is solver-call-dependent, is not an interval-held stochastic process, and should not be used as the parity reference.

## Observable contract and change of variables

An observable is now a first-class named map

```text
y = h(x),       J_h(x) = dh/dx.
```

It is distinct from a SINDy/EDMD feature dictionary. The observable chooses the modeled physical quantities; the dictionary chooses a basis in which functions of those quantities are represented.

`Observable` supports identity, a state component, stacks, batch evaluation, an analytic Jacobian, and explicit finite-difference fallback. Archive reconstruction is intentionally limited to serializable built-in descriptors unless a future callable registry is introduced.

If Koopman eigenfunctions are learned in observable coordinates, `Psi = Psi(y)`, then for a control-affine physical model

```text
x_dot = f(x) + B_x(x) u
```

the exact lifted input term is

```text
B_psi(x) = grad_y(Psi(h(x))) J_h(x) B_x(x).
```

For full-state identity observables, `J_h=I` and this reduces to the familiar `grad_x(Psi) B_x`. This chain rule is evaluated at every controller call for the exact representation. It is the state dependence that the earlier RPE investigation found to be essential.

## H-DMD, Koopman modes, and the reconstruction matrix C

For an observable `g(x)`, the Koopman expansion is

```text
g(x(t)) ~= C Psi(x(t)).
```

In a nonredundant complex representation, a retained complex coordinate implies its conjugate contribution. If column `c_j` corresponds to nonreal `psi_j`, reconstruction uses `2*Re(c_j*psi_j)`; a real coordinate uses `c_j*psi_j`. The split-real matrix therefore uses columns `[2 Re(c_j), -2 Im(c_j)]` for a complex pair.

The cleanup H-DMD path constructs DMD modes, obtains modal amplitudes from a mode pseudoinverse at its declared reference snapshot, and uses the selected observable's amplitude/mode values as `nu_n` in the Koopman expansion. Those `nu_n` values are the native H-DMD `C` coefficients.

The Korda source normalizes separate eigenfunction groups so each state component is reconstructed by summing its group. Its line equivalent to

```matlab
C = bdiag(ones(1,N_efun/n), n)
```

is a special normalization of the same observable-reconstruction idea: the first block reconstructs `x`, the second reconstructs `v`.

DynSysLab's `ObservableReconstruction` stores the declared observable, complex nonredundant coefficients, intercept, validation metrics, and provenance. `KROM.C(output="complex"|"split")` exposes it, and `KROM.reconstruct_observable` applies it.

H-DMD and SINDy eigenfunctions have independent normalization. Directly copying H-DMD `C` to newly fitted SINDy functions is only valid after scale alignment and held/common-state validation. In this run that transfer failed decisively, with relative residuals near unity. DynSysLab correctly rejected it and fitted a new scalar-`x` reconstruction for the SINDy eigenfunctions, while preserving the failed transfer diagnostics and native H-DMD `C`. This is a useful scientific result, not an error to hide.

## Three Koopman input representations

### 1. Exact, state-dependent B_psi

This is the physical default and the representation used by the successful example:

```text
B_psi(x) = grad_y(Psi(h(x))) J_h(x) B_x(x).
```

It requires differentiable eigenfunctions and the system's true input map. It remains available even when fitted approximations are attached.

### 2. Bilinear approximation

For each input channel `j`, DynSysLab fits

```text
B_psi,j(Psi) ~= M_j Psi + b_j,
```

so

```text
Psi_dot = A Psi + sum_j (M_j Psi + b_j) u_j.
```

The regression uses controlled trajectories and the continuous lifted residual `Psi_dot-A Psi`. It is the explicit version of the cleanup suite's `B1*Psi + B2` construction. The operators, offsets, time domain, sample interval, fit diagnostics, and provenance are stored in `BilinearControlMatrix`.

### 3. Constant lifted input matrix

The constant model is

```text
Psi_dot = A Psi + B_c u
```

or, in sampled form,

```text
Psi[k+1] = A_d Psi[k] + B_d u[k].
```

DynSysLab implements the Korda multi-step output regression rather than only a one-step lifted residual. With fixed `A_d` and observable reconstruction `C`, every future observed-state error contributes rows of the form

```text
sum(k=0..j-1) kron(u[k]^T, C A_d^(j-k-1)) vec(B_d).
```

Vectorization uses MATLAB/Fortran column order. The fitted discrete matrix is retained. For continuous KLQR, conversion solves

```text
B_d = (integral from 0 to dt of exp(A_c s) ds) B_c
```

by least squares; it does not confuse discrete eigenvalues or a discrete `B_d` with a continuous generator/input matrix.

## KROM and controller behavior

A KROM may now retain `exact`, `bilinear`, and `constant` simultaneously. Its default remains `exact` unless a caller explicitly promotes another representation. Both `KoopmanLQRController` and `KoopmanMPCController` accept `control_representation="exact"|"bilinear"|"constant"|"default"`, validate availability at construction, and record the selection in archive metadata.

KROM extraction retains and restricts:

- retained eigenvalues and eigenfunctions;
- selected output rows from every attached constant representation;
- the closed principal restriction of every attached bilinear representation;
- observable-reconstruction columns;
- coordinate-observable metadata.

For selected coordinates `S`, constant restriction is `B[S, :]` and bilinear
restriction is `M[:, S, S]`, `b[S, :]`. The extracted KROM is fully independent:
it neither references the source object nor evaluates removed eigenfunctions.
Because a regressed input generator need not be diagonal, the bilinear
restriction is a closed reduced approximation rather than a guarantee that the
selected parent output row is preserved pointwise. A fresh reduced regression
may be requested explicitly if higher accuracy is needed; extraction itself
never performs one silently. Exact `B_psi` and constant row restriction require
no refit.

A previously validated whole-model state decoder is invalidated by extraction and must be refitted. This prevents a reduced KROM from inheriting a decoder it no longer supports.

## Retained default Duffing SOP

The reference workflow is encoded in `DuffingKoopmanSOP` and `duffing_hdmd_sindy`:

1. Load the first 30 exact retained unit-circle ICs.
2. Integrate unforced canonical Duffing dynamics with BDF from `0` through `8`, output `dt=0.01`.
3. H-DMD uses only the scalar first component `x`.
4. H-DMD uses `window-snapshots`, matrix size `400`, sampling window `400`, rank `20`, and terminal-amplitude-energy sorting.
5. The historical global 2,000-delay setting is not silently clipped; it is infeasible for an 801-sample member, so this experiment explicitly declares 400.
6. SINDy receives the full state `(x,v)`, a degree-4 polynomial dictionary, normalization, threshold `0.05`, the Gram-eigh null-space backend, and the H-DMD spectrum.
7. H-DMD `C` is scale-aligned and validated before transfer. If invalid, SINDy's observable reconstruction is fitted and the rejection diagnostics are kept.
8. Extract the smallest KROM that passes the control experiment, retaining conjugate duals implicitly.
9. Use exact `B_psi`, dimensionless dynamics and controller coordinates, continuous RHS controller evaluation, and archived control history.

## Successful single-eigenfunction transfer-and-release example

The selected nonredundant eigenfunction is retained index `6`. Its raw full-spectrum index can change when numerically equivalent conjugate/order bookkeeping changes, so the complex eigenvalue identity—not the transient raw index—is the durable selector:

```text
continuous lambda = -0.49110199219222445 + 1.0449812136970078 i
discrete mu       =  0.9950466881082771  + 0.010398429458889748 i
```

Its conjugate is implicit, so the real KLQR realization is two-dimensional. “Single eigenfunction KROM” in this report means one nonredundant complex eigenfunction/conjugate pair, not a single real scalar coordinate.

Controller configuration:

```text
x(0)            = (0.5, 0)
target           = (-0.25, 0)
active window    = [1, 10]
release time     = 10
final time       = 20
Q                = 12500 I_2
R                = 0.01
law              = u = -K(x) (Psi(x)-Psi(target))
input map        = exact B_psi
solver           = SciPy BDF
rtol / atol      = 1e-8 / 1e-10
max internal step= 0.01
Viewer output dt = 0.0001
```

Accepted run metrics from the retained manifest:

```text
state at release                  (-0.2501115763, -0.0000019367)
final state at t=20               (-0.4837885510, -0.0120155367)
distance to left equilibrium      0.0201788058
state RMSE on t in [2,10]         0.0877554872
|u| 99th percentile on [1,10]     0.3123375429
peak |u|                          431.8161068471
```

The RMSE includes the required inter-well transit, so it is not a steady-state-only tracking statistic. The very large engagement impulse remains part of the established uncapped workflow; the user explicitly chose to retain and display around such transients rather than silently alter the target or control law.

All ten available one-pair candidates were screened. The selected pair was the best candidate that both reached the target at release and produced the correct post-release basin transition. Because it passes those gates, two nonredundant eigenfunctions are **not** the minimal KROM for this case. A two-pair attempt remains the next escalation only if a future target or stricter acceptance threshold makes every one-pair model insufficient.

The machine-readable record is `examples/duffing-minimal-exact-klqr-transfer-release/manifest.json`; the narrative record is `docs/case-studies/DUFFING_MINIMAL_EXACT_KLQR_TRANSFER_RELEASE_2026-08-27.md`.

## Code boundaries implemented

- `dynsyslib/observables.py`: state-to-observable maps and Jacobians.
- `rom/observation.py`: validated state decoders and observable `C` reconstruction.
- `rom/control_inputs.py`: typed constant and bilinear lifted input representations.
- `rom/control_regression.py`: controlled datasets, multi-step constant-`B` regression, and bilinear regression.
- `rom/koopman.py`: exact chain rule, representation coexistence/selection, `C`, reconstruction, and extraction.
- `control/koopman.py`: explicit KLQR/KMPC representation selection and provenance.
- `postprocessing/koopman/dmd.py`: scalar-observable H-DMD, terminal amplitudes, and energy sorting.
- `postprocessing/koopman/koopman_sindy.py`: full-state or general observable-coordinate SINDy with transformed derivatives.
- `postprocessing/koopman/storage.py` and `rom/storage.py`: archive round-trip for observables, `C`, and approximate input models.
- `postprocessing/koopman/workflows.py`: exact Duffing ICs, free-data SOP, controlled fixed-RK4 reproduction, and control-signal conversion.
- The dated case runner has been retired; its manifest and immutable artifacts
  retain the exact configuration. The same solve/control/publication steps use
  the supported Solver/Post-processor API and KROM/control package interfaces.

## Remaining cautions and recommended next stages

1. **Do not change generic Duffing defaults implicitly.** Promote named, tested presets and require experiment manifests to state the exact equation.
2. **Export exact MATLAB controls for parity.** A seed is insufficient because of intervening RNG calls. Compare state arrays after every fixed RK4 step before comparing learned `B`.
3. **Validate C in the destination eigenfunction normalization.** The failed H-DMD-to-SINDy transfer demonstrates why metadata alone cannot establish compatibility.
4. **Fit approximate control models without replacing truth.** Constant and bilinear maps are model reductions; exact `B_psi` remains the physical reference whenever analytic gradients and `B_x` exist.
5. **Use held-out controlled trajectories.** Report multi-step observed-state RMSE, lifted residuals, rank, condition number, and basin/target performance for all three representations.
6. **Separate reconstruction targets.** A scalar H-DMD `C` reconstructs `x`, not automatically `(x,v)`. Full-state reconstruction should use a declared full-state observable and its own validation.
7. **Preserve continuous/discrete distinctions.** Store `lambda`, `mu`, `A_c`, `A_d`, `B_c`, and `B_d` with sample-time provenance; never infer one domain from an array name.
8. **Treat two-pair KROMs as an evidence-driven escalation.** Screen one-pair models first, then pairs under the same controller/target settings, using predeclared acceptance gates.
9. **Study the initial impulse separately.** Smooth engagement, target prefiltering, or finite-horizon design can be optional experiments, but the accepted uncapped case stays the baseline.
10. **Add native acceleration only behind these contracts.** Large linear algebra or batch integration may move to optimized BLAS/LAPACK, C++, or Fortran without moving observables, provenance, archive validation, or controller policy across the boundary.

## Verification

Focused tests cover:

- observable Jacobian participation in exact `B_psi`;
- coexistence and explicit selection of all three input representations;
- recovery of a known continuous constant `B` from the multi-step discrete-output regression;
- fixed-RK4 controlled Duffing generation and the `0.5` source-command conversion;
- archive and KROM behavior exercised by the retained end-to-end run.

The final repository-wide verification completed with `301 passed, 1 skipped`.
The skip is the established optional-native-capability test when that compiled
extension is unavailable; there were no failures.

The saved example must be considered reproducible only when every artifact hash in its manifest matches. Re-running the script creates a new immutable archive lineage and updates the named manifest to that new accepted lineage; it does not mutate previous archives.
