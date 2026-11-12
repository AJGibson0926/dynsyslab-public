# Native echo-history Koopman pipeline

Date: 2026-10-02
Branch: `echo-only`

## Executive result

DynSysLab's production scientific pipeline was already predominantly native C++: dynamical systems, time integration, observables, SINDy, EDMD/HAVOK/HDMD, Koopman model evaluation, validation, control-map models, KLQR/LQR, archive I/O, the Solver workbench, and Viewer are C++. The Python-only portion was the recent echo-history experiment that pooled several trajectories into one Hankel basis, projected new echo histories through that fixed basis, fit a two-coordinate polynomial EDMD eigenfunction, selected the physical oscillatory mode, and assembled the result for Viewer.

That experimental island now has a native C++ implementation and first-class workflow registration. The Python scripts remain as research references and independent comparison tools; they are no longer required to execute the production computation.

The native viscous-RPE case reproduces the earlier Python result:

| Quantity | Native C++ result | Difference from Python |
|---|---:|---:|
| Selected continuous eigenvalue | `-0.0602447488280805 + 6.282949184160554 i` | `6.72e-12` absolute |
| First four cumulative Hankel energy fractions | `0.03107148, 0.06142573, 0.07705063, 0.09028589` | `5.52e-11` maximum absolute |
| Training trajectories | 5 independent unforced runs | same |
| Held-out trajectories | 5 independent unforced runs | same |
| Latent coordinates used by the polynomial | `xi1, xi2` from a rank-4 fixed history basis | same |

This is numerical parity, not merely a visually similar refit.

## End-to-end architecture

```text
unforced RPE archives
        |
        | scalar echo samples (stored, or derived once from an archived solve)
        v
pooled fixed Hankel basis U and singular values
        |
        | trajectory-local causal windows; no window crosses a run boundary
        v
fixed history projection  echo history -> (xi1, xi2, ...)
        |
        | same immutable U online; only the ring buffer changes
        v
pooled polynomial EDMD in selected xi coordinates
        |
        v
oscillatory-mode search and canonical polynomial conversion
        |
        v
KoopmanSINDyResult / KoopmanModel
        |
        +--> native validation on training and held-out coordinate runs
        +--> native NPZ archive and catalogue index
        +--> native Viewer heat maps, histories, spectrum, and overlays
        +--> learned approximate B_psi route for later KLQR experiments
```

### 1. Training observations

The workflow accepts multiple source archives as one pooled training stage. It enforces the echo-only experiment policy that every training and validation source is unforced and uncontrolled (`forcing_type=free`, `controller_kind=none`, `data_role=free_training`).

If a source contains an archived scalar echo channel, that channel is used directly. Existing simulated RPE archives may instead supply the state history from which the already-native `far_field_acoustic_echo` observable is evaluated once, offline. After this extraction boundary, basis fitting, projection, eigenfunction fitting, validation, and online evaluation consume echo values rather than the RPE equations of motion.

The example manifest is `examples/native/rpe-echo-history-koopman.json`.

### 2. Pooled Hankel basis

`fitPooledHankelBasis` fits a single delay basis over independent echo records. A matrix-free FFT-backed Hankel operator and deterministic randomized Rayleigh-Ritz iteration avoid materializing the large pooled matrix. RMS observation normalization, singular values, cumulative energy, seed, and source counts are archived.

Pooling means concatenating each trajectory's Hankel **columns for the linear algebra**, not concatenating raw trajectories before windowing. Every column belongs to exactly one run, so no artificial delay window bridges the end of one trajectory and the start of another.

### 3. Offline and online latent coordinates

`projectHistory` evaluates all causal windows in a complete record. `CausalHistoryProjector` is the online equivalent: it owns a fixed-length ring buffer, accepts one new echo measurement at a time, and emits coordinates only after the buffer is full. The basis is immutable during deployment.

The test suite checks sample-for-sample equality between batch and ring-buffer projection. No SVD, HAVOK refit, interpolation from `(r, rdot)`, or bubble state estimate is needed online.

### 4. Polynomial Koopman fit

The echo-history experiment uses polynomial EDMD in the latent coordinates. Snapshot pairs are formed within each coordinate run and never across trajectory boundaries. Coordinates and polynomial features are scaled for conditioning, a ridge-regularized finite-dimensional Koopman operator is fit, and its eigenvectors define polynomial eigenfunctions.

The selected polynomial is converted back to canonical raw monomials and stored in the existing `KoopmanSINDyResult` contract. Thus the durable model representation is the same kind of polynomial coefficient/exponent/eigenvalue object used by the existing SINDy path. The difference is how the coefficients are learned:

- continuous-time SINDy solves a generator-PDE regression and requires coordinate time derivatives;
- this echo-only EDMD fit enforces the sampled semigroup relation `psi(xi[k+q]) ~= mu psi(xi[k])` and requires no analytic `echo_dot`, latent derivative, state equation, or chain rule.

The archive records `coordinate_provider=causal-history`. That flag is scientifically important: the polynomial variables cannot be reconstructed from one instantaneous dynamical state.

### 5. Model evaluation and control boundary

For ordinary instantaneous coordinates, `KoopmanModel::evaluate(state)` retains its old behavior. For a causal-history model it deliberately refuses that call and requires `evaluateCoordinates(xi)`. This prevents an accidental claim that the history state is an instantaneous function of `(r, rdot)`.

An exact state-Jacobian control map is likewise unavailable for the stateful, data-driven coordinate provider. A control sensitivity still exists for a specified actuation and sampling protocol, but it must be learned from offline input/echo experiments. The existing approximate-control-map path remains usable through `controlMapCoordinates`. Consequently the next KLQR stage should regress `B_psi` (or the sampled lifted input map) from actuated echo records, then evaluate the eigenfunction online from the same ring-buffer projector. It should not invoke the RPE equations during deployment.

## Native source map

| Responsibility | Native source |
|---|---|
| Pooled Hankel fit, fixed projection, online ring buffer, polynomial EDMD, mode selection | `native/include/dynsyslab/core/HistoryKoopman.hpp`, `native/src/HistoryKoopman.cpp` |
| Canonical polynomial model and external-coordinate contract | `native/include/dynsyslab/core/KoopmanModel.hpp`, `native/src/KoopmanModel.cpp` |
| Archive persistence, including `coordinate_provider` | `native/src/KoopmanArchive.cpp` |
| Eigenfunction accuracy metrics on supplied latent runs | `native/src/KoopmanValidation.cpp` |
| End-to-end postprocess stage and policy checks | `workflow/native/include/dynsyslab/workflow/EchoHistoryKoopman.hpp`, `workflow/native/src/EchoHistoryKoopman.cpp` |
| Solver workflow registration | `workflow/native/src/Postprocess.cpp`, `workflow/native/src/Plan.cpp` |
| RPE and acoustic observable implementation | `native/src/SystemsRPE.cpp` |
| Viewer training-coordinate loader | `viewer/native/src/KromTraining.cpp` |
| Viewer eigenfunction fields, histories, spectrum, overlays | `viewer/native/src/KromWorkspace.cpp` |
| Math labels, including real Greek xi | `viewer/native/src/MathText.cpp` |

## Viewer composition correction

No new eigenfunction plotting subsystem was created. The existing native `Figure`, `Trace`, heat-map, history, domain-color, dashboard span, and KROM spectrum structures are assembled correctly for the new eigenfunctions:

- real field beside `Re(psi)` history;
- imaginary field beside `Im(psi)` history;
- one domain-color field in which hue is `arg(psi)` and brightness is normalized `|psi|`;
- `|psi|` history and `arg(psi)` history beside that single domain-color field;
- dashed black/grey training trajectories with differing opacity over the colored fields;
- `xi1` and `xi2` parsed and drawn as Greek `xi` with subscripts.

The earlier erroneous fourth history (magnitude where phase was required) and the two-plot interpretation of magnitude/phase have both been removed.

Current review images are local outputs, not Git assets:

- `viewer_outputs/rpe-viscous-echo-asymptotic-eigenfunction-field.png`
- `viewer_outputs/rpe-viscous-echo-history-native-domain.png`

## Reproduction record

Native manifest:

```text
examples/native/rpe-echo-history-koopman.json
```

Native derived archive:

```text
archives/derived/echo-history-koopman/
rpe-viscous-echo-long-training-five-ics-trajectory-0__20261002T001350Z__
echo-history-koopman__20261003T015202Z.npz
```

Representative validation results from the five independent runs:

- training one-step relative L2 by run: approximately `1.14e-6, 2.70e-6, 1.10e-6, 9.66e-7, 1.44e-6`;
- held-out one-step relative L2 by run: approximately `6.28e-6, 1.14e-6, 1.06e-6, 1.14e-6, 5.97e-6`;
- held-out aggregate one-step relative L2 from the independent Python validation report: `4.30e-6`;
- held-out anchored long-time relative L2 over the evaluated late-time segments: `1.73e-3`.

The native archive includes the fixed history basis, singular values and energy, observation scaling, EDMD operator and feature scaling, canonical polynomial coefficients and exponent table, continuous and discrete eigenvalues, source pointers, trajectory offsets, retained latent training coordinates, and training/held-out validation metadata.

## Tests and verification

New native tests cover:

- rank-two recovery and energy behavior for a known sinusoidal delay signal;
- exact agreement between offline history projection and online ring-buffer evaluation;
- recovery of a known damped rotating mode (`-0.08 + 5 i`);
- conversion to the canonical polynomial Koopman contract;
- multi-archive workflow execution and provenance;
- enforcement of unforced/uncontrolled training policy;
- state-evaluation rejection and coordinate-evaluation success for causal-history models;
- Viewer domain-color layout, magnitude and phase histories, Greek labels, and external-coordinate grids.

The compatibility pass also preserved the established Viewer aggregate layout by appending the new opacity field instead of inserting it among existing fields. The full native suite was rebuilt after the shared model structures changed, preventing stale executable/DLL ABI mismatches.

## Status and remaining scientific work

The echo-history coordinate transform, Koopman extraction, archive contract, validation, and Viewer handoff are now native C++. The production path no longer depends on Python for these steps.

What is not yet claimed:

1. Two latent coordinates are not assumed globally injective. The agreed experimental policy remains: start with two, observe failure modes, then increase to four when the evidence requires it.
2. A successful echo-only KLQR controller has not been established by this migration alone. The next prerequisite is an offline, data-only learned actuation map compatible with the fixed history projector and sampling interval.
3. The exact RPE echo formula is still convenient for producing simulated training observations. A laboratory deployment should replace that input boundary with measured echo samples while keeping the downstream history, EDMD, and evaluation pipeline unchanged.
