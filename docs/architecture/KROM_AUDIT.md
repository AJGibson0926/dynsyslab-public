# KROM closure, differential and control audits

The compiled `krom-audit` query diagnoses a frozen model without refitting it.
It complements exact [Lie analysis](LIE_ALGEBRA.md), linear-systems analysis,
and ordinary solve/post-process jobs. All numerical work lives in native services;
Python is an optional JSON client.

```powershell
dynsyslab-cli krom-audit query.json
python -m dynsyslab_client krom-audit query.json
```

```python
from dynsyslab_client.workflow import NativeWorkflow
result = NativeWorkflow().krom_audit({
    "action": "differential", "krom": "model.npz",
    "points": [1.0, 0.0, 1.2, 0.1], "input_map": "constant"
})
```

The request schema is `dynsyslab.krom-audit-query/v1`; the Python method supplies
it automatically. Results use `dynsyslab.krom-audit-result/v1`. Complex values
are `[real, imaginary]`; other matrices are flat row-major arrays. Nonfinite
summary quantities, including unbounded margins, serialize as null, never as a
finite certificate. Preserve the request alongside its result.

| Action | Inputs beyond `krom` | Meaning |
|---|---|---|
| `derivative_closure` | `extend` (default true), optional Lie `limits` | Exact spatial-derivative closure/completion of the stored polynomial bank, with a constant and dependent split rows removed. |
| `differential` | Flat physical `points`, `input_map` (default `exact`) | Analytic physical Jacobian, singular values, `L_f Psi - A Psi`, true `J_Psi g`, chosen map and its error. |
| `regress_action` | `features`, `outputs`, flat `training_features`, `training_action`, `validation_features`, `validation_action`, optional `tolerance` | Native SVD least squares for action in the affine span `[1,z]`, with separate validation error. No archive needed. |
| `local_control` | `point`, scalar `measurement`, `input_map`, `Q`, `R`, `W`, `V`, `holding_input`, `frequency_points` | Same frozen LQR/LQG design connected to learned and physical local plants; affine consistency diagnostics. |
| `trajectory` | Physical trajectory `source`, `input_map` (`constant` or `none`), optional output `stride` | Measured lift, instantaneous decoding, one-step prediction and open-loop input replay through the existing native ZOH simulator. |

`derivative_closure` treats the exact binary floating coefficients of the saved
approximation as rational numbers. It proves an identity for that polynomial
bank, not for the unknown true eigenfunctions. Derivatives are in standardized
**fitting coordinates**. Completion yields a polynomial function space; its new
basis functions are not thereby Koopman eigenfunctions. For a history provider,
this is closure in the retained history variables, not in physical state.

`differential` rejects history-only providers: a partial derivative with respect
to physical state requires an explicit instantaneous map. Its `sigma_min` is a
raw, scale-dependent quantity. Finite sampled values do not prove global
injectivity, conditioning or uniform approximation bounds. For `regress_action`,
the caller must provide genuinely separate training and validation records;
both may still be inadequate to cover the claimed operating domain.

`local_control` defaults to Q=5000, R=.01, W=1, V=1e-4, 1200 frequency points.
These are scalar diagonal design weights, not a coordinate-invariant physical
cost. The physical Jacobian includes the supplied holding input. Nonzero holding
input with a state-dependent lifted map is rejected because a frozen B alone
would omit its Jacobian term. The physical full-state gain is `K J_Psi` at the
point; the output-feedback test connects the identical compensator to the actual
local sensor Jacobian. It also returns physical equilibrium residual, lifted
generator value, actual sensor value and decoded sensor/state values. Check these
before interpreting a stable deviation model as stabilization of the desired
equilibrium. Margins are continuous frozen local quantities. They do not certify
sampling, scheduling, saturation, periodic tracking or large disturbances. The
inherited `disk_*` fields are conservative return-difference-derived gain/phase
bounds, **not** a balanced disk margin for simultaneous gain and phase variation.

`trajectory` checks physical-system compatibility and a uniform model clock.
History coordinates use the saved causal projector, with warm-up, the exact
retained coordinates and the saved sampling interval. A recorded hash identity
must match; legacy providers remain legacy evidence. `stride` only reduces the
returned arrays: it never changes acquisition, history projection or propagation.
Replay starts with the first available measured lift, then uses only the recorded
input. No future measured state or lift corrects it. This is an **input replay**,
not a new closed-loop plant simulation. `none` requires zero archived input.
The one-step residual is also normalized by the observed state change, avoiding
an artificially tiny score from predicting almost unchanged values at small dt.

## Code ownership

- `KoopmanModel::stateJacobian[Split]`: existing polynomial/scattered derivative
  backend composed with the frozen coordinate observable, in `dynsyslab_koopman`.
- `KromDifferentialAnalysis.hpp`, `KromLocalControlAnalysis.hpp`: independent
  native `dynsyslab_linear_analysis` APIs. Existing numerics own SVD, CARE,
  Lyapunov operations, discretization and loop calculations.
- `evaluateKromAuditQuery`: workflow composition, archive access and serialization.
- `dsl_workflow_krom_audit_json`: optional additive ABI-v4 symbol.
- `NativeWorkflow.krom_audit`: optional thin Python adapter.
- Batch `reconstructObservable` honors saved affine offsets, consistently with
  the online observer and pointwise model readout.
- Viewer owns every figure. No audit introduces a renderer or alternate archive.

The inverse-radius proof recognizer belongs to `dynsyslab_lie`, not to RPE or a
controller. Run `examples/native/lie-algebra/rpe-pressure-algebra.json` for the
exact recurrence certificate. The RPE report
illustrates each audit and records its limits.

The theoretical distinctions follow [Katayama–Susuki](https://arxiv.org/abs/2611.20438)
on finite Lie representations, [Strässer et al.](https://arxiv.org/abs/2312.01441)
and [SafEDMD](https://arxiv.org/abs/2402.03145) on approximation errors and robust
feedback, and [Doyle's LQG counterexample](https://doi.org/10.1109/TAC.1978.1101828)
on the absence of universal LQG robustness margins. These queries do not claim
those papers' guarantees from sampled residuals alone.
