# RPE optional production-control sequence

**Date:** 2026-08-27  
**Status:** deferred, optional add-on  
**Default:** disabled  
**Baseline retained:** `rpe-clean-isolation-klqr`

## Decision

DynSysLab will retain the established clean-isolation RPE/KLQR workflow as the
normal workflow for current experiments. It applies the moving multisine target
immediately at equilibrium and consequently produces the known large initial
control transient. That signal remains uncapped and unfiltered in the archive.
Viewer preserves the full plot and automatically adds a detail inset so the
post-transient behavior remains readable.

No smoother engagement, feasibility rejection, actuator augmentation,
coordinate whitening, KROM retraining, or control-quality modal reselection is
part of the current default. This prevents the investigation from drifting away
from established evidence before parity and mechanism questions are settled.

## Important correction: reference shaping is not controller modulation

The retained two-`tau` reduced-transient result does **not** preserve the exact
control target during startup. It replaces the first two dimensionless time
units of the target presented to the controller with a smooth C2 transition
from equilibrium, then rejoins the original multisine. The underlying
multisine after engagement is unchanged, but the complete reference trajectory
is different. This result therefore answers a target-engagement question, not
the originally intended question of how to moderate the controller while
holding the target fixed.

It remains useful evidence and is retained as an optional reference-shaping
branch. It must not be described as demonstrating that the same target can be
tracked with a smaller transient through controller modulation alone.

## Optional sequence retained for later evaluation

The following is a staged add-on, not a replacement workflow. Stage 4A and 4B
are separate research branches and must not be conflated:

1. Run a target/actuator feasibility preflight using inverse pressure, pressure
   slew, bandwidth, and energy requirements.
2. Expand training across the intended controlled domain and require held-out
   eigenfunction, gradient, `Bpsi`, observation, and tracking validation.
3. Canonicalize or whiten lifted coordinates and derive `Qpsi` and dimensionless
   `R` from physical state and actuator costs.
4. **Branch A — optional reference shaping:** engage the state reference from the actual initial condition with a smooth
   C2 path. The evidence-backed candidate is the tested quintic engagement over
   `2 tau`; it changes the supplied reference but retains
   `u=-K(Psi(x)-Psi(x_ref))`.
5. **Branch B — target-preserving controller modulation:** keep the original
   multisine reference exactly unchanged and investigate time-varying Riccati
   or gain schedules, a finite-horizon differential Riccati solution,
   coordinate-normalized physical `Q/R`, and actuator/control-rate dynamics.
   Candidate schedules must retain the recognizable law
   `u=-K(t)(Psi(x)-Psi(x_ref))`; arbitrary clipping is excluded. No candidate
   in this branch has yet been run or validated.
6. Represent real actuator bandwidth and lag as dimensionless dynamic states
   when physical specifications are available.
7. Select KROM modes using control-domain observation accuracy, `Bpsi`
   conditioning, and closed-loop robustness in addition to energy ordering.
8. Compare every optional run against the unchanged transient-producing
   baseline using separately named policies and immutable archives.

## Evidence and acceptance gate

In the isolated experiment, the two-`tau` C2 engagement reduced uncapped peak
control from `63267.48 u*` to `21.97 u*`, while second-half radius and velocity
RMSE remained approximately `0.004326` and `0.022952`. This is the strongest
tested mitigation but remains only a candidate.

There is currently no numerical evidence that a target-preserving modulation
achieves the same reduction. That branch is a documented hypothesis for later
study, and its success must be judged independently rather than inferred from
the reference-shaped result.

An optional sequence may be promoted only after it reduces peak pressure and
pressure slew without degrading established post-engagement tracking, phase,
or stability margins, and after its KROM domain and actuator assumptions pass
explicit validation. Promotion requires a separate user decision; passing the
gate does not change defaults automatically.

## Machine-readable policy boundary

`control.policies.RPE_OPTIONAL_PRODUCTION_CONTROL_SEQUENCE` records this
sequence with `enabled_by_default=False` and `opt_in_required=True`. It is
declarative and non-executable: importing it changes no controller or solve.
The current reproduction script continues to read only
`RPE_CLEAN_ISOLATION_KLQR_SETTINGS`, whose settings preserve the startup
transient (`startup_ramp=False`, no cap, no feedforward).

## Related evidence

- `docs/RPE_CLEAN_ISOLATION_KLQR_SOP.md`
- `docs/case-studies/RPE_CLEAN_ISOLATION_INITIAL_TRANSIENT_DIAGNOSIS_2026-08-27.md`
- `docs/notes/RPE_CLEAN_ISOLATION_TRANSIENT_DIAGNOSTICS_2026-08-27.json`
- `viewer_outputs/rpe-clean-isolation-klqr-latest.html`
