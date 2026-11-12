# RPE clean-isolation KLQR initial-transient diagnosis

**Date:** 2026-08-27  
**Baseline:** `rpe-clean-isolation-klqr`  
**Control law retained throughout:** `u=-K(Psi(x)-Psi(x_ref))`  
**Excluded throughout:** feedforward, clipping, saturation, command caps, and post-hoc signal filtering

**Decision status:** diagnosis retained; mitigations deferred. The unshaped,
transient-producing workflow remains the normal baseline. The recommended
production sequence is an optional add-on requiring explicit selection and is
not active in the standard workflow.

Complete numerical evidence is in
`docs/notes/RPE_CLEAN_ISOLATION_TRANSIENT_DIAGNOSTICS_2026-08-27.json`.

## Bottom line

The initial spike is primarily a reference-initialization defect, strongly
amplified by out-of-domain polynomial extrapolation and aggressive weighting in
arbitrarily normalized eigenfunction coordinates. It is not caused by
dimensionalization, sampled control, stale `Bpsi`, loss of controllability,
Viewer scaling, or pressure saturation.

A C2 reference engagement lasting two dimensionless time units reduces the
uncapped peak from `63267.5 u*` to `21.9714 u*` without changing the feedback
law or later tracking. That final peak equals the `21.9701 u*` maximum pressure
the RPE equations predict is intrinsically required once on the moving target.

This result is evidence for later consideration, not a default change. Current
experiments preserve the full initial spike, display later detail with Viewer's
automatic inset, and otherwise ignore the spike while parity and established
behavior remain the priority.

Clarification: the C2 experiment changes the reference supplied to the
controller during `0 <= tau < 2`; it does not demonstrate target-preserving
controller modulation. The original multisine is recovered after engagement,
but the full target trajectory is not identical to the baseline target.

## 1. The controller starts against an impossible state jump

At `tau=0`:

```text
plant state  = (1,0)
target state = (4,9)
state error  = (-3,-9)
```

The target never intersects equilibrium: `1.5 <= r_ref <= 6.5`. Starting a
causal second-order plant at `(1,0)` while demanding `(4,9)` immediately makes
an impulse-like command the high-gain optimum. No tuning can make that request
both instantaneous and physically moderate.

The peak occurs exactly at `tau=0`:

```text
u*(0) = -63267.48492769027
p_control(0) = -6.326748492769027 GPa
```

Moreover, `99.9635%` of the baseline integral `int u*^2 dtau` occurs during
the first single dimensionless time unit.

## 2. Polynomial extrapolation multiplies the error

Training domain:

```text
0.2327207 <= r <= 2.3333866
-27.09322 <= u <= 26.94802
```

The initial target violates the radius domain by a relative excursion of
`0.79337`; the complete target reaches `r=6.5`. Velocity is covered, but most
of the desired radius orbit is not.

The eigenfunction has training RMS `0.03191` and training peak `0.14197`, while
its extrapolated initial target value is about `6.86` in magnitude. The lifted
error becomes

```text
Psi(x)-Psi(x_ref) = (-6.14457, 3.04800).
```

Thus polynomial extrapolation turns an ordinary state mismatch into a lifted
error roughly 48 times the largest fitted training magnitude.

## 3. Aggressive cost in noncanonical lifted coordinates multiplies it again

At equilibrium:

```text
Bpsi = (-8.5873e-4, 2.09223e-2)^T
K    = (-6825.8940, 6996.4876)
||K|| = 9774.6441
closed-loop poles = -4.83348, -147.83228 per tau
```

The fast time constant is only `0.00676 tau`. The pair is controllable and its
controllability condition number is a benign `4.83`, so singular `Bpsi` is not
the cause. The spike comes from a large extrapolated lifted error multiplied by
an aggressive gain.

SINDy eigenfunctions have arbitrary multiplicative normalization. Therefore
`Q=5000 I` in raw polynomial coordinates is not a physically invariant tuning
statement. Rescaling the eigenfunction changes the effective tradeoff unless
`Q` is transformed with it.

## 4. The moving target itself requires substantial pressure

Direct inversion of the dimensionless RPE along the target gives:

```text
RMS required pressure  = 9.56618 u*
peak required pressure = 21.97005 u* = 2.19701 MPa
```

With pure error feedback, exact `Psi(x)=Psi(x_ref)` implies `u=0`. A forced
moving reference therefore cannot have mathematically zero tracking error under
this law; a small lag/error must remain to generate pressure. This explains the
small persistent later error, not the enormous startup spike.

If `2.2 MPa` exceeds the actuator's continuous capability, the target itself is
physically infeasible. No uncapped feedback tuning can change that fact.

## 5. Uncapped isolation experiments

Every case uses the same KROM, dimensionless plant, BDF solver, continuous
state-dependent `Bpsi`/CARE, and feedback law. Only the declared engagement
duration or `R` changes.

| Case | Peak `|u*|` | Control RMS | `int u*^2 dtau` | Radius RMSE, second half | Velocity RMSE, second half |
|---|---:|---:|---:|---:|---:|
| Unshaped, `R=1e-4` | 63267.48 | 632.71 | 2,002,219.6 | 0.004326 | 0.022952 |
| C2 engagement, `0.5 tau` | 96.92 | 16.11 | 2,594.4 | 0.004326 | 0.022951 |
| C2 engagement, `1 tau` | 27.18 | 9.72 | 944.7 | 0.004326 | 0.022952 |
| C2 engagement, `2 tau` | 21.97 | 8.45 | 713.4 | 0.004326 | 0.022951 |
| C2 engagement, `3 tau` | 21.97 | 8.20 | 671.9 | 0.004326 | 0.022952 |
| Unshaped, `R=1e-3` | 19474.47 | 194.94 | 190,445.3 | 0.014606 | 0.080975 |
| Unshaped, `R=1e-2` | 5634.25 | 57.56 | 17,259.3 | 0.214049 | 3.642960 |

Increasing `R` alone leaves a huge initial command and degrades tracking. At
`R=1e-2`, radius falls to `0.1473`; weakening feedback against the impossible
initial reference worsens the nonlinear excursion instead of making it safe.

The two-`tau` engagement is the shortest tested transition whose peak is no
larger than the target's ongoing inverse-dynamics requirement. It removes over
99.96% of baseline control effort while preserving later tracking.

## 6. Optional production sequence (deferred; never selected implicitly)

The following recommendations are retained as a coherent opt-in sequence for
future evaluation. They are not prerequisites for, and must not be injected
into, the current baseline reproduction. Any trial must use a distinct policy
name, archive label, and Viewer output and must compare against the unchanged
baseline.

Two future branches are intentionally distinct:

- **Reference shaping:** the tested C2 engagement modifies the startup portion
  of the reference and already has numerical evidence.
- **Target-preserving controller modulation:** keep the multisine target
  identical while studying scheduled `Q/R` or `K(t)`, a differential Riccati
  solution, coordinate-normalized costs, or explicit actuator/control-rate
  dynamics. This branch has not yet been tested. It should preserve the form
  `u=-K(t)(Psi(x)-Psi(x_ref))` and must not rely on arbitrary command caps.

### Optional stage 1: dynamics-aware reference engagement

Do not ramp or clip pressure. Construct the reference supplied to the unchanged
law so it begins at the actual IC with matching derivatives and joins the
requested orbit smoothly. The tested quintic C2 engagement is

```text
s(q)=10q^3-15q^4+6q^5, q=tau/Tengage
r_ref=1+s(r_desired-1)
u_ref=s'(r_desired-1)+s u_desired.
```

For this case, `Tengage=2 tau` is the evidence-based minimum tested. In
production, choose the shortest duration that passes plant/actuator feasibility
checks for pressure, pressure slew, bandwidth, and energy. Those specifications
shape the reference; they do not clip controller output.

### Optional stage 2: target-feasibility preflight

Before solving:

1. compare the complete target envelope and engagement path with the validated
   KROM domain;
2. evaluate plant inverse pressure along the target;
3. compare pressure, slew, bandwidth, and energy with an explicit actuator
   model/specification;
4. reject or retime infeasible targets rather than discovering infeasibility
   through a spike.

### Optional stage 3: train over the controlled domain

The free trajectory does not cover target radius. Generate an explicit campaign
covering the engagement path and at least `1.5 <= r <= 6.5`, with margin.
Validate held-out eigenfunction residuals, gradients, `Bpsi`, observability, and
tracking throughout that domain. Do not promote extrapolating KROMs merely
because one solve succeeds.

### Optional stage 4: make Q/R coordinate-invariant and physical

Canonicalize or whiten retained eigenfunction coordinates using training data.
Transform `Q` so equivalent eigenfunction rescalings produce the same
controller. Prefer deriving lifted cost from physical-state cost through a
validated decoder/Jacobian:

```text
delta x ~= Cpsi delta Psi
Qpsi = Cpsi^T Qx Cpsi.
```

Derive dimensionless `R` from actuator pressure/energy specifications using
`p_control=p0 u*`. Tune only after normalization. This preserves continuous LQR
and the exact feedback-law form.

### Optional stage 5: model actuator physics

The current plant treats pressure as instantaneously realizable. Add real
driver bandwidth, lag, pressure dynamics, and energy behavior as dimensionless
dynamical states and design LQR on that augmented system. The controller remains
state feedback, while physical pressure becomes dynamically realizable rather
than an ideal algebraic input.

### Optional stage 6: control-quality KROM selection

Select modes using held-out observability, gradient accuracy, `Bpsi`
conditioning, and closed-loop robustness across the controlled domain—not
energy order alone. The present pair is controllable at equilibrium, so mode
count is not the primary startup cause, but validated additional modes may
improve off-equilibrium robustness.

## 7. Explicit non-recommendations

- Do not clip or saturate the command to conceal startup behavior.
- Do not apply an arbitrary pressure ramp after the controller.
- Do not address the discontinuous target solely by increasing `R`.
- Do not treat successful polynomial extrapolation as validation.
- Do not add feedforward or alter the feedback law without a separate explicit
  design decision.
