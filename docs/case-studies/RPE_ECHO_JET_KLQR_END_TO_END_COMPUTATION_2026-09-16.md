# RPE echo-jet KLQR: end-to-end computation audit

Date: 2026-09-16  
System: dimensionless Rayleigh-Plesset equation, 10 micrometer air bubble in
water at 20 C  
Controlled cases: 100 kHz, 20 kPa and 40 kPa pure-sine target orbits

> **Direct answer:** yes. The target supplied to the final controller contains
> both `echo` and **autonomous echo-dot**, in that order. For every target time,
> the archived reference state `(r_ref,v_ref)` is mapped to
> `H(x_ref)=[e(x_ref),L_f0 e(x_ref)]`. Both values are then passed through the
> learned polynomial eigenfunction. KLQR penalizes the resulting two real
> lifted-coordinate errors, not the two raw observable errors independently.

This report reconstructs that statement from the saved jobs, numerical arrays,
and the native implementation. It also separates three ideas that are easy to
conflate:

1. the **two-component coordinate chart** used to fit the eigenfunction;
2. the **two-component mapped target** supplied to the eigenfunction; and
3. the **two real KLQR states**, `Re(psi)` and `Im(psi)`, on which `Q` acts.

## 1. The complete pipeline at a glance

| Stage | Input | Computation | Output used next |
|---|---|---|---|
| Broad free training | 25,001 unforced RPE states | Radius Hankel-DMD | 2,000 candidate discrete eigenvalues |
| State-coordinate fit | Same free states and candidate spectrum | Degree-four SINDy eigenfunctions in `(r,v)` | Compact parent KROM containing one nonredundant complex pair |
| Branch generation | Parent KROM plus a moderate multisine state target | Exact-map KLQR on the physical RPE | 10,001 forced states covering the intended operating branch |
| Echo-jet refit | Branch states plus the retained eigenvalue pair | Degree-four SINDy in `[e,L_f0 e]` using autonomous derivatives | Final echo-coordinate KROM |
| Target generation | Physical RPE with a known 100 kHz pressure sine | Long direct integration at 20 or 40 kPa | Reference state trajectory `(r_ref,v_ref)` |
| Target mapping | Reference state at each requested time | `(r_ref,v_ref) -> [e_ref,e1_ref] -> psi_ref -> [Re psi_ref,Im psi_ref]` | Lifted KLQR reference `z_ref` |
| Online control | Current state, mapped target, final KROM | Exact `B_z(x)`, CARE solve, `u=-K(z-z_ref)` at each RHS call | Controlled RPE trajectory and feedback pressure |
| Assessment only | Controlled pressure and known target driver | Fundamental fit and waveform comparison | Driver-recovery metrics and Viewer figures |

The two sine-target trajectories are not training data. They enter only after
the final KROM has been fixed.

## 2. Physical state, time, and actuation

The model state is

\[
  x=(r,v),\qquad r=R/R_0,\qquad v=dr/d\tau,
\]

where `R0=10 micrometers` and one model-time unit is
`T=3.0123833038 microseconds`. The plant is

\[
\begin{aligned}
  \dot r &= v,\\
  \dot v &= \frac{Z(r,v)-P u}{r},\\
  Z(r,v) &= G r^{-4}-\frac32v^2+V-P-\frac Cr-\frac{Dv}{r}.
\end{aligned}
\]

The control input is normalized by the 100 kPa ambient pressure. Therefore
`u=0.2` is 20 kPa and `u=0.4` is 40 kPa. Its control-affine direction is

\[
  g(x)=\begin{bmatrix}0\\-P/r\end{bmatrix}.
\]

This sign convention makes positive pressure compress the bubble.

## 3. First training data: excite a broad spectrum

The first dataset is one unforced transient beginning at
`x(0)=[2.3333866,0]` in dimensionless coordinates. The legacy source stores
physical time, radius, and wall velocity; post-processing converts them
explicitly to `tau`, `r`, and `v` without changing the source bytes. The
converted training record has:

- 25,001 samples over `tau=0..25`;
- uniform `Delta tau=0.001`;
- radius range `0.232721..2.333387`; and
- radial-velocity range `-27.0932..26.9480`.

This violent ring-down is useful because it excites many decaying oscillatory
components. It is not a good global chart for the final echo-polynomial model:
the mapping to `[e,L_f0 e]` folds and becomes poorly conditioned over this
large domain.

![Broad free training data and its coordinate coverage](../../viewer_outputs/rpe-echo-jet2-free-training-evidence.png)

## 4. Hankel-DMD supplies frequencies, not the final eigenfunction

Hankel-DMD is applied to the scalar radius history. A 20,000-sample window is
embedded into 2,000 delay coordinates with stride one; projected DMD retains
rank 2,000. The result is a candidate discrete spectrum `mu_j`.

The postprocessor converts each candidate to a continuous-time generator value

\[
  \lambda_j=\frac{\log\mu_j}{\Delta\tau}.
\]

Hankel-DMD therefore answers: *which exponential rates and frequencies are
present in the free ring-down?* It does not yet produce the instantaneous
echo/eigenfunction map used by control.

The candidate spectrum contains 2,000 eigenvalues. A state-coordinate SINDy
stage next fits a polynomial eigenfunction for each candidate using `(r,v)` as
its variables. Its library contains the 14 nonconstant monomials of total
degree one through four. For candidate `lambda`, it solves for a coefficient
vector in the approximate null space of

\[
  \left(\Gamma-\lambda\Theta\right)c,
  \qquad
  \Gamma=\nabla_x\Theta(x) f_0(x).
\]

An energy-ordered extraction keeps source mode 1 and its conjugate partner.
The continuous pair is

\[
  \lambda=-0.2106164458\pm4.8186345616i.
\]

The imaginary part corresponds to about 254.6 kHz in physical time. It is not
the later 100 kHz target frequency. The controller is therefore not recovering
the driver by selecting a Koopman eigenvalue at 100 kHz.

This compact state-coordinate KROM is an intermediate tool. It is used to make
better training data for the final echo-coordinate KROM.

## 5. Second training data: cover the branch to be controlled

The compact parent KROM drives the physical RPE along the moderate state target

\[
\begin{aligned}
 r_{ref}(\tau)&=1+0.09\sin(4\tau)-0.12\sin(3\tau),\\
 v_{ref}(\tau)&=0.36\cos(4\tau)-0.36\cos(3\tau).
\end{aligned}
\]

The resulting forced trajectory contains 10,001 states over `tau=0..10`, with
`r=0.7975..1.2022` and `v=-0.6444..0.7114`. The actual trajectory, rather than
the analytic target formula, becomes the second training set. It samples the
moderate branch on which later control will operate.

![Forced branch-covering training data](../../viewer_outputs/rpe-echo-jet2-branch-training-evidence.png)

Although this dataset was produced by feedback, the final SINDy refit sets
`include_input_terms=false`. At each sampled state it evaluates the autonomous
RPE vector field `f0`, not the control-modified velocity that happened to carry
the trajectory there. The forcing chooses *where the model is sampled*; it
does not redefine the autonomous Koopman generator being fitted.

## 6. What `echo_jet2` actually means

The final coordinate observable is

\[
  H(x)=\begin{bmatrix}e(x)\\e_1(x)\end{bmatrix},
  \qquad e_1=L_{f_0}e.
\]

The echo is the analytic far-field pressure-distance observable

\[
  e=R_0P_s\left(2rv^2+Z(r,v)\right).
\]

Its autonomous Lie derivative is evaluated from the state, not by finite
differencing a saved echo trace:

\[
  e_1=e_r v+e_v\frac{Z}{r},
\]

with

\[
\begin{aligned}
 e_r &=R_0P_s\left(2v^2-4Gr^{-5}+Cr^{-2}+Dvr^{-2}\right),\\
 e_v &=R_0P_s\left(4rv-3v-D/r\right).
\end{aligned}
\]

Two distinctions matter:

1. `e1` is the autonomous quantity `L_f0 e`. Under control, the total measured
   derivative is `de/dtau=L_f0 e+L_g e u`.
2. On a forced target orbit, `e1(x_ref)` is likewise not generally the time
   derivative of the target echo waveform. It is the autonomous derivative
   evaluated at the target state.

Thus “echo-dot” in this experiment is an instantaneous state observable with
derivative units in the model clock. It is not a numerical differentiator in
the feedback loop.

## 7. Fitting the final eigenfunction in the echo chart

The 10,001 branch samples are mapped to `y=H(x)`. Their two components are
standardized with

```text
center = [ 0.1472323322, -0.0894546951 ]
scale  = [ 0.5530497264,  1.9641872085 ]
```

Let `eta=(y-center)/scale`. The final library is again the 14 nonconstant
monomials of total degree at most four,

\[
  \Theta(\eta)=
  [\eta_1,\eta_2,\eta_1^2,\eta_1\eta_2,\ldots,\eta_2^4].
\]

Because the coordinate map is analytic, the training derivative is computed by
the exact push-forward

\[
  \dot y=J_H(x)f_0(x)=\begin{bmatrix}e_1(x)\\L_{f_0}^2e(x)\end{bmatrix}.
\]

The feature derivative is

\[
  \Gamma=\nabla_y\Theta(y)\dot y.
\]

For the retained `lambda`, SINDy finds the normalized sparse vector `c` that
minimizes the residual of

\[
  (\Gamma-\lambda\Theta)c\simeq0.
\]

The complex eigenfunction used by control is therefore

\[
  \psi(x)=\phi(H(x))=\Theta\!\left(\frac{H(x)-center}{scale}\right)c.
\]

The fit retains the conjugate pair in the archive, but the final controller
selects only mode 0. One nonreal mode already becomes two real controller
coordinates. Selecting its stored conjugate as a second control mode would
duplicate the same real subspace.

The accepted refit has 3.779% relative generator defect and 0.999288 alignment
with the generator target on this branch. Its median and maximum standardized
observation-Jacobian condition numbers are 1.451 and 2.786.

![Spectrum selection, normalized echo chart, learned eigenfunction, and generator closure](../../viewer_outputs/rpe-echo-jet2-spectrum-fit-evidence.png)

## 8. Exactly how the target is constructed

The 20 and 40 kPa targets are physical RPE trajectories generated beforehand
with

\[
  u_*(\tau)=a\cos(1.8927362514\tau),
  \qquad a=0.2\ \text{or}\ 0.4.
\]

The KLQR jobs identify the target as `space="observable"` and
`kind="from_state"`. They name only the target archive's time and state arrays.
At a controller evaluation time, the implementation performs

\[
  x_{ref}(\tau)
  \longrightarrow
  y_{ref}(\tau)=H(x_{ref})=
  \begin{bmatrix}e_{ref}\\e_{1,ref}\end{bmatrix}
  \longrightarrow
  \psi_{ref}=\phi(y_{ref})
  \longrightarrow
  z_{ref}=\begin{bmatrix}\Re\psi_{ref}\\\Im\psi_{ref}\end{bmatrix}.
\]

The controlled archives independently preserve this mapping. The 20 kPa result
contains a `15001 x 2` observable-target array; the 40 kPa result contains a
`30001 x 2` array. In both cases its names are exactly:

```text
far_field_acoustic_echo
far_field_acoustic_echo_dot
```

Recomputing both columns from the archived state target agrees to
`6.7e-16` for 20 kPa and `1.1e-15` for 40 kPa. This is numerical roundoff, so
there is no ambiguity about what was supplied.

Both columns materially affect the lift. As a counterfactual sensitivity test,
freezing autonomous echo-dot at its training center changes the settled target
lift by 38.0% of its centered RMS at 20 kPa and 42.8% at 40 kPa. Freezing echo
instead changes it by 112.0% and 115.8%, respectively. These are not additive
importance fractions; they show that neither coordinate is ignored by the
polynomial map.

![20 kPa observable target and lifted target](../../viewer_outputs/rpe-echo-jet2-20kpa-target-lift-diagnostic.png)

![40 kPa observable target and lifted target](../../viewer_outputs/rpe-echo-jet2-40kpa-target-lift-diagnostic.png)

The bottom-right panels show the quantities to which `Q` is actually applied:
the errors in `Re(psi)` and `Im(psi)`. The raw echo and echo-dot errors influence
the controller through the nonlinear polynomial `phi`; they do not receive
separate diagonal entries in `Q`.

## 9. Exactly how one control value is computed

For the current state, the same mappings are used:

\[
  y=H(x),\qquad \psi=\phi(y),\qquad
  z=\begin{bmatrix}\Re\psi\\\Im\psi\end{bmatrix}.
\]

For `lambda=alpha+i beta`, the one-mode split-real autonomous generator is

\[
  A=\begin{bmatrix}\alpha&-\beta\\\beta&\alpha\end{bmatrix}.
\]

The controller does not use a fitted constant input matrix. It recomputes the
exact local lifted input direction from the current physical state:

\[
  B_z(x)=\operatorname{realSplit}\!\left(
    \nabla_H\phi(H(x))\,J_H(x)\,g(x)
  \right).
\]

This chain contains all three derivatives required by the change of variables:

- `g(x)` says how pressure moves the RPE state;
- `J_H(x)` says how a state displacement moves echo and autonomous echo-dot;
- `grad_H phi` says how an echo-chart displacement moves the eigenfunction.

At every continuous BDF right-hand-side evaluation, the solver then solves

\[
  A^TS+SA-SB_zR^{-1}B_z^TS+Q=0,
\]

forms

\[
  K(x)=R^{-1}B_z(x)^TS,
\]

and returns

\[
  u(\tau)=-K(x)\left[z(x)-z_{ref}(\tau)\right].
\]

For both runs, scalar `Q=5000` becomes `5000 I_2`. The 20 kPa case uses
`R=0.0001`; the branch-stable 40 kPa case uses `R=0.006`.

The sequence is recomputed at internal adaptive trial states, not only at the
15,001 or 30,001 published output samples. Consequently, the saved control
history is a sampled record of a feedback law that was evaluated much more
often during integration.

## 10. Why the resulting pressure counts as rediscovery

The control equation contains no target pressure feedforward term. The solve
job points to `time` and `values` in the target trajectory and does not name its
stored `control_values` array. The known sine is loaded only by the later
assessment scripts for plotting and error measurement.

This has an important consequence at the initial instant. If the controlled
state equals the target state, then `z-z_ref=0`, so KLQR requests `u=0`, even
though the target orbit is about to move under its hidden sine forcing. As the
unforced controlled plant starts to separate from the moving target, feedback
error develops and the controller generates the pressure needed to stay on the
orbit. That mechanism explains both the acquisition transient and the fact that
the recovered steady pressure is close to, but not identically equal to, the
known sine.

The experiment is therefore stronger than fitting a sine to an echo record,
but narrower than general blind input reconstruction. The target state history
and exact plant/KROM derivatives are supplied; only the original driver history
is withheld.

## 11. Controlled results

All settled metrics below use `tau=100..300`.

| Metric | 20 kPa | 40 kPa |
|---|---:|---:|
| Known fundamental | 20.000 kPa | 40.000 kPa |
| Recovered fundamental | 20.215 kPa | 39.477 kPa |
| Fundamental phase error | -0.00418 deg | +0.31298 deg |
| Full-waveform driver NRMSE | 2.1158% | 18.8599% |
| Driver correlation | 0.999838 | 0.982380 |
| Radius RMSE | 0.00002199 `R0` | 0.002014 `R0` |
| Radial-velocity RMSE | 0.0003714 `R0/T` | 0.019758 `R0/T` |
| Echo NRMSE | 0.0640% | 3.8217% |
| Autonomous echo-dot NRMSE | 0.5520% | 11.5698% |
| Largest settled command | 20.468 kPa | 100.365 kPa |

The 20 kPa response is essentially locked to the target after acquisition. Its
settled lifted-coordinate RMS errors are `2.28e-5` in the real component and
`4.30e-5` in the imaginary component.

![Complete 20 kPa controlled result](../../viewer_outputs/rpe-echo-jet2-20kpa-success-detail.png)

The 40 kPa response remains on the intended long-term orbit, but the feedback
signal contains narrow phase-locked corrections. Its settled lifted RMS errors
are `1.44e-3` and `7.10e-4`. The correct 100 kHz fundamental is recovered, but
the spikes account for the larger full-waveform error and peak command.

![Complete 40 kPa controlled result](../../viewer_outputs/rpe-echo-jet2-40kpa-success-detail.png)

## 12. What this does and does not prove

The experiment establishes the following:

- `[e,L_f0 e]` is a useful local coordinate chart on the trained branch.
- A degree-four eigenfunction in that chart can approximately satisfy the
  selected Koopman generator relation.
- Mapping both current and target states through that eigenfunction provides a
  two-real-state KLQR tracking signal.
- Exact state-dependent lifted actuation can maintain held-out 20 and 40 kPa
  periodic orbits for long horizons.
- The dominant feedback pressure can recover the withheld sine driver.

It does **not** yet implement an echo-only laboratory controller:

1. The current simulated state `(r,v)` is used to evaluate `H(x)` and the exact
   `B_z(x)`. No observer reconstructs these quantities from a sensor stream.
2. The target is stored as a state trajectory and then mapped into observables.
   A laboratory target specified only as a desired echo waveform would need a
   compatible echo-dot target and a choice among any states sharing that output.
3. `L_f0 e` must be separated from the total measured derivative under input.
   A sensor implementation needs an input-aware estimator or the correction
   `L_f0e=de/dtau-L_ge\,u` with enough information to evaluate `L_g e`.
4. The controller has no pressure saturation, sensor noise, derivative filter,
   delay, parameter uncertainty, or transducer dynamics. The 20 kPa acquisition
   peak and 40 kPa correction spikes are not yet actuator-feasible.

The correct interpretation is therefore: this is a successful state-aware,
exact-map validation of an eigenfunction whose *coordinates* are echo and
autonomous echo-dot. It is strong evidence that the repaired representation is
useful, but it is not yet a complete output-feedback laboratory realization.

## 13. Implementation anchors and reproduction

The decisive native implementation points are:

| Operation | Source anchor |
|---|---|
| Map a `from_state` observable target through the KROM coordinate observable | `workflow/native/src/ControlFactory.cpp:313-321` |
| Map the resulting two coordinates through the polynomial eigenfunction | `native/src/LQR.cpp:43-55` |
| Recompute exact `grad(phi) J_H g` from the current state | `native/src/KoopmanModel.cpp:189-215` |
| Re-solve continuous LQR and apply `-K(z-z_ref)` | `native/src/LQR.cpp:56-64` |
| Fit SINDy directly in observable coordinates | `native/src/SINDy.cpp:248-267` |
| Select whether source input participates in generator derivatives | `workflow/native/src/Postprocess.cpp:259-261` |
| Archive the two-component observable target | `workflow/native/src/Solve.cpp:330-345` |

Reproduction jobs:

```text
examples/native/case-studies/rpe-echo-jet2-sindy-target-local.json
examples/native/case-studies/rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300.json
examples/native/case-studies/rpe-echo-jet2-klqr-rediscover-100khz-40kpa-r6e-3-tau300.json
```

Analysis and Viewer rendering:

```powershell
python tools/plot_rpe_echo_jet_pipeline_evidence.py
python tools/plot_rpe_echo_jet_target_lift_diagnostics.py
python tools/plot_rpe_echo_jet_low_amplitude_detail.py
python tools/plot_rpe_echo_jet_success_detail.py
```

Machine-readable target-lift checks are retained in:

```text
viewer_outputs/rpe-echo-jet2-target-lift-diagnostics.json
```
