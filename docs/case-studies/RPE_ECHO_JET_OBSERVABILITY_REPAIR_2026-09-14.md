# RPE echo-jet observability repair

Date: 2026-09-14  
System: dimensionless Rayleigh-Plesset equation, air in water  
Observable coordinates: `e = h(x)` and `e1 = L_f0 h(x)`

## Verdict

The earlier echo-plus-echo-dot failure does **not** show that echo history lacks the
state information needed for control. It shows that one global polynomial
Koopman fit, trained on a violent free transient, did not produce usable
eigenfunctions.

On the moderate multisine target branch, the map

\[
  \mathcal O_2(x) = \begin{bmatrix}h(x)\\L_{f_0}h(x)\end{bmatrix}
\]

is locally full rank and well conditioned. A target-local refit in these same
two coordinates reduces the generator defect from approximately 100% to 3.78%
and restores stable KLQR tracking. The repaired controller reproduces the
inverse-dynamics driver with correlation 0.99994 and 1.06% NRMSE.

This is a **local** result. The same coordinate map has severe near-self-intersections
and near-rank loss over the original large transient, and the present exact-map KLQR
still reads the simulated physical state when it evaluates the lifted input
map. A physical echo-only implementation still needs an input-aware observer or
delay estimator.

## What the 2026-09-11 attempt actually tested

The failed experiment combined five distinct hypotheses:

1. the two-coordinate observation map contains state information;
2. a single degree-four polynomial in those coordinates represents Koopman
   eigenfunctions over the entire training domain;
3. eigenvalues extracted from another observable remain suitable;
4. the target's mapped autonomous derivative is the right second reference;
5. KLQR with an exact state-evaluated input map is numerically and physically
   usable.

Only the first hypothesis is about observability. Failure of the combined
pipeline cannot isolate it.

The coordinate named `echo_dot` was specifically the autonomous Lie derivative
`L_f0 h`. Under actuation, a laboratory time derivative is instead

\[
  \dot e = L_{f_0}h(x) + L_g h(x)u.
\]

The old run disabled input correction and mapped its target to
`[h(x_ref), L_f0 h(x_ref)]`. It therefore did not use an echo-history estimate
of the controlled output derivative and did not implement a measurement-only
observer. Its exact lifted input map also evaluated the full physical state.

## Why the global fit failed

| Diagnostic | Original global transient | Target-local trajectory |
|---|---:|---:|
| Samples | 25,001 | 10,001 |
| Radius range | 0.233 to 2.333 | 0.797 to 1.202 |
| Velocity range | -27.09 to 26.95 | -0.644 to 0.711 |
| Echo range | -4.52 to 379.83 | -0.560 to 1.615 |
| `L_f0 echo` range | -121,645 to 119,749 | -4.663 to 4.238 |
| Raw degree-four feature condition | `1.02e18` | — |
| Standardized degree-four feature condition | `7.26e5` | `71.0` |
| Median standardized observation-Jacobian condition | `346.6` | `1.45` |
| Maximum standardized observation-Jacobian condition | `2.51e7` | `2.79` |
| Generator residual / generator RMS | `1.00008` | `0.03779` |
| Generator/target alignment | `0.0161` | `0.99929` |

The old training set contains multiple physical branches at nearly the same
echo jet. With a temporal exclusion of 0.2 tau, 3,914 of 25,001 samples had a
nonlocal neighbor within 0.01 standardized echo-jet units while the states were
more than 0.5 standardized units apart. No such collision was detected on the
target curve using the same thresholds.

Normalizing the original global fit changes its coefficient units and reduces
the absolute residual from 335.3 to 40.6, but the scale-invariant defect remains
approximately 100% and alignment remains 0.016. Normalization alone therefore
does not solve the global-domain problem.

## Repaired experiment

The repair keeps the same physical coordinates and prior successful KROM
spectrum, but changes the identification domain:

- training states come from the moderate target-following trajectory produced
  by the state-basis baseline;
- autonomous generator derivatives are evaluated at those states, so the known
  training input is excluded from `f0`;
- the two echo-jet coordinates are standardized before the degree-four fit;
- KLQR uses one representative of the conjugate eigenvalue pair. Including both
  separately creates a redundant four-real-state realization and caused the
  continuous Riccati solve to fail for this refit.

Nominal tracking from `[1, 0]` gives:

| Metric | Repaired echo-jet KLQR | Old global echo-jet KLQR | State-basis baseline |
|---|---:|---:|---:|
| Radius RMSE | `0.000344` | `5.6251` | `0.001352` |
| Velocity RMSE | `0.002839` | `1.1104` | `0.004732` |
| Echo NRMSE | `0.375%` | `0.760%` | `1.563%` |
| Autonomous echo-dot NRMSE | `0.893%` | `221.0%` | `1.771%` |
| Control RMS | `0.46293` | `1.2200` | `0.45548` |
| Control peak | `1.43056` | `6.3155` | `1.40211` |

Against the inverse-dynamics input required by the prescribed multisine state
trajectory, the repaired KLQR input has 0.00492 RMSE, 1.06% NRMSE, and 0.99994
correlation over the full run. After tau 2 the driver NRMSE is 0.65%.

The repair is not merely hiding the old runaway: radius remains between 0.7953
and 1.2047 instead of growing to 10.24 by tau 10.

## Displaced-start check

Starting at `[1.1, 0]` instead of the target's `[1, 0]`, the original high-gain
weights reacquire the trajectory by tau 1.545 under the criterion
`|r-r_ref| <= 0.005` and `|v-v_ref| <= 0.02` thereafter. However, the
unconstrained controller requests a dimensionless peak input of 2,647, which is
not physically acceptable.

Raising the input penalty from `R=1e-4` to `R=10` reduces the peak to 8.09 but
makes acquisition slow and phase dependent. This is a control-constraint issue,
not evidence that the echo jet lacks local state information. A hard-constrained
KMPC or saturation-aware design is needed for acquisition outside the nominal
orbit.

## Does measured echo history contain the information?

Locally, yes—provided the applied input is known and the trajectory remains on
the target branch.

For the repaired nominal trajectory, the known-input measurement map

\[
  x \mapsto \left(h(x), L_{f_0}h(x)+L_gh(x)u\right)
\]

has standardized minimum singular value 0.311, median condition 1.406, and
maximum condition 2.720. A finite-difference derivative of the simulated echo
agrees with the analytic total derivative to relative RMS error `4.62e-6` at
the archived sampling interval. The applied-input term therefore does not
destroy local invertibility in this experiment.

The autonomous second target differs from the actual derivative of target echo
by 0.1211 RMS, or 5.99% of total derivative RMS. That difference is exactly the
known `L_g h u_ff` contribution. A real observer must account for it rather than
treating raw echo-dot as `L_f0 h`.

## Required architecture for a physical implementation

1. Measure echo and estimate its derivative with a bandwidth/noise-aware filter,
   or retain a short echo-delay window.
2. Include the applied input and its timing in the observer. Locally solve
   `[e, de/dtau] = [h(x), L_f0 h(x)+L_g h(x)u]` by continuation from the previous
   estimate, or use an EKF/MHE formulated in the physical state.
3. Monitor the smallest singular value and innovation. Do not use the local
   inverse after the estimate leaves its validated branch.
4. Evaluate the existing state KROM from the estimated state, or evaluate the
   repaired local echo-jet KROM after recovering its autonomous second
   coordinate.
5. Use a constrained controller for acquisition; reserve local KLQR for the
   verified neighborhood.
6. Retrain on independent forced trajectories that cover a tube around the
   target with multiple initial conditions, then validate on held-out forcing,
   noise, derivative filtering, model mismatch, and actuator saturation.

A longer delay vector can help resolve global branches, but Takens-style delay
embedding is not equivalent to using only one derivative. With periodic forcing,
the forcing phase and input history are part of the effective state. Delay
reconstruction, sampled-operator identification, polynomial Koopman closure,
and closed-loop controllability must remain separate tests.

## Artifacts

- Metrics: `viewer_outputs/rpe-echo-jet2-observability-repair.json`
- Figure: `viewer_outputs/rpe-echo-jet2-observability-repair.png`
- Global normalization probe:
  `examples/native/case-studies/rpe-echo-jet2-sindy-prior-r-spectrum-normalized.json`
- Target-local fit:
  `examples/native/case-studies/rpe-echo-jet2-sindy-target-local.json`
- Nominal repaired control:
  `examples/native/case-studies/rpe-echo-jet2-control-target-local.json`
- Displaced-start probes:
  `examples/native/case-studies/rpe-echo-jet2-control-target-local-displaced.json`
  and `rpe-echo-jet2-control-target-local-displaced-tuned.json`

## Method provenance

The interpretation follows the DynSys Class peer records for Takens delay
coordinates, sampled-dynamics identifiability, and observability Gramians. In
particular, it keeps separate the implications
`measurement -> state reconstruction -> sampled operator -> continuous
generator -> feedback design`; success or failure at one arrow is not evidence
for all the others.
