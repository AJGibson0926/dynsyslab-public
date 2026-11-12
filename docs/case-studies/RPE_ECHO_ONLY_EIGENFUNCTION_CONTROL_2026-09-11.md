# Echo-only SINDy eigenfunctions and RPE control

Date: 2026-09-11. Native DynSysLab on main; MATLAB was not modified.

## Result

The requested prior-r-spectrum experiment is complete. Echo-only polynomial
functions fitted with the previously successful r-H-DMD eigenvalues track the
mapped echo target closely, but they do not reproduce the successful state
trajectory. Radius reaches 10.279 times its equilibrium value. Independent
500-column and 2,000-column echo-H-DMD spectra produce essentially the same behavior.

All three echo-basis fits and their control comparisons completed successfully
to tau = 10. Successful integration does not imply successful state control.

![Control comparison](../../viewer_outputs/rpe-echo-only-control-comparison.png)

All errors below use the same 5 <= tau <= 10 interval. Echo NRMSE is echo RMSE
divided by the RMS of the supplied echo target, expressed as a percentage.

| Eigenfunction coordinates / spectrum | Echo NRMSE | Radius RMSE | Velocity RMSE | Maximum radius | RMS input | Peak absolute input |
|---|---:|---:|---:|---:|---:|---:|
| State / prior r-H-DMD baseline | 1.5713% | 0.001213 | 0.004121 | 1.2022 | 0.4555 | 1.4021 |
| Echo / prior r-H-DMD pair | 0.022030% | 7.4424 | 0.8476 | 10.2785 | 1.2380 | 7.6324 |
| Echo / echo-H-DMD, 500 columns | 0.022034% | 7.4434 | 0.8475 | 10.2794 | 1.2385 | 7.6878 |
| Echo / echo-H-DMD, 2,000 columns | 0.022030% | 7.4423 | 0.8476 | 10.2784 | 1.2380 | 7.6304 |

Radius and velocity are dimensionless, and input is p_control/p0. Input statistics
and maximum radius use the whole interval 0 <= tau <= 10. The three echo-controlled
trajectories almost overlap in the figure. The fresh state-basis baseline closely
reproduces the previously successful exact-map reference.

## Source, fit and spectrum selection

The source is the main RPE reference training trajectory: 25,001 samples over
0 <= tau <= 25, with delta tau = 0.001 and initial state [2.333386598601735, 0].
Its original SHA-256 is

e4529d55dd549e12c66dbb78c9521c6bc273aa8bb135637319568a2edfdc9eba.

The legacy archive contains physical seconds, radius and velocity, despite its
dimensionless integration-form label. Native import now converts using its saved
R0 = 1e-5 m and T = 3.0123833037972497e-6 s. The archived parameter overrides,
including kappa = 4/3, are preserved; the source archive is unchanged.

For the explicitly requested comparison, the generator eigenvalues are retained
exactly from the old compact KROM:

\[
\lambda=-0.21061644578341265\pm4.818634561556464i.
\]

For the separate echo-H-DMD run, the Hankel shape is 20,000 by 500, stride 1,
full retained rank and SVD tolerance 1e-10. From that spectrum, select the
positive-imaginary eigenvalue nearest the old successful eigenvalue and its
conjugate. This gives -3.6554268101 +/- 6.3863922908 i. The selection rule is
spectral proximity, not post hoc optimization of control outcomes.

The 20,000 by 2,000 echo-H-DMD matrix retains rank 2,000. Applying the same
selection rule gives -0.10497352123 +/- 4.99950227211 i, again modes 0 and 1.
The full spectrum was fitted before extracting this compact pair.

All three fits use only instantaneous echo e as their polynomial variable:

\[
\phi(e)=c_1e+c_2e^2+c_3e^3+c_4e^4.
\]

There is no constant term, delay coordinate, radius coordinate or velocity
coordinate in this polynomial library. SINDy uses threshold 1e-6, one iteration,
no coordinate normalization, conjugate reuse, and analytical observable
Jacobians composed with the known RPE RHS. The extracted archives retain two
state columns for provenance, one training-coordinate column, and one exponent
column. Their evaluation functions depend only on e.

For the positive-imaginary member of the prior-r-spectrum fit:

| Power | Real coefficient | Imaginary coefficient |
|---|---:|---:|
| e | 0.9999300818439264 | -6.659453798795293e-5 |
| e^2 | -0.011824749326952573 | 3.1138342552581863e-7 |
| e^3 | 4.7875049641960246e-5 | -6.0701815364294594e-15 |
| e^4 | -6.189648623158826e-8 | -1.8792920888247764e-12 |

The other member has conjugate coefficients.

## Exact meaning of the supplied target

This experiment uses DynSysLab's existing state-defined echo observable, in Pa m:

\[
h_0(r,v)=R_0P_s r\left(2v^2+\frac{Z(r,v)}r\right),\qquad
Z=p_g r^{-3\kappa}-\frac32v^2+p_v-p_0-\frac{c}{r}-\frac{\nu v}{r}.
\]

Here v = dr/dtau, Z/r is the zero-input acceleration, and the pressure groups,
c and nu are the native model's dimensionless coefficients. The expression is
retained verbatim in the KROM so that its identity, coefficients and derivatives
remain reproducible after reload.

The previously successful state reference is

\[
r_{\rm ref}=1+0.09\sin(4\tau)+0.12\sin(3\tau+\pi),\qquad
v_{\rm ref}=0.36\cos(4\tau)+0.36\cos(3\tau+\pi).
\]

The actual supplied scalar target is h0(r_ref, v_ref), through target.space =
observable and target.kind = from_state. It is not a radius target mislabeled
as an echo target. Every saved scalar target sample agrees exactly with
re-evaluation of the retained expression on the state reference (maximum
floating-point difference 0 in all three echo runs).

The name “echo” here refers to this existing state observable. It does not
substitute acceleration from the controlled RHS. If the same expression were
evaluated with that acceleration, its value would acquire the direct-input term
R0 * ambient_pressure * pressure_input_sign * u, equal to -u Pa m in this study.
The measured errors in the table concern h0, not that input-dependent quantity.

All compared controllers start at [1,0] and use KLQR with the exact lifted input map,
Q = 5000 I, R = 1e-4, continuous-RHS control evaluation, and BDF integration
with relative/absolute tolerance 1e-9 and maximum step 0.001. There is no control
cap, feedforward, or target shaping. Q acts in real split eigenfunction
coordinates, so the same numerical Q does not imply the same physical cost
when the basis changes.

## Why good echo tracking is insufficient

The fitted functions have large generator defects on their own training data.
For dpsi = J_phi(h0) J_h0 f0, define the scale-independent diagnostic

\[
D=\frac{\|\dot\psi-\lambda\psi\|_2}
        {\|\dot\psi\|_2+\|\lambda\psi\|_2}.
\]

D is 0.96850 for the prior-r-spectrum echo fit and 0.95277 for the selected
500-column echo-spectrum fit. The corresponding generator residual RMS values
are 407.93 and 408.48. The 2,000-column echo-spectrum fit gives D = 0.96742
and residual RMS = 407.97. These numerical fits are poor Koopman eigenfunction
approximations on this trajectory, even though feedback tracks the output.

The expression itself supplies a stronger structural check than comparing
similar-looking curves. With r = 1, both v = 0 and v = 2 nu =
0.24147361152683358 give h0 = 0 algebraically. Compiled evaluations give values
within 1.4e-16 of zero, but their zero-input time derivatives are approximately
0 and -1.0461073343. Thus h0 does not admit a globally single-valued autonomous
evolution law on the full state space. This does not prove that every local or
specialized eigenfunction fit must fail; it explains why a delayed H-DMD
spectrum need not admit an accurate instantaneous-echo-only representation.

The echo-controlled radius first exceeds the training maximum near tau = 1.81;
about 81.9% of saved samples then lie above that radius. The experiment therefore
also involves substantial state-domain extrapolation. It shows output tracking
with poor internal-state behavior, not successful reproduction of the original
state-control result.

The eigenfunctions are echo-only, but the tested exact input map is
B_psi(x) = J_phi(h0(x)) J_h0(x) B(x). Its evaluation still uses the plant state.
This is not a demonstration of an echo-sensor-only controller. An approximate
input map would need to be fitted afresh in this echo basis; a map trained for
the old state-coordinate KROM cannot be transplanted unchanged.

## Implementation, reproduction and inspection

The [native API guide](../architecture/NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md)
describes symbolic spatial derivatives, repeated Lie derivatives using the
native equations of motion, serialized coordinate expressions, generic
observable SINDy coordinates, and KLQR/KMPC observable targets. These operations
are available through the C++ workflow, C ABI, CLI, and optional thin Python
client. The solver app exposes observable targets and the postprocessor exposes
independent H-DMD/SINDy observable selections. MATLAB remains untouched.

Reusable jobs are in examples/native/case-studies:

- rpe-echo-sindy-prior-spectrum.json and rpe-echo-control-prior-spectrum.json.
- rpe-echo-hdmd500-sindy.json and rpe-echo-control-hdmd500.json.
- rpe-echo-hdmd2000-sindy.json and rpe-echo-control-hdmd2000.json.
- rpe-echo-extract-hdmd500.json and rpe-echo-extract-hdmd2000.json select the
  comparison pairs from the full fitted models.
- rpe-echo-state-baseline.json.

Run a job with dynsyslab-cli run JOB --archive-root archives. Pipeline runs
publish fresh archive paths; point subsequent control jobs at the new compact
KROM if rerunning. The checked-in control jobs reproduce this report's saved
models. The [machine-readable comparison](../../viewer_outputs/rpe-echo-only-control-comparison.json)
contains exact source/model/solution paths, expression definitions, parameters,
fit diagnostics and metrics. The source archives are immutable.

The existing “RPE successful minimal-KROM control study” Viewer Mask remains
active and includes the new compact models and controlled trajectories. The
installed background Viewer API verified the actual 10,001-point echo target,
with no initial-condition markers on its time series. Its
[echo plot](../../viewer_outputs/rpe-echo-control-viewer.png) and
[KROM workspace](../../viewer_outputs/rpe-echo-krom-viewer.png) were rendered and
visually checked without displaying application windows.

Validation: 77 native CTest tests passed, isolated package smoke passed, and
three UI parity tests passed. A later Windows long-path publication fix passed
archive, observable-pipeline and KROM-postprocess regressions and the real
compact-KROM extraction that originally exposed the issue. Installed runtime
diagnostics report ok. Builds, commands and Viewer checks ran in the background.


## Can the pipeline use measurements alone?

Yes in principle: record t, e, edot and the applied input u. H-DMD needs only e.
An echo-only polynomial fit can use Gamma_ij = theta_j'(e_i) * edot_i, without
state values or an RPE RHS call. Controlled data then allow a regression such as

\[
\phi'(e)\dot e-\lambda\phi(e)\approx\widehat B_\psi(e)u.
\]

Constant, polynomial or lifted bilinear libraries may be used for the input
map. The free training trajectory alone cannot identify input authority. Inputs
must sufficiently excite the relevant directions and be distinguishable from
autonomous dynamics. Because the present eigenfunction defects are large,
drift-model error must be fitted or assessed separately rather than all being
assigned to the input map. See the original
[SINDYc paper](https://arxiv.org/abs/1605.06682).

For this particular h0, the same-echo counterexample above also produces
opposite input gains: L_g h0 is +0.1207368058 at [1,0] and -0.1207368058 at
[1,0.2414736115]. Multiplication by the same phi'(0) cannot eliminate this
conflict for the fitted functions. Thus a globally exact single-valued B_psi(e)
is impossible across both states. Local or distribution-dependent regression
is still possible and requires independent controlled validation.

Using edot only as derivative data leaves phi a single-variable function.
Alternatively, phi(e,edot) and measurement-dependent input maps can use both
measurements as coordinates, without requiring radius or velocity. Their
sufficiency must be checked: two measurements do not guarantee a globally
injective state description. A continuous-time two-coordinate fit needs a
derivative of edot; a discrete-time fit or measurement delays can avoid explicit
second-derivative estimation. The measured edot also depends on the applied
input during control, so input timing and possible direct dependence on u must
be handled explicitly. Delays can include known input history. See
[measurement-based Koopman control with delays](https://arxiv.org/abs/1804.05291).

The current native regression still requires dynamical-state arrays to form
lifts (native/src/ControlRegression.cpp, residualData), and the present SINDy
runs used known-RHS derivatives. A strict measurement-only import, regression
and feedback path has not been implemented in this study. Such a path should
accept only measurement arrays and their timing, prohibit state/RHS access in
learning and control, and reserve any hidden simulator state for independent
validation. The prior r-H-DMD eigenvalue can remain a supplied comparison
constant; the entirely echo-derived arm must obtain its spectrum from echo.
