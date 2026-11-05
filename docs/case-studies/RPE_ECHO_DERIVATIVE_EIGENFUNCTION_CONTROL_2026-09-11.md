# RPE control with echo and autonomous echo-dot eigenfunctions

Date: 2026-09-11. Native DynSysLab on main. MATLAB was not modified.

## Result

The native pipeline now fits functions of echo and autonomous echo-dot, with analytic second and spatial derivatives. Both requested spectrum comparisons completed to tau = 10. Input-dependent derivative corrections were disabled in this study, as requested. No H-DMD computation was repeated.

Adding echo-dot reduced the generator residual but did not restore the intended state motion. Radius still reached approximately 10.24 equilibrium radii. Echo tracking was less accurate than in the preceding scalar-echo fit; the new derivative coordinate also failed to track its mapped target.

![Control comparison](../../viewer_outputs/rpe-echo-jet2-control-comparison.png)

Errors use 5 <= tau <= 10. Echo NRMSE is RMSE divided by the RMS of the supplied echo target. Maximum radius and input statistics use 0 <= tau <= 10.

| Eigenfunction coordinates / spectrum | Echo NRMSE | Radius RMSE | Velocity RMSE | Maximum radius | RMS input | Peak absolute input |
|---|---:|---:|---:|---:|---:|---:|
| State / prior r spectrum, reference controller | 1.5713% | 0.001213 | 0.004121 | 1.2022 | 0.4555 | 1.4021 |
| Echo / saved echo spectrum, previous run | 0.022030% | 7.4423 | 0.8476 | 10.2784 | 1.2380 | 7.6304 |
| Echo / prior r spectrum, previous run | 0.022030% | 7.4424 | 0.8476 | 10.2785 | 1.2380 | 7.6324 |
| Echo + autonomous echo-dot / saved echo spectrum | 0.65407% | 7.3992 | 0.8512 | 10.2415 | 1.2199 | 6.2961 |
| Echo + autonomous echo-dot / prior r spectrum | 0.65806% | 7.3995 | 0.8512 | 10.2418 | 1.2200 | 6.3155 |

Radius r=R/R0 and velocity v=dr/dtau are dimensionless; input is p_control/p0. The autonomous echo-dot RMSE is 3.5482 Pa m for both new runs, compared with target RMS 1.8262 Pa m (about 194.3% normalized error). Echo RMSE is 0.0034980 and 0.0035193 Pa m, respectively. The two new trajectories nearly overlap. This is a negative result for these fitted functions and settings, not a general impossibility result for two-observable control.

## Native derivative and fit definitions

Write the existing RPE echo as e=h0(x), x=(r,v), and the zero-input vector field as f0. The new native observables are

\[
 e_1=L_{f_0}h_0,\qquad e_2=L_{f_0}^{2}h_0,
 \qquad \frac{\partial e_1}{\partial r},\quad
 \frac{\partial e_1}{\partial v}.
\]

The catalog aliases are echo_dot, echo_ddot, echo_dot_dr and echo_dot_dv. The echo_jet2 preset is exactly [e,e1]. They use the existing analytic expression graph and equations of motion; their spatial Jacobians are analytic. Derivatives use native time tau: physical-time first and second derivatives divide by T and T^2. The spatial derivatives use native r and v.

The eigenfunctions have the form

\[
 \phi(e,e_1)=\sum_{1\le a+b\le4}c_{ab}e^a e_1^b,
 \qquad
 \Gamma_j=\frac{\partial\theta_j}{\partial e}e_1+
          \frac{\partial\theta_j}{\partial e_1}e_2.
\]

The library has 14 mixed monomials and no constant. No radius or velocity was added as a polynomial coordinate. Both fits use the same full free training trajectory, raw coordinates (normalize=false), threshold 1e-6, one sparse-fit iteration and conjugate reuse. These retain the preceding scalar fit's numerical settings while enlarging the library. The archives retain training coordinates, their time derivatives, and their full spatial Jacobians; compact-KROM extraction preserves this evidence. The identity dy0/dtau=y1 holds to 1.26e-14 relative maximum discrepancy in the saved training arrays.

## Saved spectra and control target

The main training archive contains 25,001 points over tau in [0,25], step 0.001, with initial state [2.333386598601735,0]. Its SHA-256 remains e4529d55dd549e12c66dbb78c9521c6bc273aa8bb135637319568a2edfdc9eba. Native import uses the archived physical-coordinate conversion, R0=1e-5 m, T=3.0123833037972497e-6 s and kappa=4/3.

The saved echo H-DMD calculation used a 20,000 by 2,000 Hankel matrix. Its retained continuous-time pair, read from the existing compact KROM, is approximately

\[
 \lambda=-0.1049735212\pm4.9995022721i.
\]

The comparison retains the previously successful r-H-DMD pair

\[
 \lambda=-0.2106164458\pm4.8186345616i.
\]

Each new postprocess job contains only a SINDy stage, with spectrum_archive referencing an existing model. The exact stored eigenvalues and source hashes are in the machine-readable report.

The same canonical state reference is

\[
 r_{\rm ref}=1+0.09\sin(4\tau)+0.12\sin(3\tau+\pi),\qquad
 v_{\rm ref}=0.36\cos(4\tau)+0.36\cos(3\tau+\pi).
\]

The controller target is [h0(x_ref), e1(x_ref)]. Its first component is exactly the previous echo target. Its second component is the autonomous derivative observable evaluated on the reference state; it is not the time derivative of the arbitrary reference echo waveform. Saved targets were verified against these definitions (maximum discrepancies 0 and 8.9e-16).

Control starts at [1,0] and runs to tau=10 with output step 0.001, BDF tolerances 1e-9, maximum step 0.001 and continuous-RHS control evaluation. The existing KLQR framework uses Q=5000 I and R=1e-4, with the exact lifted input map and no added caps or feedforward. The same numerical Q in a different fitted basis is not an identical physical tracking cost.

## Optional input-dependent terms

For x'=f0(x)+g(x)u and b=L_g h0, the implemented correction is

\[
 e'=e_1+bu,
\]
\[
 e''=e_2+(L_g e_1+L_{f_0}b)u+(L_gb)u^2+b\dot u.
\]

The C++ Observable::inputFlowDerivative API and versioned observable query support these terms and their spatial Jacobians. Query include_input_terms=false ignores supplied input and rates; true includes them. A missing input rate means held-constant input. For corrected derivatives, query observable="echo" with flow_order=1 or 2. The named catalog derivatives always retain their autonomous definitions.

SINDy's include_input_terms setting controls whether archived input enters the RHS used to construct Gamma. Its default remains true for compatibility. The app exposes “Include known input in generator derivatives.” Both study jobs explicitly set it false; the derivative-query default is also false. The echo_jet2 coordinate definition stays autonomous in either case.

Turning off these corrections does not disable the physical actuator or the lifted input map. This experiment uses

\[
 B_\psi(x)=J_\phi(e,e_1)J_{[e,e_1]}(x)g(x).
\]

Consequently this is not a strictly measurement-only B_psi regression or controller: analytic derivatives, target mapping and the exact input map still use the known RPE model/state. The accepted autonomous approximation applies to the echo-derived coordinates. No external control framework or paper was used.

## Fit quality and interpretation

| Spectrum | Previous scalar generator residual RMS | New generator residual RMS | Previous normalized defect | New normalized defect |
|---|---:|---:|---:|---:|
| Saved echo | 407.9667 | 335.4333 | 0.96742 | 0.96347 |
| Prior r | 407.9308 | 335.3284 | 0.96850 | 0.96471 |

For the positive-imaginary mode, normalized defect is

\[
 \frac{\|\Gamma c-\lambda\Theta c\|_2}
 {\|\Gamma c\|_2+\|\lambda\Theta c\|_2}.
\]

The approximately 18% raw residual reduction is accompanied by only a small improvement in normalized defect, which remains large. Raw residuals also depend on eigenfunction scaling; full coefficients are retained in the JSON report. Training echo ranges from -4.52 to 379.83, while e1 spans roughly -121645 to 119749. Raw degree-four monomials therefore have widely different scales. This study did not tune scaling, degree or gains. It establishes that adding the autonomous derivative under these matched settings did not fix control; it does not establish whether a better conditioned fit would do so.

## Reproduction and validation

The [machine-readable report](../../viewer_outputs/rpe-echo-jet2-control-comparison.json) contains all model/solution paths, original source hashes, coefficients, fit diagnostics, settings and metrics. The [PDF figure](../../viewer_outputs/rpe-echo-jet2-control-comparison.pdf) is suitable for export. Original archives were preserved; new archive publication refuses replacement and adds a unique suffix when fits finish in the same second.

Portable native jobs are in examples/native/case-studies:

- rpe-echo-jet2-sindy-echo-spectrum.json
- rpe-echo-jet2-sindy-prior-r-spectrum.json
- rpe-echo-jet2-control-echo-spectrum.json
- rpe-echo-jet2-control-prior-r-spectrum.json

The SINDy jobs reuse archived spectra. The control examples point at this report's completed KROMs; use the newly emitted model path if reproducing a new fit. examples/native/rpe-echo-second-derivative-input-off.json and the corresponding -on.json demonstrate the query switch. See the [native observable guide](../architecture/NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md) for API details.

All 78 native CTest tests passed, including analytic derivative identities, independent derivative checks, input switching, two-coordinate fitting, target identity, archive round-trip and non-overwriting publication. The isolated package smoke gate and three UI parity tests passed. The updated native runtime was installed and its doctor reports ok; the installed input switch was verified in both positions. All four source/model hashes are unchanged.

The existing “RPE successful minimal-KROM control study” Mask remains active and now includes the two new KROMs and two controlled solutions (42 visible tree items, 1,622 hidden). Installed Viewer inspection verified the four new derivative observables, each with a time series followed by a heat map; only echo and echo-dot carry the supplied 10,001-point targets. Observable time series and heat maps have no initial-condition or equilibrium markers. The KROM workspace has a single vertically scrolling arrangement with the eigenvalue plot on the right. The [derivative target plot](../../viewer_outputs/rpe-echo-jet2-derivative-viewer.png) and [KROM workspace](../../viewer_outputs/rpe-echo-jet2-krom-viewer.png) were rendered and visually checked. Both report zero visible windows. All computation and rendering ran in the background.
