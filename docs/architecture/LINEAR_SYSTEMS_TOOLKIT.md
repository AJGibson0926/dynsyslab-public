# Classical linear-systems toolkit

Every linear model DynSysLab produces is an ordinary LTI realization once a
measurement and an input are chosen: a K-ROM in real split coordinates
(`A = generatorSplit`, `C` = named reconstructions, `B` = a named control map),
an H-DMD model, a DMDc/EDMDc model, or the linear part of a SINDy fit. This
toolkit applies standard linear-systems theory to all of them through one
`StateSpaceModel` value type. It complements, and does not replace, the
earlier [observability/controllability diagnostics](JOINT_HDMD_AND_LINEAR_ANALYSIS.md)
and the KLQE/KLQR runtime controllers.

## Where it lives

| Layer | Files | Role |
|---|---|---|
| Dense kernels (`dynsyslab_numerics`, `EigenNumerics*.cpp`) | `StateSpace.hpp`, `BalancedTruncation.hpp`, `SystemStructure.hpp`, `ControlDesign.hpp`, `SystemNorms.hpp` | algorithms on row-major `std::vector<double>` matrices |
| Classical frequency domain (`dynsyslab_numerics`, Eigen-free) | `TransferFunction.hpp`, `LoopMargins.hpp` | channel transfer functions and Bode data; loop robustness of LQR, Kalman, LQG and LTR loops |
| Adapters (`dynsyslab_linear_analysis`, Eigen-free) | `KoopmanStateSpace.hpp`, `analyzeStructure` (Kalman report plus PBH check) | K-ROM to `StateSpaceModel`; structure verdicts |
| Estimation and control sidecars (`dynsyslab_control`) | `ReducedOrderEstimation.hpp`, `ReducedOrderControl.hpp` | balanced-truncation LQE and LQR for K-ROMs |
| Workflow | `linear_systems` post-processing method (`LinearSystemsPostprocess.cpp`, `LinearSystemsClassical.cpp`) | KROM or A/B/C/D archive in, one immutable `.npz` out |
| Viewer | `LinearSystemsPayload.cpp`, `LinearSystemsLoopPayload.cpp` | `hankel_spectrum`, `gramian_spectrum`, `frequency_response`, `pole_zero`, `bode`, `nyquist`, `loop_bode`, `sensitivity` |

The library has no dependency on Koopman types except `KoopmanStateSpace.hpp`,
which converts a K-ROM into a `StateSpaceModel`.

`ControlDesign.hpp` is the single gain path. KLQR's `solveContinuousLQR` calls
`designLQR`, and KLQE's `solveContinuousLQE` calls `designLQE` (CARE by
duality). Only `ControlDesign.cpp` calls the numerics Riccati solver. Finite-horizon
KLQR still calls `numerics::solveFiniteHorizonLqr` directly, because the toolkit
has no finite-horizon design.

Eigen stays inside `dynsyslab_numerics` (see [NATIVE_COMPILE_GRAPH.md](NATIVE_COMPILE_GRAPH.md)).
No public header exposes Eigen types. The toolkit shares numerics' kernels rather
than duplicating them: one Lyapunov/Stein solver, one finite-horizon gramian, the
numerics zero-order hold, and the numerics CARE.

## Capabilities

**Model and simulation** (`StateSpace.hpp`)
- Continuous and discrete models (`TimeDomain`), validation, poles, stability.
- Discretization: zero-order hold, Tustin (exact inverse `continuousFromTustin`),
  forward Euler. ZOH simulation, impulse and step responses.
- Parallel connection (`G1 ± G2`, used for error systems) and similarity transforms.
- Stable/unstable additive split by a real matrix-sign spectral projection.
  Modes with `Re λ ≥ -margin` (`|μ| ≥ 1 - margin`) count as unstable, so marginal
  Koopman modes (equilibrium offsets, λ ≈ 0) are kept intact.

**Gramians** (`StateSpace.hpp`)
- Continuous Lyapunov `A X + X Aᵀ + Q = 0` and discrete Stein `A X Aᵀ − X + Q = 0`
  by Bartels–Stewart on the complex Schur form, with residual checks.
- Infinite-horizon controllability/observability gramians (stable models only).
- Finite-horizon gramians: Van Loan block exponential (continuous, any stability)
  and finite sums (discrete).
- Empirical gramians (Lall–Marsden–Glavaski): impulse responses and ±perturbed
  initial states averaged over perturbation scales; `snapshotGramian` builds
  `Σ w_k x_k x_kᵀ` from recorded snapshots, which is the data-driven route for
  nonlinear plants.

**Balancing and reduction** (`BalancedTruncation.hpp`)
- Hankel singular values by the square-root algorithm (SVD of `L_oᵀ L_c`).
- Balanced realization; directions with negligible HSV are dropped (numerically minimal).
- Balanced truncation and singular-perturbation residualization. Residualization
  keeps the DC gain exactly; truncation keeps `D`. Both satisfy
  `σ_{r+1} ≤ ‖G − G_r‖∞ ≤ 2 Σ_{i>r} σ_i`, and both are reported.
- Order selection by explicit order, by error tolerance, or numerically minimal.
- Unstable plants: the unstable/marginal part is kept exactly and only the stable
  part is balanced (`stablePart` and `reducedStablePart` give the exact error system).
- Every reduction returns `projection` `P` (`z_r = P x`) and `lift` `L` (`x ≈ L z_r`, `P L = I`).
- Data-driven realization: ERA (eigensystem realization from Markov parameters),
  OKID (Markov parameters from arbitrary input/output records through an ARX
  observer), and balanced POD (Rowley) from primal/adjoint snapshots.

**Structure** (`SystemStructure.hpp`)
- Controllability/observability staircases (orthonormal Krylov bases).
- Four-part Kalman decomposition with a structure residual, minimal realization.
- `analyzeStructure`: controllable/observable dimensions, uncontrollable and
  unobservable modes, stabilizable/detectable, and the PBH (`Pbh.hpp`) verdicts
  and smallest margins as an independent check.

**Design** (`ControlDesign.hpp`)
- Pole placement: Ackermann (single input, repeated poles allowed) and
  Kautsky–Nichols–Van Dooren robust eigenvector assignment (multi-input).
  Observer placement by duality. Every placement is verified.
- LQR: continuous CARE (the suite's Hamiltonian solver, reused) and discrete DARE
  (structure-preserving doubling).
- LQE: continuous and discrete steady-state Kalman filters (predictor gain and,
  for discrete models, filter gain).
- LQG compensator and closed loop; the closed-loop poles equal the union of the
  regulator and observer poles (separation principle, tested).

**Input/output** (`SystemNorms.hpp`)
- Frequency response (Schur-based, O(n²) per frequency) and singular values.
- H2 norm; H∞ norm by the Boyd–Balakrishnan/Bruinsma–Steinbuch Hamiltonian
  iteration (discrete models via the bilinear map).
- Transmission zeros of square systems (Rosenbrock pencil, QZ). Non-square
  systems are rejected explicitly.
- DC gain, modal damping table, SISO gain/phase margins.

**Transfer functions** (`TransferFunction.hpp`)
- `transferFunction(model, output, input)`: poles `eig(A)`, invariant zeros of the
  channel, gain `k` in `G = k Π(s−z)/Π(s−p)`, numerator and monic denominator
  polynomials, relative degree, and near-coincident pole/zero pairs
  (`cancellations`: uncontrollable or unobservable modes of the realization, which
  are reported rather than removed). `channel` extracts one SISO channel.
- `bode`: magnitude in dB and phase unwrapped along the grid. `unwrappedPhaseDeg`.

**Loop robustness** (`LoopMargins.hpp`, SISO loops, negative feedback)
- Loops: `stateFeedbackInputLoop` `K(sI−A)⁻¹B` (LQR at the input),
  `estimatorOutputLoop` `C(sI−A)⁻¹L` (Kalman filter at the output),
  `outputFeedbackInputLoop` `−C_K(s)G(s)` and `outputFeedbackOutputLoop`
  `−G(s)C_K(s)` for any compensator `u = C_K(s) y`, for example
  `LqgDesign::compensator`.
- `loopRobustness`: classical gain margins both ways (`gainMarginUpper`,
  `gainMarginLower`) and phase margin (angular distance from −1, 0–180°) with
  crossovers; the disk measures from `α = min|1 + L(jω)| = 1/‖S‖∞` (tolerated input
  gains `[1/(1+α), 1/(1−α)]`, disk phase margin `2 asin(α/2)`); `‖S‖∞`, `‖T‖∞`;
  nominal closed-loop stability; and `L(jω)` on the grid. The grid extends upward
  until a high-gain loop has rolled off, and includes `L(0)` and `L(∞)`.
- `loopTransferRecoveryCovariance`: Doyle–Stein recovery `W + q² B Bᵀ`. As `q²`
  grows, the LQG input loop tends to the LQR loop for minimum-phase plants; with a
  right-half-plane zero it does not (check `transferFunction(...).zeros`).
- Guarantees used as tests: an LQR loop and a Kalman loop (scalar `R`, `V`) keep
  `α ≥ 1` (gain margin `[1/2, ∞)`, phase margin ≥ 60°); an LQG loop has no
  guaranteed margin (Doyle 1978).

**Empirical (nonlinear) observability** (`EmpiricalObservability.hpp`, `dynsyslab_linear_analysis`)
- `empiricalObservability(system, observable, x0, t0, options)`: from `x0`, perturb each
  state by ±e_i, simulate the nonlinear ODE over a horizon with the recorded input
  (fixed-step RK4, input sampled at every stage), and integrate the output sensitivities,
  `W_ij = ∫ dy_i·dy_j dt`. For a smooth system and small e this is the finite-horizon
  observability gramian of the linearization along the trajectory; on a linear plant it
  equals `finiteObservabilityGramian` (tested). Also returns eigenvalues/vectors, the
  condition number, the instantaneous output Jacobian, and whether every perturbed run
  stayed valid.
- `scaledObservability(result, scale)`: `S W S` re-decomposed, so states with different
  units are comparable (scale = orbit RMS in the workflow).
- `estimationBoundStd(result, V)`: Cramér–Rao bound `sqrt(V (W⁻¹)_ii)` on estimating the
  state from the output window under white measurement noise of intensity V. It is local
  (linearization) and says nothing about global ambiguity (distinct states with similar
  output histories).

## Balanced-truncation sidecars for estimation and control

`designReducedLQE` / `designReducedKoopmanLQE` balance `(A, G, C)`, where `G` is
the process-noise input (identity unless a control map is named). They design
the LQE on the reduced model and lift the gain back as `L_full = L · L_r`. They
also report `eig(A − L_full C)` of the full model.

`designReducedLQR` / `designReducedKoopmanLQR` balance `(A, B, C_w)` with
performance output `C_w` (identity unless measurements are named). They design
the LQR on the reduced model with `Q_r = Lᵀ Q L` and lift the gain back as
`K_full = K_r P`. The acceptance check is `eig(A − B K_full)` on the full model.

The established full-order KLQE and KLQR are unchanged unless a stage asks
for a reduction.

### Stage option `balanced_truncation`

The KLQE stages (solve-job `estimator`, output feedback, offline `klqe`) and the
KLQR controller (`controller.settings`) accept one option:

| Value | Meaning |
|---|---|
| absent, `null`, `false` | full-order design (default) |
| `true` | numerically minimal order |
| integer `r` | reduced order `r`, counting retained unstable/marginal modes |
| `{"order", "error_tolerance", "method", "stability_margin"}` | as in `linear_systems`; `method` is `truncation` or `residualization` |

- **KLQE:** `reduceKoopmanLQEDesign` balances `(A, I, C_y)`, designs the LQE on
  the reduced model with `P W Pᵀ`, and writes the lifted gain, the full-KROM
  observer poles and the covariance `L P_r Lᵀ` into the design. The observer
  still runs at full order on the full lifted state, so archives keep their
  shapes. Metadata gains `balanced_truncation` (offline) or
  `estimator_balanced_truncation` (online), with the Hankel singular values and
  bounds. A lifted observer that is not asymptotically stable is rejected.
- **KLQR:** `KoopmanLQRController::setBalancedReduction` makes every gain solve
  balance `(A, B, I)`, solve on the reduced model with `Lᵀ Q L`, and apply
  `K_r P`. State-dependent `B_ψ` re-solves this whenever the input map changes,
  as the full KLQR does. A lifted closed loop that is not asymptotically stable
  throws. Archives gain `klqr_reduced_order`, `klqr_reduction_error_bound`,
  `klqr_gain` and `klqr_closed_loop_poles`, and evidence gains
  `controller_balanced_truncation`. KMPC rejects the option.

## `linear_systems` workflow method

Source: a K-ROM archive, or any archive with arrays `A` (and optionally `B`, `C`, `D`).

| Setting | Meaning |
|---|---|
| `measurements` / `measurement` | K-ROM reconstruction names stacked as `C` |
| `control_map` | K-ROM control map for `B` (constant maps directly; state-dependent maps need `operating_state` or `operating_lifted`) |
| `balancing_input` | `control_map` (default with a map), `process_noise` (`B = I`, estimation view), `none` |
| `selected_modes` | K-ROM mode subset (indices into the archived ordering) |
| `time_domain`, `sample_time` | `continuous` (default) or `discrete` (stored discrete input map, else ZOH) |
| `reduction_order`, `reduction_error_tolerance`, `reduction_method` | balanced reduction (`truncation` or `residualization`) |
| `stability_margin` | stable/unstable split margin (default `1e-8`) |
| `structure_tolerance` | Krylov staircase relative tolerance (default `1e-10`) |
| `frequency_points` | logarithmic grid size (default 400) |
| `lqr_state_weight`, `lqr_input_weight` | scalar Q, R; adds full and reduced (lifted) LQR |
| `lqe_process_weight`, `lqe_measurement_weight` | scalar W, V; adds full and reduced (lifted) LQE |
| both LQR and LQE weights | adds the coordinated output-feedback (LQG) design: LQE feeding LQR |
| `ltr_weights` | list of `q²`: loop-transfer-recovery LQG designs with `W + q² B Bᵀ`, one input loop each |
| `loop_frequency_points` | loop-robustness grid size (default 2000) |
| `readouts` | K-ROM reconstructions (e.g. `["state"]`) to analyze through the measurement: full observability gramian SVD and directions, each direction's readout signature, the fraction of each readout in the unobservable subspace, and with LQE weights the Kalman-predicted readout standard deviations and the covariance anisotropy inside the observable subspace |

The output schema is `dynsyslab.linear-systems/v1` (method `linear-systems`).
Arrays: `system_{A,B,C,D,poles}`, `modal_damping` (Re, Im, ωn, ζ),
`hankel_singular_values`, `hankel_error_bound[r] = 2 Σ_{i≥r} σ_i`,
`reduced_{A,B,C,D,poles}`, `reduction_projection`, `reduction_lift`,
`{controllability,observability}_gramian[_eigenvalues]` (stable part),
`frequencies`, `{full,reduced}_singular_values`, `{full,reduced}_phase` (SISO),
`transmission_zeros`, and optional `lqr_*`, `reduced_lqr_lifted_*`, `lqe_*`,
`reduced_lqe_lifted_*`, and with both weights `lqg_closed_loop_poles`,
`lqg_compensator_{A,B,C}`, `reduced_lqg_compensator_{A,B,C}` (designed on the
reduced model) and `reduced_lqg_full_closed_loop_poles` (that compensator closed
around the full model). Complex arrays are stored as `(count, 2)` real/imaginary.
Metadata records the structure report (with PBH verdicts and margins), the
reduction summary, the actual `reduction_error_hinf`, and the norms.

Classical frequency domain (always, when the model has inputs and outputs):
`transfer_function_poles`, `transfer_function_denominator`,
`transfer_function_zeros__y{i}_u{j}`, `transfer_function_numerator__y{i}_u{j}`, and
`bode_frequencies`, `bode_magnitude_db`, `bode_phase_deg` (frequency × channel);
metadata `transfer_functions` (gain, relative degree, cancellations per channel) and
`bode_channels`. Loop robustness (single input; LQR and/or LQE weights): for each
loop `lqr_input`, `kalman_output` (single output), `lqg_input`, `lqg_output`
(single output) and `ltr_{k}_input`, arrays `loop__<name>__frequencies`,
`loop__<name>__response` (complex) and `loop__<name>__return_difference`
(`|1 + L|`); metadata `loop_names` and `loop_robustness.<name>` with every margin
(null = infinite). The loops use the same split realization and scalar weights as
the runtime KLQR (`Q`, `R`) and KLQE (`W`, `V`) designs. Stable-part
norms are added when unstable modes make the full norms infinite.

### `empirical_observability` workflow method

Source: a solution archive of an ODE system. Settings: `observable` (default
`far_field_acoustic_echo`), `start`, `stop`, `points` (sample times spread over the window
on the output grid), `horizon`, `samples`, `substeps`, `relative_perturbation` (times the
state scale; default 1e-5), `state_scale` (default: orbit RMS over the window),
`noise_intensity` (default 1e-4). The run's applied input is replayed: the controller-clock
`applied_input` when present, otherwise `control_values`. Schema
`dynsyslab.empirical-observability/v1`; arrays `sample_times`, `sample_states`, `gramians`,
`eigenvalues`, `eigenvectors`, `condition_number`, `scaled_eigenvalues`,
`scaled_eigenvectors`, `scaled_condition_number`, `estimation_bound_std`,
`relative_estimation_bound`, `output_jacobian`, `valid`, `state_scale`; metadata medians,
the weakest scaled direction and invalid points. Plot with
`DynSysLabViewerCli time-series --time-key sample_times`.

## Numerical cautions

- **Split doubling.** A K-ROM that keeps both members of a conjugate pair as
  separate complex modes gives four real states for one oscillation. Every pole
  appears twice, and a single output observes at most two of the four states.
  The structure report and the HSVs expose this. `makeRealModalRealization`
  (`EstimationRealization.hpp`) is the remedy.
- **Near-repeated eigenvalues.** The Krylov staircase decides rank against an
  absolute tolerance. When eigenvalues nearly repeat (H-DMD conjugate copies
  that agree only to rounding), it can count numerically distinct but practically
  unobservable directions. The PBH margins (`pbh.min_observability_margin`) and
  the Hankel singular values are the more reliable evidence. `pbh_agrees: false`
  flags the case.
- **Infinite gramians and H2/H∞** need asymptotic stability. Marginal and growing
  Koopman modes are excluded from balancing and kept exactly.
- **H∞** returns infinity for unstable models. It is an L∞ norm only for stable ones.

## Data roles

The toolkit analyzes a model; it reads no trajectories (`trainingDataRole` is
`any`). Report each analyzed model with the categories in
[CONTROL_DATA_ROLES.md](CONTROL_DATA_ROLES.md), because those are inherited from
the K-ROM's eigenvalues, eigenfunctions and `B_ψ`.

## Tests

- `dynsyslab_linear_systems_tests`: analytic Lyapunov/Stein solutions,
  finite/empirical gramian convergence, first-order HSV `|bc|/(2a)`, diagonal
  balanced gramians, both balanced error bounds, residualization DC gain,
  unstable-mode retention, ERA/OKID recovery, balanced POD, Kalman dimensions,
  minimal realization transfer equality, PBH agreement, split-doubled pairs,
  Ackermann/observer/KNV placement, DARE golden ratio, CARE, discrete Kalman
  filter, LQG separation, H2/H∞ closed forms (first order and resonant peak
  `1/(2ζ√(1−ζ²))`), discrete H∞, zeros of `(s+2)/((s+1)(s+3))`, damping, phase
  margin of `1/(s(s+1))` (51.83°), gain margin of `2/(s+1)³` (4), Tustin
  warping and round trip, ZOH.
- `dynsyslab_reduced_order_tests`: K-ROM realization entries and DC gain, sampled
  transition, refusal of state-dependent maps without an operating point, and
  lifted reduced KLQR/KLQE, `reduceKoopmanLQEDesign`, the KLQR controller with
  `setBalancedReduction`, and full-order balancing reproducing the KLQE/KLQR gains.
- `dynsyslab_linear_systems_postprocess_tests`: workflow archive contents, the
  unstable part retained, the actual error inside the bounds, the LQG
  output-feedback arrays, and the `balanced_truncation` parser.
- `dynsyslab_pipeline_boundary_tests`: offline `klqe` with `balanced_truncation`
  records the reduction and still estimates the full lifted state.
- `dynsyslab_loop_robustness_tests`: LQR and Kalman `α = 1` (reaching it at
  roll-off, also for a high-gain loop), loop algebra `−C_K G` at input and output,
  Doyle's counterexample (LQR keeps `α ≥ 1`, LQG drops below 0.1), LTR converging to
  the LQR loop on a minimum-phase plant, both gain margins and the sensitivity
  peaks of `2/(s+1)³`, transfer function of `(s+2)/((s+1)(s+3))`, a reported
  cancellation, Bode magnitude and unwrapped phase.
- `dynsyslab_empirical_observability_tests`: the empirical gramian of a linear Duffing
  plant equals the finite-horizon gramian, descending positive eigenvalues, the output
  Jacobian, both Cramér–Rao bounds, scaling `S W S`, a nonlinear case.
- `dynsyslab_empirical_observability_postprocess_tests`: a forced solve analyzed by the
  workflow method (input replayed, exact gramian at every point, window checks).
- `dynsyslab_linear_systems_payload_tests`: all Viewer views, including the
  output-feedback pole figure, plant Bode, Nyquist (input and output loops,
  clipping, palette), loop Bode and sensitivity.

## Demonstration

See `examples/native/case-studies/rpe-echo-linear-systems-20261115/README.md`.
