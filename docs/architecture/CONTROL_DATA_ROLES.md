# Data roles of the RPE Koopman controllers (traced 2026-11-15)

**Goal.** An echo-only controller that can be set up in a lab with nothing but acoustic sensors. Offline use of radius data is allowed for now, for target computation and possibly H-DMD, but every method must declare it. This page records, method by method, which information each controller uses offline and online. It was traced through the code (`Postprocess.cpp` SINDy, `CoordinateSindy.cpp`, `InputMap.cpp`, `ControlFactory.cpp`, `HistoryControl.cpp`, `EstimatorFactory.cpp`, `EstimatedControl.cpp`, `ControlTargetLift.cpp`) and through the manifests of the runs reported in `docs/reports/`.

## Categories

| Tag | Meaning |
|---|---|
| **OFF-r** | Offline radius or state record (r, and r' where used) |
| **OFF-RPE** | Offline use of the RPE equations: right-hand side, input field g(x), or the echo formula |
| **OFF-e** | Offline echo records only (the echo e(t), its history ξ, and applied inputs) |
| **ON-state** | Online access to the true state x = (r, r') |
| **ON-RPE** | Online evaluation of RPE expressions: g(x), e(x), ė(x), or Jacobians |
| **ON-e** | Online measurement of the echo (and ė) only |
| **TGT-offline** | Target lifted into the controller's coordinates offline (`control_target`) |

## Method table

| | ψ(r, r'), exact B_ψ | ψ(e, e'), exact B_ψ | ψ(ξ), regressed B_ψ | KLQE → KLQR, ψ(e, e') |
|---|---|---|---|---|
| Eigenvalues λ | **OFF-r** + **OFF-RPE** | OFF-r + OFF-RPE | **OFF-r** + OFF-RPE | OFF-r + OFF-RPE |
| Eigenfunctions ψ | OFF-r + OFF-RPE | OFF-r + OFF-RPE | **OFF-e** | OFF-r + OFF-RPE |
| B_ψ | OFF-RPE (analytic) | OFF-RPE (analytic) | **OFF-e** (`held_input`) | OFF-r + OFF-RPE (`exact_map`) |
| Target | TGT-offline | TGT-offline (OFF-RPE) | TGT-offline (OFF-RPE) | TGT-offline (OFF-RPE) |
| Online state | **ON-state** | **ON-state** | no | no |
| Online RPE | **ON-RPE** | **ON-RPE** | no | no |
| Online echo only | no | no | **ON-e** | **ON-e** (e and ė) |
| Steady drive error (mild / strong / hard) | 0.01 / 0.01 / 0.2% | 0.00 / 0.03 / 0.19% | 8.7 / 6.1 / 37% | 0.7 / 0.3 / 3.1% |

How each entry arises:

- **Eigenvalues (every method).** H-DMD runs on the radius (`methods_1`: `hdmd`, `observable: 0`). A state SINDy with `gamma_method: system_rhs` (the RPE right-hand side) follows, and `krom_extract` keeps the oscillatory pair. Every eigenfunction fit, ξ included, reads λ from that archive (`spectrum_archive`).
- **Eigenfunctions in (r, r') and (e, e').** `system_rhs` computes ψ̇ targets from the RPE right-hand side at recorded states. For (e, e'), the coordinates are also evaluated from those states with the RPE echo formula.
- **Eigenfunctions in ξ.** `runCoordinateSindy` uses finite differences of the echo history and refuses a physical right-hand side.
- **B_ψ.**
  - **exact:** ∇ψ·g(x) with g = (0, −Eu/r), evaluated online at the true state. For (e, e') the chain rule goes through the RPE Jacobian of (e, ė).
  - **held_input (ξ):** ψ from the frozen history of the echo along forced records, plus the applied input. Measurement only.
  - **exact_map (the KLQE runs):** a least-squares constant fitted to the analytic map ∇ψ·g sampled at the forced records' **states**. It is not an input–output regression.
- **Online paths.**
  - **Exact-B controllers:** `controller->value(time, state)`.
  - **ξ:** `historyControlLaw` → `measurement.evaluateReal(state)` (the simulated echo sensor).
  - **KLQE:** the sensors (echo, echo_dot) feed ψ(e, ė) (`measurement_transform: eigenfunctions`) and then the estimator.
- **Targets (since this change).** The `control_target` post-processing step evaluates the state target in the controller's coordinates offline, and the controller reads a `trajectory` target. The offline lifting uses the radius waveform plus the RPE echo formula, for (e, e') and ξ.

## ξ with an echo-only spectrum (2026-11-15)

The ξ spectrum no longer needs the radius. Running the pipeline with `methods=["hdmd"]` and `method_settings.hdmd.observables=["far_field_acoustic_echo"]` makes H-DMD read the echo. DynSysLab inserts the coordinate change that evaluates the echo on the training record's states, and there is no state SINDy and no RPE right-hand side. `artefact` then keeps the dominant pair.

- **Eigenvalues:** λ = −0.060276 ± 6.282856i from the echo, against −0.060058 ± 6.282857i from the radius. The linearization at equilibrium gives −0.06037 ± 6.28290i.
- **Tracking:** unchanged at 8.69% (mild), 6.14% (strong) and 36.98% (hard).

ξ's tags then become: eigenvalues **OFF-e**, eigenfunctions **OFF-e**, B_ψ **OFF-e**, online **ON-e**. The target is still TGT-offline (OFF-r + OFF-RPE), which is allowed. The training and forced echoes themselves come from the simulated sensor; see the caveats below.

Two fixes made this route work inside pipelines:

- **Observable channels on coordinate archives.** H-DMD's `observables` list now works when its input is a coordinate archive such as the identity `coords` stage. The channels are evaluated on the archive's stored `source_states` (`CoordinateChangeKinds.cpp`).
- **Singular stage references.** `${stage.archive}` names the method's result when the only other archive is that automatic coordinate change (`Pipeline.cpp`).

## KLQE → KLQR with echo-only training (2026-11-15)

Every offline step of the KLQE run now reads echo data only. The pipeline (`20kpa.json`, training r₀ = 0.3747) runs with:

- `methods=["hdmd"]`, `method_settings.hdmd.observables=["far_field_acoustic_echo"]`: λ from echo H-DMD.
- `artefact={"mode_selection":{"rule":"oscillatory_pairs","count":1,"max_real_part":0}}`: on this record the echo H-DMD ranks a growing pair first (+0.051 ± 6.18i), so the rule now skips growing modes. The kept pair is −0.0615 ± 5.941i; the radius-based K-ROM had −0.016 ± 5.943i.
- `echo_sindy.gamma_method=finite_difference`: ψ(e, e') from finite differences of the measured (e, ė), with no RPE right-hand side.
- `input_map` `fit: controlled_residual` on the three forced records: B_ψ by input–output regression (central differences of ψ along the forced (e, ė) records against the applied input), giving B_ψ = −0.213 + 0.772i. The model-based `exact_map` gives 0.094 + 0.842i on the same K-ROM.

| Steady drive error (offline targets, τ 100–300) | mild | strong | hard |
|---|---|---|---|
| Before: radius λ, `system_rhs` ψ, `exact_map` B_ψ | 0.7% | 0.3% | 3.1% |
| Echo λ, finite-difference ψ, `controlled_residual` B_ψ | **0.72%** | **0.40%** | **2.83%** |
| Echo λ, finite-difference ψ, `exact_map` B_ψ (reference) | 0.96% | 1.14% | 2.65% |

KLQE tags then become: eigenvalues **OFF-e**, eigenfunctions **OFF-e**, B_ψ **OFF-e**, online **ON-e** (e and ė). The target is TGT-offline (OFF-r + OFF-RPE), which is allowed.

**Still state dependent: the online ė.** The second KLQE sensor, `far_field_acoustic_echo_dot`, is the formula ∇e·f(x) evaluated at the simulator state, not a differentiated signal. The estimator option `measurement_derivatives: [0]` (with `measurement_derivative_method` `backward` or `backward2`) replaces it with a backward difference of the sampled echo. Results with that option:

- **Open loop** (±20 kPa pulses): the measured ψ changes by only 0.3–1.9%.
- **Closed loop at the tracking gain (R = 0.01):** it diverges within 0.1 τ. The true derivative of e along the forced flow contains the drive, (∂e/∂ṙ)(−Eu/r)·u. That term is a direct input feedthrough into ψ(e, ė), which the trained coordinates (u = 0) do not contain.
- **Gentler gains:**
  - R = 1 runs, at 7.0% on mild.
  - A constant feedthrough D, regressed from a forced record's measured ψ alone (second differences at pulse switches, R² = 0.51), makes mild work at R = 0.1 (1.2%). Strong and hard still diverge.
  - Drive caps (±300 kPa, ±1000 kPa) make it chatter.
- **Next steps:** a state-dependent (for example bilinear in ψ) feedthrough, or a sensor model that includes the drive term in e itself (caveat 1 below).

## Verification: offline targets reproduce the live-target results

All twelve best runs were repeated with `control_target` lifting the target offline (300 τ; scores on τ 100–300):

| Steady drive error (live → offline) | mild | strong | hard |
|---|---|---|---|
| ψ(r, r'), exact | 0.01 → 0.01% | 0.01 → 0.02% | 0.20 → 0.19% |
| ψ(e, e'), exact | 0.00 → 0.00% | 0.03 → 0.03% | 0.19 → 0.20% |
| ψ(ξ), regressed | 8.69 → 8.69% | 6.14 → 6.14% | 36.98 → 36.98% |
| KLQE → KLQR | 0.7 → 0.7% | 0.3 → 0.3% | 3.1 → 3.1% |

The `continuous_rhs` controllers (the exact-B_ψ rows) evaluate the target between samples, so they need a fine offline step: `target_step` = 1e-4. At 0.001 the linear interpolation of the target costs 0.36% (r, r') and 1.55% (e, e') on the mild target. Sampled controllers (ξ, KLQE) read the target exactly on their 0.001 τ clock.

## Remaining caveats for a lab setup

1. **The simulated echo sensor is a state formula.** e = S·r(2r'² + momentum), which uses the unforced acceleration. It omits the drive's direct term in R̈, and ė is a flow derivative along the unforced drift, not a differentiated signal. A hydrophone would record the drive term as well.
2. **Control and estimation run at 0.001 τ (3 ns), with no sensor noise.**
3. **Settings were selected by scoring against the exact inverse drive.** That covers R, training windows and training starting radius, and the drive is computed from the RPE model.
4. **Training and forced records are RPE simulations.** That is fine as a stand-in for experiments; the table above says which signals of them each method reads.

## Moving toward echo-only

These are the steps still needed per method:

- **ξ:**
  - λ from echo H-DMD: done (section above);
  - an echo target (TGT-offline from a measured or specified echo waveform), if radius-based targets are ever to be dropped.
- **KLQE → KLQR:**
  - λ from echo data, ψ(e, e') by finite differences, B_ψ by input–output regression: done (section above);
  - ė online from the echo signal rather than a state formula: the option exists, but it diverges at tracking gains (input feedthrough);
  - the same echo target.
- **All methods:** a sensor model that includes the drive term, numerically differentiated ė, noise, and a realistic sample rate.
