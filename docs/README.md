# DynSysLab documentation map

- [Mathematics toolkit](architecture/MATHEMATICS_TOOLKIT.md) - exact arithmetic, p-adic finite computations, graphs, Fourier calculus, resonances and distribution pairings.
- [Multiwell systems](MULTIWELL_SYSTEMS.md) - native triple-well chain and three-leaf Hamiltonian clover, observables and presets.

- Nonlinear pendulum derivative-closure proof - exact harmonic obstruction, native certificates and comparison with Duffing/Van der Pol.

- [Oscillator observables and spectral dictionaries](architecture/OSCILLATOR_OBSERVABLES_DICTIONARIES.md) - formulas, Colbrook sources, CLI/API and ownership.

- [KROM closure and physical-control audits](architecture/KROM_AUDIT.md) - differential defects, causal replay and learned versus physical local margins.

- [Lie algebra and Koopman representation analysis](architecture/LIE_ALGEBRA.md) - central research question, exact algorithms, public APIs and control-paper links.

- [Spectral KROM and native Koopman analysis](architecture/SPECTRAL_KROM.md) - methods, typed measures, observable references, Post-Processor and API.
- [Identification options](architecture/IDENTIFICATION.md) - pooled H-DMD and SINDy, domain SINDy, SINDy fit options, weak-form and weighted vector-field SINDy, delay taps, IC tables, spectrum comparison.
- [Colbrook native adaptation provenance](provenance/COLBROOK_NATIVE_ADAPTATION.md) - pinned references and redistribution notices.

- [Coordinate transforms and dynamics SINDy](architecture/COORDINATE_DYNAMICS_SINDY.md) - shared modules, vector-field fits and validation.
- [Lagrangian coherent structures and spectral submanifolds](architecture/LCS_AND_SSM.md) - FTLE fields, ridges and variational LCS; SSMs, backbone curves and data-driven SSMs (`ftle`, `ssm`).
- [Field models](architecture/FIELD_MODELS.md) - data-driven vector fields beyond one dictionary: compact radial, local tangent-chart, paired amplitude/phase, structured second-order; waveform fit and shooting calibration (`field_model`).
- [Trajectory kernels, orbit geometry and model libraries](architecture/TRAJECTORY_KERNELS.md) - cross-orbit injectivity, affine propagation, prehistory, clock audits, high-order differences, SINDy library families.
- [Receiver arrays](architecture/RECEIVER_ARRAYS.md) - retarded multi-receiver echoes, neighbour and linear state readouts, ring errors, sensor-design experiments.

- [Native pipeline boundaries and runtime APIs](architecture/PIPELINE_BOUNDARY_REFACTOR_2026-10-31.md)
- [Experiments](architecture/EXPERIMENTS.md) - one submitted deterministic study: stages, references, parameters, extends/overrides, provenance manifest.

Start with [Start here](START_HERE.md). The categories below distinguish
operating guides, architecture contracts, scientific workflow references, and
dated evidence. Case studies are reproducibility records, not default policy.

## Use DynSysLab

- [Overview](OVERVIEW.md) — the suite at a glance: modules, experiments and diagrams.
- [Install DynSysLab](INSTALL.md) — verified prebuilt Windows release, Codex setup, and first solve.
- [Data packs](DATA_PACKS.md) - version-pinned job files that regenerate published data and figures offline.
- [Archive management](ARCHIVE_MANAGEMENT.md) - archive root, compression, shared time grids, audit, metadata edits, storage locations and relocation.
- [Viewer reports and export](architecture/VIEWER_REPORTS_AND_EXPORT.md) - SVG/HTML export, multi-page reports, experiment reports, result documents and research figures.
- [Start here](START_HERE.md) — five-minute human and AI-agent entry point.
- [Native build and SDK](BUILD_NATIVE.md) — compiled installation, release
  gates, incremental builds, and extension development.
- [Optional Python client](PYTHON_CLIENT.md) — dependency-free CLI and workflow
  binding over the installed compiled engine.
- [Solver and Post-processor guide](SOLVER_POSTPROCESSOR.md) — native solve,
  campaign, and independent post-processing workflows.
- [Tutorial](TUTORIAL.md) — installation, systems, solves, campaigns,
  analysis, visualization, and reproducibility.
- [AI context and maintainer guide](AI_CONTEXT.md) — package contracts, data
  model, archive rules, and development checklist.
- [Agent harness](AGENT_HARNESS.md) — durable operating boundaries and
  protected-data rules.
  for scientific archives, Viewer products, screenshots, and report assets
  that remain local rather than being committed to Git.

## Suite and UI architecture

- [Viewer panel composer](architecture/VIEWER_PANEL_COMPOSER.md) — native plot assembly.

- [Native-core migration](architecture/NATIVE_CORE_MIGRATION.md) — canonical
  C++ ownership, optional Python bindings, staged parity gates, and progress.
- [Native compilation graph](architecture/NATIVE_COMPILE_GRAPH.md) — private
  Eigen kernels, runtime-sized matrix contracts, incremental-build rules, and
  measured before/after performance.
- [Native UI parity gate](architecture/NATIVE_UI_PARITY.md) — complete
  retained-workflow coverage across the compiled apps and release evidence.
- [Suite architecture remediation](architecture/DYNSYSLAB_SUITE_ARCHITECTURE_REMEDIATION_2026-08-31.md) — implemented refactors, before/after measurements, compatibility invariants, and final verification.
- [Script lifecycle policy](architecture/SCRIPT_LIFECYCLE_POLICY.md) — boundary
  between reusable suite features and retired one-off experiment code.
- [Solver/Post-processor C++ workbench](architecture/SOLVER_POSTPROCESSOR_CPP_WORKBENCH.md) — shared archives, versioned jobs, and panel design.
- [Viewer C++ workbench](architecture/VIEWER_CPP_WORKBENCH.md) — native Viewer information model and performance boundaries.
- [Viewer CLI](architecture/VIEWER_CLI.md) — headless figures from archives and arrays: views, windows, spectra, comparisons.
- [Integral action and PID](architecture/INTEGRAL_AND_PID_CONTROL.md) — general PID and LQI components with KLQR adapters.
- [TimeSeries and Fourier contract](architecture/TIME_SERIES_AND_FOURIER.md) —
  native sample ownership, zero-copy solver views, lazy spectral analysis, and
  archive/Viewer boundaries.
- [Bubble models and observables](architecture/BUBBLE_MODELS_AND_OBSERVABLES_2026-08-30.md) — RPE, Marmottant, and Reddy-Szeri model contracts.
- [Classical linear-systems toolkit](architecture/LINEAR_SYSTEMS_TOOLKIT.md) - gramians, Hankel singular values, balanced truncation/residualization, Kalman decomposition, pole placement, LQR/LQE/LQG, H2/Hinf norms, zeros, ERA/OKID; the `linear_systems` method and Viewer views.
- [Mode selection for echo-only KLQE/KLQR](architecture/MODE_SELECTION.md): choosing the eigenfunction subset by echo/state readouts, observability, controllability and B_psi quality
- [Joint H-DMD and linear diagnostics](architecture/JOINT_HDMD_AND_LINEAR_ANALYSIS.md) � named reconstructions, matrix/KROM Gramians, SVD directions, and estimator subset criteria.
- [Observable Hankel matrices](architecture/OBSERVABLE_HANKEL.md) covers shared windows, tapped channels, observable derivatives, POD coordinates, and DMD modes.
- [H-DMD scalar observables](architecture/HDMD_SCALAR_OBSERVABLES.md) — arbitrary scalar-observable Hankel inputs.
- [Native observable expressions and targets](architecture/NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md) — analytic derivatives, echo-only fits and scalar control targets.
- [Observable-coordinate SINDy and KLQR](architecture/OBSERVABLE_COORDINATE_SINDY_AND_KLQR.md) — change-of-variables and Jacobian ownership.
- [Koopman error evaluation](architecture/KOOPMAN_ERROR_EVALUATION_FRAMEWORK.md) — free/test trajectory and control-map error contracts.
- [Eigenfunction and KROM API](architecture/EIGENFUNCTION_MODEL_API.md) — algorithm-independent polynomials, explicit coordinates/spectra, separate evidence and compatibility migration.
- [KLQE first implementation](architecture/KLQE_FIRST_PASS_IMPLEMENTATION_2026-08-30.md) — implemented estimator interfaces and gates.
- [KLQE MATLAB trace and design](architecture/KLQE_KOOPMAN_ESTIMATION_FIRST_PASS_2026-08-30.md) — evidence and first design pass.
- [Mermaid architecture source](DYNSYSLAB_ARCHITECTURE.mmd) — compact visual reference;
  historical PDF exports are listed in the local asset pointers.

## Koopman and post-processing workflows

- [H-DMD to SINDy pipeline](HDMD_SINDY_PIPELINE.md) — matrix, spectrum,
  conversion, performance, and persistence contract.
- [Composable Koopman spectral pipeline](KOOPMAN_COMPOSABLE_SPECTRAL_PIPELINE.md) — independent eigenvalue/eigenfunction stages, gradients, observation maps, and control regressions.
- [KROM observation maps and training domains](KROM_OBSERVATION_AND_TRAINING_DOMAIN.md) — reconstruction validation and extrapolation gates.
- [Duffing Koopman architecture and control report](DUFFING_KOOPMAN_ARCHITECTURE_AND_CONTROL_REPORT_2026-08-27.md) — observables, reconstruction, exact/bilinear/constant input models, and retained SOP.

## Operations, security, and retained notes

- [Chat snapshot policy](CHAT_SNAPSHOT_POLICY.md)
- [Shared Antenna integration](ANTENNA_INTEGRATION.md)
- [Antenna security boundary](ANTENNA_SECURITY.md)

studies and are intentionally not listed as narrative documents.

- [Analytic trajectories and phase-space context](ANALYTIC_SOLUTIONS_AND_PHASE_CONTEXT.md)
- [Periodic spectral analysis and estimation](architecture/PERIODIC_SPECTRAL_ANALYSIS.md)
- [Repeatable periodic observable spectral comparisons](architecture/PERIODIC_SPECTRAL_COMPARISON.md)
- [Causal packet features and nonlinear periodic identification](architecture/PERIODIC_FEATURE_IDENTIFICATION.md)
- [Analytic tori and Hamiltonian perturbation tools](architecture/ANALYTIC_TORI_AND_PERTURBATIONS.md)
- [Orbit-family, damped-orbit, Melnikov and Gramian tools](architecture/ORBIT_FAMILY_TOOLS.md)
- [Koopman resolvent tools](architecture/KOOPMAN_RESOLVENT_TOOLS.md) - Laplace and inverse transforms, pole/residue models with Jordan blocks, GLA, pseudospectra, harmonic and Volterra response, orbit-shell boundary values, frequency-domain Gramians, Prony identification, continuation and contraction checks.
- [Spectral diagnostics tools](architecture/SPECTRAL_DIAGNOSTICS_TOOLS.md) - eigenvalue-lattice candidates, rational-kernel shell packets, and Reeb graphs of sampled scalar functions.
- [Certified bilinear control tools](architecture/CERTIFIED_BILINEAR_CONTROL_TOOLS.md) - dense LMI/SDP solver with an independent verifier, bilinear K-ROM fit, residual sector, and robust state-feedback design.
- [Toda lattice and Lax API](architecture/TODA_LATTICE.md)
- [Inverse-scattering building blocks](architecture/INVERSE_SCATTERING.md)
- [Mathematics toolkit](architecture/MATHEMATICS_TOOLKIT.md)
- [System-owned research metadata and bibliographies](SYSTEM_RESEARCH_METADATA.md)

- [Gravity-free double-pendulum analytic reference](GRAVITY_FREE_DOUBLE_PENDULUM.md)

- [Periodically forced phase context](PERIODIC_PHASE_CONTEXT.md) — forcing clocks, return-map shooting, section manifolds, cache identity and Viewer reconstruction.
- [Poincare sections and bifurcation sidecar](POINCARE_AND_BIFURCATION_SIDECAR.md) — affine first returns, projected stability, pseudo-arclength branches, and candidate event diagnostics.

- [Manifold surface continuation](MANIFOLD_SURFACES.md)
- [Hénon–Heiles analytic symmetry-axis orbits](HENON_HEILES_ANALYTIC.md)

- [Measured-output KLQE and KLQR composition](architecture/KOOPMAN_OUTPUT_FEEDBACK.md) - separate observer/controller runtimes, sampled solver wiring, applied-input feedback, and diagnostic archives.

Records of past studies (reports, case studies, notes) are kept in the private development repository, listed in `docs/STUDY_RECORDS.md` there. Every shipped example generates its own data.
