# DynSysLab documentation map

- [Native pipeline boundaries and runtime APIs](architecture/PIPELINE_BOUNDARY_REFACTOR_2026-10-31.md)

Start with [Start here](START_HERE.md). The categories below distinguish
operating guides, architecture contracts, scientific workflow references, and
dated evidence. Case studies are reproducibility records, not default policy.

## Use DynSysLab

- [Install DynSysLab](INSTALL.md) — verified prebuilt Windows release, Codex setup, and first solve.
- [Local asset pointers](LOCAL_ASSET_POINTERS.json) — hashes and relative locations of local-only report media.
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
- [Local artifact pointers](LOCAL_ARTIFACT_POINTERS.md) — locations and policy
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
- [Suite architecture audit](architecture/DYNSYSLAB_SUITE_ARCHITECTURE_AUDIT_2026-08-31.md) — measured policy conformance, dependency risks, usability findings, and roadmap.
- [Suite architecture remediation](architecture/DYNSYSLAB_SUITE_ARCHITECTURE_REMEDIATION_2026-08-31.md) — implemented refactors, before/after measurements, compatibility invariants, and final verification.
- [Script lifecycle policy](architecture/SCRIPT_LIFECYCLE_POLICY.md) — boundary
  between reusable suite features and retired one-off experiment code.
- [Solver/Post-processor C++ workbench](architecture/SOLVER_POSTPROCESSOR_CPP_WORKBENCH.md) — shared archives, versioned jobs, and panel design.
- [Viewer C++ workbench](architecture/VIEWER_CPP_WORKBENCH.md) — native Viewer information model and performance boundaries.
- [TimeSeries and Fourier contract](architecture/TIME_SERIES_AND_FOURIER.md) —
  native sample ownership, zero-copy solver views, lazy spectral analysis, and
  archive/Viewer boundaries.
- [Bubble models and observables](architecture/BUBBLE_MODELS_AND_OBSERVABLES_2026-08-30.md) — RPE, Marmottant, and Reddy-Szeri model contracts.
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

- [Input-map regression requirements](reports/BPSI_REGRESSION_REQUIREMENTS_AND_ECHO_HISTORY_EXTENSION_2026-10-05.md) — state, analytic-observable, and history-coordinate evidence requirements.
- [Native echo-history pipeline](reports/ECHO_HISTORY_NATIVE_PIPELINE_MIGRATION_2026-10-02.md) — native computation and migration evidence.
- [Sampled echo-only eigenfunctions](reports/ECHO_ONLY_SAMPLED_KOOPMAN_EIGENFUNCTION_RESULTS_2026-10-01.md) — conservative-system experiment and accuracy diagnostics.

- [High-resolution viscous ring-down](reports/RPE_VISCOUS_HIRES_RINGDOWN_2026-11-04.md) — full-sampling Viewer plots and native echo-history Hankel coordinates.
- [Echo-only causal Hankel–Koopman control roadmap](reports/ECHO_ONLY_HANKEL_KOOPMAN_CONTROL_ROADMAP_2026-10-01.md) — approved model-free research plan for causal echo-history coordinates, sparse sampled eigenfunctions, statistical input maps, discrete control, and laboratory validation.
- [H-DMD to SINDy pipeline](HDMD_SINDY_PIPELINE.md) — matrix, spectrum,
  conversion, performance, and persistence contract.
- [Composable Koopman spectral pipeline](KOOPMAN_COMPOSABLE_SPECTRAL_PIPELINE.md) — independent eigenvalue/eigenfunction stages, gradients, observation maps, and control regressions.
- [KROM observation maps and training domains](KROM_OBSERVATION_AND_TRAINING_DOMAIN.md) — reconstruction validation and extrapolation gates.
- [RPE clean-isolation KLQR SOP](RPE_CLEAN_ISOLATION_KLQR_SOP.md) — retained standard H-DMD/SINDy/KLQR workflow.
- [Duffing Koopman architecture and control report](DUFFING_KOOPMAN_ARCHITECTURE_AND_CONTROL_REPORT_2026-08-27.md) — observables, reconstruction, exact/bilinear/constant input models, and retained SOP.

## RPE case studies

- [Echo-history input regression and KLQR experiments](reports/ECHO_HISTORY_BPSI_EXPERIMENT_2026-10-06.html) — native measurement-only control-map fits, free-trained eigenfunctions, 104 control attempts, mild-target partial success and nonlinear-target failures, with Viewer plots.
- [Echo-jet KLQR end-to-end computation audit](case-studies/RPE_ECHO_JET_KLQR_END_TO_END_COMPUTATION_2026-09-16.md) — exact two-component target mapping, training lineage, eigenfunction fit, online Riccati computation, and 20/40 kPa evidence.
- [Echo-jet KLQR control-signal derivation](case-studies/RPE_ECHO_JET_CONTROL_SIGNAL_DERIVATION_2026-09-15.md) — consolidated derivation, immutable lineage, control-law audit, and 20/40 kPa comparison.
- [Echo-jet sine-driver rediscovery](case-studies/RPE_ECHO_JET_SINE_DRIVER_REDISCOVERY_LONG_TERM_2026-09-14.md) — full training-to-control pipeline and long-horizon 20 kPa result.
- [Nonlinear echo-jet sine stress test](case-studies/RPE_ECHO_JET_NONLINEAR_SINE_STRESS_2026-09-14.md) — target search, gain selection, and long-horizon 40 kPa result.
- [Echo-jet observability repair](case-studies/RPE_ECHO_JET_OBSERVABILITY_REPAIR_2026-09-14.md) — separates local observability from global polynomial closure, repairs the target-local fit, and verifies KLQR driver recovery.
- [Echo and autonomous echo-dot control](case-studies/RPE_ECHO_DERIVATIVE_EIGENFUNCTION_CONTROL_2026-09-11.md) — analytic derivative coordinates, optional input corrections, saved spectra and control comparison.
- [Echo-only eigenfunction control](case-studies/RPE_ECHO_ONLY_EIGENFUNCTION_CONTROL_2026-09-11.md) — prior-r and echo spectra, scalar targets and state behavior.
- [Echo-only KLQR driver rediscovery](case-studies/RPE_ECHO_ONLY_DRIVER_REDISCOVERY_2026-09-14.md) — a one-tone RPE echo target, matched state baseline, and a negative instantaneous-echo identifiability result.
- [Echo observability and long-horizon stability](case-studies/RPE_ECHO_OBSERVABILITY_AND_TAU40_STABILITY_2026-09-14.md) — nonlinear rank/conditioning diagnostics, laboratory radius sensing, and the four-times-longer multi-sine KLQR run.

- [Native Duffing/RPE KLQR parity replay](case-studies/NATIVE_KLQR_PARITY_REPLAY_2026-09-05.md)
- [Minimal-KROM sampled-KLQR stability](case-studies/RPE_MINIMAL_KROM_KLQR_STABILITY_CASE_STUDY_2026-08-26.md)
- [Clean-isolation initial transient](case-studies/RPE_CLEAN_ISOLATION_INITIAL_TRANSIENT_DIAGNOSIS_2026-08-27.md)
- [Optional production-control sequence](case-studies/RPE_OPTIONAL_PRODUCTION_CONTROL_SEQUENCE_2026-08-27.md)
- [Prior-control clean-pathway rerun](case-studies/RPE_PRIOR_CONTROL_CASES_CLEAN_RERUN_2026-08-27.md)
- [Offset-1.5 minimal-KROM failure](case-studies/RPE_OFFSET1P5_MINIMAL_KROM_FAILURE_CASE_STUDY_2026-08-27.md)
- [Designed control-identification data](case-studies/RPE_DESIGNED_CONTROL_IDENTIFICATION_DATA_2026-08-29.md)
- [Exact versus approximate control maps](case-studies/RPE_APPROXIMATE_BPSI_KLQR_COMPARISON_2026-08-29.md)
- [RPE versus Duffing approximate-control-map study](case-studies/RPE_VS_DUFFING_APPROXIMATE_BPSI_KLQR_CASE_STUDY_2026-08-30.md)

## Duffing case studies

- [Successful minimal exact-KLQR transfer/release](case-studies/DUFFING_MINIMAL_EXACT_KLQR_TRANSFER_RELEASE_2026-08-27.md)
- [Korda–Mezić control-matrix candidates](case-studies/KORDA_MEZIC_DUFFING_CONTROL_MATRIX_CANDIDATES_2026-08-27.md)
- [Optimized-spectrum/SINDy control regression](case-studies/DUFFING_OPTIMIZED_SINDY_CONTROL_REGRESSION_2026-08-28.md)
- [`duffrom2` MATLAB/Python parity](case-studies/DUFFING_DUFFROM2_MATLAB_PYTHON_PARITY_2026-08-28.md)
- [Minimal independent control-map restriction](case-studies/DUFFING_MINIMAL_INDEPENDENT_CONTROL_MAPS_2026-08-29.md)
- [Minimal fresh control-map refit](case-studies/DUFFING_MINIMAL_FRESH_CONTROL_REFIT_2026-08-29.md)
- [Two-mode fresh-refit KLQR](case-studies/DUFFING_TWO_MODE_FRESH_REFIT_KLQR_2026-08-29.md)
- [Bounded bilinear KROM search](case-studies/DUFFING_BILINEAR_KROM_SEARCH_2026-08-29.md)

## Operations, security, and retained notes

- [Chat snapshot policy](CHAT_SNAPSHOT_POLICY.md)
- [Shared Antenna integration](ANTENNA_INTEGRATION.md)
- [Antenna security boundary](ANTENNA_SECURITY.md)
- [RPE minimal-KROM reproduction notes](notes/RPE_MINIMAL_KROM_KLQR_STABILITY_REPRO_NOTES_2026-08-26.md)
- [RPE MATLAB parity checkpoint](notes/RPE_MATLAB_PARITY_CHECKPOINT_2026-08-26.md)

Machine-readable JSON diagnostics in `docs/notes/` accompany the relevant case
studies and are intentionally not listed as narrative documents.

- [Analytic trajectories and phase-space context](ANALYTIC_SOLUTIONS_AND_PHASE_CONTEXT.md)
- [System-owned research metadata and bibliographies](SYSTEM_RESEARCH_METADATA.md)

- [Gravity-free double-pendulum analytic reference](GRAVITY_FREE_DOUBLE_PENDULUM.md)

- [Periodically forced phase context](PERIODIC_PHASE_CONTEXT.md) — forcing clocks, return-map shooting, section manifolds, cache identity and Viewer reconstruction.
- [Poincare sections and bifurcation sidecar](POINCARE_AND_BIFURCATION_SIDECAR.md) — affine first returns, projected stability, pseudo-arclength branches, and candidate event diagnostics.

- [Manifold surface continuation](MANIFOLD_SURFACES.md)
- [Hénon–Heiles analytic symmetry-axis orbits](HENON_HEILES_ANALYTIC.md)
