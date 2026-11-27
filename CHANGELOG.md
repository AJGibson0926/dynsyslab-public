# Changelog

## 0.10.0 — 2026-11-26

- **Eigenpairs, spectrum and pseudospectrum are core types.** New `core::Eigenpair` (an
  `Eigenvalue`, its `Eigenfunction` and an optional `EigenpairResidual {value, definition, scope}`)
  and `core::Pseudospectrum` (grid points z with the minimum residual; `levelSet(ε)`); `core::Spectrum`
  now holds the point part with each pair's residual, a continuous part of `SpectralMeasure`
  components and an optional pseudospectrum, and `KoopmanModel::eigenpair(j)` joins them. The krom
  stores them once on its spectrum (`eigenpair_residual`, `pseudospectrum_*`,
  `continuous_spectrum_json`) and a spectral fit's invariance, wave packets, stochastic residuals and
  RKHS adjoints as `spectral_diagnostics_json`; the separate `spectral_data_json` block, the copy in
  the embedded operator contract, SINDy's `training_residual_rms` and `EigenfunctionFitDiagnostics`'
  residual RMS are gone. ResDMD and the `koopman_spectral` family attach one-step residuals and the
  `residual_grid` pseudospectrum (Colbrook & Townsend 2024), SINDy its generator residual of the
  published eigenfunction (Kaiser, Kutz & Brunton 2021), history EDMD its successor-normalized
  pair score, H-DMD (single and pooled) the one-step residual of its evaluated representation
  along the training trajectory (Arbabi & Mezić 2017), replacing `eigenfunction_koopman_residual`;
  DMD, EDMD and hand-given expansions leave residuals empty. `krom_extract` keeps the selected pairs'
  residuals, and `KoopmanModel::pairs()` iterates the eigenpairs in mode order without copying. The Viewer's K-ROM spectrum shows each pair's residual and a
  lattice pseudospectrum as a heat map. Kroms written before this keep loading, without residuals.

- **SINDy given eigenfunctions fits the generator null space again.** Since 3dd2122 SINDy anchored
  its fit to the eigenfunctions' sample values by default and took their scale. For H-DMD -> SINDy
  on a decaying trajectory those values cover only the last window, near rest. The anchored
  polynomial was then dense and no eigenfunction, with a residual of order its values, and it took
  H-DMD's arbitrary scale. That broke KLQR in rpe-viscous-klqr-multisine. The default is now `fit:
  "nullspace"` (Kaiser, Kutz & Brunton 2021); `anchored` stays opt-in. A mode that disagrees with
  its references beyond `agreement_tolerance` keeps the unit-RMS normalization, recorded in
  `agreement.calibration`.

- **Side-car model families, pulses and input-map audits are native.** `core/FieldModels.hpp`
  holds compact radial, local tangent-chart, paired amplitude/phase, second-order, pair-network and
  damped-interior fields as JSON definitions, the filtered-library waveform fit, the waveform
  objective and fixed-support shooting calibration. The pooled stage `field_model` (validation
  rollouts, `calibrate`, replay `controls`) uses them (docs/architecture/FIELD_MODELS.md). New solve
  input `pulse_branches` and initial conditions `archive_states`; stages `pulse_response` (paired
  causal sensitivities) and `control_map_audit` (path representation error, design rank and
  condition; core `ControlMapAudit.hpp`); `input_map fit: "symmetric_pulse"` reads pulse-branch
  pairs. The named RPE identification policy is the experiment
  examples/native/experiments/rpe-control-identification. Deleted: the side-car compact_field,
  local_field/local_rhs, paired, structured, waveform_fit, waveform_objective, calibration,
  sparse_stlsq, compat_rpe and audit modules; the echo_pulse programs; and
  scripts/generate_rpe_control_identification_data.py.
- **Remaining identification and eigenfunction rows are native.** `koopman_validation` gains
  `segments` and the 2-D occupied-cell `binned_defect`; `krom_refit` gains `orient` and
  `coefficient_summary`; `sindy_dynamics` gains `derivatives: central_high_order` and `monomials`
  (Laurent) libraries; `svd`/`pca`/`pod` gain `amplitude_queries` (state-interpolated POD amplitudes
  and reconstruction); `control_error_compare` adds `absolute_maximum`. New stages `delay_pairs` and
  `delay_consistency`, and core `extendLockedBasis`, `historyEndpointDecoder` and `TangentReadout`.
  The rpe-corrected-echo `echo_rate_fit` program is deleted; its job is a postprocess job.
- **Second pass: K-ROM, solver and queue ledgers closed.** Multi-input `klqr`/`kmpc` through
  `actuation: "additive_state"` with an n x m `B` (core `KoopmanModel::setAdditiveInputMatrix`). Solver:
  native `rayleigh-plesset-dimensional-baseline`, `history.controller_call_log`,
  `history.keep_partial_history_on_failure`, PDE `communication: "halo_exchange"` (optional
  `DYNSYSLAB_WITH_MPI` transport) and adaptive repartitioning. Queue: checkpoint-before-deadline with
  requeue and continuation, `history_policy`, `queue.reuse: "extend"`, `queue.staging` with digests,
  `campaign restructure` and `campaign review`; `checkpoint_profile: recoverable` checkpoints inside
  one integration. Retired: the MATLAB `.mat` `U_source` path and the Python Duffing bilinear
  campaign (`operations/duffing_bilinear`, `scripts/search_duffing_bilinear_krom.py`), replaced by
  `krom_control_search`; the Antenna link and the policy `audit` stay thin Python by decision.
- **Remaining identification and eigenfunction rows are native.** `koopman_validation` gains
  `segments` and the 2-D occupied-cell `binned_defect`; `krom_refit` gains `orient` and
  `coefficient_summary`; `sindy_dynamics` gains `derivatives: central_high_order` and `monomials`
  (Laurent) libraries; `svd`/`pca`/`pod` gain `amplitude_queries` (state-interpolated POD amplitudes
  and reconstruction); `control_error_compare` adds `absolute_maximum`. New stages `delay_pairs` and
  `delay_consistency`, and core `extendLockedBasis`, `historyEndpointDecoder` and `TangentReadout`.
  The rpe-corrected-echo `echo_rate_fit` program is deleted; its job is a postprocess job.
- **Shadow KLQE→KLQR study and KLQE tuning sweep are experiments.** `control_score` gains
  `estimation` scores: per-component estimate error, bias and relative error, innovation relative
  size and whiteness at time lags, and the force error of a shadow controller. The reference signal
  is optional with them, and `run_values` stores a per-run abscissa. Experiment sweeps read fields of
  object-valued points (`${sweep.setup.krom}`). `viewer_figure` gains `fit_y_to` (model-scale pages),
  `time_key` and `scale`. The Viewer adds the `solenoid_suspension` and `cycle_clock` research figures
  and animates phase covers over time. The rpe-echo-mode-selection shadow and tuning Python scripts
  and their generated jobs are removed.
- **Case-study programs become experiments.** The forced-EDMDc, forced-EDMDc validation, manual
  and retrained PDF eigenfunction and joint H-DMD/KLQE case studies are experiment jobs
  (`examples/native/case-studies/*/experiment.json`); their C++ programs are deleted. New for them:
  `dmdc`/`edmdc` `pool_sources` and `validation_sources`, `krom_refit` `gauge_reference`,
  `linear_analysis` `klqe-assessment` `role: validation`, `select` `stands_for`, and
  `state_readout` mode `affine`.
- **Viewer reports, export and result documents are native (phase H).** `DynSysLabViewerCli render`
  writes SVG (`--output plot.svg`) or a self-contained interactive HTML page of every page
  (`--output plot.html`: hover, zoom, pan, page tabs). `DynSysLabViewerCli report SPEC.json` renders
  multi-page `dynsyslab.report/v1` reports (paragraphs, tables, Viewer figures, images) to HTML or PDF,
  and `dynsyslab-cli experiment report MANIFEST` builds one for an experiment: stage graph, method
  tables, results and Viewer plots. The Viewer opens periodic spectral comparison results, phase-cover
  geometry and `dynsyslab.research-figure/v1` documents (finite-cycle spectra, cocycles, coordinate
  reconstruction, clock changes) directly. New: analytic query `phase_cover` with its clock-ratio
  check; pipeline methods `bifurcation_points` (post-transient maxima) and
  `observable_spectrum_comparison`; `figure_dataset` operations norm, difference, diff and psi;
  `viewer_figure` curves with their own source and clock; Koopman-mode participation, observer poles
  and the estimation-error norm in the observability view; spectra that open on their dominant
  band; per-system default views and cameras; automatic archive refresh; campaign type labels
  (docs/architecture/VIEWER_REPORTS_AND_EXPORT.md). Removed: the Python Viewer research and
  presentation modules, notebook packaging, `tools/experiment_report.py`,
  `postprocessing/formal_report.py` and `postprocessing/koopman/phase_covers.py`.
- **Experiment sweeps and validation-only selection; receiver arrays; trajectory kernels.**
  Experiments gain `sweep` stages (a block of stages repeated per point of its variables) and
  `select` stages (the best point by a held-out score, with a selection record), and
  `${id.archives}` lists. New methods: `orbit_geometry` (cross-orbit crossings, peak closure and
  state ambiguity of coordinate charts), `receiver_signals`, `state_readout` and `ring_errors`
  (retarded multi-receiver echoes, neighbour and linear delay readouts, per-ring errors; the
  receiver side-car is removed), `affine_dynamics` (A, c and exact propagation), `prehistory`,
  `clock_audit` and `run_invariant` (a conserved polynomial from run identity). Derivative method
  `central_high_order`, `sindy_dynamics` library families (bounded, rational, localized, gated,
  paired_gated, harmonic) and `dmdc`/`edmdc` windows from record metadata
  (docs/architecture/TRAJECTORY_KERNELS.md, RECEIVER_ARRAYS.md, EXPERIMENTS.md).
- **Identification methods are native (phase D).** `hdmd` and `sindy` pool several trajectories
  (`additional_sources`) with no window or difference across a boundary; SINDy trains on a
  state-space `domain` grid or `domain_points`, and takes `rescale`, `sparsification_policy`
  (`single-absolute`), `directional_derivative` (`finite-difference`), any `dictionary` (radial,
  trigonometric, expression terms) and `readout_transfer` (keep a parent krom's C behind a validation
  gate, else refit). `sindy_dynamics` adds `ridge`, `sample_weights`, `speed_scale` and weak-form
  `weak_durations`; `dynamics_validation` scores rollouts on their own clock with peak matching
  (`peak_prominence`) and a `kinematics` chain-rule check. New:
  coordinate kind `delay_taps` (physical-time lags, quintic delay jets), IC `mode: "file"` with a
  SHA-256 check, HAVOK `derivative_order`, per-step `svd_driver` and `analysis_coordinates`,
  eigenvalue summaries (`growth_rate`, `angular_frequency`, `cyclic_frequency`, `quality_factor`)
  and empirical spectral weights. `examples/native/duffing-hdmd-sindy-sop.json` is the Duffing
  H-DMD to SINDy procedure (docs/architecture/IDENTIFICATION.md).
- **K-ROM control and estimation tools are native.** Classical in-loop estimators (`lqe`,
  `luenberger`, `static_map`) and KLQE design reuse (`design_archive`); multi-input `lqr`/`mpc`/
  `pid`, SDRE with user factorisations, PCHIP, `expression` and lifted targets, and the KLQR/KMPC
  `target_domain_policy` with an archived domain report. Post-processing gains mode rules
  `state_observability`, `nearest_frequency` and `validation_score`; validation windows, binned
  defects and the psi history; a held-out quality gate for `observable_reconstruction`;
  symmetric-pulse B_psi estimation; p95 error fields; and the methods `control_error_compare`,
  `krom_rollout`, `pod_galerkin`, `krom_refit` (fixed-multiplier refit, gauge, pruning) and
  `reconstruction_transfer`. `control_score` adds tone/harmonic fits, settling time, low-pass
  drive error, estimate RMSE, innovation whiteness and run differences; `linear_analysis` adds
  `modal-gramian` and `klqe-assessment`; `krom-audit` adds `control_map_path`. The K-ROM control
  search takes an explicit schedule, closure-weighted ordering, map variants on path, arc-length
  and tube samples, and `resume_from`. `dynsyslab-cli krom-training-plan` designs fingerprinted
  training-domain plans from targets and probes (RPE preset included). Presets for the RPE
  clean-isolation KLQR and training screening are tracked under examples/native/presets. The
  Python K-ROM training planners, their screening and `rom.pod_galerkin` are removed
  (docs/architecture/KROM_TOOLS.md).
- **Phase B finished.**
  - The `coordinate_transform` kind `resample` resamples any post-processing source (PCHIP,
    linear, cubic or nearest) and publishes a record that names its source.
  - Mackey-Glass and lattice solves write checkpoint archives (`history.checkpoint_interval`).
  - The Python copies now covered natively are deleted: `postprocessing/resampling.py`,
    `postprocessing/koopman/resampling.py`, `solver/checkpoints.py`, and the controlled-training
    generators in `postprocessing/koopman/workflows.py`.
  - Their users now run native jobs through `operations/native_controlled_training.py`.
  - Python fixed-lag methods point non-uniform sources to the native resample kind.
- **Native Lagrangian coherent structures and spectral submanifolds.** The `ftle` method computes
  forward and backward FTLE fields from the flow-map Cauchy-Green tensor on a grid of initial
  conditions, their second-derivative ridges and variational LCS (strainlines), as a grid field
  Viewer draws. The `ssm` method computes spectral submanifolds and their reduced dynamics by the
  parameterization method about an equilibrium (normal-form or graph style, exact Taylor
  coefficients by power-series arithmetic), or fits them from trajectories (SSMLearn style), and
  publishes the reduced model as a `dynamical_model` with its backbone curve as a field
  (docs/architecture/LCS_AND_SSM.md). The placeholder `postprocessing/lcs` and `postprocessing/ssm`
  MATLAB-toolbox adapters are removed.
- **The native job queue schedules like a batch system.** Jobs take a `queue` block: `retry_limit`,
  `depends_on` with `afterok`/`afterany`/`afternotok`/`aftercorr`, `resources`, `wall_time_seconds`,
  `hold` (side tray), `reuse: "identical"` and spec claims (no two runs of one scientific spec at
  once). `execution.max_threads`/`max_memory_mib` and `--max-threads`/`--max-memory-mib` set
  execution limits. Campaign jobs submit as job arrays chunked by `systems_per_job` with
  `max_concurrent`; `queue serve` runs them with pool admission (`--max-running`, `--pool-threads`,
  `--pool-memory-mib`), a pause watchdog that quarantines unresponsive workers, and an optional
  `--command-template` for SSH or scheduler submission. `queue pause|resume|tray|release|remove|
  watchdog|resources` complete the queue. Campaigns gain linspace and `driver.*` axes, quality,
  execution and checkpoint profiles, cycle-based time (`transient_cycles`, `retained_cycles`,
  `samples_per_cycle`, Lyapunov cycle profiles), a campaign registry and `campaign
  list|status|validate|repair|profiles`. Progress files carry a heartbeat and remaining-time
  estimate; derived archives record stage timing; post-processing plans report Hankel and SINDy
  geometry and refuse jobs over their memory budget before allocating. `dynsyslab-cli benchmark`
  and `benchmark gate` replace the Python `performance` harness; `template solve|postprocess|
  estimate` emits starter jobs (docs/architecture/JOB_QUEUE.md).
- **Native solver: phase B of the port.** New ODE methods: RK23 (Bogacki-Shampine), DOP853,
  Radau IIA, CVODE Adams, and LSODA-style Adams/BDF switching. They run in one integrator driver
  that also handles:
  - arbitrary output times (`time.output_times`);
  - state guards (`solver.state_bounds`, `max_state_magnitude`) and a per-trajectory step budget;
  - `on_failure: "reject"`;
  - post-solve resampling with monotone PCHIP or nearest (`time.sample_times`, keeping `raw_*`
    arrays);
  - checkpoints, pause and resume (`history.checkpoint_interval`, `pause_file`, `resume_jobs` in
    `<manifest>.outcome.json`);
  - continuation into linked child archives (IC mode `continue_from`);
  - dimensional integration of RPE and Marmottant (`time.coordinate_system: "dimensional"`);
  - vector, additive-state actuation (`actuation: "additive_state"`).

  `ScalarControlPolicy` is now `ControlPolicy`. Also new:
  - designed inputs `tabulated`, `command_tape`, `telegraph_multisine` and `symmetric_pulse`
    (archived `input_design` with `pair_id` and `polarity`);
  - IC modes `stratified`, `target_tube` and `concatenate`, with `repeat` and `clip`;
  - the Marmottant `shell_law` sigmoid and Gompertz;
  - preset drivers, references and processing metadata;
  - method-of-steps delay rollout;
  - the adaptive method of lines for PDEs;
  - extension entry points for delay, field and lattice models.

  See `docs/architecture/NATIVE_SOLVER.md`.
- **A K-ROM is one class and one data object.** `core::KoopmanModel` is a spectrum and its
  eigenfunctions (sample values, origin and named representations, the first evaluated), optionally
  with input maps B and readouts C; it absorbs `core::Eigenfunctions`, which is gone. Every fitting
  method (H-DMD, SINDy, DMD, EDMD, history EDMD, the `koopman_spectral` family including ResDMD,
  boundary and hand-given eigenfunctions, polynomial fits, mode selection) publishes one data-object
  kind, `krom` (`objects/krom/`, schema `dynsyslab.krom/v1`), with its readouts where it fits them.
  Steps that add to a K-ROM (`observable_reconstruction`, `input_map`, `control_regression`,
  `history_control_regression`, `krom_extract`, SINDy on a K-ROM, polynomial fits) publish a child
  krom that references its parent (path, SHA-256, modes) and stores only what it adds;
  `loadKoopmanModel` reads any krom, following parents. A K-ROM without B says so when control needs
  one. The kinds `eigenfunctions`, `reconstructions` and `koopman_model` and the `reconstructions`
  step setting are gone; regenerate old archives (docs/architecture/DATA_OBJECTS.md).

- **Archive management is native.** Archives can be written with LZMA, deflate or stored ZIP
  entries (`--compression`, or `DYNSYSLAB_COMPRESSION`), and `dynsyslab-cli archive recompress`
  re-encodes existing ones losslessly. Archives that keep their time mesh in a shared grid
  (`time_grid_ref`) read natively, and `archive share-time-grids` writes such grids. `archive
  edit-metadata` writes a new archive with patched metadata and records the edit. `archive audit`
  gains `--campaign` scope and `--report DIR` (JSON and Markdown files); `archive campaigns` lists
  every campaign from its manifests and solution metadata. `dynsyslab-cli storage` registers
  storage locations, copies with SHA-256 verification and relocates a campaign. The CLI honours
  `DYNSYSLAB_ARCHIVE_DIR` (docs/ARCHIVE_MANAGEMENT.md). The Python `archiver` audit, deduplication
  and CLI, `recompress_archive`, `rewrite_archive_metadata`, `deduplicate_solution_time_grids` and
  `recompress_derived` are removed; `archiver.storage` is a thin client over the CLI.
- **Measured data and external control maps are usable again.** `measurement_origin: "external"`
  lets laboratory data with no recorded system drive a physical K-ROM in estimation, and
  `control_map_origin: "external"` attaches a control map that records no K-ROM identity. Both are
  recorded in the results. Control regressions again accept inputs held per sample or per interval
  whatever the record's alignment flag says.
- **Type aliases removed.** Code names the real types: `numerics::ContinuousLqrResult` (was
  `RegulatorDesign`), `numerics::FiniteHorizonLqrResult` (was `FiniteHorizonLQRSolution`),
  `HankelLeftBasis` (was `PooledHankelBasis`), the `Eigenfunction*` validation types (were
  `Koopman*` aliases), and `ObservabilityMethod`/`ObservabilityOptions` for controllability, the
  dual problem (were `Controllability*` aliases).
- **One in-memory time-series container.** The workflow's `SourceData` now holds a
  `core::TimeSeriesRecord` (series, inputs, channels, origin) instead of its own sample arrays,
  and `core::ControlledData` is removed: control regressions take a `TimeSeriesRecord` whose
  `inputs` are per sample or per interval. The `timeSeriesRecord`/`sourceView` converters are gone.
- **The Python application layer is retired; its capabilities are native.** Every capability of
  the legacy Python packages, the MATLAB parity scripts and the side-car experiments was ported into
  the C++ core and job pipeline first (see `docs/architecture/native-port/`, every row COVERED), and
  then the Python copies were removed: `postprocessing`, `rom`, `control`, `estimation`, `dynsyslib`,
  `solver`, `job_manager`, `campaigner`, `archiver`, `archives`, `spectral`, `infrastructure`,
  `analysis_data`, `performance`, the Python Viewer, the Tk `gui`, `app.py`, `hierarchy.py`,
  `time_series.py`, `provenance.py`, the Python job executor, `scripts/`, the computational and
  presentation scripts in `tools/`, and `tools/experiments/`. Python is now the installable thin
  client `dynsyslab_client`, the `python -m dynsyslab` front door (native CLI and workflow
  adapters, `ui`/`viewer` launch the built executables), the audit (which forbids NumPy/SciPy in
  client packages), the test runner, the MCP harness and the diagram and maintenance tools. The
  pinned Colbrook upstream manifest moved to `docs/provenance/colbrook/`.
- **One schema key.** Every versioned document, native or Python, names its contract under `schema`; archive and
  manifest metadata no longer write `format`, K-ROMs and DMDc/EDMDc archives no longer write
  `model_schema`, and solves no longer write `time_series_format`. A method that names its own
  contract (for example `dynsyslab.koopman-model/v3`) keeps it; other post-processing archives are
  `dynsyslab.derived/v1`. Archives written before this change are not read as K-ROMs.
- **Last compatibility switches retired.** `unverified_transform_reason`,
  `unverified_legacy_reason` on estimation sources, and `allow_unverified_coordinate_contract` on
  periodic spectral comparisons are gone: a K-ROM without a measurement-map identity, a physical
  K-ROM measured from a source with no physical system, and an unreadable coordinate contract are
  refused. Display coordinates no longer special-case Python-reference archives.
- **Legacy paths retired.** The unverified-legacy binding (`allowUnverifiedLegacy`,
  `unverified_legacy_reason` on K-ROM loads) is gone: an archive that records no physical system
  binds to the system its caller names, with status `requested-system`. An external control map
  that records no K-ROM identity is refused (`unverified_control_map_reason` is gone). The Viewer
  and the archive index no longer read `dynsyslab.hierarchy/v1` manifests.
- **Duplicate types merged.** One LQR result (`RegulatorDesign`, the numerics backend's
  `ContinuousLqrResult`) replaces `ContinuousLQRSolution`; `FiniteHorizonLQRSolution` is the
  backend's finite-horizon result. The real-split bilinear form of a K-ROM's bilinear map is a
  `BilinearKrom` (A, B0, Bᵢ) instead of `ContinuousBilinearSplit`. Pooled history bases are
  `HankelLeftBasis` (`PooledHankelBasis` is now its alias; `delays` → `rows`, `trainingWindows` →
  `columns`, `totalHankelEnergy` → `totalEnergy`, `trainingTrajectories` → `trajectories`). Archive
  keys are unchanged.
- **A K-ROM references its eigenfunctions.** A `koopman_model` (`dynsyslab.koopman-model/v3`) names
  an eigenfunctions object (path, SHA-256, modes, representation) and stores only its own readouts,
  input maps, spectral diagnostics and dataset references; nothing of the eigenfunctions or their
  spectrum is copied. Every step that makes a K-ROM publishes the eigenfunctions object first, so
  those steps now publish two archives. Older K-ROM archives are refused with a message to
  regenerate them, and the legacy loaders (`loadPolynomialKoopmanModel`,
  `savePolynomialKoopmanModel`, `loadKoopmanEigenfunctionModel`) are gone. An extracted or
  regressed K-ROM no longer copies its parent's training states, eigenfunction values or spectrum.
- **One spectrum encoding.** Every archive that stores eigenvalues writes the same keys
  (`continuous_eigenvalues`, `discrete_eigenvalues`, `sample_interval`, `logarithm_branches`,
  `eigenvalue_domain`; `SpectrumArchive.hpp`), and readers use them only; the duplicate
  `eigenvalues` array is gone. `DMDResult`, `HankelDMDResult`, `EDMDResult` and `DmdControlResult`
  hold a `Spectrum` (`result.spectrum`) in place of separate eigenvalue vectors. `Spectrum` moved to
  the foundation library. The Viewer reads a K-ROM's spectrum through its eigenfunctions reference
  and classifies K-ROMs by their declared kind only.
- **One fit result.** `KoopmanSINDyResult` is gone. SINDy, EDMD, DMDc, history EDMD and manual
  coefficients (`polynomialEigenfunctions`) all return `EigenfunctionFitResult`. Removed with it:
  `KoopmanModel::eigenfunctions()`, `evaluatePolynomialEigenfunction`,
  `gradientPolynomialEigenfunction`, `DmdControlResult::expansion()` and
  `EdmdControlResult::expansion()` (use `eigenfunctionSet()`), and the SINDy option
  `materializeTrainingValues`.
- **Data-object kinds are declared, never inferred.** Every writer names its object's kind
  (`publishObject`, `stampDataObject`, `withDataObject`); the name-to-kind table is gone, and every
  archive is stamped when it is written. An archive's `method` is its step's name (`hdmd`, `sindy`),
  the same as `data_object.produced_by.method`. Readers check kinds and contents instead of
  method names, and the Viewer catalogue classifies by kind. Fixed: the DMDc/EDMDc fits are
  `dynamical_model`s, `specrkhs` output is a `report`, history transforms are `coordinate_map`s, the
  mode-selection K-ROM is a `koopman_model`, the validated K-ROM is filed and registered, and H-DMD's
  saved Hankel matrix (one builder with `hankel_matrix`) can be reread with `hankel_archive`.
  `DataObject` and `TimeSeriesRecord` moved to the foundation library. Old K-ROM layouts are no
  longer recognized by the catalogue.
- **SINDy adds representations to existing eigenfunctions.** SINDy accepts eigenvalues, a spectrum,
  eigenfunctions, a K-ROM, or training data alone. Given eigenfunctions (H-DMD → SINDy), it fits a
  polynomial anchored to their values, calibrates its scale, records the agreement per mode, and
  publishes the same eigenfunctions with the polynomial added: one object, three representations
  (`polynomial`, the default; `interpolation`; `time-delay`). It extends H-DMD's object by
  reference, copying nothing. Given training data alone, it also finds the eigenvalues (Kaiser,
  Kutz & Brunton 2021, Algorithm 1). Representations now have unique names. `krom_extract` takes
  `eigenfunction_representation`, and a time-delay K-ROM is refused. Earlier SINDy archives still
  read. H-DMD → SINDy experiments now produce extensions, so their SINDy stage results change.
- **One eigenfunction, all its representations.** An `Observable` can hold several named
  representations (`withRepresentation`, `as`, `representationNames`); the first is its default.
  `Eigenfunctions` (formerly `EigenfunctionSet`, the in-memory eigenfunctions object) holds named
  `FunctionRepresentation`s directly, and `eigenfunctions.eigenfunction(j)` returns mode j with every
  representation. The wrappers `EigenfunctionRepresentation`, `SharedRepresentation` and
  `WindowMapRepresentation` are gone: H-DMD's window map is an ordinary expansion on the window
  domain (`windowMapFromModes`). `readEigenfunctionSet`/`writeEigenfunctionSet` are now
  `readEigenfunctions`/`writeEigenfunctions`. Archives are unchanged.
- **Stage reuse checks the pipeline.** An experiment reuses a stage only when the same pipeline
  computed it: the manifest records each stage's `pipeline`, a SHA-256 over the running DynSysLab
  binaries. After a rebuild, rerunning an experiment recomputes it and shows the change. Before
  this, an identical job reused results from older code.
- **Staged studies are experiments.** The pipeline is DynSysLab's modules and data objects; a job
  that chains stages into a study is an *experiment* that uses it. The job kind `"pipeline"` is now
  `"experiment"` (the old kind is refused with a message to rename it). Runs record
  `archives/campaigns/experiments/<label>__<utc>.experiment.json` (`dynsyslab.experiment/v1`), and
  stage reuse reads only those records, so the first run of an existing experiment recomputes its
  stages. The examples moved to `examples/native/experiments/`, the guide to
  `docs/architecture/EXPERIMENTS.md` and the report tool to `tools/experiment_report.py`.
- **Hankel windows are causal by default.** `window_alignment` now defaults to `"end"` for the
  legacy keys (`matrix_size`, `sampling_window`) as well as the `window` object: a window is
  labelled by its newest sample, so an eigenfunction's value at time t uses only history. Jobs that
  relied on the old legacy default, `"start"`, change their time labels; give `"window_alignment":
  "start"` to keep them.
- The Viewer's phase-plane view names a window-map eigenfunction "time-delay". Several runs of
  one observable overlay on one phase plane (`--compare-source`). Each keeps its domain colouring,
  brightness uses one modulus scale for all runs, and the legend names no runs. Training runs are
  drawn bolder than test runs, from the trajectory's recorded data role.
- Every colour plot has a vertical colour bar right of the plot (time, scalar values, or phase from
  −π to π). Phase portraits mark each initial condition with an open circle labelled with its
  coordinates, and verified equilibria with filled circles; crowded labels move to a clear row.

## 0.9.0 — 2026-11-25

A maintenance release that brings the public tree up to date before a larger release. It contains
no data packs. **Archive formats changed**: eigenfunction archives written by 0.8 are refused with a
message to regenerate them.

**Archives hold data objects.** Each archive a processing step writes is one data object (see
`docs/architecture/DATA_OBJECTS.md`). It carries a `data_object` envelope: its kind, and the
method, settings and inputs (with SHA-256) that produced it. Archives are filed by kind under
`objects/<kind>/`, not by method under `derived/`.

- H-DMD publishes two objects. The eigenfunctions object holds the spectrum and stand-alone
  eigenfunctions; each eigenfunction has its own window-map row and the shared interpolation. The
  reconstructions object holds readouts written over those eigenfunctions. H-DMD copies none of its
  input.
- SINDy, DMD, EDMD and `eigenfunction_expansion` publish eigenfunctions objects. SINDy copies no
  training data.
- `krom_extract` and `control_regression` assemble a K-ROM from eigenfunctions and the
  reconstructions named by `reconstructions`.
- Pipelines can reference `${stage.object:KIND}`. `${stage.archive}` names a method's primary
  object.

**One interpolation framework for every observable.** `core::interpolateSamples` builds an
interpolation from the samples' full states, found through the series' origin and mapped through
the domain's state map (the full state by default). Each observable keeps its own interpolant. The
domain belongs to the eigenfunction; the K-ROM loader never re-domains.

**Observables along time series.** `Observable::evaluate(TimeSeries)` returns a time series, and
eigenfunctions inherit it. A domain can be a series window, so a window-map eigenfunction evaluates
new windows without the method that made it. `coordinate_transform` with `kind: observable` accepts
an `eigenfunction`.

**Coordinate changes store a record, not copies.** Reading one replays it from its digest-checked
source. `dynsyslab-cli materialize` writes a self-contained copy.

**Viewer.**
- The "Observable on the phase plane" view shows an observable's values along a trajectory,
  domain-coloured.
- Viewer shows composed objects, such as eigenfunctions with their reconstructions.
- K-ROM training states are found through references instead of copies.
- Also new: domain fields with domain colouring, 3-D slices, and an animation player (toolbar and
  `DynSysLabViewerCli animate`).

**Also new.**
- The system eigenfunction library.
- Domain grids and coordinate maps: pull-back, push-forward, round trip and fold set.
- Time-series crop, resample, spectrum and window histories.
- Redesigned architecture diagrams.
- I/O-delay control models.

## 0.8.1 — 2026-11-23

First published build of the 0.8 line; the 0.8.0 tag failed its release gate and was never released.

- Data packs keep exact bytes on Windows checkouts, so their manifests validate.
- Public documentation no longer links to private-only reports.

## 0.8.0 — 2026-11-23

**Data packs: version-pinned job files that regenerate published results offline.** See
`docs/DATA_PACKS.md`.

- Job kinds: `job_group` (child job files, dependency order, symbolic archive links, resume ledger,
  reference-counted cleanup, free-space floor) and `viewer_figure` (headless Viewer CLI render).
- Post-processing: `krom_recenter`, `krom_readout`, `reference_signal`, `control_score`, `figure_dataset`.
- `dynsyslab-cli datapack inspect|validate|plan|run|verify|seal|export|trace|rehash`.
- Archives record the resolved job, version, source commit and build identity (`run_provenance`).
- **Changes results:** the Hankel basis start matrix now uses a compiler-independent normal stream
  (`mt19937_64+as241-inverse-cdf/v1`; Wichura, Applied Statistics 37(3), 1988), so history-coordinate
  (xi) results can differ from earlier releases.

**Koopman output feedback (KLQG) away from the RPE.** Upgrades made during the 2026-11-20 KLQG study
of classical systems, so a K-ROM can be identified around an operating point and used by KLQE -> KLQR:

- `dmdc` / `edmdc`: `time_window` and `state_window` (`[[low, high], ...]` per snapshot component)
  keep only snapshot pairs whose two snapshots lie inside both windows, so a fit can be restricted to
  the neighbourhood of a target (for example an unstable equilibrium). `available_pairs` is recorded
  beside `pairs`. A window that keeps fewer than two pairs is an error.
- `edmdc`: `export_krom: true` also publishes the eigenfunction bank as a standard
  `koopman-model/v2` K-ROM (`derived/koopman-edmdc-krom`), with the fitted (windowed) training states,
  their eigenfunction values and a per-mode `mode_consistency` (relative one-step residual including
  the modal input term, in the spirit of ResDMD: Colbrook & Townsend, Comm. Pure Appl. Math. 77, 2024). `krom_extract`, `input_map`, KLQE and KLQR then use it directly.
- `krom_extract`: `mode_selection: {"rule": "principal", "target_state": [...]}` selects the principal
  eigenfunctions of an operating point: modes whose gradients at the target are independent, the first
  pick being the most consistent growing mode. Products and harmonics of selected eigenfunctions (zero
  gradient at the target) are never selected; at an unstable target they make the lifted (A, B)
  unstabilizable for a constant input map. `max_residual` (default 0.3) sets the consistency filter.
  Principal eigenfunctions follow Mauroy & Mezić (IEEE TAC 61(11), 2016) and Kvalheim & Revzen
  (Physica D 425, 2021).
- `input_map`: `state_window` restricts the B_psi regression to forced samples inside a state box;
  each run of consecutive inside samples is its own record.

## 0.7.0 — 2026-11-19

**DynSysLab Hub: a plain chat session greets you.**

- `dynsyslab-cli hub`: the DynSysLab greeting. On the first login after an install (tracked per
  version in `%LOCALAPPDATA%\DynSysLab\hub.json`) it includes the welcome and an overview of the
  suite and its architecture; every time it ends with "What would you like to do?" and numbered
  options. `--format json`; `--peek` does not record a login.
- `dynsyslab-cli mcp`: a Model Context Protocol server (stdio) for Claude Code, Codex and other
  agents. Its instructions open the Hub at the start of a DynSysLab conversation; tools: `hub`,
  `quickstart`, `list_systems`, `doctor`, `validate_job`, `plan_job`, `run_job`, `fourier`,
  `render_figure` (returns the PNG); the `hub` prompt is the same greeting as a slash command.
- `dynsyslab-cli agent-setup [--remove]`: registers the MCP server with Claude Code
  (`claude mcp add --scope user`, including the copy bundled with the Claude desktop app) and Codex
  (`codex mcp add`, or `~/.codex/config.toml`), adds a Claude Code SessionStart hook
  (`dynsyslab-cli hub --hook`: greets on the first login after an install and in the DynSysLab
  workspace, silent elsewhere), and writes `AGENTS.md`, `CLAUDE.md` and `.dynsyslab-workspace`
  into the workspace (`Documents\DynSysLab`). `--remove` undoes all of it and keeps entries it did
  not write.
- `install.ps1` runs `agent-setup` after installing (`-NoAgentSetup` skips it), then `quickstart`.
  `release-install.json` adds `workspace` and `hub_command`.

**Data policy: DynSysLab ships no solver output.** The release and the public repository contain
settings, parameters, schemas and instructions only; every shipped example generates its own
data and runs with an empty archive folder; numerical results live only in `archives/`.

- Converted to self-generating pipelines (same parameter values, data solved by a first stage):
  `examples/native/rpe-echo-history-koopman.json` (training and validation ring-downs),
  `pipelines/rpe-klqe/echo-to-state.json` (training and hold-out),
  `pipelines/rpe-coordinate-dynamics/study.json` (training and validation initial conditions),
  `pipelines/rpe-echo-jet-klqr-rediscovery/20kpa.json` (the 20 kPa, 100 kHz target orbit),
  `duffing-hdmd-sindy-postprocess.json` and `spectral-koopman-job.json` (their sources).
- Private development repository only (left out of the public export and the release with
  `export-ignore` in `.gitattributes`, which packaging reads too): `examples/native/case-studies/`
  (records of past studies, including `tuning-results.jsonl` and the fitted drive in
  `forced/hard.json`), `examples/duffing-minimal-exact-klqr-transfer-release/`, the mode-selection
  pools `pipelines/rpe-echo-mode-selection/{echo-pool,hdmd-pool,range-pool}.json` (they need
  derived K-ROMs and closed-loop runs), `docs/reports/`, `docs/notes/` (including the controller
  diagnostics JSON), `docs/case-studies/`, the local asset pointers, the RPE isolation SOP, the
  August architecture audit, `tools/analyze_rpe_echo_jet_nonlinear_sine.py`, two `scripts/` that
  read archived runs and the `rpe-successful-control` Viewer focus profile. The docs index lists
  them in the private `docs/STUDY_RECORDS.md`.
- Removed `docs/case-studies/assets/viewer-insets-2026-11-15/*.png`; `.gitignore` covers
  `/docs/case-studies/**/assets/` and `/docs/**/assets/*.png`.
- Release gate `tools/Test-ReleaseContents.ps1` (repository and package, in both workflows) fails on
  `archives/(solutions|derived|campaigns|hierarchies)/` paths, `.npz .npy .mat .h5 .parquet`,
  PNG/JPG outside the application icons, `.jsonl`, and text naming an archived run.
- The release install test validates and plans every shipped example on a fresh install with an
  empty archive folder and runs every one that finishes within 5 minutes
  (`examples/native/long-running.json` lists the longer RPE campaigns).

**CLI output paths.**

- CLI output paths: `run --progress progress.jsonl` (a bare file name) no longer fails with
  `create_directories: ... ""`. `--progress` and `--archive-root` resolve against the current
  directory when parsed; missing parent folders are created; output paths are checked after
  job validation and before any stage runs. Path failures are `dynsyslab.cli-error/v1` JSON
  naming the option and the resolved path (every CLI error is now that JSON). Progress and
  manifest writers, `coordinate-inverse` outputs and `DynSysLabViewerCli --output` share one
  helper (`core::prepareOutputFile`).

## 0.6.2 — 2026-11-19

- `dynsyslab-cli quickstart`: the first run as one command. It prints the welcome and
  overview, runs the bundled examples (Duffing double well, Lorenz 1963), reports each run
  (samples, solver steps, dominant frequency, final state, archive) and Viewer CLI figures as
  Markdown (`--format json` for tools), and ends with things to try next. No window opens.
- Installer: runs `quickstart` at the end (`-NoQuickstart` skips it) and no longer addresses
  agents in its output; `release-install.json` records `version`, `cli`, `welcome_document`,
  `archive_root` and `quickstart_command`.
- Release notes start with an "Install and run" procedure, generated from
  `.github/release-notes-template.md` by the release workflow.
- Docs: `AGENT_WELCOME.md` is plain user-facing content; INSTALL describes install and run as
  three commands. The release install test checks the first-run output, that no desktop
  application starts, repeated runs and a long working directory.

## 0.6.1 — 2026-11-19

- CLI: without `--archive-root`, results go to the archive folder the installer configured
  (`runtime.json`), or `<checkout>/archives` from a source build, never the working
  directory; `doctor` reports the folder and its source. `<command> --help` prints usage.
- Archives: runs in the same second get distinct names, manifests and batch identities
  (`...Z-2.npz`) instead of replacing one another; archive and manifest paths past the
  Windows 260-character limit work; a failed publication leaves no result without its
  manifest.
- Installer: never starts the applications; it ends with DynSysLab's welcome and overview,
  which agents present in the chat after installing.
- Source: personal absolute paths removed from tracked case-study files.

## 0.6.0 — 2026-11-19

- Classical linear-systems toolkit: gramians, balanced truncation, structure, LQG,
  transfer functions and Bode data, loop margins (classical, disk, sensitivity) for LQR,
  Kalman, LQG and loop-transfer-recovery loops; Viewer Bode, Nyquist, loop Bode and
  sensitivity views.
- K-ROM observability measures: empirical (nonlinear) observability gramians along
  trajectories, readout observability, obsv in discrete modal coordinates, and
  per-eigenfunction observability with Cramér–Rao and KLQE bounds.
- Estimation and control: Estimator module with echo-only KLQE, `estimated_state`
  feedback for any controller, shadow KLQE → KLQR runs (`applied_controller`), PID and
  LQI components with KLQR adapters, KLQE process covariance shaped from K-ROM residuals
  or the input direction, balanced-truncation KLQE/KLQR.
- Modelling: `mode_selection` for eigenfunction subsets, smooth polynomial fits of
  scattered eigenfunctions, offline control-target lifting, noise-robust numerical
  differentiation, deterministic pipeline jobs with plug-and-play inputs.
- Solver: `on_failure: "truncate"` archives the computed part of a failed run.
- Viewer and its CLI: archive figures without scripts, comparisons, time windows,
  spectra, stem plots with open-circle tips, integer index ticks, 3-D stem plots, mode weights.
- Agent front end: a first-run welcome and overview (`docs/AGENT_WELCOME.md`) for
  Claude Code, Codex and other agents; release now also attaches the Python client wheel.

## 0.5.2 — 2026-11-04

- Restored Python 3.10 compatibility in the optional reference HTML renderer.
- Made archive-discovery tests compare resolved Windows paths, including hosted
  runner short-path aliases.
- Replaced a machine-specific report link with its declared local-asset pointer,
  preserving the source-only public distribution.

## 0.5.1 — 2026-11-04

- Fixed the extension test build dependency so clean native release builds include
  the runtime-loaded test DLL instead of relying on a previous local build.
- Publish release downloads only after the ZIP, checksum, and installer have all
  uploaded, and avoid duplicating the suite build for the same tagged commit.

## 0.5.0 — 2026-11-04

- Added native DMDc and fourth-order-default EDMDc, manual polynomial
  eigenfunction construction, explicit continuous/discrete spectra, and
  separate training/validation dataset references.
- Established immutable echo-history transforms and reusable basis/projection
  stages, coordinate-driven control and estimation, model/map identities, and
  physical-parameter compatibility checks while retaining the feedback laws.
- Added true interactive 3-D Viewer plots, repaired mathematical axis labels,
  and improved eigenfunction domain/time-history composition.
- Added standard Duffing force and pendulum torque observables.
- Added a checksum-verified prebuilt Windows installer, release packaging,
  bundled first-run examples, and instructions discoverable by repository name.
- Kept generated scientific plots, PDF exports, and design-source graphics
  local, recording pointers and hashes instead of including them in the current
  source checkout. Historical Git commits remain unchanged.

- Added a Python-free, read-only archive audit that reports repeated plans and
  byte-verified exact solution duplicates, and separated the operations DLL from
  the scientific core link graph so archive/job edits do not invalidate Viewer
  or workflow builds.
- Standardized the compiled job queue on versioned JSON responses and added
  explicit cancellation, making submit/list/status/recover/run-next suitable for
  native UI and automation clients without console-text parsing.
- Routed the optional Python compute CLI through the ABI-versioned compiled
  workflow service by default, retaining the former Python executor only as an
  explicit `python-reference` parity backend.
- Added Python-free workflow planning for ODE, delay, lattice, and PDE jobs,
  including hierarchy counts and storage estimates, and made native planning
  the desktop and optional Python CLI default.
- Made the public Python compute API delegate to the native workflow and
  renamed the former implementation as an explicit reference-only backend.
- Standardized DDE, lattice, PDE, and ODE solution metadata so every campaign
  trajectory carries its campaign ID, kind, point index, hierarchy path, and
  controller configuration directly.

This project follows [Semantic Versioning](https://semver.org/).

## 0.10.0 — 2026-11-26

- Consolidated the Solver's physical parameters, dimensionless groups, and
  derived scales into one divided, vertically scrollable parameter book, and
  promoted **Open Viewer** to a context-free global workbench action. Parameter
  values now edit exclusively in place and round-trip decimal, `e`, rendered
  superscript, and TeX-style scientific notation. The parameter book now uses a
  client-drawn blue-and-black ornamental scrollbar instead of Windows chrome.
- Moved controller editing out of the System page. Uncontrolled systems now
  expose a compact **+ add control** action; controlled systems expose separate
  System and Control folder tabs connected to a framed editor. Solver/execution
  settings now remain in an independent, job-owned lower region while Batch,
  System, Cluster, and Trajectory selection replaces only the upper editor.
  The lower region remains enabled and unchanged at every tree level, and its
  execution contract is broadcast to every System only when the job is
  serialized for planning, execution, queueing, or saving. Viewer now
  routes data-tree selections only to native payload views,
  so selecting archived data can no longer launch Visual Studio through a
  Windows file association.
- Corrected the folder-tab chrome so its frame is painted behind the scientific
  controls instead of occupying an overlapping child window. Parameter and
  Control pages remain unobstructed, and tabs use the requested rounded-leading,
  diagonal-trailing orientation.
- Aligned Viewer's completed-solve navigator with the Solver's exact Batch →
  System → Cluster → Trajectory hierarchy. Campaign/free/forced identity is now
  metadata rather than an intervening tree level, live archives replace stale
  bundled snapshots at launch, and legacy hierarchy-v1 NPZ records remain
  readable without migration.
- Added compiled Solver UI surfaces for the durable native job queue and archive
  maintenance. Jobs can be submitted, inspected, run, cancelled, and recovered
  without Python; archive audit and deduplication retain scoped dry runs,
  confirmation, reference-safety checks, and recoverable moves.
- Centralized queued execution in a reusable native workflow service shared by
  the desktop app and CLI, with throttled cancellation polling so adaptive RHS
  and controller evaluation do not acquire a database query per step.
- Added a machine-readable native UI parity ledger and release gate covering
  Solver, operations, Post-processor, and Viewer workflows, with source and
  regression evidence for every retained human capability.
- Added paged native Viewer dashboards for registered scalar-observable time
  histories and two-state phase-plane heat maps, plus initial/final lattice
  configuration views. Corrected KLQE routing so state time/phase requests are
  not masked by the same archive's observability matrices.

- Added a native system coordinate-scaling contract and enforced it for
  prescribed and built-in forcing. RPE, Marmottant, and Reddy--Szeri waveform/pulse
  pressure, Hz, and time settings now convert to the same nondimensional
  coordinates used by their dynamics; multitone waveforms and explicit
  model-coordinate overrides are supported. Solution archives expose control
  units, physical scales, time scales, and controller-setting provenance for
  Viewer.

- Added native observability analysis through the classical block matrix,
  finite discrete and continuous Gramians, and a stable infinite-horizon
  Lyapunov solve. Derived archives retain the Gramian, ordered observable
  directions and strengths, rank, conditioning, effective rank, and method
  provenance for Viewer and KLQE workflows.

- Added native classical state-space LQR, sampled local SDRE-style LQR, PID,
  and finite-horizon MPC beside the existing KLQR/KMPC family. Added exact
  sampled controller clocks for fixed and BDF integration, while preserving
  continuous adaptive-RHS `B_Psi(x)` and Riccati recomputation as the KLQR
  default. Controller-family/clock provenance is archived and exposed in the
  compiled CLI and Solver UI.

- Added Python-free native Mackey--Glass DDE and 2-D/3-D Ising plus 2-D Potts
  lattice workflows, each in an independently compiled DLL with canonical
  archives, system-family metadata, progress, and installed-package gates.
- Added an independently compiled native structured-field PDE engine and 1-D
  heat-equation workflow, including stability/resource gates, field ensembles,
  canonical archives, campaigns, and native Viewer space-time plots.
- Added deterministic multidimensional PDE partition plans, boundary-aware
  ghost regions, bounded local worker execution, and serial/partitioned parity
  gates to the independently compiled PDE module.
- Generalized native cartesian campaigns across ODE, DDE, and lattice solves,
  including system-parameter and indexed-array axes, monotonic aggregate
  progress, and correct free-versus-forced Viewer modifiers.
- Corrected the retained Python Potts energy observable to count both positive
  lattice bond directions, matching its Metropolis Hamiltonian and the native
  implementation.
- Completed the suite architecture remediation: all maintained Python and
  native sources are at or below 500 lines, static package cycles and measured
  boundary violations are zero, and the strict audit enforces a 70% outward
  public-documentation floor.
- Split Viewer rendering, JobManager services, Duffing bilinear search, KROM
  capabilities, and the remaining preferred-size scientific modules behind
  compatible façades and versioned contracts.
- Added neutral analysis-data, spectral, persistence, controlled-data,
  coordination, and estimation seams to remove reverse package dependencies.
- Generalized KROM training-domain planning and screening across registered
  systems while retaining RPE as a convenience preset.
- Added a low-overhead native Solver progress bar driven by existing throttled
  heartbeat and hierarchy signals.
- Added a shared native UI-support library and Help/About diagnostics for the
  exact project, archive, UI, and schema versions.
- Expanded the `dynsyslab` front door to solve, post-process, archive, jobs,
  estimation, and generic KROM-training workflows.
- Retired one-off parity/rerun scripts whose evidence is already retained, and
  documented the lifecycle boundary for reusable research tools.
- Added Python-version, installed-wheel, sdist, strict-audit, and native-UI CI
  gates; fixed the wheel to include the compute workbench and exclude runtime
  archive records.

- Added first-class state observables with analytic/explicit numerical
  Jacobians, plus validated Koopman observable reconstruction matrices `C`.
- KROMs can now retain exact state-dependent, fitted bilinear, and fitted
  constant control-input representations simultaneously; KLQR/KMPC selects
  and archives the representation explicitly.
- Added typed controlled-trajectory ensembles, Korda--Mezic multi-step output
  regression for a constant lifted `B`, continuous/discrete conversion, and
  continuous bilinear lifted-residual regression.
- Added the retained Duffing scalar-H-DMD/full-state-SINDy SOP, exact reference
  IC table, controlled fixed-RK4 reproduction with explicit source-input
  scaling, and the accepted minimal exact-KLQR transfer/release example.
- Added explicit automatic, sampled, and continuous-RHS controller evaluation
  policies. Continuous mode rejects stateful plug-ins and separates its full
  internal evaluation tape from the accepted-state control series.
- Replaced the floating sample clock with integer ticks and made controller
  call logs record actual evaluations rather than every held-command reuse.
- Added a reproducible, bounded feedback-only RPE/KLQR policy diagnostic CLI.
- Corrected the default RPE input to physical external acoustic pressure, where
  positive pressure compresses the bubble, across dimensional, dimensionless,
  Jacobian, KROM, and native-adapter paths.
- Added an explicit legacy input convention for old archive reconstruction and
  provenance-preserving analytic-KROM control-system rebinding.
- Added fast grouped polynomial KROM observability selection, immutable CLI
  publication, and coefficient-evaluated subset values that avoid loading full
  eigenfunction matrices.
- Added KROM training-domain reports and default error gates for out-of-domain
  Koopman LQR/MPC state or target evaluation.
- Added first-class affine state-observation maps with explicit real/imaginary
  lifted coordinates, held-out validation thresholds, immutable archive
  publication, extraction invalidation, and independent reconstruction
  capability gates.
- Added vectorized polynomial observation-map fitting and a dry-run-first RPE
  KROM training-domain planner derived from canonical target envelopes, with
  dimensionless/physical boundary provenance and a mandatory screening gate.
- Added a non-persistent canonical ODE-kernel boundary for bounded diagnostics,
  target-path-aware RPE campaign screening, content-addressed screening
  records, filtered candidate reporting, and a hard publication gate.
- Added provenance-preserving KROM extraction by retained index or exact
  eigenvalue/eigenfunction, with default conjugate-dual expansion, independent
  spectral copies, constant-input row restriction, and independent closed
  bilinear principal restriction over the selected eigenfunction coordinates.
  Extracted KROMs never retain or evaluate through a parent KROM.
- Added executable KROM archive loading and immutable subset reservation with
  dictionary, coordinate, preprocessing, and source-lineage preservation.
- Koopman LQR/MPC now lift canonical state-space targets internally while
  archiving a physical-state display target; dimensionless solves convert the
  controller time/state boundary explicitly.
- Added compact Koopman-controller archive metadata that avoids traversing the
  cyclic KROM/dual graph, and made controller sample periods consistently use
  the public physical-time solver boundary.
- Added optional, provenance-recorded physical input saturation to KLQR so
  long controlled RPE evaluations can enforce explicit pressure bounds.
- Added immutable canonical harmonic RPE control targets with analytic radius
  and radial-velocity references.
- Added a small-amplitude, equilibrium-start multisine RPE target for safe
  long-horizon controller evaluation without target-induced bubble collapse.
- Enforced physical equilibrium `(R0, 0)` as the initial condition for every
  controlled RPE solve, rejecting request- and system-level overrides.
- Controlled solves now archive target histories and exact controller-call
  logs alongside the existing solution-aligned control history.
- Viewer recognizes controlled solution archives, overlays the target as a
  dashed red state/phase trace, and displays actual control histories in
  separate time-series panels.
- Viewer now displays archived RPE pressure histories as `p_c/P0` when the
  pressure scale is available, keeps the physical Pa archive untouched, and
  uses full-width horizontal legends/control panels to avoid squeezed plots.
- Added an explicit `auto`/`python`/`native` numerical execution policy and
  backend provenance to solver results.
- Added a separate optional `dynsyslab-native` C++17 wheel with pybind11
  bindings for single and batched Duffing/van der Pol Euler and RK4 kernels.
- Added reference/native numerical parity tests, native capability discovery,
  clear no-fallback behavior for explicit native requests, and cross-platform
  cibuildwheel CI.
- Added a versioned benchmark CLI, runtime/SIMD fingerprints, and same-runtime
  performance admission/regression gates.
- Retained the compiler-free Python/SciPy wheel and reference execution path.
- Extracted one shared result lifecycle so scalar and native-batched hierarchy
  solves preserve identical interpolation, archives, Viewer notifications,
  postprocessors, telemetry, and failure handling.
- Added single-call native cluster integration while retaining one canonical
  archive per trajectory and fixed hierarchy manifests to honor configured
  archive roots.
- Added Job resource-to-thread-pool limits to prevent BLAS/OpenMP
  oversubscription in serial and multi-process execution.
- Added autonomous KROM rollout and reconstruction with a vectorized NumPy
  reference; retained the slower comparison C++ kernel as explicit-only.
- Added an optional official SUNDIALS/CVODE adaptive backend, version and
  algorithm provenance, and cross-platform optional-backend CI.
- Added an explicit native 2-D Ising kernel with recorded PRNG semantics; auto
  retains the established NumPy stream because stochastic parity is not
  bitwise.
- Added within-run relative performance admission gates.
- Added a coordinate-explicit dimensionless Rayleigh--Plesset Euler/RK4 kernel
  to native 0.2.0, including natural-period/pressure-driver scaling, Python and
  end-to-end parity tests, model capability provenance, and a representative
  25,000-step performance gate.
- Consolidated pooled-trajectory Hankel DMD onto the shared SVD/rank/resource
  kernel, exposed the same SVD driver and tolerance controls, and recorded
  explicit no-cross-boundary shifted-pair provenance.
- Fixed autonomous KROM construction so unactuated plug-in systems can evolve
  and reconstruct without being forced to provide a control input matrix.
- Added an explicit, provenance-recorded H-DMD SVD thread limit and structured
  BLAS thread-pool benchmark fingerprints; the default continues to inherit
  caller/Job resources because optimal threading is hardware- and shape-specific.
- Added a separate reduced-eigendecomposition thread limit after exact
  2,000-mode operator benchmarks showed a different thread optimum from SVD.
- Added native adapter/binary version-skew gating so `auto` falls back instead
  of dispatching dimensionless RPE into an older 0.1 extension; legacy Duffing
  and van der Pol capability remains recognized.
- Added consistent `python -m` entry points for performance, post-processing,
  solver diagnostics, archive management, and the desktop GUI, backed by
  subprocess CLI smoke tests.
- Converted the public `solver` package into a lazy façade so importing request
  types or backend diagnostics no longer loads SciPy integration engines,
  archive sessions, hierarchy execution, DDE, PDE, and lattice pathways.
- Separated archive capability gateways from the archive-management façade and
  made SciPy interpolation lazy, so Viewer, solver, campaign, job, and
  post-processing imports no longer pull in auditing/storage machinery merely
  to obtain read/write authority.
- Converted the top-level post-processing package into a lazy façade: archive
  primitives import immediately, while PCA/POD, Koopman workflow orchestration,
  and SciPy resampling load only when requested.
- Converted `archiver` exports to lazy resolution so importing a read/write
  gateway no longer initializes audit or storage-management implementations.
- Converted the Viewer package façade to lazy exports so routing/API discovery
  no longer initializes plotting, sessions, NumPy arrays, or archive access.
- Converted the umbrella `dynsyslab` module to lazy subsystem/value exports,
  preserving its convenience API without coupling every import to all engines.
- Made the optional SUNDIALS extra's declared nanobind runtime dependency
  explicit and repaired the local environment so package consistency checks
  pass.
- Fixed fixed-step PDE Euler/RK4 dispatch and duplicate times at adaptive
  checkpoint segment boundaries.
- Extracted post-processing algorithm execution and derived-publication rules
  into closed Koopman and decomposition handler registries, leaving the
  processor responsible only for workflow ordering, source/coordinate policy,
  lineage, and timing.
- Added typed post-processing source resolution, analysis preparation, and
  per-trajectory workflow context; repeated attachments reuse one prepared
  coordinate view, invalid resampling controls fail explicitly, cyclic source
  graphs are rejected, and spectrum/archive lineage advances atomically.
- Added the conventional `PostprocessingPipeline` public API while retaining
  lowercase `processor` as a compatibility alias, and routed read-only planning
  through the same coordinate/resampling preparation contract as execution.
- Replaced mutable process-wide standard workflow defaults with the immutable
  named `TRAJECTORY_HDMD_SINDY_POLICY`; preset instances receive independent
  setting copies and propagate the policy name through descriptions, run
  summaries, and result attachment provenance.
- Added a typed complete-DMD linear-algebra backend boundary. The admitted
  SciPy/LAPACK path records its actual BLAS provider; an optional explicit-only
  PyTorch CUDA path keeps transfer, full `gesvd`, projection, eigendecomposition,
  and mode construction on-device with memory preflight and no silent fallback.
- Added backend discovery to diagnostics and `dynsyslab-postprocess backends`,
  `--linalg-backend` pipeline selection, and paired full-spectrum SciPy/CUDA
  benchmark cases for parity and performance admission.
- Made the standard H-DMD-to-SINDy policy publish coefficient-first SINDy
  archives. Polynomial eigenfunction sample arrays are now exact virtual arrays
  reconstructed lazily. Viewer shows the full spectrum and every retained KROM
  point while reconstructing only five dashboard modes by default; explicit
  display selection never truncates the model. Repeated Plotly coordinate arrays
  are stored once in the page payload. Generic
  archive calls keep their materialized-array default for compatibility, while
  the CLI exposes `--sindy-storage full` as an explicit override.
- Added archive-size planning and virtual-array contract tests. On the complete
  25,001-sample, 2,000-eigenfunction RPE result, the new payload measured
  459,218 bytes and 0.036 seconds to publish versus the prior 419,158,070-byte,
  24.30-second materialized result.
- Removed the same redundant matrix from the standard SINDy algorithm's live
  result: coefficient-first runs now defer training-value evaluation, avoiding
  a 25,001-by-2,000 complex allocation (about 763 MiB). Callable polynomial
  eigenfunctions and KROM behavior are unchanged, and full storage reconstructs
  the exact matrix when explicitly selected.
- Added a selectable lazy-eigenfunction reconstruction benchmark and
  `dynsyslab-benchmark --case NAME` / `--list-cases`, so focused backend and
  Viewer regressions no longer require running the entire benchmark suite.

## 0.1.0 — 2026-08-09

- Established the first version-controlled dynsyslab release.
- Added a durable agent-harness record and chat-snapshot policy.
- Added local historical chat archives with optional S3 mirroring.
