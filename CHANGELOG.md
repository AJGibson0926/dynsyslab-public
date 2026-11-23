# Changelog

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

## Unreleased

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
