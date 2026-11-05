# Using the DynSysLab Solver and Post-processor

## Launch the native workbench

```powershell
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\DynSysLab.exe"
```

The Start-menu and taskbar shortcuts point to this stable installed executable.
The workbench opens maximized so its dense scientific controls can use the
available display without manual window setup. Its DPI-aware layout separates
the proposed hierarchy from the selected node's System or Control page;
Post-processor similarly separates stage selection from ordinary and advanced
settings.
Use the in-app **Refresh** action to re-index systems, presets, KROMs, and
archives. Developer rebuilds are installed to the same path, so the shortcuts
remain valid without a Python launcher. The native installer records the shared
project and archive roots in `%LOCALAPPDATA%\DynSysLab\runtime.json`; both apps
read that small configuration directly, and explicit command-line roots still
take precedence.
The pinned executable and its application windows use the shared DynSysLab
phase-portrait icon from `native_ui/assets`; Viewer uses the same product mark
from its own executable. The taskbar groups each application by its executable
path, without a separate shell registration that can become stale.

The left side is mode-specific. **Post-processor** keeps the same searchable,
system-centric archive tree used by Viewer. **Solver** replaces it with an
editable proposed-solve tree: Job → Batch → System → Cluster → Trajectory.
The right side shows metadata/relationships and the latest plan or execution
result. Solver and Post-processor are separate central tabs. Drag either
vertical divider to resize the active tree, central work area, and right-side
context panel; the dividers enforce useful minimum widths and resize live.

**Save job…** validates the active Solver or Post-processor definition and
writes its versioned JSON contract to disk. **Open job…** accepts the same file
and offers either a read-only plan/validation pass or immediate execution.
This is the same contract accepted by `dynsyslab-cli`; saved workflows therefore
remain reusable without Python or manual form reconstruction.

**Native jobs…** opens the durable SQLite-backed queue inside the app. It can
submit the active panel definition, inspect exact queued manifests and status,
run the next job asynchronously, cancel a queued or running job, and recover a
stale worker after an explicit age. Queue execution uses the same workflow
service as **Run and archive**; cancellation is checked at throttled numerical
safe points so it does not add a database query to every RHS evaluation.

**Archives…** exposes read-only archive auditing and recoverable exact-duplicate
cleanup. The dialog scopes by age, section, system, and derived records. Preview
uses the full reference-safety gate; **Apply** requires confirmation and moves
only verified, unreferenced duplicates into archive-local recovery storage. It
does not permanently delete files or rewrite provenance.

**Import data…** is the compiled front door for the retained tagged DynBatch
CSV format. It previews and validates the hierarchy before publication, keeps
the source file read-only, records its SHA-256 identity and source line ranges,
and writes canonical solution, hierarchy, and import records into the shared
archive. Legacy dynamics names are preserved by default. Mapping them to an
executable native system, preset, and parameter set is explicit because a name
alone is not proof that the equations agree. Controlled imports retain both the
full sample-aligned control history used by Viewer and the initializer-free
interval inputs used by control-map regression. Newly imported records appear
in Solver and Viewer after the catalog refresh.

The same service is available to automation without becoming a CLI-only human
workflow:

```powershell
dynsyslab-cli import legacy-csv data.csv --archive-root archives --label experiment
```

**Open Viewer** launches the native Viewer on the shared archive tree. Viewer provides
header-level **Copy plot** and **Save PNG/PDF** actions for every native plot view.
Viewer's **Post-process** action makes the reverse transition: it opens this
compiled application on the Post-processor tab with the selected archive
already populated as the source.
Viewer's tree contains only objects with native scientific display routes.
Legacy HTML/PDF reports and generic JSON records stay on disk but are not sent
through Windows file associations, so tree selection cannot launch Visual
Studio or another unrelated application.

**Help / About** reports the native UI and JSON contract versions together with
the exact project and archive roots provided at launch. It can open the source
checkout's Start Here guide when that guide is available.

The footer progress bar is indeterminate only while a backend cannot expose a
safe intermediate fraction. Native ODE, DDE, PDE, and lattice execution publish
throttled progress, and ensembles advance at trajectory archive boundaries.
The UI reads one small atomic status file four times per second; no numerical
state is copied into the UI process.

## Solver panel

The default draft is **New Job** with RPE selected. The **Job name** above the
Solver tree names the complete job, not an individual System. Use the buttons
below the tree to add Batches, Systems, Clusters, or Trajectories, or remove the
selected node. The upper middle workspace follows the selected tree level:

- **Batch** summarizes its Systems, Clusters, and Trajectories.
- **System** owns equations, presets, physical parameters, control, target, and
  KROM-training configuration. Coordinates, time, solver, and campaign policy
  belong to the complete proposed Job.
- **Cluster** owns its initial-condition generator, seed, and distribution. Its
  right-hand plot previews all proposed member ICs in a square phase plane for
  2-D systems or a projected 3-D view for 3-D systems. Plot ranges fit the
  generated points automatically.
- **Trajectory** shows its exact generated IC and owning Cluster.

The middle workspace has two independent vertical regions. The upper editor is
replaced by Batch, System, Cluster, or Trajectory content as selection changes.
The lower **Solver and execution** region remains visible and enabled at every
level. It contains the one job-level coordinate, integration-clock, numerical,
campaign, and controller-evaluation policy. Tree selection and System/Control
tabs never reload, disable, move, or repaint this region. When a job is planned,
run, queued, or saved, that execution contract is copied into every System solve
record so the existing `solve`/`solve_set` schema and native engines remain
compatible. The final job action bar remains fixed below both regions and
operates on the complete job. Cluster and Trajectory views do not expose System
parameter editing.

An uncontrolled System editor begins with a compact **+ add control** action and
contains no empty controller block. Adding control creates connected folder-label
**System** and **Control** tabs and selects Control. Their active tab opens into
the fully framed upper editor body. Choosing `none` in the Control page's
controller selector removes that page and restores the compact uncontrolled
editor.

The **System** page contains the system and preset selectors and three compact
scientific grids: physical/model parameters, dimensional quantities derived
during nondimensionalization, and genuinely dimensionless groups. Symbols use
conventional
mathematical abbreviations and subscripts instead of internal snake-case
identifiers. Units are omitted entirely from the dimensionless grid. All three
views come from the same compiled system resolution API used by execution and
refresh immediately after preset or override changes.

The **Control** page changes its entire configuration surface with the selected
control law. Prescribed waveform and random-pulse inputs are open-loop and
therefore have no target. Feedback laws expose a separate target workspace
only when the law can track a reference. Koopman feedback adds compatible KROM
and `B_Psi` selectors; classical LQR exposes `A`, `B`, `Q`, and `R`. Controller
evaluation is a per-System Control setting: **continuous RHS** is the concrete
default, while choosing **sampled** reveals the associated **Δt control** field.
Uncontrolled and open-loop Systems ignore those fields.

Editable values follow a read-mostly interaction: they render as mathematical
text until selected. Clicking a value temporarily creates an editor in that cell;
Enter or clicking elsewhere validates and commits the change, while Escape
discards it. The physical-parameter grid and numeric System-panel fields share
this interaction. Derived quantities remain read-only.

Scientific rows use the compact `symbol = value unit` form and default to four
significant figures. The **long** checkbox at the lower right of the parameter
section expands all three grids to round-trip double precision. It changes only
presentation; archived and in-memory numerical values always retain full
precision. The single parameter book scrolls vertically with DynSysLab's
blue-and-black ornamental rail; it does not use Windows scrollbar chrome or
per-family horizontal scrollbars.

Scientific canvases are painted as complete off-screen frames and presented in
one copy, so parameter scrolling does not expose an intermediate cleared
window. System/Control page changes use descendant compositing and suppress
redundant visibility repaints; the independent Solver and execution region is
not erased during the transition. Dropdown-list faces use the same
blue-and-black ornamental treatment as the scientific scrollbar while
retaining native popup, keyboard, focus, and accessibility behavior.

**Plan**, **Run and archive**, **Save job…**, and **Native jobs…** operate on
the whole proposed tree. They encode `kind: "solve_set"` in the existing
`dynsyslab.compute-job/v1` envelope. Every proposed Batch must contain at least
one configured System. Execution publishes all Trajectory archives plus one
merged hierarchy manifest per Batch, so the completed archive tree mirrors the
proposal without duplicating numerical arrays. Older single `solve` job files
remain loadable and appear as a one-Batch, one-System tree.

The proposal is configuration, not scientific data. Ordinary editing and
planning keep it in application memory and do not write the archive. A user may
explicitly save a reusable job definition outside the archive root or submit it
to the durable operational queue; neither appears in Viewer/Post-processor's
scientific data tree. The archive scanner rejects compute-job definitions and
the queue's `job_manager` state. Only **Run and archive** (or execution of an
explicitly queued job) publishes solved Trajectories and their completed
hierarchy.

1. Select a dynamical system and either **model default** or a named preset.
2. Inspect the single parameter book in the System page. It presents physical
   parameters first, derived physical scales second, and dimensionless groups
   last. Click a physical-parameter value to edit it in place; explicit changes
   become overrides while preset defaults remain authoritative. One compiled
   resolver serves ODE, DDE, PDE, and lattice systems.
3. The System's concrete coordinate convention is selected from its declared
   default when the System is chosen; the UI does not expose a deferred
   `automatic` value. Choose its time interval and output resolution directly
   below the parameter book. These values belong to that System draft,
   so Systems in one proposed job may retain different dimensional or
   dimensionless coordinates and different sampling intervals.
4. Configure job-wide numerical execution directly in the persistent **Solver
   and execution** region. Method, initial/fixed step, adaptive tolerances,
   maximum internal step, maximum internal-step count, and backend are inline
   controls. There is no modal solver
   dialog or hidden JSON-editing step. Common target states, Q/R diagonal
   penalties, and MPC limits accept ordinary comma-separated values; JSON
   remains optional interchange rather than a prerequisite for normal operation.
   **Build campaign…** sits below the proposed job tree, above its hierarchy
   construction buttons. PDE partitioning is not a job-wide solver setting: its
   shape, ghost width, and worker count appear on the System page only when that
   System is a PDE.
5. Select a Cluster and choose its IC mode:
   - `default`: repeat the model/controller default.
   - `explicit`: one state, history, field, or lattice per trajectory.
   - `uniform`: independent uniform draws from per-state `[low, high]` ranges.
   - `normal`: clipped normal draws centered in those ranges; each range spans
     six standard deviations.
   - `sphere`: an isotropic shell or annulus about a supplied center.
   - `ellipsoid`: an axis-scaled shell or annulus about a supplied center.

   These modes retain the useful continuous-state concepts from the earlier
   MATLAB research suite without copying its implementation. Poisson sampling
   is intentionally omitted because it is a count distribution rather than a
   useful general distribution for continuous state coordinates.

   For lattice systems, `default` and `uniform` generate discrete sites from
   the recorded seed. `explicit` accepts one nested lattice, or one outer row
   per trajectory. DDE initial values are constant history functions. The
   heat-equation PDE accepts a default sinusoidal field, an explicit field per
   trajectory, or one shared/per-cell random range.

   ODE Cluster generation is provided by the compiled workflow API, not by UI
   code. The UI sends the proposed solve specification to
   `resolveOdeInitialConditions`, and execution calls that same function.
   Consequently the seed, membership, and plotted points are exactly those the
   eventual solve will use. Generation and preview write no archives. The same
   API is available without the UI:

   ```powershell
   dynsyslab-cli generate-ics solve.json
   ```

   It returns `dynsyslab.initial-conditions/v1` JSON and performs no
   integration. Optional language adapters can call the stable
   `dsl_workflow_generate_initial_conditions_json` C ABI instead of launching
   the CLI. **Edit ICs** is family-aware. ODEs use named state components, DDEs edit the
   constant history value, and PDEs edit each grid value with its physical
   coordinate. Lattices use a compact pattern/custom-site builder so a 3-D
   field does not create tens of thousands of desktop controls; uniform random
   fields stay seed-driven. Explicit mode edits each trajectory independently;
   uniform and normal continuous modes edit componentwise sampling bounds while
   leaving deterministic seeded sampling in the compiled solver.
   Controlled-system defaults resolve through
   `controlledInitialState()` (for nondimensional RPE, `(1,0)`), while free
   defaults continue to use the system's normal initial state. The adjacent
   JSON is an inspectable advanced representation of what these native editors
   serialize, not a prerequisite for normal setup.
5. Add or remove Trajectory and Cluster nodes in the tree. Trajectories are
   divided deterministically into contiguous, balanced Cluster membership
   ranges. Each Cluster may use a different generator and seed.
   For a KROM-guided expansion, switch to Post-processing and use **Expand KROM
   training domain…** instead of
   assembling the batch manually. The native dialog reads the selected KROM's
   retained training states, shows their current axis-aligned envelope, accepts
   a required envelope, grid density, margin, and mandatory anchor states, and
   previews the exact number of proposed initial conditions. Applying the plan
   switches to Solver and appends a new Batch containing the source System, one
   Cluster, the exact Trajectory/IC set, full resolved parameters, coordinates,
   free-training archive role, and integration horizon. Before any solution is published, each IC is integrated in an
   independent short native screening solve with optional componentwise lower,
   upper, and absolute-state limits. The submitted IC grid is re-derived from
   the recorded plan, so later manual edits cannot silently detach the solve
   from its provenance. A rejected candidate blocks the batch before solution
   publication; accepted solutions retain the source KROM and their screening
   diagnostics. Screening progress advances only between trajectories and
   occupies the first ten percent of the ordinary Solver progress bar.
6. Choose free dynamics, prescribed waveforms/random pulses, PID, classical
   LQR, state-dependent LQR (`sdre`), finite-horizon MPC, KLQR, or KMPC. The
   panel is controller-specific rather than a shared form with irrelevant
   fields:

   - **Waveform** is a predefined open-loop input. Its left side previews the
     planned signal; its right side edits component frequencies, amplitudes,
     phases, offset, coordinate convention, and active window. It never creates
     a target object or target-history array.
   - **LQR** exposes only `A`, `B`, `Q`, and `R`. Explicit row-major `A` and `B`
     are accepted; leaving them absent preserves the system-derived
     linearization and input matrix.
   - **KLQR/KMPC** expose the KROM, `B_Psi` method, `Q`, `R`, and the applicable
     predictive settings.
   - **Feedback target workspace** selects a fixed point, an archived trajectory,
     a kinematic waveform, or a compatible canonical preset. It is absent for
     prescribed inputs.

   Every Riccati controller accepts scalar, diagonal-vector, or full-matrix Q
   and R. KLQR/KMPC expose the KROM and `exact`, `constant`, or `bilinear`
   control map as first-class controls. The KROM selector is filtered to the
   active system (plus legacy models without an archived system label), and
   RPE-only canonical targets appear only for RPE. MPC/KMPC likewise expose their horizon
   and optional scalar/vector control limits. **Controller evaluation** chooses
   automatic, continuous adaptive-RHS evaluation, or an exact sampled clock;
   a separate sample-time field overrides the output step. Less common
   controller-specific options remain available in the advanced JSON box.
   **Configure…** provides controller-specific native fields: physical or
   model-coordinate multisine amplitudes/frequencies/phases/offset and active
   windows; random-pulse amplitude, duration, seed, and window; PID gains; LQR
   linearization state; and optional Koopman source-mode/external-map details.
   The dialog clears mutually exclusive physical/model keys and preserves
   unrelated advanced settings.
7. For a feedback law, configure the target in its dedicated side workspace.
8. Optionally supply a cartesian campaign over system parameters, controller
   settings/targets, time, solver settings, or IC fields. Every point is
   archived through the same hierarchy and grouped by one campaign identity;
   each point may still contain a clustered IC ensemble. **Build campaign**
   provides a native axis editor with system-parameter and common forcing/time/
   solver/IC path suggestions. Values may be entered as a JSON array or a
   comma-separated list; arbitrary workflow paths remain available.
   In Viewer, Batch and Cluster nodes can be opened directly as ensemble time,
   phase-space, and (when present) control-history plots. A parameter-preserving
   Campaign also provides terminal-response, amplitude, and asymptotic
   tail/bifurcation summaries. Physical-
   parameter changes remain separated by equations fingerprint even when they
   were launched by one campaign.
9. Use **Plan / validate** before **Run and archive**.
   For unattended or ordered work, use **Native jobs…** to queue this same
   validated definition and run it through the persistent native scheduler.

The Post-processor's **Configure standard H-DMD → SINDy** dialog exposes the
canonical Koopman chain without requiring JSON editing. H-DMD accepts any
scalar observable declared by the selected source system, with an editable
fallback for legacy/custom names, and exposes matrix size, sampling window,
stride, absolute SVD tolerance, optional rank cap, and terminal-amplitude
energy sorting. Rank zero is the standard full-spectrum policy. SINDy exposes
observable coordinates, reconstruction observables, polynomial degree,
sparsity threshold, iteration count, exact-RHS versus finite-difference
generator derivatives, normalization, conjugate-fit reuse, and optional
training-value materialization. The standard path stores coefficient-first
eigenfunctions and performs `lambda=log(mu)/delta_t` inside the compiled
workflow service. Unrelated advanced JSON keys survive structured edits.

**Optimize eigenvalues** is an independent, spectrum-only native stage. Its
dialog selects a seed spectrum, one registered scalar or vector observable,
additional compatible trajectory archives, optional mode groups, complex-plane
bounds, L-BFGS memory, iteration limit, tolerance, and conjugate-pair policy.
The compiled variable-projection service eliminates each trajectory's boundary
amplitudes analytically while optimizing the discrete-time eigenvalues. Its
archive retains both the seed and optimized spectra, the exact conversion
`lambda=log(mu)/delta_t`, per-observable residuals, iteration/evaluation counts,
convergence status, all trajectory sources, and the seed-spectrum relation.
Viewer provides a dedicated diagnostic dashboard, and every plot in that view
uses the ordinary clipboard and PNG export actions. The stage deliberately
constructs no eigenfunctions; place it before SINDy or fitted boundary
functions in the ordered job to select the next representation explicitly.

An H-DMD result is also a self-contained scattered-interpolation KROM. It
stores the aligned physical training states, modal/eigenfunction values,
continuous and discrete spectra, and the H-DMD reconstruction row for the
selected scalar observable. Native evaluation shares one k-d-tree neighbor
search across all modes. Its eigenfunction descriptor records
`scattered-interpolation` with a finite-difference gradient backend, so the
same extraction, exact `B_Psi`, control-regression, KLQR/KMPC, validation, and
Viewer services used by analytic SINDy KROMs can consume it. No Python
interpolator or parent archive is required at runtime.

**Boundary eigenfunctions** is a separate construction stage. It accepts an
immutable discrete/continuous spectrum produced by H-DMD, eigenvalue
optimization, or another supplied source; one or more free trajectories on a
common sample grid; and a scalar or multi-output registered observable. The
stage regresses only trajectory-specific boundary amplitudes. Mode groups map
the supplied eigenvalues to observable outputs, and the resulting scattered
eigenfunctions carry a reconstruction `C` map whose assigned modal weights are
one by construction. The native dialog exposes the spectrum, observable,
additional trajectories, mode groups, interpolation neighbors, finite-
difference step, and SVD tolerance. It cannot mutate the input spectrum.

KROM extraction and control-map regression also have structured native
editors. **Select modes…** reads the source KROM spectrum, supports extended
multi-selection and conjugate-partner completion, shows inherited ancestral
mode IDs, and submits local zero-based indices to the independent extraction
service. **Configure…** beside control regression filters candidate KROMs by
system and exposes the regression algorithms independently. The constant map
may use the Korda--Mezić multistep observable-output fit, a continuous lifted-
residual fit, or a pointwise exact-map mean. The multistep method requires a
named reconstruction `C`, retains its raw discrete `B_d` and sample interval,
and converts it to the continuous `B_c` used by KLQR. The bilinear map may use
controlled lifted residuals or a pointwise exact-map projection, with optional
conjugate-complete features and state subsampling. A native multi-file picker
adds any number of compatible trajectory archives to one pooled fit; the
workflow rejects sources whose resolved equations fingerprint differs. The
same editor exposes
ridge regularization, optional mode selection, 2-D error-map bounds, and grid
resolution. Each run produces a self-contained KROM plus native Viewer heat
maps against exact `B_Psi`; methods and time-domain conversion are explicit in
the archive metadata.

**KROM control search…** in Post-processing opens the bounded control-candidate
campaign without requiring a hand-written job file. Select a parent KROM, add one or more controlled
trajectory archives for pooled control-map regression, and enter the initial,
controlled-target, and post-release reference states. The dialog exposes the
integration method and tolerances, continuous-RHS or sampled controller clock,
KLQR Q/R weights, reconstruction `C`, constant and bilinear regression methods,
error-grid resolution, success tolerances, and per-case timeout. Q and R accept
a scalar or comma-separated diagonal. The campaign can be named, planned,
saved, queued, or run entirely inside the compiled app.

Candidate screening is read-only and does not consume an attempt. It ranks
system-agnostic sub-KROMs using target separation, exact input authority,
controllability across the initial/midpoint/target states, spectral spread, and
mode-count cost. At most 25 retained candidates may proceed through the costly
sequence: independent KROM extraction, pooled control-map regression, and
otherwise identical exact, constant, and bilinear KLQR solves. This is a hard
workflow limit, not merely a dialog suggestion. Every completed attempt is
published immediately to a durable ledger; the final campaign archive retains
mode indices, screening evidence, release/final-state errors, settled RMSE,
peak control, control-map error summaries, success decisions, and links to all
independent result archives. A cooperative timeout bounds each control case;
the numerical integrators check cancellation at safe points.

Viewer recognizes the final archive as a **Control campaign** and displays a
six-panel comparison of release error, final error, peak control, mean
`B_Psi` error, screening score, and retained mode count. The ordinary native
copy, PNG, and PDF actions apply to this dashboard. Selecting a linked result
archive opens its detailed controlled trajectory or control-map heat maps.

**Validate KROM** consumes a free trajectory and an independently selected KROM
for the same system. Its native editor marks the dataset as `training` or
`test` and accepts explicit semigroup horizons or chooses them from the output
sampling and duration. It never refits the eigenfunctions. The resulting
archive contains per-mode generator, one-step, anchored long-term, and
finite-horizon errors, together with the fraction and maximum excursion of
samples outside the KROM's recorded training bounds. Viewer opens this as a
dedicated validation dashboard and preserves a direct lineage link to the KROM.

**Configure decompositions…** provides independent native settings for DMD,
HAVOK, EDMD, SVD, PCA, POD, bagged PCA, robust PCA, and Hankel-POD. Each stage
owns its own rank and centering policy, so one method cannot silently change
another. HAVOK accepts a state index or any scalar observable published by the
selected system. EDMD exposes polynomial and RBF dictionaries, regularization,
and explicit RBF centers.

The **Observability / KLQE Configure…** editor filters KROMs by the selected
trajectory's system and discovers the KROM's named reconstruction maps (its
observable-specific `C` matrices). It exposes optional mode selection,
process/measurement covariance scalars or diagonals, detectability policy, all
four observability constructions, numerical rank controls, and the KLQE
initial estimate.

The controller names preserve two distinct coordinate spaces:

| Kind | State used by the feedback law | Gain update | Automatic clock |
|---|---|---|---|
| `lqr` | dynamical state `x` | one CARE at `linearization_state` (or the first IC) | continuous RHS |
| `sdre` | dynamical state `x` | local Jacobian/input map and CARE at each update | sampled |
| `mpc` | dynamical state `x` | local finite-horizon Riccati recursion | sampled |
| `pid` | dynamical-state tracking error | stateful P/I/D update | sampled |
| `klqr` | split-real Koopman state `Psi(x)` | recompute `B_Psi(x)`, CARE, and `K(x)` | continuous RHS |
| `kmpc` | split-real Koopman state `Psi(x)` | recompute the selected `B_Psi` and finite-horizon recursion | sampled |

`sdre` is the suite's successive-local-linearization/SDRE-style state-space
mode: the native workflow uses the system Jacobian as the local `A(x)`. A
future system may expose an exact state-dependent coefficient factorization
without changing the controller API. `klqr` is not an alias for `lqr` or
`sdre`; it retains the established Koopman law
`u=-K(x)(Psi(x)-Psi_target)` and can use the objectively true, constant-fit,
or bilinear-fit `B_Psi` representation.

The **Controller evaluation** control writes `solver.control_evaluation` as
`automatic`, `continuous_rhs`, or
`sampled`. Sampled mode uses exact integer ticks and a right-continuous
zero-order hold; set `solver.control_sample_time` or
`controller.settings.sample_time`. Automatic mode keeps KLQR continuous so
adaptive RK45/BDF trial states recompute `B_Psi(x)` and the Riccati gain.
PID, state-dependent state-space LQR, MPC, and KMPC require a sampled clock.
An explicit sampled clock may also be selected for KLQR experiments.

Classical and Koopman controllers return the actual scalar system input `u`;
the dynamical system applies its own state-dependent input matrix. The archive
records `controller_family`, `control_evaluation`, `control_sample_time`, the
dynamics coordinate system, and the corresponding dimensionful or
dimensionless control units.

The Solver always publishes through the canonical archive lifecycle. It saves
the generated input history alongside every forced solution and saves a target
history only for a feedback solve that actually has a target. Prescribed
waveforms and random pulses do not publish a fictitious target.
Lattice archives additionally retain exact integer `lattice_initial` and
`lattice_final` snapshots beside their order-parameter/energy time series.
PDE archives retain the spatial coordinate array and complete field history;
Viewer renders them as a space-time heat map plus a final spatial profile.
For explicit PDE methods, the PDE System may own an optional local mesh
partitioning contract without changing the archive contract:

```json
"system": {
  "pde_partitioning": {
    "enabled": true,
    "partition_shape": [4],
    "ghost_width": 1,
    "local_workers": 2
  }
}
```

The partition shape has one positive extent per spatial axis. `local_workers`
is bounded by the shared execution thread limit. Every RK stage constructs
boundary-aware ghost regions and synchronizes before the next stage; metadata
records the execution mode, partition shape, halo width, and worker count.
MPI and adaptive refinement are not selected by this block. Older jobs with
`solver.parallel` remain readable and are migrated to this System-owned form.

## Post-processor panel

Select a trajectory, controlled dataset, or KROM in the archive tree; its path
is copied into the source field. You may also paste an archive path. Select
ordered stages and edit their settings:

- H-DMD accepts any eligible scalar system observable, not just state 0.
- Eigenvalue optimization accepts a seed spectrum and one or more compatible
  trajectories, and publishes a spectrum-only result with convergence evidence.
- SINDy consumes the spectrum from the preceding H-DMD/DMD stage or an
  explicitly supplied spectrum archive.
- HAVOK has independent delay, rank, and scalar-observable settings.
- DMD, POD, PCA, and SVD have independent settings and can run independently.
  The same panel also exposes
  seeded bagged PCA, principal-component-pursuit robust PCA, and Hankel POD.
  Their native archives retain reproducibility/convergence diagnostics and
  open directly as singular spectra, component maps, temporal scores, or
  low-rank/sparse traces in Viewer.
- **DMDc** jointly identifies a discrete transition and input map from recorded
  snapshots and controls. **Configure DMDc…** edits its independent settings;
  selecting it clears the default H-DMD/SINDy stages. Its core and compiled
  workflow also work without either desktop application (details below).
- **EDMDc** is the separate dictionary-lifted version. **Configure EDMDc…**
  defaults to the SINDy total-degree-4 polynomial library, including mixed
  monomials, a constant, and frozen coordinate normalization. Both DMDc and
  EDMDc can be selected for a side-by-side fit to the same source.
- `observability` and `klqe` can analyze the classical block matrix, a
  finite discrete Gramian, a finite continuous Gramian, or the continuous
  infinite-horizon Gramian. The last form is accepted only for a strictly
  stable lifted generator.
- **Extract KROM** accepts zero-based `selected_modes`. It copies the selected
  eigenvalues/eigenfunctions, compatible reconstruction columns, constant-map
  rows, and bilinear principal blocks into a new immutable archive. The new
  KROM embeds its fully resolved dynamical-system configuration and is loaded
  without consulting the parent at runtime. Selecting this operation clears
  the default H-DMD/SINDy boxes because it consumes a KROM rather than a
  trajectory. It also carries the reduced discrete spectrum and aligned
  training states, time, coordinates, and sampled eigenfunction values when
  available. Control-regression output inherits the same inspection evidence
  from its source KROM, so Viewer can render either result without reopening a
  parent KROM as a runtime dependency.
- **Validate KROM** evaluates existing eigenfunctions on the selected free
  trajectory. It distinguishes training-data fidelity from later test or
  extrapolation fidelity and publishes one-step, generator, anchored, and
  finite-horizon error evidence without altering the KROM.

Observability settings use `observability_method` with one of `matrix`,
`discrete-gramian`, `continuous-gramian`, or
`continuous-infinite-gramian`. `observability_blocks` is the integer horizon
for the matrix and discrete forms. `observability_horizon` is the model-time
horizon for the continuous finite form, and `observability_sample_time`
converts a continuous generator to a discrete transition for the discrete
form. The derived archive always retains the Gramian, ordered observable
directions, strengths, numerical rank, tolerance, effective rank, and
condition number; the raw block matrix is retained when that method is used.

**Run selected stages** never reruns the source dynamics. Published derived
archives immediately enter the shared catalog and are available in Viewer.
Their compact metadata carries direct source/KROM relationships; Viewer turns
those into clickable lineage links without inflating the numerical arrays.

### DMD with control: native computational component

DMDc is implemented in `native/include/dynsyslab/core/DmdControl.hpp` and
`native/src/DmdControl.cpp`, inside the independently compiled spectral module.
It reuses the compiled SVD, products, eigensystem and least-squares services;
it does not add Eigen templates to the workflow or UI. The thin
`workflow/native/src/DmdControlPostprocess.cpp` adapts recorded arrays to that
API and publishes the usual immutable derived NPZ. The UI only edits the same
JSON method settings. No existing controller, SINDy fit, or B-psi regression
implementation is replaced.

The implementation follows the unknown-input-matrix formulation in
[Proctor, Brunton and Kutz, *Dynamic Mode Decomposition with Control* (2016)](https://doi.org/10.1137/15M1013857).
It fits the sampled model

$$q_{k+1}\approx A_q q_k+B_q u_k.$$

Here **q is the supplied snapshot vector**, not automatically the physical
state. It can contain states, measured observables, or already computed latent
coordinates. DMDc does not construct an echo-history buffer, perform another
Hankel decomposition, differentiate measurements, or evaluate equations of
motion. Keep the existing frozen echo-history transform outside this module.
Using a latent array does not authorize the fit to look at r or r-dot.

Let X and Y contain paired snapshots as columns and U contain their interval
inputs. The algorithm computes a truncated SVD of Omega = [X; U], solves the
joint least-squares problem Y approximately equals [A_q B_q] Omega, and obtains
an output basis P from a separate SVD of Y. It stores a reduced model

$$z_k=P^Tq_k,\qquad z_{k+1}\approx Az_k+Bu_k,\qquad q_k\approx Pz_k.$$

The saved A and B are in **z coordinates**. Recover the represented full-space
maps as P A P^T and P B when needed; an unnecessary dense n-by-n A is not
archived. This model has no affine intercept and performs no centering or
automatic variable scaling. Use consistently scaled snapshots/inputs, and
include a constant feature explicitly if an affine model is intended.

#### What eigenfunctions does it provide?

Write A W = W diag(mu). The **right** eigenvectors W give reconstruction modes,
not eigenfunction coefficients. When W is sufficiently well-conditioned, the
core forms L = W^-1 P^T and the linear functions

$$\psi_j(q)=L_jq,\qquad
\psi(q_{k+1})\approx\operatorname{diag}(\mu)\psi(q_k)+W^{-1}Bu_k.$$

They are exact eigenfunctions of the identified autonomous reduced linear
model, and only candidate approximations to the real zero-input Koopman
eigenfunctions. If q = (xi_1, xi_2), these functions are linear polynomials in
those two latent variables. Output rank r gives at most r such coordinates;
two real input coordinates often yield one conjugate pair, not two independent
complex oscillations or a large nonlinear eigenfunction bank.
If each xi is a fixed linear projection of an echo
window, psi is consequently a linear function of that window. They are not
automatically nonlinear SINDy polynomial eigenfunctions. Select the separate
**EDMDc** stage below for a dictionary lift; ordinary DMDc remains unchanged.

While u is nonzero, the additive modal input term remains; pure exponential
behavior is expected only for the zero-input fitted model. Fitting A jointly
from driven data does **not** certify its spectrum against free trajectories
or replace our existing unforced eigenfunction-training bank.

If cond(W) exceeds 1e6, the archive retains A, B and the spectrum but omits the
unreliable modal coordinate/eigenfunction arrays and records a warning. Optional
continuous rates use the principal log(mu)/dt; they are omitted when a multiplier
is exactly zero. They are branch-dependent sampled rates, and do not turn B
or W^-1 B into a continuous-time B-psi matrix. No direct KLQR handoff is claimed.

#### Inputs, settings and headless use

The source must have uniform, finite, increasing `time[N]`, a selected real
`snapshot_array[N,n]` (a scalar `[N]` is also accepted), and recorded controls:

- Prefer `control_interval_values[N-1,m]`; row k drives q[k] to q[k+1].
- Otherwise accept `control_values[N,m]` under that same start-of-interval
  convention; the final row is unused. Scalar input `[N]` / `[N-1]` works too.
- Missing inputs, all-zero inputs, wrong shapes and non-finite fitting data
  are rejected. When both input arrays exist, the interval array takes priority.
- An instantaneous sample of a rapidly varying driver is not automatically
  equivalent to a held interval input. Resolve the sampling or supply a
  physically justified interval input before interpreting B.

The stage uses the chosen arrays in their **stored units**, with no automatic
physical-state nondimensionalization. This avoids treating two latent columns
as bubble radius/velocity. Source metadata, selected array names, input
convention, dt, settings and coordinate policy remain in the derived record.

```json
{
  "schema": "dynsyslab.compute-job/v1",
  "kind": "postprocess",
  "label": "DMDc measured coordinates",
  "sources": ["C:/data/controlled-measurements.npz"],
  "methods": [{
    "name": "dmdc",
    "settings": {
      "snapshot_array": "values",
      "augmented_rank": 0,
      "rank": 0,
      "svd_tolerance": 0,
      "save_eigenfunction_values": true
    }
  }],
  "archive": {"enabled": true}
}
```

Run it with `dynsyslab-cli run dmdc-job.json --archive-root C:/data/archives`.
The same job is accepted by `executeJob` and `dsl_workflow_run_json`; there is
no Python requirement or UI-only execution path. Zero ranks select numerical
rank; a positive `svd_tolerance` is an absolute singular-value cutoff, while
zero selects epsilon times matrix size times its largest singular value.
Input/output ranks are independent: the augmented rank must leave capacity
for both snapshot and input directions. Requested ranks above numerical rank
are capped and the actual retained ranks are recorded.

Multiple job sources produce independent models. For a pooled fit, the public
C++ API accepts explicit sample-major X, Y and U pairs, allowing callers to
concatenate valid pairs from several trajectories without introducing false
transitions between their endpoints.

```cpp
#include <dynsyslab/core/DmdControl.hpp>
using namespace dynsyslab::core;
// X,Y: pairCount x snapshotDimension; U: pairCount x inputDimension.
auto model = dynamicModeDecompositionWithControl(
    X, Y, U, pairCount, snapshotDimension, inputDimension, dt, options);
auto prediction = predictDmdControl(model, newSnapshots, newInputs, sampleCount);
auto psi = evaluateDmdControlEigenfunctions(model, newSnapshots, sampleCount);
auto polynomialPsi = model.polynomialEigenfunctions(); // shared degree-one representation
```

#### Saved results and interpretation

The `derived/koopman-dmdc/` record uses `dynsyslab.derived/v1` with
`model_schema: dynsyslab.dmdc/v1`. Its main arrays are:

| Array | Meaning |
|---|---|
| `operator`, `control_matrix` | Reduced discrete A and B |
| `projection_basis` | P, mapping z back to supplied q coordinates |
| `eigenvalues` | Discrete multipliers mu of A |
| `modes`, `exact_modes` | Projected P W and DMDc A-hat P W modes, respectively |
| `eigenfunction_coefficients` | Rows L, evaluated directly against q |
| `modal_input_matrix` | Discrete modal forcing map W^-1 B |
| `eigenfunction_values`, `time` | Optional sampled psi and its aligned times |
| `augmented_singular_values`, `output_singular_values` | Two SVD spectra |

Source snapshots and control histories are referenced, not duplicated. Viewer
classifies the result under Koopman spectra and explicitly labels it **not a
control-ready KROM**; existing polynomial KROM loaders are not repurposed.

Diagnostics report ||P(A P^T X + B U)-Y||_F / ||Y||_F, the output projection
error, retained condition numbers, numerical and retained ranks, and warnings.
These are **training diagnostics**, not evidence of held-out control accuracy.
Row-rank deficiency in [X; U] means full A/B separation is not identifiable;
independent excitation and representative trajectories are still required.
Full row rank is only a necessary algebraic condition, not proof of causal
identification, noise robustness, or Markov closure.

For echo-history experiments, validate on separately held-out free and driven
trajectories, including multi-step prediction after a pulse has ended.
Jointly fitting A/B can expose a different useful approximation, but it cannot
magically restore missing information or remove the forced-history closure
failure found in the pulse diagnostics. This addition does not change the
existing continuous-time controller or make claims of successful bubble control.

Native tests cover exact A/B recovery, complex left-eigenfunction identities,
forced modal evolution, free exponential evolution, rank truncation, pooled
pairs, scalar/multiple inputs, defective eigenbases, collinearity, input
alignment, invalid data, measurement-only archives and immutable publication.

### EDMD with control: shared polynomial representation

The separate `edmdc` stage implements the dictionary-lifted joint regression
described in [Korda and Mezić, *Linear predictors for nonlinear dynamical
systems: Koopman operator meets model predictive control* (2018), Sections 3–4](https://arxiv.org/abs/1611.03537).
It reuses `DmdControl` for the numerical regression and
`polynomialExponentTable` for the exact SINDy monomial ordering. Its public
computational API is `native/include/dynsyslab/core/EdmdControl.hpp`;
`native/src/EdmdControl.cpp` contains the lifting adapter, not a second SVD,
eigensolver, or controller implementation.

For supplied coordinates q, let s[j]=(q[j]-center[j])/scale[j]. The default
dictionary contains every monomial of total degree 0 through 4, in SINDy's
graded order. In two variables it has 15 features:

$$\Theta(q)=[1,s_1,s_2,s_1^2,s_1s_2,s_2^2,\ldots,s_2^4]^T.$$

Degree is total degree, not degree 4 independently on every axis. Four
coordinates give 70 features. Degree, feature count, retained model rank,
number of eigenfunctions, and echo-history length remain different quantities.
For safety this dense implementation permits degree 1–32 and at most 4096
features including the constant. Compress a long history before polynomial
lifting; do not polynomial-expand thousands of raw delays.

Normalization uses the mean and population standard deviation of X (the
current snapshots in the fitting pairs), with scales <=1e-15 replaced by 1.
Y and all later observations use those frozen values. `normalize:false` uses
raw variables, center=0, scale=1. This is statistical normalization, not
physical nondimensionalization. No equations of motion or hidden states are
consulted, and no delay/Hankel transform is constructed by this stage.

The model and its candidate eigenfunctions are

$$z=P^T\Theta(q),\qquad z_{k+1}\approx Az_k+Bu_k,$$
$$\psi_j(q)=\sum_\alpha L_{j\alpha}
 \prod_i\left(\frac{q_i-center_i}{scale_i}\right)^{\alpha_i}.$$

L=W^-1 P^T as for DMDc, but now its columns index monomials. Modes reconstruct
dictionary features, not directly physical coordinates. `predictEdmdControl`
decodes the predicted degree-one feature entries with the fixed affine inverse
normalization. To propagate the lifted linear model without re-lifting a
decoded prediction, use the stored A/B/P or `predictDmdControl` on lifted data.

This is a dense discrete-time least-squares fit, **not** the SINDy sparsity or
generator-PDE fitting objective. It selects its spectrum from A; the default
dictionary is shared, not the optimization method. The fitted coefficients
export through `EdmdControlResult::eigenfunctionSet()` to the shared immutable
`PolynomialEigenfunctionSet`, with explicit coordinates and a discrete spectrum.
Zero multipliers are retained without inventing finite continuous rates.
The legacy `polynomialEigenfunctions()` DTO adapter remains available for
nonzero multipliers. Neither export invents SINDy residuals or training metrics.
See [Eigenfunction and KROM API](architecture/EIGENFUNCTION_MODEL_API.md) for
manual construction, naming migration, optional physical-system adapters and
`dynsyslab-cli krom-info` inspection.

```json
{
  "schema": "dynsyslab.compute-job/v1",
  "kind": "postprocess",
  "label": "Compare DMDc and polynomial EDMDc",
  "sources": ["C:/data/controlled-measurements.npz"],
  "methods": [
    {"name":"dmdc", "settings":{"snapshot_array":"measured_latents"}},
    {"name":"edmdc", "settings":{
      "snapshot_array":"measured_latents",
      "dictionary":"polynomial", "polynomial_degree":4,
      "include_constant":true, "normalize":true,
      "rank":0, "augmented_rank":0, "svd_tolerance":0,
      "save_eigenfunction_values":true
    }}
  ],
  "archive":{"enabled":true}
}
```

The input alignment, sample-grid requirements, C ABI and CLI entry points are
the same as DMDc. No differentiation is needed. Each job source fits a separate
model; explicit paired C++ arrays support pooled fits without false endpoint
transitions. Other dictionary families are rejected rather than ignored.

Output is an immutable `derived/koopman-edmdc/` NPZ with
`model_schema:dynsyslab.edmdc/v1`. It retains the DMDc array names plus
`polynomial_exponents[features,coordinates]`, `state_center[coordinates]`,
and `state_scale[coordinates]`. `eigenfunction_coefficients[rank,features]`
uses that exact dictionary. The archive records the effective degree,
constant/normalization settings, input convention, source pointer, dt, and
separate lifted-space and decoded-snapshot one-step training errors. It does
not duplicate raw snapshots, control histories, or a large lifted training
matrix. Viewer identifies it as a polynomial input/output spectral model,
not automatically as a continuous-control-ready KROM.
`loadEdmdControlModel(path)` in `EdmdControlArchive.hpp` reloads the fitted
dictionary and model for prediction, evaluation and shared-polynomial export
without accessing its training source or a physical system.

**Limitations:** constant additive B after nonlinear lifting is approximate;
input-state products and nonlinear input terms need not lie in this model
class. Good snapshot error does not imply good lifted error, free eigenfunction
accuracy, multistep prediction, or successful control. Validate independently
on held-out free and forced trajectories. Joint EDMDc fitting does not replace
our unforced SINDy training policy. No controller or continuous B-psi is changed.

### Dataset lineage and autonomous validation

KROMs retain lightweight `dataset_references` using the versioned
`dynsyslab.krom-datasets/v1` contract. Each reference has `path`, an independent
`role` (`training`, `validation`, or `test`), and `dynamics` (`forced`, `unforced`,
or `unknown` for unclassified legacy records). `input_archive` optionally points
to separately recorded inputs. Save/load and mode extraction preserve these
references; neither data arrays nor forcing classifications are invented.
Forced training and unforced validation are therefore compatible in one KROM.
In-memory fits to anonymous arrays must be given their source references by the
caller with `KoopmanModel::addDatasetReference`.

For a state-coordinate EDMDc KROM, use the read-only reuse-first planner:

```powershell
dynsyslab-cli krom-validation-plan C:/data/krom.npz --archive-root C:/data/archives
```

The shared C++ API is `workflow/KromValidationPlan.hpp`. It searches existing
links and the archive catalog, supplements a stale catalog from `solutions/`,
matches resolved physical parameters/equations, and excludes this KROM's own
training sources. Matching free records become validation pointers and an
ordinary `koopman_validation` job. If none match, it returns `simulation-offer`
with the system parameters and missing IC/duration choices. It never simulates
without a separate approved solve. Native-solution EDMDc fits and extraction of
EDMDc-derived KROMs invoke this lookup. Arbitrary measured/latent snapshots need
their frozen coordinate provider; they are not silently reinterpreted as states.

`koopman_validation` checks unforced metadata and recorded inputs, full physical
parameter agreement, and coordinate compatibility before evaluating the frozen
model. Its default validation step is the KROM's fitted step. A finer source
grid is selected at an exact integer stride, without interpolation or rewriting
the dataset. Sampled-origin models require the fitted step or an integer multiple.
Optional `start`, `stop`, and `validation_delta_t` select the evaluation interval.
Forced test pointers remain valid provenance, but the autonomous defect must not
be applied to forced data as if its input were zero.

The existing shared `validateKoopmanSemigroup` computes one-step, anchored
long-term, finite-horizon, and finite-difference generator errors. A supplemental
`one_step_change_relative_l2` divides the one-step residual norm by the norm of
the observed change, rather than the signal norm; `one_step_change_defined`
marks a zero denominator. Generator errors are diagnostics of the saved
continuous logarithm branch, not proof of a continuous eigenfunction from
sampled data. Default validation publication creates an immutable,
coefficient-identical KROM revision containing the dataset link and compact
training bounds. The result's `validated_krom` points to that revision; the
original is not overwritten. Set `retain_validation_pointer:false` only when
the supplied immutable KROM already contains the intended references.

### Manual eigenfunctions and discrete/continuous eigenvalues

Fitted and manually supplied polynomials use the same `KoopmanSINDyResult`
and `KoopmanModel`, evaluators, gradients, reconstruction/control-map setters,
archive writer/reader and validation services. The historical SINDy name does
not require the coefficients to have come from SINDy.

```cpp
#include <dynsyslab/core/EdmdControl.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
using namespace dynsyslab::core;

// psi(q1,q2) = 2 + q1 + i*q2; lambda = -0.2 + 3i.
auto psi = KoopmanSINDyResult::fromCoefficients(
    2, {0,0, 1,0, 0,1}, {{2,0},{1,0},{0,1}}, {{-0.2,3}});
psi.sampleTime = 0.01;
auto mu = psi.discreteEigenvalues();      // exp(lambda * saved dt)
auto mu2 = psi.discreteEigenvalues(0.02);  // same lambda, different dt

// Alternatively provide sampled multipliers, an explicit dt, and optionally
// one integer logarithm branch per mode (default zero = principal branch).
auto sampled = KoopmanSINDyResult::fromDiscreteCoefficients(
    2, psi.exponents, psi.coefficients, mu, 0.01);

// A data-driven fit exports the identical polynomial class:
auto fit = extendedDynamicModeDecompositionWithControl(X,Y,U,pairs,n,m,dt);
auto fittedPsi = fit.polynomialEigenfunctions();

// For history coordinates, preserve their identity before KROM construction:
fittedPsi.coordinateProvider = "causal-history";
fittedPsi.coordinateObservableNames = {"xi_1","xi_2"}; // must match n
KoopmanModel model(fittedPsi, systemConfiguration);
auto value = model.evaluateCoordinates(latentValues.data());
```

For direct system-state variables, the convenience factory
`KoopmanModel::fromPolynomialCoefficients(system,dimension,exponents,coefficients,
continuousEigenvalues,center,scale)` constructs the same object in one call.
For analytic observables, set the existing coordinate-observable metadata on
the polynomial result before constructing `KoopmanModel`. No automatic physical
state interpretation is applied to an explicitly declared history provider.

The coefficient matrix is mode-major; the exponent table is feature-major.
Arbitrary unique monomials, constants, and complex coefficients are supported.
Empty center/scale mean raw variables. Shape mismatches, nonfinite values,
nonpositive scales, duplicate monomials and identically zero polynomials are
rejected. Manual specification asserts a candidate, not its accuracy: no
training samples or residuals are fabricated. Use the existing validation
service to test the asserted eigenvalue relation.

Continuous rates are the canonical polynomial storage. Sampled-origin rates
retain dt, origin and logarithm branches; archive round trips and mode
extraction preserve them. Uniform-grid SINDy fits record their sampling dt.
When an older archive has no dt, pass it explicitly. This is the sampling
interval in the trajectory's time units, not an adaptive solver's internal
step or a guessed system parameter.

$$\mu_j=e^{\lambda_j\Delta t},\qquad
\lambda_j=\frac{\log\mu_j+2\pi i\,b_j}{\Delta t}.$$

Sampling alone cannot determine the integers b_j (frequency aliasing). Zero
multipliers remain valid DMDc/EDMDc results but cannot be exported to this
finite-continuous-rate KROM representation; conversion rejects them explicitly.
**Converting the spectrum does not convert the input matrix.** An explicit
continuous control-map identification/conversion is still required before
using an EDMDc model with continuous KLQR.

## CLI and API

```powershell
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" systems
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" generate-ics solve.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" validate solve.json
& "$env:LOCALAPPDATA\Programs\DynSysLab\bin\dynsyslab-cli.exe" run solve.json
```

For example, a Python-free Potts job uses the same schema as an ODE solve:

```json
{
  "schema": "dynsyslab.compute-job/v1",
  "kind": "solve",
  "label": "potts-reference",
  "system": {
    "name": "potts",
    "family": "lattice",
    "preset": "q3_reference",
    "parameters": {"size": 32, "states": 3, "temperature": 1.0}
  },
  "time": {"start": 0, "stop": 100, "output_step": 1},
  "initial_conditions": {
    "mode": "default", "count": 1, "clusters": 1, "seed": 7,
    "values": [], "ranges": []
  },
  "controller": {
    "kind": "none", "settings": {},
    "target": {"kind": "constant", "state": [0, 0]}
  },
  "solver": {"method": "metropolis"},
  "archive": {"enabled": true, "data_role": "free_training"}
}
```

Optional Python callers use thin bindings over the same native services. No
UI-only numerical pathway exists.

Campaign paths use dotted JSON components and may index existing arrays, for
example `controller.settings.frequencies.0`. ODE, DDE, PDE, and lattice parameter
campaigns share this dispatcher. A family still rejects a controller it does
not define; a campaign never changes that scientific capability boundary.
