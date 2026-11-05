# DynSysLab Solver and Post-processor Workbench

## Purpose

The Solver and Post-processor is a native C++ scientific-workflow UI over the
canonical DynSysLab C++ services. It is not a second solver,
post-processing suite, archive, or object model. Its panels construct versioned
JSON jobs and invoke the same compiled API/CLI gates used by scripts and tests.

The application has two independent central panels:

- **Solver** configures systems, physical presets and overrides, time and
  coordinate policy, initial-condition ensembles and clusters, controllers,
  targets, solver options, and cartesian campaigns.
- **Post-processor** selects an existing trajectory or KROM and applies an
  ordered set of spectral, decomposition, regression, estimation, or KROM
  extraction stages without rerunning the dynamics.

The panels share the metadata inspector but intentionally use different left
navigation modes. Post-processor shares Viewer's archive tree. Solver owns a
lightweight proposed-job model and tree whose nodes are Batch, System, Cluster,
and Trajectory. Successful jobs refresh the same catalog file Viewer consumes,
and the global **Open Viewer** action launches Viewer on the shared archive tree
without binding it to the currently selected solve or archive node.
Viewer also opens this workbench directly on the Post-processor tab with its
selected source, so the two applications form one bidirectional workflow.

The workbench and Viewer consume the same native visual-language layer. Shared
UI support owns fonts, colors, DPI conversion, header and toolbar metrics,
standard spacing, owner-drawn controls, ornamental dropdowns and scrollbars,
scientific-number formatting, mathematical labels, and inline value editing.
Workflow panels may compose these primitives differently, but neither
application defines an independent scale or formatting vocabulary.

The workbench can validate and save either active panel as a versioned job file,
then reopen any valid solve or post-processing job for planning or execution.
The file is passed unchanged through the canonical workflow dispatcher; the
desktop app and compiled CLI do not maintain separate persistence formats.

The native queue is exposed by a compact modal workbench rather than a third
scientific panel. A separate `dynsyslab_queue` workflow DLL owns validation,
durable submission, atomic claims, throttled cancellation checks, progress,
stale-worker recovery, and execution. Both the UI and compiled CLI call this
service; the dialog owns only presentation and user intent. This keeps SQLite
and scheduler policy out of event handlers while avoiding any Python worker.

Archive maintenance follows the same boundary. Its native dialog returns a
typed scope and requested action; the asynchronous workbench backend invokes
the operations DLL and renders the versioned JSON report. Read-only audit,
deduplication preview, and confirmed recoverable application therefore share
the compiled CLI's reference-safety rules without embedding filesystem policy
in Win32 callbacks.

System parameters are resolved by the canonical native system factory rather
than duplicated in the UI registry. The inline parameter book presents every
value and unit for the selected preset, records clicked-and-edited values as
explicit overrides, lists the resulting dimensionless groups, and serializes
only changed values into the job. Changing systems clears stale overrides;
changing presets re-resolves the same system while preserving deliberate
overrides.

The System panel presents physical/model parameters, dimensionless groups, and
derived physical scales as divided families inside one vertically scrollable
parameter book. It never creates per-family horizontal scrolling. Compact
mathematical rows retain four significant figures by default, expand through
the **long** option, and make only the resolved physical/model values editable
in place. Every parameter is measured as an independent mathematical display;
its value follows the equals sign directly, and intrinsic-width displays flow
across each family without table columns or reserved alignment slots.
Numeric editing accepts decimal and `e` notation as well as mathematical forms
such as `1 x 10^{-3}`, `1 × 10⁻³`, and `1*10^(-3)`; committing or clicking
outside always returns a valid value to its mathematical display form.
The parameter book uses the shared client-drawn DynSysLab scrollbar rather than
Windows scrollbar chrome. Its black rail, blue thumb, diamond caps, and line
filigree preserve the application's visual language while retaining wheel,
track-page, and drag interaction.

## Presentation and repaint boundary

Custom scientific surfaces use a persistent 32-bit DIB back buffer. A widget
renders its complete new state off screen and publishes only the invalid region
with one `BitBlt`; background erase is suppressed. The workbench additionally
uses Win32 descendant compositing because its editable controls must remain
real child windows. Visibility updates avoid redundant `ShowWindow` calls and
page transitions invalidate only after the new projection has been assembled.
This gives scrolling and System/Control page changes the atomic-frame behavior
of a real-time renderer without replacing native editing, focus, popup, or
accessibility semantics.

`GdiBackBuffer` is a small shared presentation primitive rather than a second UI
framework. `PlotCanvas` already follows the equivalent GDI+ bitmap pattern.
Ornamental dropdown-list faces use that same buffer and subclass the native
ComboBox only for painting; selection, focus, keyboard navigation, and the
popup list remain owned by Windows. A future Direct2D/DirectWrite backend may
replace the rasterizer behind custom surfaces if profiling justifies it, but a
graphics API change alone would not fix flicker caused by sequential child
window composition.

An uncontrolled System page has a small **+ add control** action above its
scientific content and no empty controller section. Adding control creates
peer **System** and **Control** folder tabs and opens Control. A rounded leading
shoulder, diagonally cut trailing edge, omitted active-tab baseline, and shared
body frame make the selected tab visually continuous with the upper editor. System owns
the dynamical-system definition, parameters, and KROM training-domain setup.
Control exclusively owns the controller family, controller-specific
configuration, target (when the family has one), and KROM/control-map/Q/R
policy. Choosing `none` removes the Control page and returns to the single
uncontrolled System page. The tabs change only the editor projection; the same
versioned System draft remains the sole serialization source.

Numerical execution is a separate persistent region below the node-specific
editor. Batch, System, Cluster, and Trajectory selection replaces only the upper
region. Coordinates, integration clock, solver policy, campaign configuration,
and controller evaluation clock remain below it as one job-owned contract. The
lower region is never reloaded or disabled by upper-node selection or by the
System/Control tabs. At the plan/run/save boundary, the workbench broadcasts
that contract into each child `solve` record; this preserves schema compatibility
while keeping UI ownership unambiguous. Uncontrolled Systems ignore the
controller-clock fields.

Preset/default resolution is a public compiled workflow service shared by the
planner and UI. It dispatches ODE, DDE, PDE, and lattice families through their
own scientific modules, so parameter editing and campaign-axis discovery do
not assume every system is an `ODESystem`.

Initial-condition editing follows the same family boundary. Named ODE states,
constant DDE histories, and coordinate-labelled PDE fields share a compact
vector editor; lattice systems use a pattern/custom-site editor whose memory
cost is independent of the number of Win32 controls. Default and random lattice
fields remain reproducibly seed-driven by the native lattice engine.

## Architectural boundaries

```text
Native Solver/Post-processor panels
             |
             | dynsyslab.compute-job/v1 JSON
             v
compiled compute service / native CLI
       |                         |
       v                         v
solver + hierarchy + control    postprocessing pipeline
       |                         |
       +------------+------------+
                    v
             canonical archives
                    |
                    v
       shared Viewer catalog and data tree
          |                       |
          v                       v
       Viewer UI        Post-processor archive tree

Solver proposed-job tree --solve_set--> same compiled compute service
```

Both native executables link the `dynsyslab_ui_support` static library, which
owns the shared `Theme`, `Json`, `Process`, `Model`, `NavigationPane`,
`InspectorPane`, and Help/About modules. The apps therefore compile visual
primitives and catalog semantics once per build while keeping their specialized
workbench controllers independent.

Help/About reports the native support and JSON schema versions plus the exact
project and canonical archive roots resolved by the executable. No Python
launcher or interpreter is involved.

## Scientific hierarchy

The Solver's editable root is one named Job. `SolveJobModel` owns only proposed
hierarchy and per-System job drafts; it contains no numerical solver logic.
`SolveJobPane` owns Win32 tree presentation and versioned serialization.
`SolverPanel` contains two deliberately separated surfaces: a node-specific
single-System editor and a job-specific execution editor. The workbench stores
the active upper form into its System draft before selecting another System,
then reloads only the new System-owned state. Cluster and Trajectory selection
resolves to the owning System without touching the job-level execution state.

This model is deliberately non-archival. Editing and `planJob` make no archive
writes. Explicitly saved job files are kept outside the archive root, while
durable queued definitions live in operational `job_manager` state that the
catalog excludes. Scientific Batch/System/Cluster/Trajectory records are
created only by the numerical execution boundary after trajectories exist.

A normal solve is archived as:

```text
Batch
  System
    Cluster
      Trajectory -> canonical solution archive
```

An ensemble may contain many clusters and trajectories. A control sweep is:

```text
Batch
  System (sweep coordinate 0, controller 0)
    Cluster
      Trajectory -> solution + target + control history
  System (sweep coordinate 1, controller 1)
    ...
```

`solver.solve_hierarchy` accepts one controller per `System`. Batch, System,
and solution metadata retain campaign ID, manifest, sweep coordinates,
controller configuration, initial conditions, physical parameters, and
coordinate-system policy. Each trajectory remains the sole owner of its
numerical arrays; hierarchy archives retain references rather than duplicates.

## Versioned process contracts

- `dynsyslab.compute-registry/v1`: systems, presets, editable parameters,
  observables, targets, KROM archives, controllers, and processing methods.
- `dynsyslab.compute-job/v1`: one `solve`, composite `solve_set`, or
  post-processing request.
- `dynsyslab.compute-result/v1`: a plan or completed archive result.

`compute_workbench.contracts.validate_job` is the common UI/CLI/API gate.
`compute_workbench.execution.plan_job` performs no numerical work.
`compute_workbench.execution.execute_job` is the only workbench execution
entry point. Numerical algorithms remain in their owning packages.

## Coordinate and controller policy

The job chooses `automatic`, `dimensionless`, or `dimensional` coordinates.
Dimensionless UI time and explicit/generated ICs are converted through the
system's declared transforms before entering the solver request; the solver
then applies its established coordinate gate. Controllers are bound by the
solver to the same active coordinates as the dynamics. Waveform and random
pulse settings are specified in physical input/time units and use system
transformations for a dimensionless solve. RPE, Marmottant, and Reddy--Szeri
therefore interpret pressure amplitudes in Pa, frequencies in Hz, and pulse or
active-window durations in seconds, while their archived control histories
remain in the model convention \(p_\mathrm{control}/p_0\). Abstract models
without a declared physical scale use the same fields as model-coordinate
values. Explicit `frequencies`/`model_amplitudes`,
`model_pulse_duration`, and `model_active_time_window` fields bypass
physical conversion. Waveforms accept parallel `frequencies_hz`,
`amplitudes`, and `phases_rad` arrays for multitone excitation.

Every solution records model control units, the scale back to physical input,
the dimensional input units, time scale, availability flags, and whether the
controller settings were physical, model-native, or mixed. Viewer can
therefore label a control history without inferring units from the system
name. KLQR/KMPC load an executable KROM,
validate Q/R dimensions, select an available control representation, and use
the KROM/system coordinate boundary already enforced by `control.koopman`.

The controller family, KROM, Koopman map representation, Q/R weights,
continuous-versus-sampled evaluation policy, control sample time, MPC horizon,
and MPC limits are native structured controls. Their values overwrite any
same-named key in the optional advanced JSON box, preventing contradictory
duplicate settings. This makes continuous adaptive-RHS KLQR and each exact,
constant, or bilinear `B_Psi` path visible rather than implicit.
The KROM picker filters by the selected system's identity, while retaining
unlabeled legacy KROMs for compatibility; canonical RPE targets are suppressed
for non-RPE systems. The execution boundary remains the final authority and
rejects any incompatible legacy selection.

Controller-family details have their own structured dialog. Prescribed
waveforms expose vector amplitudes, frequencies, phases, offset, active window,
and an explicit physical-versus-model coordinate choice; random pulses expose
amplitude, duration, seed, and window. PID gains, LQR linearization state,
common active windows, and optional Koopman source-mode/external-map fields are
also editable without hand-writing JSON. Physical and model keys are mutually
exclusive at serialization, preserving the system-owned nondimensionalization
boundary. Q/R, KROM, map representation, controller clock, and MPC policy stay
visible on the main panel.

The standard H-DMD → SINDy chain likewise has a native structured editor.
H-DMD scalar-observable suggestions are resolved from the selected archive's
system catalog, but the field remains editable for a custom or legacy
observable. Matrix geometry, stride, absolute SVD tolerance, optional rank cap,
and energy sorting are explicit; rank zero retains the full spectrum. SINDy's
observable-coordinate change of variables, reconstruction observables,
polynomial library, sparse-fit policy, exact-RHS/finite-difference derivative,
normalization, conjugate reuse, and training-value materialization policy are
equally visible. The editor merges into, rather than replaces, advanced method
JSON. Numerical execution still belongs only to the compiled post-processing
service, including discrete-to-continuous eigenvalue conversion at the
spectrum-to-eigenfunction boundary.

KROM extraction uses a native spectrum selector rather than requiring manual
mode-index entry. It displays the source KROM's local indices and eigenvalues,
annotates ancestral IDs when present, supports extended selection, and can add
available conjugate partners. The resulting job still enters the single
compiled extraction service, which copies compatible reconstruction, evidence,
and approximate-control blocks into a fully independent archive.

Control-map regression has a system-filtered KROM selector and structured
ridge, selected-mode, conjugate-feature, phase-space-bound, and error-grid
controls. It accepts only a selected archive with a real control history. One
compiled run fits both constant and bilinear representations and archives the
absolute and relative error maps used by Viewer's native heat-map display.

The other decomposition methods also have a typed native surface. DMD, HAVOK,
EDMD, SVD, PCA, POD, bagged PCA, robust PCA, and Hankel-POD retain separate
settings objects rather than sharing a mutable catch-all object. HAVOK's
observable picker is populated from the selected system's scalar-observable
catalog, and EDMD exposes both polynomial and explicit RBF dictionaries.

The observability/KLQE editor resolves candidate KROMs by dynamical-system
identity and reads each selected KROM's `reconstruction_names` field to
populate measurement and state-reconstruction choices. It makes mode
selection, scalar or diagonal process/measurement covariance, detectability,
classical/discrete/continuous observability construction, horizon, sampling,
rank tolerance, and initial estimator state explicit. Numerical work and
archive publication still occur solely in the compiled estimation service.

## Control campaigns

The native contract supports cartesian axes over `system.parameters.*`,
`system.preset`, `controller.settings.*`, `controller.target.*`, `time.*`,
`solver.*`, and `initial_conditions.*`. Array elements use numeric path
components. For example:

```json
{
  "enabled": true,
  "mode": "cartesian",
  "axes": [
    {"path": "controller.settings.frequencies_hz.0", "values": [100000, 150000]},
    {"path": "controller.settings.amplitudes.0", "values": [25000, 50000]}
  ]
}
```

ODE initial conditions have a separate state-aware native editor. It resolves
the selected preset and explicit parameter overrides, chooses the ordinary or
controlled default seed according to the active controller, and edits either
one full vector per requested trajectory or one bounded interval per state.
Random values are never drawn in the UI; the seed and sampling policy remain
inside the canonical native executor. This preserves the rule that equilibrium
is the control default, not a forced restriction, and that free systems retain
their own defaults.

This becomes four point-level system runs, each with its own waveform
controller and clustered trajectory batch. Their shared campaign identity and
manifest let Viewer group them under one physical realization and forcing
campaign without merging their numerical arrays.

The compiled campaign builder edits this contract directly. It discovers every
resolved parameter name for the selected system/preset and suggests ordinary
forcing, target, time, solver, and IC paths. Each axis accepts an arbitrary
workflow path and typed JSON values, so the editor is system-agnostic and does
not embed an RPE-only planning branch. Advanced keys are preserved when the
structured campaign fields are updated.

The same dispatcher supports free parameter campaigns, such as
`system.parameters.temperature` for Ising/Potts or `system.parameters.beta`
for Mackey--Glass. Such trajectories retain the `free` modifier; campaign
membership alone never mislabels them as forced or advertises a nonexistent
control history. Child-engine progress is mapped monotonically into the global
campaign fraction rather than restarting the UI bar at every point.

## Modularity policy

The source-root `AGENTS.md` makes the lightweight-module policy persistent.
Every first-pass compute workbench source is below 500 lines. The hard maximum
for future source files is 1,000 lines; new code should normally split before
500. Public operations require documented API and CLI gates, and UI code must
not import or reproduce numerical algorithms.
