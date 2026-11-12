# Native C++ Viewer Workbench

## Purpose

The native Viewer workbench is the compiled desktop application for traversing
and examining DynSysLab's related scientific objects. It is not a second
numerical framework: it consumes the canonical archives and public native
services and makes their analysis/display contracts quickly accessible from
one persistent interface.

The implementation is Windows-native C++17 using Win32, common controls, and
GDI+. It has no web runtime, plotting-library dependency, Python interpreter,
or Python subprocess. Native archive and payload services read canonical NPZ
objects directly while preserving nondimensionalization metadata,
control-history interpretation, target overlays, Koopman semantics, and
archive compatibility.

Koopman validation archives have a dedicated dashboard: per-mode generator,
one-step, and anchored long-term error comparisons plus a mode-by-horizon
semigroup-error heat map. Catalog lineage points from the validation result to
the evaluated KROM, and training/test role plus extrapolation coverage remain
visible metadata rather than being inferred from the plot.

The KROM dashboard is representation-neutral. Analytic-polynomial SINDy KROMs
and scattered-interpolation H-DMD/boundary KROMs use the same spectrum and
selected-mode eigenfunction workspace. Viewer asks the native KROM service to
evaluate each representation; it does not reimplement polynomial evaluation
or scattered interpolation in the windowing layer.

Viewer uses one draggable vertical divider between a wide context column and
the plotting canvas. The selected-object Inspector, archive tree, and
search/refresh row are stacked in that order in the left column; the Inspector
no longer consumes a separate right-side lane or repeats the tree's selected
object name. The plot canvas therefore retains its former width,
moves right by exactly the added context-column width, and reaches the right
window margin. Viewer and the Solver/Post-processor use the same small splitter
control and constraint model so resizing behavior remains consistent.

Viewer and the Solver/Post-processor also consume one native visual-language
layer. Font families, colors, DPI conversion, header and toolbar metrics,
standard margins and gaps, owner-drawn buttons, ornamental dropdowns and
scrollbars, scientific-number formatting, mathematical labels, and inline
editing patterns belong to shared UI support. An application may arrange those
primitives for its workflow, but must not introduce private pixel scales or an
incompatible formatting convention.

## Information model

The workbench treats a file as a scientific object with identity, type,
available views, properties, and typed relationships. The important graph is:

```text
training trajectory ──trains──> KROM ──contains──> eigenfunctions
       │                           │                    │
       │                           │                    ├─ representation
       │                           │                    ├─ free/test errors
       │                           │                    ├─ trajectory values
       │                           │                    └─ phase-space field
       │                           │
controlled dataset ──regresses──> control map          │
       │                           │                    │
       │                           ├─ exact B_psi       │
       │                           ├─ constant fit      │
       │                           ├─ bilinear fit      │
       │                           └─ error heatmaps    │
       │                                                │
control target ──tested by──> control attempt <─────────┘
                                  │
                                  ├─ actual trajectory
                                  ├─ dashed target trajectory
                                  ├─ independent control history
                                  └─ tracking/control diagnostics
```

Only relationships present in provenance are shown. Missing provenance is not
guessed from filenames. The native catalog resolves direct `source_archive`,
`source_krom`, and controller-KROM links into clickable Inspector records and
generates their reverse links, so navigation works in either direction.
Hierarchy manifests retain their explicit nested nodes; derived archives remain
independently selectable. Viewer adds a derived, non-scientific navigation
hierarchy for the full catalog:

```text
Data
└─ Dynamical system                  Bubble · Rayleigh-Plesset
   └─ Physical realization           Air in water · 20 C
      └─ Major lane                  Training data / Models & assets /
         │                           Control & estimation / Other data
         └─ Study or role            Free / Forced / KROM models /
            │                        Bψ identification / KLQR controls / …
            └─ Source KROM           KROM #1 · H-DMD 2000
               └─ Control target      Target · EC
                  └─ Comparison       Bilinear fit / pulse data / …
                     └─ Record        Attempt #1 / Fit #1 / Training set #1
```

Only levels supported by the record type and provenance appear. For example,
ordinary completed solves retain `Run #n → Parameter set #n → IC set #n →
IC #n`, while multiple Bψ datasets for one KROM share that KROM branch and split
by evidence type. Search campaigns without retained KROM provenance remain
honest about that limitation but group their outputs as `Search attempt #n →
Representation comparison / Error fields`. KLQR attempts made with one KROM
share its branch and then split by target and control-map strategy. If a required KROM relationship is
absent, Viewer says `KROM source not recorded`; it never invents one from a
filename.

Terminal labels are deliberately limited to a short noun and sequence number.
Model methods, forcing modes, evidence types, targets, and similar context live
in their branches; for example, `KROM models → H-DMD 2000 → KROM #1` and
`Periodic forcing → Run #1`. A source-KROM branch may combine its compact
identity and method because it provides context for all of the descendants.

Tree names deliberately exclude solver settings, control implementation detail,
timestamps, and KROM source paths. The original archive label, source path,
properties, and typed relations remain available in the Inspector and are
included in filtering. Presentation groups are not selectable payloads and do
not modify archives, hierarchy manifests, or scientific identities.

Any source record connected by archive provenance to a KROM extraction or fit
is classified under `Training data` and named `Training set #n`. Explicit
forcing metadata separates these records into `Free` and `Forced` subgroups;
legacy role or label evidence is used only when that metadata is absent. KROMs,
spectra, validation results, Bψ fits, and related reduced-model artifacts live
under `Models & assets`. Controller and estimator trials live under `Control &
estimation`. These lanes begin below the system and physical-realization levels;
the existing category and system organization above them is unchanged. Viewer
opens the complete catalog and initially selects a
fully contextualized training set, preferring the newest peer when several are
available.

Under the presentation layers, completed solves preserve exactly the same
structural hierarchy as the Solver's proposed-job tree. Viewer builds it only
from committed scientific archive records:

```text
Run #n
└─ Parameter set #n
   └─ IC set #n
      └─ IC #n
```

The System node owns the exact physical parameters, dimensionless groups,
coordinate convention, preset, and equation fingerprint. Campaign, free,
forced, KLQR-test, and Bψ-identification status remain tags and properties; they
do not insert extra levels into the solve hierarchy. Batch, System, and Cluster
nodes expose bounded aggregate plots, while Trajectory nodes expose the actual
finished solution. The catalog reads both current `dynsyslab.hierarchy/v2`
manifests and legacy `dynsyslab.hierarchy/v1` NPZ manifests without rewriting
the archives. Draft and queued Solver jobs remain excluded.

KROMs, spectra, control-map fits, error fields, and estimates remain
independently selectable records beneath study groups. Their provenance links
still point to the archived solve objects that produced or used them.

Eigenvalue-bearing data is not automatically called a KROM. A KROM must retain
an executable analytic or sampled eigenfunction representation plus KROM
metadata. Eigenvalue-only H-DMD output is cataloged under **Koopman spectra**,
and standalone statistical control fits remain under **Control maps**.

## Workbench layout

The application has four stable regions:

1. **Catalog and filter** — one compact `Data` tree organized by dynamical
   system, physical realization, major lane, study, and only the provenance-backed
   subgroups useful to that study. Scientific leaves use short nouns and stable
   ordinal labels such as `KROM #2`, `Attempt #1`, and `IC #4`; their detailed
   archive names remain searchable. Nonselectable group nodes replace redundant
   role badges. Three icon-only controls above the tree sort only terminal data
   records by date, trajectory count, or total archive size; hover tooltips name
   each control. Structural system, physics, study, KROM, target, and comparison
   branches retain their inherent catalog order. Date is the default and places
   the newest leaves first within each final group. Count and size are likewise
   descending.
2. **View strip** — buttons appear only for views supported by the selected
   object. Additional plot pages use previous/next controls, the mouse wheel over the margins,
   or the page selector on ordinary dashboards. KROM and Observables workspaces
   scroll vertically through all their panels. **Compare**
   offers every other archive object supporting the current view and builds a
   synchronized, page-aligned comparison without eagerly loading the catalog's
   numerical arrays.
3. **Plot workspace** — canonical views begin with up to four native figure
   cards per page. **+ Panel** opens an asynchronous family/quantity selector;
   the six-dot card grip reorders, `×` removes, and the lower-right grip changes
   row/column spans. Custom compositions are session-local and selecting a
   canonical view resets them without touching an archive. The grid adapts
   to panel count and window dimensions. Phase, spectrum, and heat-map axes use
   square drawing regions. Scalar heat maps scale each coordinate independently
   to fill that region, so a physically anisotropic domain remains readable;
   complex-plane spectra retain equal numerical axis units so circles and modal
   geometry remain undistorted. Ordinary trajectory plots omit redundant
   `value versus time`, `Control history`, and `Phase space` card titles, expand
   the axes into the reclaimed space, and use the larger body font for tick and
   axis labels. Specialized titles, including color-linked KROM eigenfunction
   titles, remain visible. **Axes** selects any panel on the visible page
   and applies explicit numeric X/Y bounds; **Use auto** restores bounds derived
   from the archived values. These are display-only overrides and never mutate
   the scientific archive.
   Within each plot, the pointer snaps to the nearest rendered curve and shows a
   transient sample marker. Clicking pins that marker with its interpolated
   `(x, y)` coordinates; dragging the pinned marker moves it along the same
   curve. The plot toolbar has two mutually exclusive drag tools: **Move plots**
   (the default) rearranges panels when dragging a plot background; **Pan curves**
   pans both displayed axes. The selected tool stays active when changing views
   and pages. Plot backgrounds keep a steady move cursor, while `×` uses the
   pointer cursor in either mode. Card borders and headers always rearrange
   panels, and the resize corner always resizes them. Axis limits and pinned
   probes follow their panel when rearranging or removing other panels; switching
   tools does not reset them. Probe overlays are
   session-only and excluded from clipboard/disk exports. Legends are transparent
   and render TeX-style Greek commands, dotted variables, subscripts, and
   superscripts rather than exposing markup as plain text.
4. **Inspector** — tags, compact properties, source provenance, supported
   views, and clickable relations, stacked above the archive tree in the same
   wide left context column. It disappears below the responsive width
   threshold so plots remain usable in a narrow window.

### True three-dimensional plots

Figures whose payload kind is `three_d` use Viewer's native three-dimensional
camera and renderer. They are not flattened with a fixed isometric shift. The
camera preserves the relative scale of all three data coordinates, supports
perspective and orthographic projection, sorts trajectory segments by depth,
and draws a bounded three-axis frame, an orientation triad, time coloring, an
IC marker, and live azimuth/elevation/zoom state. Existing `three_d` payloads
upgrade automatically; `z_label` is optional for legacy payloads and defaults
to `z`.

The plot-local icon gizmos provide rotation, movement, zoom, projection, and
home/reset actions with 40 logical-pixel targets and hover names. Direct input
uses the following equivalent controls:

- left drag rotates; Shift constrains the drag to its dominant rotation axis;
- Ctrl+left drag or middle drag moves the camera frame;
- Alt+left drag or right drag zooms, and the wheel zooms around the pointer;
- arrow keys rotate, Ctrl+arrow keys move, and Alt+arrow keys zoom;
- `+`/`-` zoom, `P` toggles perspective/orthographic projection, and `Home` or
  `0` restores the deterministic camera.

The canvas is keyboard-focusable and retains each panel's camera in its
session state, including workspace/tab restoration and visible-page clipboard
or disk export. Camera operations alter only display state and never the
archived samples. The independent `Plot3DCamera` API owns projection and hit
testing, and `DynSysLabViewerCli verify-interactions --json` is the compiled
headless behavior gate.

Every native plot view exposes **Copy plot** and **Save PNG** controls in the
header. They export the currently visible plot page, including multi-panel
dashboards, target overlays, detail insets, legends, and page-specific content.
Copy places a bitmap on the Windows clipboard; Save writes a PNG through a
normal overwrite-confirming file dialog. Export never modifies an archive or
asks the scientific backend to regenerate data.

Pairwise comparison overlays compatible curve panels and prefixes every legend
entry with its archive-object label; the second object's ordinary curves are
dashed. Non-overlayable scalar fields such as control-map error heat maps are
placed side by side. Named KROM pages stay synchronized by eigenfunction page,
and the resulting comparison can use the same page selector, axis controls,
clipboard export, and PNG export as an ordinary view.

The design borrows the useful *structural* lessons from the Starloom Session
Manager—dense status, docked context, live refresh, and small independent
components—but uses an original scientific-workbench layout and visual system.
Both native applications embed the same multi-resolution DynSysLab
phase-portrait icon. The icon resource and source artwork live in `native_ui`,
while each executable keeps its own stable Windows application identity.

## View coverage

The native workbench uses the same `PlotRequest` routes as Viewer:

| Scientific object | Native views |
|---|---|
| Free/training trajectory | overview, separate state time series, phase space, Fourier spectra, combined time/frequency |
| Controlled trajectory/attempt | all trajectory views, dashed-red target overlay, separately plotted control history |
| ODE trajectory observables | scrolling scalar/vector/complex output time series, each followed by a full-panel phase map for two-state systems; exact state-identity deduplication |
| KROM | fixed discrete spectrum on the right; selected eigenfunction time panels then phase-space maps on the left; real/imaginary, magnitude/phase, or domain coloring |
| Control-map evaluation | native relative-error heat maps, paged when many fits are present |
| Controlled/regression dataset | state histories, phase coverage, control excitation |
| Campaign | navigable manifest hierarchy; plots are requested from member trajectories |
| Observability analysis | singular values and observable-direction heat map |
| KLQE result | observability panels plus estimated/true state time and phase views when retained |
| PCA/POD/SVD family | singular values, variance, component/right-vector heat maps, temporal scores/vectors, robust low-rank/sparse traces |
| PDE result | space-time field heat map and final spatial profile |
| Lattice result | square initial/final configuration maps; 3-D systems use the central slice, alongside ordinary observable histories |

Legacy HTML/PDF artifacts and generic JSON operation records remain available
on disk but are not scientific payloads in the native data tree. Data-tree
selection never calls a Windows file association: every advertised view is
rendered by Viewer itself. This prevents an HTML, JSON, or archive path from
being handed to Visual Studio (or any other unrelated registered application).
An explicitly supplied stale catalog entry containing only the old `open` view
is safely rejected in the plot canvas rather than launched externally.

Observable pages evaluate the selected archive's states through the canonical
compiled system observable catalog. They do not store a second copy of the
observable or use display-only formulas, so RPE echo/energy/surface-tension,
Duffing energy, and future registered observables stay consistent with H-DMD,
SINDy, and control coordinate changes.

Selecting a KROM in the data tree opens its split workspace, including legacy
catalog entries that advertise only a spectrum. The right-hand spectrum keeps
a dashed unit circle behind the eigenvalues and stays fixed while the left panels scroll.
Click selects one retained mode; Ctrl+click adds/removes modes. The wheel zooms
clustered eigenvalues and double-click resets the spectrum. All and None
operate on the retained spectrum. The left side uses the existing draggable
panel layout. All selected training time histories precede their heat maps.
Round radio buttons choose Real / Imaginary, Magnitude / Phase, or Domain
coloring. The first two choices produce two components per eigenfunction.
Domain coloring produces one magnitude-versus-time curve colored by argument,
and one complex phase-space map. Hue encodes argument cyclically; brightness
is sqrt(magnitude / maximum) on the map. Time curves retain full hue brightness
and encode magnitude by their height. Scalar maps use jet. Heat maps have no IC or equilibrium markers.

The public C++ entry point is buildKromWorkspace(path, KromWorkspaceOptions).
ViewerSession::workspace accepts those same options. Background CLI examples:

    DynSysLabViewerCli inspect --id KROM_ID --modes 0,1 --complex-display real-imaginary
    DynSysLabViewerCli render --id KROM_ID --modes all --complex-display domain-coloring --scroll-y 0 --output krom.png

Mode indices are local retained indices; mathematical labels preserve archived
ancestral identities when available. Selection only changes the display.
Permanent extraction remains a Post-processor operation.

Native KoopmanModel evaluators handle polynomial and scattered models. The
native evaluator-only import loadKoopmanEigenfunctionModel preserves explicit
legacy conjugate modes while leaving control imports' independent-mode
realization unchanged. Legacy system parameters are resolved only from
explicit training provenance, validated against retained native states.
Training timestamps are retained or recovered from aligned source trajectories
with their coordinate conversion. Time resets and archived offsets keep
trajectories separate. Missing timestamps use the archived time step with a
warning, or sample indices if unavailable. Two-state systems get a 101×101
map on the same automatic state-phase domain as their training data; higher
dimensional systems retain time histories without inventing a phase slice.
No archive or model is modified by inspecting, selecting, or rendering modes.

Trajectory spectra use the native core `TimeSeries` Fourier sidecar rather than
a Viewer-owned DFT. Catalog construction loads no trajectory arrays. Selecting
**Fourier** or **Time + frequency** materializes only that archive's raw time
and value arrays, applies any evidence-backed display-coordinate conversion,
and requests an O(N log N) transform. Viewer reduces plotted bins only after
the full transform and retains local peaks. Derived frequency, coefficient,
amplitude, power, and phase arrays live only in process memory and are never
published to the primary archive.

## Scientific and UI boundaries

### Native catalog and payload services

Archive discovery is deliberately outside the display loop. Archiver owns
`archives/management/archive_index.sqlite`; every published NPZ has a compact
row containing its archive-relative pointer, byte size, nanosecond modification
stamp, SHA-256 content digest, array names, and a versioned compact metadata
projection. Successful native solve and post-processing publication update just
their new records. The
`dynsyslab-cli archive index` and optional `python -m archiver.cli index`
commands perform an explicit reconciliation for externally copied or legacy
files. Reconciliation uses path/size/mtime as its cheap test and only reopens and
rehashes new or changed records. Missing pointers are removed transactionally.

Viewer and the Solver/Post-processor trust this Archiver-owned database at
startup; opening either application does not scan or reconcile the archive
backend. A generation-keyed presentation snapshot avoids reconstructing an
unchanged catalogue on every launch. Archive publication advances the generation
only when a pointer or its compact metadata actually changes, invalidating that
snapshot. The next consumer then rebuilds it once. An unchanged publication
preserves the snapshot. The `↻` control is the explicit request for the longer
disk reconciliation; if it finds no changes, the retained snapshot remains
valid. If an installation has no central index at all, startup may construct it
once as a recovery path. Numerical arrays are still opened only when a view is
requested.

`Model`, `HierarchyCatalog`, and `ArchiveCatalog`

- Read hierarchy manifests plus NPZ member names and `metadata_json` directly.
- Use the lightweight NPZ index API, which never inflates trajectory, field,
  eigenfunction, or error-map arrays during catalog construction.
- Assign stable path-derived identities, classify objects by archive capability,
  and resolve explicit source/KROM/controller lineage into bidirectional
  clickable relations.
- Preserve the Batch → System → Cluster → Trajectory hierarchy and attach
  derived objects without mutating or coalescing archives.
- Canonicalize system identifiers and common aliases for presentation, so
  `RPE`, `rayleigh-plesset`, and `Rayleigh-Plesset` share one system node while
  their original archive metadata remains unchanged.

`NativePayload`

- Inflates arrays only after the user requests a particular object/view.
- Produces bounded native payloads for trajectories, control histories and
  targets, PDE fields, spectra, control-map errors, observability/KLQE data,
  canonical KROM dashboards, and the complete PCA/POD family.
- Orders combined trajectory displays as control history, state time series,
  then phase space. When no control history exists, state time series still
  precede phase space.
- Fills each phase panel; preserves square spectrum/error geometry and automatic detail insets.
- Converts legacy canonical dimensional RPE/Marmottant `TimeSeries` arrays to
  dimensionless display coordinates when their metadata says the integration
  was dimensionless; native archives explicitly record their stored coordinate
  system, preventing both false labels and double conversion.

### C++ application

`Backend`

- Resolves project-relative sources.
- Calls the native catalog and payload services in-process.
- Caches each object/view payload and invalidates it when its source changes.

`Model` and `Json`

- Parse the two versioned contracts with no third-party JSON dependency.
- Convert UTF-8 metadata to native wide strings at the boundary.

`NavigationPane`

- Leads with the existing category and system organization, then physical
  realizations, major lanes, compact studies, KROMs, targets, comparisons, and
  numbered-record groups.
- Uses icon-only date, trajectory-count, and total-size sort controls with
  explanatory hover tooltips; only terminal records use the selected descending
  order, while conceptual branches retain their inherent order.
- Nests derived KROMs under their explicit provenance parents at arbitrary
  depth; top-level KROM counts therefore describe scientific model lineages.
- Filters compact and original archive labels, sources, properties, systems,
  roles, tags, and descendants without rebuilding scientific data, retaining
  the ancestor path of every match.

`PlotCanvas`

- Draws lines, scatter data, phase-colored complex eigenfunctions, projected
  3-D traces, spectra, and heat maps.
- Preserves dashed target styling, full-panel phase geometry, independent control
  panels, automatic smaller-scale detail insets, and IC/EC labels centered
  immediately beside their reference points (switching sides at the plot edge).
- Pages large reports instead of constructing one large bitmap.
- Supports session-local add/remove, drag reordering, and responsive row/column
  spans through the separate `PanelWorkspace` packing model.

`PanelCatalog` and `PanelWorkspace`

- Discover and deduplicate the selected object's independently renderable plot
  panels without placing archive logic in the dialog.
- Decode a native archive once for an explicit add-panel request, then reuse the
  in-memory record for its time, state-space, Fourier, observable, and other
  available payloads.
- Implement stable first-fit displacement and adaptive grid sizing as a small,
  Win32-independent C++ API with focused tests.

`InspectorPane`

- Draws compact metadata and clickable related-object records.
- Maintains its own scroll state and sends relation navigation to the workbench.

`Workbench`

- Owns layout, selection, asynchronous work, status, and view dispatch.
- Filters obsolete external-only `open` views and routes every selectable tree
  object through the native payload pipeline; tree activation never shells an
  archive or report path through Windows.
- Allows only one asynchronous payload request at a time, avoiding duplicate
  archive inflation while keeping the Windows message loop responsive.
- Opens the selected archive directly in the compiled Post-processor.
- Exposes Help/About with the shared UI/schema versions and the launcher's exact
  project and canonical archive roots.

`native_ui`

- Builds the shared model, JSON, process, theme, navigation, inspector, and
  Help/About sources into one `dynsyslab_ui_support` static library.
- Lets Viewer and the Solver/Post-processor share UI contracts without linking
  their specialized workbench controllers or duplicating compiled support code.

No implementation file in this pass exceeds 500 lines.

## Performance policy

- Catalog creation reads ZIP member names and inflates only `metadata_json`.
  Large archive arrays remain compressed on disk until a view is requested.
- Plot payloads are per-object and per-view, display-bounded, persisted outside
  the repository under `%LOCALAPPDATA%\DynSysLab\Viewer`, and reused until the
  source changes.
- Native rendering draws vector geometry directly with GDI+; it does not create
  intermediate image files.
- Ordinary reports with more than four figures are paginated. KROM and Observables
  use vertical scrolling and render only panels intersecting the viewport.
- Scientific generation remains asynchronous, so the native message loop stays
  interactive while Viewer prepares a payload.

## CLI and application API

Viewer opens the complete catalog by default and does not load a persistent
Mask. It can still apply a temporary, non-destructive archive allowlist when
`--archive-focus <manifest.json>` is supplied explicitly. The
`dynsyslab.viewer.archive-focus/v1` manifest accepts exact project-relative or
absolute archive sources. Required Batch → System → Cluster ancestors remain
visible, while links to masked objects are removed. Missing or explicitly inactive
manifests permit the full catalog. An active manifest with no matching sources
keeps the Mask label and an empty tree; malformed or unreadable existing manifests
report an error without publishing an unfiltered catalog.
This is only a Viewer presentation mask; scientific archives and provenance are
not moved or modified.

An explicitly applied Mask appears in the window title. Startup diagnostics record its
manifest path, label, applied state, and hidden-item count. `--smoke-test` loads
the actual catalog and Mask without opening a window and writes those diagnostics
to the selected cache; it no longer stops after constructing the workbench.

The installed application is `bin/DynSysLabViewer.exe`. It accepts explicit
`--project-root`, `--archive-root`, `--cache`, and optional `--source`
arguments. The main `DynSysLab.exe` launches it with those values, and Viewer
uses the same values when opening Post-processor. A no-argument installed
launch restores the shared roots from `%LOCALAPPDATA%\DynSysLab\runtime.json`;
source-tree builds discover their repository roots directly. The compiled
`dynsyslab-cli.exe` owns automation for workflow planning/execution; Viewer
does not invent a separate scientific CLI.

Plot interaction has a separate, headless test surface so cursor and gesture
contracts do not require Windows UI automation. `PlotInteractionSession` in
`dynsyslab_viewer/PlotInteractionModel.hpp` accepts an arbitrary `Figure` and
`PlotTransform`, then exposes deterministic `move`, `press`, `drag`, `release`,
and state-inspection calls. The independent `dynsyslab_viewer_interaction`
library also exposes curve interpolation, probe restoration, plot panning, and
TeX-style label parsing; it has no GDI rendering dependency.

The installed console companion exercises the same API and emits stable JSON:

```powershell
DynSysLabViewerCli.exe verify-interactions --json
DynSysLabViewerCli.exe parse-math '$\psi_{12}^{*}$'
```

The first command verifies curve-pointer selection, click-to-pin, constrained
probe dragging, background-move selection, two-axis panning, and math-span
classification. It exits nonzero on any contract failure. The second command
returns the parsed normal/subscript/superscript/fraction spans for a supplied legend
label. These checks are also registered with CTest.

`dynsyslab_plot_cursor_tests` additionally drives a hidden native `PlotCanvas`
through `WM_SETCURSOR`/`WM_MOUSEMOVE` sequences and checks the cursor between
messages, including plot backgrounds, curves, card grips, margins, and captured
panning. The canvas consumes client cursor requests so Windows cannot restore
the class arrow between the cursor updates owned by mouse-move handling.
The same test clicks the native `PlotDragToolbar` toggle buttons, verifies their
exclusive selection and background-drag behavior, and checks that close buttons,
border drags, and panel axis limits remain correct in both modes. Passing an
output directory to the test executable also renders compact and wide toolbar
previews for visual inspection.

## Build

```powershell
cmake --preset windows-fast-release
cmake --build --preset windows-fast-release --parallel
```

The executable links the packaged DynSysLab native DLLs and Windows UI
libraries. It never locates or launches Python.

## Verification

The native catalog/payload tests and optional Python compatibility tests verify:

- catalog classification for trajectories, controlled attempts, KROMs,
  control-map errors, POD data, and campaigns, plus exclusion of non-renderable
  reports and generic records from the native scientific tree;
- relationship resolution;
- display-bounded generic state/phase/control payloads;
- catalog-cache reuse;
- physical-parameter isolation, forcing-campaign summaries, mixed-formulation
  grouping, and compact free-ensemble discovery;
- compact system/physics/lane/study presentation, provenance-based training,
  model-asset, control/estimation grouping, free/forced training roles,
  numbered IC/EC terminology, and original-label searchability;
- leaf-only newest/count/size ordering with stable structural branches,
  natural numeric ties, and de-duplicated subtree-size accounting;
- the icon-only sort toolbar, mouse-over tooltips, and default date selection;
- CLI contracts; and
- the sub-500-line module boundary.

The release executable is also visually inspected against real RPE controlled
solution and spectrum archives. This checks native target overlays, separated
control history, full-panel phase geometry, inspector readability, UTF-8 labels,
and asynchronous view switching.

## Remaining native UI extensions

- User-authored relationship annotations do not yet have a catalog editing UI;
  the lineage graph is generated from immutable archive provenance.
- Comparison currently operates on two archive objects at a time; a persistent
  comparison tray for three or more objects remains a future extension.
- Publication export writes the visible page as PNG or a self-contained
  landscape PDF, or copies it as a bitmap. The PDF embeds the rendered page so
  it preserves native plot geometry without requiring Python or a print driver;
  editable vector SVG remains a future extension.
- The immutable archive graph is generated from retained provenance. Editing
  scientific relationships in place is intentionally unsupported.

These constraints preserve the current clean scientific gates while leaving
clear extensions for the next design review.

## Plot presentation and coordinate controls

The workspace toolbar places **Dimensionless** and **Dimensional** radio buttons
beside **Move plots** and **Pan curves**. The coordinate selection transforms
Viewer-owned samples, targets, initial-state/equilibrium references, manual axis
limits, spectral frequencies/amplitudes, and phase time colors together. It does
not modify archives. A mode is enabled only when every panel supplies its current
coordinate identity or an explicit alternate scale; no physical scale is guessed.
RPE state/time scales come from archived radius and time metadata; control pressure
uses the separately recorded input scale, not the equation pressure scale.
Added panels adopt the current workspace coordinates when possible. Comparison
panels overlay only compatible coordinate maps.

Axes and legends share the native mathematical renderer, supporting LaTeX
fractions (including nested fractions), Greek commands, and subscripts/superscripts.
For example `\frac{dr}{d\tau}` is laid out with a fraction rule. Coordinate words
are removed from axis labels; physical unit symbols remain. `DynSysLabViewerCli
parse-math` exposes the same expression tree with numerator/denominator children.

Trajectory views retain all published samples up to 100,000 rather than reducing
RPE records to 2,400. This matches the chat figure's method: connect the saved
samples directly, with antialiased lines and a jet-colored segment collection
using each segment's starting time. Curves use 1.8-point strokes, time-colored
phase curves use 1.25 points, and targets use 1.5 points. Target dashes use four
stroke widths of ink followed by six stroke widths of space, exposing the
underlying trajectory. Legends and detail insets share the same pen styling.
Trace payloads may also supply an opacity in the closed interval from zero to
one; Viewer applies it consistently to ordinary, time-colored, domain-colored,
legend, probe, and 3-D series. The default remains fully opaque for older
payloads.
The native canvas uses a 2× supersampled surface (within the pixel budget),
rounded joins, and high-quality downsampling. No adaptive solver nodes, spline
interpolation, or newly solved dynamics are substituted for the published data.
Nonfinite samples break curves; drawing is clipped to
its axes. Phase portraits fill the available rectangular panel. Their solution
is colored by actual archived time, with one labeled time scale per panel,
including comparisons. Spectra and scalar-field plots keep their scientific
aspect policies.

Initial conditions (IC) use red open circles only on state phase portraits
(including 3D state trajectories). Time-series plots, including controls and
observables, and frequency-domain plots show curves and optional targets without
IC or equilibrium-condition (EC) dots. Known ECs use red filled circles only on
2D state phase portraits; observable plots never show ECs. Heat maps show
neither IC nor EC markers. Coincident IC and EC points on a state phase portrait
show a filled center inside a larger open ring and use the `IC / EC` annotation. The same
visibility policy governs drawing, inspection/export, and automatic axis bounds,
including older payloads that still contain hidden reference records.
Scalar heat maps and time colorbars use the standard 256-color Matplotlib jet
palette. An unknown EC is not inferred from an IC.

`dynsyslab_plot_presentation_tests` checks physical scaling (including pressure),
round trips, native round radio-button styles and exclusivity, manual axes, sample retention, mathematical
fractions, phase time colors, and use of the published sample series. Optional arguments `<output-directory>
<trajectory.npz>` render native dashboard and phase PNGs for visual inspection.
Payload cache names include the presentation revision so old cached plain-text
figures do not bypass the updated builder.

Scroll redraws reuse a 64 MiB least-recently-used cache of plot panels, rendered
with the same 2× supersampling. Only the panel images and interactive overlays
move while scrolling. A low-priority timer prepares nearby offscreen panels
between interactions, one at a time, within the same memory budget. Data, coordinates, scaling, mode selection, panel size,
fonts, and axis changes refresh the affected images. Wheel input respects the
Windows scroll distance and accumulates into a short 16 ms timer animation;
scrollbar dragging and CLI scroll positions remain direct.

The Viewer executable embeds its own spectrum icon from
`native_ui/assets/dynsyslab-app-icon-v6-spectrum.png`, packaged as
`dynsyslab-viewer.ico`; the solver keeps its existing icon. Viewer retains its
separate `DynSysLab.Viewer` taskbar application identity.

## Background navigation and plot inspection

Add `--benchmark-scroll 60` to a CLI `render` call to report cold-frame time,
scroll redraw median/p95, panel render count, cache hits, and cached bytes.

Use the installed `DynSysLabViewerCli` for automated Viewer work; do not launch
or manipulate desktop windows for verification. The CLI shares `Backend`, the
optional explicit Mask, payload builders, coordinate conversion, and `PlotCanvas` with
Workbench. Rendering uses an invisible, non-activating surface and writes a PNG;
it never opens Viewer, a console window, a file dialog, or the clipboard.

```powershell
DynSysLabViewerCli tree --search rpe-reference
DynSysLabViewerCli tree --parent <id>
DynSysLabViewerCli inspect --id <id> --view phase
DynSysLabViewerCli panels --id <id>
DynSysLabViewerCli render --id <id> --view phase --coordinates dimensionless --output phase.png --width 1100 --height 580
```

Use `--source <archive-path>` instead of `--id` for an exact archive selection.
Selection is restricted to the visible, masked data tree. `tree` returns stable
IDs, parent IDs, category/system grouping, views, relations, the Mask name, and
visible/hidden counts. `--parent` returns immediate children; `--search` filters
records. `inspect` returns mathematical labels, sample counts and bounds,
references, and coordinate maps; `--samples` includes all plotted arrays in the
round-trippable `dynsyslab.viewer.payload/v1` contract. JSON is printed to stdout
or written with `--output`. Errors are versioned JSON on stderr with nonzero exit.

`render --figure N` fills the workspace with one figure; `--page N` selects a
normal dashboard page. Both use zero-based indices. PNG exports include the
coordinate/drag toolbar and report actual axis bounds and page count in JSON.
The full inspected payload contains every figure, regardless of page.

Runtime overrides: `--project`, `--archives`, `--cache`, `--catalog`, `--focus`.
The full catalog is the default; `--focus` is an explicit temporary override.
No command reconciles or writes the scientific archive arrays. The public C++ API is `ViewerSession`, `inspectPlot`, and
`renderPlotPng` in `dynsyslab_viewer/ViewerSession.hpp`, linked through the
`dynsyslab_viewer_services` CMake target.

`dynsyslab_viewer_cli_background_contract` tests tree traversal, explicit-Mask
parent retention and hidden-item rejection, source selection, coordinate transforms,
payload export/reload, PNG output, and invalid selection handling without showing
windows. Presentation/interaction tests also use hidden test surfaces.

For Windows automation, the launcher must also remain hidden. Use direct process
execution with `shell: false`, `windowsHide: true`, and captured stdout/stderr;
do not start a visible PowerShell host just to hide its child. Supply the normal
Windows profile paths (`LOCALAPPDATA`, `APPDATA`, `USERPROFILE`) if the automation
runtime provides a reduced environment. File inspection can use filesystem APIs
without launching a shell.

`DynSysLabViewerCli windows` reads the actual running Viewer title, executable,
process ID, and tree-node count without activation. `reload <process-id>` queues
that Viewer's normal catalog refresh without restoring or activating its window;
the JSON response acknowledges the request, not its completion. An unresponsive
tree reports a null node count rather than a false zero.

`DynSysLabViewer --background-check --background-refresh-count 3 --cache <test-cache>`
creates the complete desktop workbench without ever showing its top-level window,
loads its catalog and initial plot, exercises repeated refreshes, and paints its
inspector and plot into memory. It writes `workbench-check.json` and exits with a
nonzero status on failure. `--background-refresh-count` is optional (0–20).
The same project/archive/catalog/Mask overrides as ordinary Viewer startup apply.
Unlike `--smoke-test`, this checks the actual tree and desktop event loop. Refresh
uses the Archiver's normal index reconciliation; use a fixture archive for tests.

An explicitly supplied active Mask with no matching records displays an empty
tree with its name. Malformed or unreadable existing Mask files report an error
and do not publish an unfiltered catalog. A missing preference on first launch,
or an explicitly inactive Mask, permits the full tree. Catalog replacement clears
the inspector's borrowed item pointer before freeing old records; tree deletion
can synchronously repaint sibling controls. The Mask title and tree diagnostics
are updated before initial plot loading. The desktop regression test reproduces
that repaint during tree replacement while every window remains hidden.

### Overlay time series from saved arrays

Use the native CLI directly to plot several estimates against a reference in a chat or report. A single NPZ supplies a strictly increasing time vector and one real array per curve. Each curve has shape N or N × C on that clock. Viewer produces a labeled panel for each component using its standard PlotCanvas renderer.

```powershell
DynSysLabViewerCli time-series --source run.npz --output comparison.png --curve "corrected=Monopole plus near field" --color "#DE9B00" --line-width 4 --curve "actual=RPE solution" --color "#D8242A" --line-width 2.7 --line-style solid --curve "leading=Leading monopole" --color "#004DB8" --line-width 1.65 --line-style time-dash --axis "Radius r" --axis "Radial velocity" --start 0.3 --stop 7.2
```

Repeat `--curve` in drawing order. Each curve can take a hex `--color` and a `--line-style` of `solid`, `dash` or `time-dash`. The last spaces dashes horizontally so overlapping estimates remain distinguishable even during fast oscillations. Use `--axis` once per component or accept the default labels. Optional arguments include `--time-key`, `--time-label`, `--title`, `--width`, `--height` and zero-based `--page` for more than two components. The command needs only the NPZ and output destination. For dense PNGs it retains endpoints and both extrema in chronological order within each subpixel time bin. The public C++ function `dsl::timeSeriesOverlay` returns a Viewer payload and keeps every sample by default. Pass the resulting payload to `dsl::renderPlotPng` for a background export.

## Main and Observables workspace tabs

The plot workspace begins with **Main** and **+ observables**, using the same
folder-tab drawing component as the Solver/Post-processor. Clicking **+ observables**
opens **Observables**. Main retains the current canonical plots and custom panel
composition. Switching tabs preserves each workspace's payload, page or scroll position, manual axes,
and pinned probes; the coordinate and drag controls continue to operate on the
active workspace. Selecting another archive rebuilds the workspace from that
selection, and comparisons retain both selected sources for observable evaluation.

Observables shows every registered system observable. For a two-state dynamical
system, each output's time series is immediately followed by its phase-space heat
map. The complete list scrolls vertically, with no page selector or left/right buttons. Higher-dimensional systems retain time
series only. Vector outputs receive one panel per component; complex
outputs receive real and imaginary panels. The catalog's public evaluator receives
states in the system's native coordinates, including conversion of historical
physical RPE archives. Time uses the same archived samples and scaling as Main,
retaining up to 100,000 samples per trajectory. Published target and true-state
overlays are evaluated too. Known output units come from the observable descriptor;
no scaling for an arbitrary derived observable is invented. Hierarchy selections
reuse Main's trajectory discovery and 96-trajectory limit.

Heat maps evaluate the same public observable catalog on a 101 × 101 grid over
Main's phase-space bounds, including manual axes. Returning from Main regenerates
the maps when those bounds or coordinates changed. The shared
`automaticAxisLimits` function supplies both canvas limits and the map domain;
heat maps add no second margin. Scalar fields use jet and a labeled value
range. State-observable values and grid axes follow the coordinate selector;
complex real/imaginary parts and vector components each receive their own map.
Compared/grouped sources use a common domain; fields with the same catalog identity,
system/equations/parameters, output index, component and coordinates are deduplicated,
and distinct parameter-dependent fields remain adjacent to their time series.

The shared native API is `Backend::loadObservables(items, optional<PhaseSpaceDomain>)`
and `ViewerSession::workspace`. The `observable_time` route remains time-only,
`observable_heatmap` supplies fields, and `observable_workspace` supplies pairs.
Background inspection and rendering select the same workspace without opening UI:

```powershell
DynSysLabViewerCli inspect --id <id> --workspace observables
DynSysLabViewerCli render --id <id> --workspace observables --scroll-y 0 --output observables.png
DynSysLabViewerCli inspect --id <id> --workspace observables --compare-id <second-id>
```

`--workspace main` is the default. Plot inspection reports the jet palette
and red open/filled reference-marker styles. `--background-workspace-check`
exercises the actual Main/Observables controls, asynchronous evaluation, cached
return, and Main restoration while the complete workbench remains hidden.
The observable regression gate compares every real/vector/complex output and
target overlay against the public system API, checks sample and coordinate
preservation, heat-map evaluator values, pair ordering, exact automatic/manual
phase bounds, and marker pixels and shared palette values. Plot inspection includes
axis limits and scalar value ranges.

## White panels, continuous scrolling, and value scales

Plot panels and insets use white backgrounds with dark mathematical labels, ticks,
grid lines and legends. The surrounding workbench stays dark. Ordinary curve colors
are adjusted for contrast; the exact jet palette and complex domain coloring are retained.
Drawing lives in PlotFigureDrawing; canvas layout and scrolling are separate modules.

KROMs have one **KROM** workspace tab. All selected eigenfunction time series, then
all phase-space maps, are arranged in a vertically scrolling grid on the left. The
spectrum on the right remains fixed. Hover inside a plot and roll the wheel up
to zoom in or down to zoom out around the pointer. Double-click the plot to
reset its view. Roll the wheel over a panel margin or use the vertical scrollbar
to scroll between panels. The ordinary **Observables**
tab uses the same scrolling layout, keeping each time series immediately before its
map. Offscreen panels are culled from painting, not removed from the payload.

The native Observable API exposes stateComponentIndex(output), populated exclusively
by stateComponent and identity factories and propagated through stack with exact indices.
Viewer removes a vector component only when a scalar observable has the same proven
projection, coordinate identity and units. Custom evaluators with matching names or
coincident values remain distinct. Figure observable_identity records the system equation
fingerprint, catalog function, output index and real/imaginary part. Composition never
uses curve similarity or sampled field equality as a function-identity test.

Round **Auto scale**, **Linear**, **Log**, and **Symmetric log** radios control values:
time-series ordinates and scalar heatmap color normalization. Time, state phase-space
coordinates, the spectrum, and complex domain-color hue retain their definitions.
Auto uses the 5th/95th percentiles of nonzero magnitudes and selects a nonlinear scale
when that central range spans at least three decades. Positive data uses log; data
including zero or negative values uses symmetric log. An explicit Log request also
falls back to symmetric log if necessary to preserve nonpositive samples. A zero-only
or constant series remains linear in Auto. Symmetric log is linear inside a data-relative
threshold and logarithmic outside it. Thresholds and decisions transform with units.
Each nonlinear axis/colorbar identifies its effective scale. Inspection reports the
requested scale, effective kind, target, threshold, and reason. Samples are never rewritten.
Axis padding, curves, probing, dragging and pan transforms share the native PlotScale API.

Public C++ calls are setValueScale(payload, requested), resolveValueScale(figure), and
PlotRenderOptions::scrollY. CLI inspect/render accept --scale auto|linear|log|symlog;
render accepts --scroll-y pixels and reports scroll_y, scroll_maximum and visible_figures.
--figure still renders one full-size plot; --page remains available for ordinary dashboards.
The scale/scroll test covers wide positive and signed data, constants/outliers, unit
invariance, exact identity/index propagation, raw-sample preservation, logarithmic probes
and pans, white pixels, round radios and scroll restoration. The hidden workbench gates
verify absent paging controls, the single KROM tab, and access to the final plot.
