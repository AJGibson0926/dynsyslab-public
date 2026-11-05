# Native Viewer panel composer

## Contract

Viewer's plot workspace is a session-local composition over immutable
scientific archives. A user can add, remove, reorder, and resize plot panels;
none of those display choices alter archive arrays, metadata, catalog records,
or provenance. Selecting a normal Viewer view again discards the composition
and reconstructs the canonical dashboard.

The composer supports every native panel that the selected archive object can
produce, including:

- control histories and state-component time series;
- state and control Fourier-amplitude spectra;
- two- and three-dimensional state-space plots;
- scalar-observable histories and phase-plane fields;
- KROM spectra/eigenfunction panels, control-map error fields,
  observability diagnostics, lattice views, PDE fields, and campaign summaries.

The **+ Panel** command discovers these choices asynchronously. Its modal
selector first chooses a scientific plot family and then the exact quantity,
such as `u0 Fourier spectrum`, a state-component history, or a phase portrait.
The selected panel is appended without retaining the unselected numerical
payloads.

## Interaction

Every plot card has three compact affordances:

- the six-dot upper-left grip drags the panel into a new ordering;
- the upper-right `×` removes it;
- the filigreed lower-right corner changes its column and row span.

The toolbar below the view/page selectors selects what dragging an empty plot
background does. **Move plots** is selected initially and rearranges panels using
the same ordering as a card-border drag. **Pan curves** moves the data within its
axes by changing the displayed X/Y limits. Exactly one button remains selected;
the choice lasts for the Viewer session, including across view/page changes.
Curves retain their hover/click/drag probe interaction in either mode. Card
borders always rearrange panels, and the `×` close button always uses the pointer
cursor. Resizing, reordering, and closing other panels preserve a panel's axis
limits and pinned probe; no gesture changes the scientific samples.

Dragging the resize corner one or more grid cells creates partial-width,
full-row, partial-height, or full-column panels. The stable first-fit packer
reserves the enlarged rectangle and displaces later panels into the next
available cells. There is no overlap and no absolute pixel position to become
stale after a window resize.

The grid recomputes its column count from panel count and available aspect
ratio, subject to a minimum useful scientific-panel width. A narrow Viewer
collapses toward one column; a wide Viewer uses up to four. Explicit spans are
clamped to the current grid. This gives a responsive dashboard without hiding
the axes in undersized tiles.

## Boundaries

`PanelWorkspace`

- Public, Win32-independent C++ model API.
- Selects adaptive column counts and packs arbitrary row/column spans.
- Provides stable reorder semantics suitable for any future desktop front end.

`PanelCatalog`

- Converts the selected catalog item's supported native routes into deduplicated
  panel choices.
- Separates plot-family/quantity semantics from the dialog and canvas.
- Uses one NPZ decode for all native choices; independent payload generation
  reuses that immutable in-memory record.

`PanelPickerDialog`

- Owns only modal selection and the shared DynSysLab visual language.
- Does not open archives or perform Fourier/observable calculations.

`PlotCanvasInteraction`

- Translates hit regions and mouse gestures into composition operations.
- Keeps reorder/resize state out of scientific payload generation and drawing.

`PlotDragToolbar`

- Provides a compact pair of native toggle buttons using the shared theme and
  DPI metrics, with explanatory text when the available width permits.
- Sets `PlotCanvas::setBackgroundDragMode` to `rearrange_panels` or `pan_curves`.
- Finishes an active gesture on a mode change and preserves panel display state.

`PlotCanvas` and `PlotPanelChrome`

- Draw scientific cards and their compact composition affordances.
- Derive all card geometry from the packer's placements and double-buffer the
  result, so drag/reflow does not leave resize artifacts.

`Backend::loadPanelChoices`

- Performs archive work off the Windows message loop.
- Reuses one decoded NPZ record and tolerates an unavailable optional view
  without invalidating the remaining choices.

These are display APIs over immutable scientific archives. The compiled
`DynSysLabViewerCli` exposes `tree`, `inspect`, `panels`, and `render` for
background navigation and visual inspection using the same services and canvas.
See `VIEWER_CPP_WORKBENCH.md` for commands and versioned JSON contracts.

## Fourier behavior

`control_fourier` is a first-class native payload route. It reads the saved
control `TimeSeries`, requests the same lazy O(N log N) Fourier sidecar as
state spectra, and labels both time frequency and control amplitude from the
archive's coordinate metadata. Control values never receive state-coordinate
multipliers. The route remains read-only and does not persist derived spectra.

## Performance and lifecycle

- Catalog startup does not enumerate or calculate panels.
- The selector is populated only after **+ Panel** is requested.
- Native NPZ choices inflate the selected archive once, not once per view.
- Fourier and observable work stays on the existing asynchronous Viewer task.
- Display traces remain bounded by the native payload policies.
- Closing the dialog releases all unselected panel payloads.
- Loading another canonical view releases the entire custom composition.

The temporary archive focus mask is applied before panel discovery, exactly as
it is for ordinary Viewer navigation. The composer cannot reveal or mutate a
masked archive object.

## Verification

Native tests cover adaptive column selection, full-row and tall-column packing,
collision displacement, stable reorder, route filtering, panel deduplication,
control-spectrum coordinates, and catalog discovery. The complete existing
Viewer/native test suite continues to cover archive immutability, focus masks,
payload semantics, export, and hierarchy behavior.

The toolbar also selects **Dimensionless** or **Dimensional** coordinates through
standard round radio buttons. The drag-mode controls remain push-style toggles. Existing manual ranges and added panels follow that selection when the
archive provides the required conversion. These controls replace the drag-hint
sentence. Phase plots stretch across their panel, with actual-time coloring,
initial-state markers, and known RPE equilibrium references. Axis and legend
mathematics use the same fraction/Greek/script renderer.
