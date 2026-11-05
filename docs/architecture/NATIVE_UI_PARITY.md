# Native UI parity release gate

Native UI parity is a release requirement for the compiled DynSysLab product.
It is not satisfied merely because a native CLI command exists. Every retained
human workflow must have an operable surface in either `DynSysLab.exe` or
`DynSysLabViewer.exe`, and that surface must call the same native service used
by the CLI/API.

The causal history stages (`history_basis`, `history_projection`, and
`echo_history_koopman`) currently use Solver's **Open job…** interface: open a
post-processing JSON job, choose **No** to plan or **Yes** to execute, and inspect
the ordinary progress/output pane. They do not yet have dedicated field editors.
The installed `resources/examples/rpe-echo-history-koopman.json` is a template;
replace its local archive paths with your own. DMDc and EDMDc have their own
configuration dialogs in the Post-processor. The ledger records these distinct
surfaces rather than claiming every method has a specialized dialog.

The machine-readable ledger is
[`NATIVE_UI_PARITY_CONTRACT.json`](NATIVE_UI_PARITY_CONTRACT.json). Its static
gate verifies that every declared capability:

1. is marked `implemented` rather than CLI-only or planned;
2. names a C++ service owner and a C++ desktop surface;
3. points to source and regression evidence that still exists; and
4. contains UI tokens that remain present in the implementation.

The ledger covers the complete retained workflow in three groups:

- Solver: systems and presets, parameters, coordinates, initial-condition
  ensembles, KROM-guided training-domain expansion and safety screening,
  controllers and targets, evaluation clocks, campaigns, resource
  planning, execution progress/cancellation, saved jobs, durable queues,
  archive maintenance, tagged legacy trajectory import, and Viewer handoff.
- Post-processor: decomposition, H-DMD, trajectory eigenvalue optimization,
  SINDy and fitted-boundary eigenfunction construction from independently
  supplied spectra, EDMD/HAVOK/DMD, KROM
  extraction, free-trajectory KROM validation, control-map regression/error
  grids, observability, and KLQE.
- Viewer: completed Batch → System → Cluster → Trajectory archive hierarchy and
  lineage; native Batch, System, and Cluster aggregate views; trajectories and control histories; scalar
  observables; all-component time series; phase and three-state projections;
  native Fourier spectra, spectra/eigenfunctions, eigenvalue-optimization diagnostics,
  KROM validation errors, control-map
  errors, decompositions, estimation, non-ODE data, comparison, axis controls,
  and plot export.

Batch, System, and Cluster nodes are directly plottable as trajectory ensembles
rather than serving as inert folders. Campaign identity and response views are
properties of the completed Batch instead of an extra hierarchy level, and
free/forced status remains a tag. New hierarchy manifests
store solution references relative to themselves, while the Viewer retains a
backward-compatible reader for legacy v1 NPZ hierarchy records and older
project-relative paths. Tentative Solver jobs are never admitted to this tree.
The compact view selector is populated from the complete view list of the
selected object; it has no fixed nine-view ceiling for rich KROMs or estimates.

Free-trajectory KROM validation is a first-class retained workflow. The native
Post-processor evaluates each existing eigenfunction without refitting it and
archives generator, one-step, anchored long-term, and selectable finite-horizon
semigroup defects. Training-domain bounds travel with the KROM, so validation
against later test trajectories records extrapolation coverage explicitly.
Viewer renders the error comparisons and the mode-by-horizon heat map through a
dedicated archive type and links the result back to the evaluated KROM.

Trajectory eigenvalue optimization is likewise end-to-end native. The
Post-processor owns a structured configuration dialog and the Viewer owns a
diagnostic view comparing the initial and optimized unitless spectra, residuals,
and numerical effort. Merely exposing the optimizer through the CLI would fail
the parity ledger because `eigenvalue_optimization` is included in the CLI/UI
coverage equality check.

KROM training-domain expansion is likewise a first-class native Koopman
post-processing entry point, not a specialist Python planner. Its dialog
derives the current domain from archived KROM training states, constructs a
system-agnostic Cartesian expansion with anchors, and appends an exact
explicit-IC Batch to the ordinary Solver. A reusable core service supports
arbitrary state dimension.
The workflow re-derives the grid at validation time and independently screens
every IC before archive publication; source KROM, domain policy, and per-
trajectory screening evidence remain in the resulting native metadata.

Tagged DynBatch CSV import is also a native archive operation with a dedicated
desktop dialog. Preview/validation, optional explicit native-model mapping,
control-input scaling, canonical hierarchy publication, and immutable source
provenance all use the same compiled service as the CLI. A scientist never has
to invoke a Python importer or manually construct archive metadata.

Bounded KROM control-candidate search is a first-class cross-app workflow. The
Solver dialog exposes source selection, target/release semantics, solver and
controller clocks, Q/R tuning, both control-map regressions, screening limits,
success criteria, plan/save/queue/run actions, and a hard 25-attempt ceiling.
The workflow records each attempt durably and Viewer renders the resulting
exact/constant/bilinear comparisons. Raw job JSON and archive editing are not
part of the human path.

The gate deliberately excludes proposed net-new features such as free-form
annotations or three-or-more-object comparison trays. Those remain product
extensions, not missing parity with a retained workflow. Visible Viewer pages
can be copied to the clipboard or saved as PNG or self-contained landscape PDF
without an external reporting runtime.

Run the gate from a development checkout with:

```powershell
python -m unittest tester.tests.test_native_ui_parity_contract -v
```

Python is used only by this development-time source audit. Neither desktop app
nor any scientific workflow imports or launches Python at runtime.

The isolated package smoke gate also starts both compiled desktop entry points
in headless initialization mode. This exercises their DLL loading, runtime-path
discovery, native registry construction, and shared archive-catalog startup from
the installed tree. It catches a package that contains nominal UI executables
but cannot actually initialize without the development checkout.
