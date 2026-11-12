# DynSysLab native interface system

This file records the approved Solver and Viewer presentation direction. It is
the implementation companion to `docs/reports/ui-refactor-prd-2026-09-14/`.

## Intent

- Human: a scientist configuring numerical experiments or inspecting archived
  dynamics, often with dense parameter and provenance context.
- Primary tasks: configure, validate, and run in Solver; orient, inspect, and
  compare in Viewer.
- Feel: a well-lit computational laboratory—precise, calm, compact, and native
  to Windows rather than terminal-like or dashboard-generic.

## Domain and signature

- Domain: phase portraits, trajectories, spectra, experiments, operators,
  observables, provenance, and control loops.
- Signature: a restrained phase-trajectory line marks active context in tabs,
  selected tree rows, progress, and plot-selection affordances.
- Avoid: a dark terminal shell, equal-weight toolbar actions, generic card
  dashboards, decorative gradients, and color without scientific meaning.

## Tokens

- Lab paper `#F4F8FB`: application canvas.
- Instrument mist `#E8F1F6`: navigation and inset working surfaces.
- White `#FFFFFF`: raised working surfaces and plot continuity.
- Graphite ink `#13212B`: primary text.
- Supporting slate `#526879`: secondary text.
- Hairline blue-gray `#C8D8E2`: structural borders.
- Trace selection `#D6F4F0`: selected scientific context.
- Oscilloscope teal `#0B7D75`: active scientific context.
- Focus cobalt `#0B6FBF`: primary actions and keyboard focus.
- Reference red `#C73E45`: IC/EC, destructive, and error cues.
- Warning amber `#A34F00`: warning state.

All colors are consumed through semantic constants in `Theme.hpp`; screen code
must not introduce a private palette. Scientific plot palettes remain governed
by their existing contracts.

## Type, depth, and spacing

- Display/title: Bahnschrift SemiBold where Windows provides it, with native
  substitution allowed.
- UI text: Segoe UI Variable Text with clear weight and color hierarchy.
- Type scale: 9 pt caption, 10 pt body, 11 pt heading, 19 pt title.
- Depth strategy: surface-color shifts and quiet hairline borders. Shadows are
  reserved for native floating surfaces.
- Spacing: 4 px base unit; 8 px micro, 12–16 px component/pane spacing.
- Controls: 34–36 logical px high; icon actions prefer 40 px and never fall
  below the native 24 px WCAG 2.2 target floor.

## Layout patterns

- One focal task per view; one visually primary action.
- Solver: Job explorer | staged editor | Review rail, culminating in persistent
  Plan / validate and Run & archive actions.
- Viewer: tree-first context rail | selected-object context and scientific
  workspace | compact Details/Lineage drawer.
- At constrained widths, contextual rails collapse before the scientific
  workspace becomes unusable. Primary actions remain reachable.

## Interaction and accessibility

- Every icon-only action has a tooltip and Windows UI Automation name.
- Keyboard focus uses a visible 2–3 px cobalt indicator and is never encoded by
  color alone.
- Native controls keep native keyboard semantics. Dialogs restore focus.
- Dynamic status is exposed through native status/progress controls.
- Verify 100%, 150%, and 200% DPI, keyboard-only workflows, contrast, target
  geometry, and UI Automation Name/Role/Value before release.
