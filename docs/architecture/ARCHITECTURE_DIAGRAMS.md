# Architecture diagrams

`tools/architecture_diagrams` draws pipeline and module diagrams from a small JSON spec, so every
agent and author produces the same kind of figure. It needs only the Python standard library.

```powershell
python tools/architecture_diagrams new module-input-map --title "Module: input_map"
python tools/architecture_diagrams check docs/diagrams/module-input-map.diagram.json
python tools/architecture_diagrams render docs/diagrams/module-input-map.diagram.json --png
python tools/architecture_diagrams all --png
python tools/architecture_diagrams icons
```

- `new` scaffolds a spec from a template, either `module` (the default) or `overview`. The
  scaffold renders as-is, so you can edit from a working diagram.
- `check` validates specs and reports the first error with its block or edge.
- `render` writes the SVG next to the spec, and `--png` writes a PNG next to that. `render SPEC`
  can be shortened to `SPEC`.
- `all` re-renders every `docs/diagrams/*.diagram.json`. Run it after changing the renderer or an
  icon.

PNG export uses a headless Edge or Chrome, if one is installed.

To change a diagram, edit its JSON and run `render`: add a block to the `children` of the block
that should contain it, and connect it with an entry in `edges`. Block ids are the handles for
edges and waypoints.

## Spec: `dynsyslab.architecture-diagram/v1`

```json
{"schema": "dynsyslab.architecture-diagram/v1", "title": "...",
 "root": {"id": "pipeline", "kind": "group", "layout": "row", "children": [...]},
 "edges": [{"from": "hdmd", "to": "spectrum"}]}
```

A **block** has `id`, `kind`, `label`, and optionally `icon`, `notes`, `layout`, `gap`, `align`
and `children`.

- `kind` is `module` (a processing stage), `data` (something a module produces or consumes) or
  `group` (invisible layout only).
- Blocks nest. A module inside a module is drawn as a block inside a block, with the container
  shaded one step darker per level. For example, H-DMD and SINDy sit inside Post-Processing.
- `layout` is `row` or `column`, and `align` is `start`, `center` or `end`.
- `notes` attach the settings a job used, as a list of strings or a `{key: value}` object. They
  are shown verbatim in a dashed card.
- `flags` list audit findings, such as a hidden dependency or a mismatch with the intended design.
  They are drawn in a red card under the notes.
- In labels and edge labels, a single character after `_` (`B_ψ`, `ψ_j`) or `_{...}` is a
  subscript. Identifiers such as `krom_extract` stay literal.

An **edge** joins block sides: `from_side` defaults to `right` and `to_side` to `left`. It is
routed orthogonally with rounded corners. To route an edge around other blocks, give `route`
waypoints. Each waypoint may set `x`, `y` or both, as a number or as `block.anchor±offset`, where
the anchor is `left`, `right`, `top`, `bottom`, `cx` or `cy`. Set `style` to `dashed` for an
optional or estimated signal, or `flagged` (red, dashed) for a data path the audit questions, and
use `label` to name the edge.

## Two scales

- **Overview** (`pipeline-overview`): one block per module, with data blocks between them. Keep it
  at this scale.
- **Module** (`module-<name>`): one module's internal steps, in execution order, together with its
  input and output archives and the next module that reads them. Notes cite the source file and
  the job whose settings they repeat.

## Current diagrams (`docs/diagrams/`)

| Diagram | Shows |
|---|---|
| `library-layers` | The native libraries and what each owns, left to right in dependency order: foundation (functions, domains, time series, systems), eigenfunctions, Koopman and observation, workflow jobs, Viewer. |
| `object-model` | Observable, Eigenfunction, Eigenvalue and Spectrum; FunctionExpansion and FunctionInterpolation behind FunctionRepresentation; FunctionDomain with its one optional state map. |
| `pipeline-overview` | The RPE control run on (r, ṙ) eigenfunctions (`examples/native/experiments/rpe-viscous-klqr-multisine/base-r.json`). |
| `module-domain-field` | The `domain_field` and `animation` jobs, the grid-field archive, Viewer and the HTML player. |
| `module-time-series` | TimeSeries operations and the window history built from crops. |
| `module-hdmd`, `module-sindy`, `module-krom-extract`, `module-controller`, `module-estimator` | Those modules' internal steps. |

`docs/DYNSYSLAB_ARCHITECTURE.mmd` is the older Mermaid sketch, kept for the PDF build; these
diagrams supersede it.

## Icons

| Icon | Used for |
|---|---|
| `monitor` | Solver, Post-Processing, or any compute stage |
| `matrix` | H-DMD: a Hankel matrix with constant anti-diagonals |
| `library` | SINDy: a shelf of candidate functions |
| `spectrum` | Koopman spectrum: eigenvalues around a dashed unit circle |
| `heatmap` | Eigenfunction: Re ψ over a phase plane |
| `time_series` | Solver output: trajectories from several initial conditions |
| `krom` | K-ROM: diagonal Λ, input column B_ψ and observable rows C |
| `controller` | Feedback loop: summing junction, gain K, plant and return path |
| `estimator` | Magnifier recovering a smooth state from scattered measurements |

To add an icon, write one function in `icons.py` that draws into a 100 × 100 box and register it
in `ICONS`.

## Conventions

- Keep one spec per diagram under `docs/diagrams/`, named `<name>.diagram.json`, and commit both
  the spec and its SVG.
- Notes should repeat settings from a real job file, never invented values.
- Diagrams describe the suite. Study results, runs and figures stay in the private research
  folders.
