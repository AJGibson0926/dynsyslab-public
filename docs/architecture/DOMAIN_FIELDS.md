# Functions on their domains and coordinate maps

A function (an observable or an eigenfunction) lives on a domain: the coordinates q it takes as
input (`FunctionDomain`, see [Spectrum and eigenfunction objects](EIGENFUNCTION_OBJECTS.md#domains)).
This guide covers sampling a function on its domain, moving functions between state spaces through a
coordinate map, the Koopman action of the flow on a function, and how Viewer draws all of these.

## Core objects

| Object | Header | What it is |
|---|---|---|
| `DomainGrid` | `core/DomainGrid.hpp` | A rectangular grid on some axes of a `FunctionDomain`. Coordinates that are not gridded stay at `anchor`, so a grid is also a slice of a larger domain. The first axis varies fastest. |
| `GridField` | `core/DomainGrid.hpp` | Complex values at every grid point (points × outputs) and a status per point: `valid`, `undefined` (evaluation failed or was not finite), `conflict` or `noPreimage`. Push-forwards also keep the fibre `spread` and sample `counts`. |
| `evaluateOnGrid` | `core/DomainGrid.hpp` | The one grid service: an `Observable` (grid in its state space), a `FunctionRepresentation` (grid in its own q) or any point function. Viewer's observable and K-ROM heat maps use it. |
| `dynsyslab.grid-field/v1` | `core/GridFieldArchive.hpp` | The archive of one or more `GridField` frames on one grid, written by the `domain_field` job and read by Viewer. |

## Coordinate maps (`core/CoordinateMaps.hpp`)

An instantaneous coordinate transformation is an `Observable` h: X → Q of the state; the
`coordinate-transform` job applies such a map to time series and pushes rates through Dh. History
(Hankel) coordinates are measurement maps, not maps of the state, so they transform time series
but not domains. On functions a map acts two ways:

- **Pull-back** `pullBack(G, h)` = G ∘ h, with the exact chain-rule Jacobian. It is defined wherever
  h(x) lies in G's domain. A domain's state map is exactly this pull-back.
- **Push-forward** `pushForward(g, h, source, target)` = E[g | h = q], the average of g over the
  source samples that land in each target cell. It is a function of q only where g is constant on
  each fibre h⁻¹(q). The fibre spread max |g − mean| measures how far it fails: above the tolerance
  (relative to max |g|) the cell is a `conflict`, and cells no source sample reaches have
  `noPreimage`. (Conditional expectation; Lasota and Mackey, *Chaos, Fractals, and Noise*, 1994.)
- **Round trip** `roundTripError` = |h^* h_* g − g| on the source grid: zero exactly where g factors
  through h, and large where h loses information about g. This marks where a function pushed to,
  say, (e, e′) cannot be pulled back cleanly to (r, r′).
- **Fold set** `mapVolumeFactor(h, grid)` = √det(Dh Dhᵀ), the local volume factor (area formula;
  Evans and Gariepy, 2015). It vanishes where Dh loses rank, where sheets of h meet.
- **Flow map and Koopman action.** `flowMap(system, t)` is the flow Φᵗ of the unforced system
  (classical RK4, Hairer, Nørsett and Wanner 1993) as a map X → X. Its pull-back is the Koopman
  operator, Uᵗ g = g ∘ Φᵗ (Koopman 1931, *PNAS* 17). `koopmanAction(g, system, grid, dt, steps)`
  samples g ∘ Φ^{k dt} for k = 0…steps; it shows the flow acting on the function rather than on
  the state, and Uᵗ ψ = e^{λt} ψ for an eigenfunction.

The objects stay layered: `FunctionDomain`, `Observable`, `DomainGrid` and the expansions are in the
foundation library, the coordinate maps in the observation library, the job in the workflow and
the drawing in Viewer. Nothing below depends on anything above it.

## The `domain_field` job

```json
{"name": "domain_field", "settings": {
   "function": {"expression": "x_1*cos(x_2)"},
   "domain": {"lower": [-1, -1], "upper": [1, 1], "points": 101},
   "operation": "koopman_action", "time_step": 0.1, "steps": 8}}
```

- `function`: `{"observable": name}`, `{"observables": [...]}`, `{"expression": text or [texts],
  "variables": [...]}` (parameters are constants), `{"library": id, "reference_state": [...]}` or
  `{"krom": path, "modes": [...]}`. A K-ROM whose coordinates cannot be computed from the state is
  sampled on its own q.
- `domain`: bounds, points, `axes` to grid and `anchor` for the rest. The default is the source
  trajectory's bounds with 4% padding, on every coordinate up to three.
- `operation`: `evaluate` (default), `koopman_action`, `push_forward`, `pull_back` (the function is
  then written in the map's output coordinates), `round_trip` or `fold_set`.
- `map` (a coordinate map, given like a function) and `target` (the grid on Q, by default the
  bounds of h over the domain) for the map operations; `tolerance` for conflicts.

## Viewer: plotting a function on its domain

Opening a grid-field archive offers **Function on its domain**:

- **2-D domain.** Heat maps. A complex scalar defaults to domain coloring (hue = phase,
  brightness = magnitude; Wegert, *Visual Complex Functions*, 2012), with a magnitude panel beside
  it. The toolbar switches to real/imaginary or magnitude/phase in place, without re-evaluating the
  grid; the same switch now works on observable phase maps of complex observables.
- **Status.** Points that are not valid are blank in the value panels. A separate status panel
  shows undefined, conflict and no-preimage points, and push-forwards add a fibre-spread panel.
- **Frames.** Each frame, such as each step of the Koopman action, is its own page.
- **3-D domain.** Five slices across the third axis as 2-D panels, plus a 3-D view of the same
  slices as translucent stacked planes (farthest drawn first, invalid points transparent) that the
  3-D camera can orbit.
- **1-D domain.** Curves (real and imaginary parts for complex values).
- **Four or more dimensions.** A grid gridding two or three axes with the rest at `anchor` is a
  slice; the job grids at most three axes by default.

Observable phase maps and K-ROM heat maps evaluate through the same `evaluateOnGrid` service.

## Animations

An animation is a Viewer payload whose pages are frames. The plot toolbar then shows first, back,
play/pause, forward and stop controls and plays the frames on a timer; animations play on opening
by default.

- **Frame series.** A grid-field archive with several frames (the Koopman action) is an animation.
  Its real-valued panels share one colour range across all frames, so colours compare between
  frames. `"animation": {"frame_seconds": 0.25, "loop": true, "autoplay": true}` in the
  `domain_field` settings sets the timing.
- **Composed animations.** The `animation` job strings any views of any archives into one
  animation (`dynsyslab.animation/v1`):

  ```json
  {"kind": "postprocess", "sources": ["${action.archives[0]}", "${solve.archives[0]}"],
   "methods": [{"name": "animation", "settings": {
     "title": "Pendulum", "frame_seconds": 0.25, "columns": 2,
     "panels": [{"view": "domain_field"}, {"view": "phase"}]}}]}
  ```

  Each source is a panel drawn in its view. Frame k shows page k of every panel; a panel has
  either one page (held in every frame) or the full frame count, so panels stay in place.
- **HTML player.** `DynSysLabViewerCli animate --source ARCHIVE --view animation --output
  player.html [--every n] [--width W --height H]` renders every frame with the Viewer canvas into
  one self-contained HTML page with the same controls (first, back, play/pause, forward, stop), a
  frame slider and a speed choice, for chats, reports and browsers.
- **Catalog.** Animation and grid-field archives are classified by their declared `format`, so a
  name such as "koopman-action" no longer makes Viewer open them as K-ROMs.
- **Example.** `examples/native/animations/koopman-action-demos.json` animates the Koopman action
  on a state component for the undamped pendulum, Duffing oscillator and Rayleigh-Plesset bubble
  beside their phase portraits, and pushes the bubble radius forward into the far-field echo
  coordinates (e, e′).
