# Viewer reports, export and result documents

The native Viewer presents what the pipeline computed. It arranges supplied numbers into figures;
it fits, integrates and estimates nothing. This page covers vector and interactive export,
multi-page reports, experiment reports, JSON result documents, and the smaller presentation
features added with them.

## Export: SVG and interactive HTML

`DynSysLabViewerCli render ... --output plot.svg` writes the selected page as SVG 1.1. The
traces, value scales (linear, log, symlog), focus ranges, heat maps (viridis) and three-state
views (orthographic projection at the figure's camera) match the on-screen canvas.
`--output plot.html` writes every page of the payload into one self-contained HTML file. The page
holds the plotted samples and a short script: hover readout, wheel zoom, drag pan, double-click
reset and page tabs. It opens offline in any browser. PNG output is unchanged (`--output plot.png`).

Notebook packaging (standalone fragments, iframes, executed notebooks) is retired: an HTML export
can be embedded in any notebook or page as it is.

## Multi-page reports

`DynSysLabViewerCli report SPEC.json --output report.html|report.pdf` renders a
`dynsyslab.report/v1` spec:

```json
{"schema":"dynsyslab.report/v1","title":"...","subtitle":"...","footer":"...",
 "pages":[{"title":"...","blocks":[
   {"kind":"paragraph","text":"...","style":"muted|bad"},
   {"kind":"table","title":"...","columns":["name","value"],"rows":[["a","1"]]},
   {"kind":"figure","source":"run.npz","view":"time","panel":"label","views":["krom_dashboard","time"],
    "page":0,"width":1400,"height":800,"coordinates":"dimensional","start":0,"stop":10,
    "scale":"log","title":"...","caption":"..."},
   {"kind":"image","path":"plot.png","caption":"..."},
   {"kind":"svg","markup":"<svg>...</svg>"}]}]}
```

A figure block names an archive (relative paths resolve against the spec) and a route or panel,
exactly as `render` takes them. `views` picks the first route the archive offers. A JSON result
document (below) can be a figure source too. HTML output is one file, with figures as inline SVG
and images embedded, one section per report page, and print page breaks. PDF output is landscape
US Letter, one PDF page per report page, continued onto further pages when a page overflows. Text
is in the standard Helvetica fonts and figures are rendered by the Viewer canvas (ISO 32000-1:2008).
`svg` blocks are HTML only.

## Experiment reports

`dynsyslab-cli experiment report MANIFEST.experiment.json [--output report.html|report.pdf]`
turns an experiment manifest (`dynsyslab.experiment/v1`) into a report spec and renders it.
The report contains the dependency graph of the stages (edges from `${stage.archive}` references),
then one section per stage:

- the graph with that stage highlighted;
- its method settings;
- its results (arrays and shapes, leading eigenvalues, scalar metadata);
- its Viewer plot, from the route preferred for its role.

Failed, reused and not-run stages are marked. The spec stays beside the output as
`<output>.report.json`.

## JSON result documents

The Viewer opens these query results directly. Use `render --source file.json`, a report figure
block, or a catalog entry.

| Document | Pages |
|---|---|
| periodic-analysis `spectral_compare` / `spectral_sweep` result (`"method":"periodic-spectral-comparison"`) | one page per observable, family and packet centre: raw spectral density (symlog) and RMS packet amplitude across orbits; `--view parameter` uses orbit parameters as the abscissa |
| `dynsyslab.phase-cover-geometry/v1` (analytic query `phase_cover`) | one panel per tower level: the sampled orbit coloured by mapped base phase, with its compatible phase points (`display.base_phase_turns`, `display.digits`); with `display.elapsed_base_periods` an animation, one frame per time, the residues carrying as the solenoid flow advances |
| `dynsyslab.research-figure/v1` | the research figures below |

Research figures (`"figure"`):

- **Finite-cycle spectra:** `cycle_spectrum_sequence`, `return_sections`,
  `prediction_comparison`, `regularized_measures`.
- **Cocycles:** `multiplier_branch`, `perturbation`, `truncation_amplitude`, `cyclic_partition`.
- **Coordinate reconstruction:** `coordinate_reconstruction`, `character_mass`, `error_comparison`.
- **Clock changes:** `clock_character`.
- **Odometer and solenoid schematics:** `solenoid_suspension` (the cut-open suspension of the
  adding machine and the Cantor refinement of one fiber), `cycle_clock` (finite quotient clocks of a
  divisibility chain, one animation frame per step, exact integer residues).

The record fields are listed in `viewer/native/src/ResearchFigurePayload.cpp`.

### Phase-cover query

`dynsyslab-cli analytic query.json` with `"operation":"phase_cover"`, a `system`,
`initial_states` (one orbit per tower level), `prime`, `samples` and
`relative_frequency_tolerance` samples each orbit uniformly in model time. It then checks the clock
ratios of the circle p-cover tower, `|omega_j p^j / omega_0 - 1| <= tolerance`. A failed check is
an error, never a drawn cover. The result is a finite sampled covering. It is not a conjugacy or an
infinite-limit certificate. The inverse limit of the tower is the p-adic solenoid (Vietoris 1927;
van Dantzig 1930).

## Pipeline figure data

- **`figure_dataset` operations** (`"operation"` per array):
  - `norm`: Euclidean norm across the selected columns.
  - `difference`: minus another spec, `"minus":{...}`.
  - `diff`: successive differences; with `"per_time":true`, divided by the clock steps.
  - `psi`: K-ROM eigenfunctions along the source's rows, `"krom"`, `"modes"`,
    `"part":"real|imag|abs|phase"`. Rows are physical states or the K-ROM's own coordinates; a
    causal-history K-ROM, ψ(ξ), takes coordinates only, for example ξ lifted from a run or its
    target by `control_target`. `"domain":"state|coordinates"` overrides the choice.
- **`viewer_figure` curve sources:**
  - Each curve may name its own `"source"` archive, `"clock"` and `"interpolation"`.
  - The figure's clock is `"figure.clock":{"source","array"}`.
  - `"figure.fit_y_to":[curve indices]` fits the ordinate to those curves only (a model-scale page);
    `"figure.time_key"` plots against another array (for example per-run values of a score table);
    `"figure.scale"` is linear, log or symlog.
  - Presentation: `"figure.legend":"auto|top"`, `"figure.spectrum":true` with
    `"frequency_scale"`, `"figure.stack_axes":true`, and `"axis":"right"` on a curve, which the
    `right_axis.label` titles (docs/architecture/VIEWER_CLI.md).
  - Every curve is resampled onto that clock (linear by default, NaN outside its span) before
    rendering.
- **`bifurcation_points`:** a peak-amplitude bifurcation diagram.
  - Each run contributes the post-transient local maxima of one component (`tail_fraction`,
    default 0.25), placed at its parameter; a run with no peak in its tail (settled on an
    equilibrium) contributes its final value.
  - The parameter is a dotted metadata path (`"parameter":"system_parameters.driver.amplitude"`),
    or set explicitly per run (`"runs":[{"archive","parameter"}]`).
  - The Viewer shows it as the `bifurcation` view. Campaign batches use the same maxima in
    `campaign_bifurcation`, with the campaign's `bifurcation_tail_fraction`. (Strogatz,
    *Nonlinear Dynamics and Chaos*, 2015, sec. 10.2; Parlitz & Lauterborn, Phys. Rev. A 36, 1428,
    1987.)
- **`observable_spectrum_comparison`:** H-DMD spectra of different scalar observables of the
  same states (Arbabi & Mezic 2017).
  - The sources are hdmd archives fitted with identical settings; differing settings are rejected.
  - The output holds the spectra and the symmetric nearest-neighbour RMS distance matrices, in
    the discrete and continuous domains (Borgefors 1988).
  - Viewer view `observable_spectrum_comparison`: the unit circle with the spectra, the continuous
    spectra, and the distance matrix.

## Presentation details

- **Observability:**
  - The `observability` view adds the Koopman-mode participation heat map. Entry [g, j] is the
    squared participation of nonredundant mode g in observability direction j; each column sums
    to one (Perez-Arriaga, Verghese & Schweppe 1982).
  - For KLQE archives it also adds the observer poles and the state-estimation error norm
    `||x - x_hat||`.
- **Focused spectra:** Fourier panels open on the band holding 99.5% of the amplitude, and at
  least 1.25 times the peak frequency. Resetting the axes shows the full range.
- **System defaults:**
  - Trajectory archives open on the system's recommended view.
  - Lorenz, Rossler, Chua and Hindmarsh-Rose open on the three-state view, with an oblique camera
    (azimuth -46.8, elevation 19.9 degrees).
- **Automatic refresh:** the Viewer watches the archive index generation every 5 s. When new
  results are indexed, it reloads the catalog in the background and reopens the selected archive.
- **Campaign types:** campaign batches show a `campaign type` property: bifurcation, lyapunov,
  trajectory ensemble, parameter sweep or generic.
