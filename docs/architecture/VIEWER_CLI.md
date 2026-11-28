# Viewer CLI: headless plots from archives and arrays

`DynSysLabViewerCli` renders every Viewer plot without opening a window. Use it for
all result figures: time series, phase portraits, Fourier spectra, control histories,
estimator records and comparisons between runs. It draws with the same `PlotCanvas`,
payload builders, coordinate conversion and value scales as the desktop Viewer, so a
CLI figure and the on-screen plot are the same plot. No plotting script is needed.

There are two entry points:

| Data | Command |
|---|---|
| An archive the Viewer catalogs (solutions, control runs, KROMs, post-processing) | `render` / `inspect` |
| Arrays you saved yourself in an `.npz` (one clock, any curves) | `time-series` |

Both write a PNG and print a `dynsyslab.viewer.render/v1` JSON report (page count,
axis bounds). `render` and `time-series` with `--output plot.svg` write vector graphics, and with
`--output plot.html` a self-contained interactive page of every page. Both exports, like the PNG
and the Viewer window, draw the right ordinate, stacked axes and the legend placement, and the
HTML and window readouts give right-axis curves in their own units. `report SPEC.json` renders multi-page HTML or PDF
reports and JSON result documents render directly
([Viewer reports and export](VIEWER_REPORTS_AND_EXPORT.md)). Errors are `dynsyslab.viewer.error/v1` JSON on stderr with a nonzero exit.
A relative `--output` (including a bare `plot.png`) resolves against the current directory and
missing folders are created; a path that cannot be created is an error with `option` and `path`.

## Archive views

### Find what an archive offers

```powershell
DynSysLabViewerCli panels --source archives\solutions\rpe\RUN.npz
```

Each entry has a `route` (`time`, `phase`, `fourier`, `control_fourier`, `estimator`,
`observable_time`, ...) and a `label` (`r versus time`, `Control history`,
`r Fourier spectrum`, `Phase space`, ...). Select archives with `--source <path>` or
`--id <catalog id>` (`tree --search text` lists ids).

### Select a view

| Option | Meaning |
|---|---|
| `--view route` | every panel of that route (`--view time` gives the control history and each state) |
| `--view "label"` | one panel by its `panels` label, or a route by its catalog name (`"Phase space"`, `"Frequency spectra"`) |
| `--view route --panel "label"` | one panel of a route, when a label exists in several routes |

An unknown view lists the archive's routes and panel labels in the error.

### Routes for trajectory and control archives

| Route | Panels |
|---|---|
| `time` | `Control history` (when the run is forced), one `<state> versus time` per state; targets dashed, the online estimator's state estimate (`estimate`) gold dashed |
| `phase` | `Phase space`, colored by time, with target and estimate |
| `fourier` | one-sided amplitude spectrum of each state (mean removed, Hann window) |
| `control_fourier` | the same for each control input |
| `estimator` | `Measurement and estimator prediction` (the measured output against the estimator's predicted measurement) and `Innovation` |
| `time_frequency`, `three_d`, `observable_time`, `observable_heatmap` | as in the desktop Viewer |

**Control is never filtered or resampled.** When a sampled controller archives the
input it applied at every controller step (`applied_input` on `estimator_time`),
`Control history` and `control_fourier` use it, so chatter above the output grid's
Nyquist frequency is shown. Otherwise they use `control_values` on the output grid.

### Window, coordinates, scales, ranges

| Option | Meaning |
|---|---|
| `--start t --stop t` | time window in the archive's stored clock units (the `time` array, e.g. τ). Applies to every route: time series, phase space, **and spectra**, which are computed over the window only. Arrays on other clocks (e.g. `estimator_time`) are windowed on their own clock. |
| `--coordinates dimensional\|dimensionless` | display coordinates. Dimensional spectra are in Hz when the archive records seconds per time unit (`time_scale_to_seconds`); dimensional RPE control is in Pa. |
| `--scale auto\|linear\|log\|symlog` | ordinate scale; `log` on spectra shows the full dynamic range |
| `--x-range low:high`, `--y-range low:high` | fixed axis ranges in display units. A missing range is fitted to the data; the ordinate is fitted to the samples inside `--x-range`. |
| `--title text` | replaces the title of a single panel, or prefixes every panel title |
| `--width`, `--height`, `--page`, `--figure` | image size, dashboard page, single figure (zero-based) |

### Compare runs

Overlay further archives on the selected one. Repeat the options in order; the
selected archive is the first run.

| Option | Meaning |
|---|---|
| `--compare-source path`, `--compare-id id` | another run (repeatable) |
| `--label text` | one display label per run, selected run first |
| `--color #RRGGBB` | one color per run, selected run first |

Figures at the same position are overlaid when their kinds and axes agree, otherwise
placed side by side. With `--color`, each run's own curves (solution, control,
spectrum) take its color and stay solid; its estimate keeps its dash. Targets identical in several
runs are drawn once in dark grey: solid for the most widely shared target, dashed for
the others, labelled with the runs that share them. Without `--color`, the second and
later runs are dashed.

### Routes for Koopman spectra

| Route | Panels |
|---|---|
| `koopman_spectrum` | a krom's eigenvalues on the unit circle over its pseudospectrum (log10 of the minimum residual, when the fit's `residual_grid` is a rectangular lattice), then each selected mode's eigenfunction; `--figure 0` renders the eigenvalue panel alone |
| `mode_weights` | stems of each observable's modal weight on the mu and lambda planes |
| `spectral_measures` | a krom whose spectrum holds bands or measures: per observable output, the atoms \|c_j\|² as 3-D stems at the multipliers on the unit circle with each unit-circle band's density drawn around the circle, and every band's density against its parameter with the atoms as stems (bands of one observable overlay: a system's exact reference band beside Rigged DMD's); one-point bands (packet centres) are markers |

A krom that only refers to its parent (a child krom) is recognised from the archive index summary
as well as from its full metadata.

### Routes for linear-systems archives

| Route | Panels |
|---|---|
| `bode` | magnitude (dB) and unwrapped phase of every plant channel, with gain and relative degree |
| `nyquist` | `L(jω)` of the archived loops at the plant input (LQR blue, LQG red, LTR sweep gold) and output (Kalman gold, LQG red dashed), mirror image dotted, −1 and the `|1 + L| = 1` circle; `|L| > 10` is not drawn, so zoom with `--x-range`/`--y-range` |
| `loop_bode` | magnitude and phase of the same loops |
| `sensitivity` | `|S| = 1/|1 + L|` and `|T| = |L/(1 + L)|` in dB per loop |
| `hankel_spectrum`, `gramian_spectrum`, `frequency_response`, `pole_zero` | as before; `pole_zero` also marks channel zeros of non-square models |

```powershell
DynSysLabViewerCli render --source archives\derived\linear-systems\RUN.npz --view nyquist --x-range -3:1.5 --y-range -2.5:2.5 --output nyquist.png
DynSysLabViewerCli render --source archives\derived\linear-systems\RUN.npz --view bode --output bode.png
```

### Recipes

Echo-only KLQE→KLQR run (`E`) against full-state feedback with the same K-ROM (`F`) and
a reference run (`R`), DynSysLab colors red, blue and gold:

```powershell
$runs = @('--compare-source', $F, '--compare-source', $R,
  '--label', 'echo-only KLQE→KLQR', '--label', 'full-state, same K-ROM', '--label', 'reference',
  '--color', '#004DB8', '--color', '#D8242A', '--color', '#DE9B00', '--width', '1400', '--height', '760')

# radius and its KLQE estimate, start-up and steady state
DynSysLabViewerCli render --source $E --view "r versus time" --start 0 --stop 30 @runs --output radius_startup.png
DynSysLabViewerCli render --source $E --view "r versus time" --start 85 --stop 100 @runs --output radius_steady.png
# raw applied drive (controller clock, unfiltered)
DynSysLabViewerCli render --source $E --view "Control history" --start 85 --stop 100 @runs --output drive_steady.png
# spectra over the steady window, Hz, log amplitude, first 1 MHz
DynSysLabViewerCli render --source $E --view fourier --panel "r Fourier spectrum" --start 40 --stop 100 --coordinates dimensional --scale log --x-range 0:1e6 @runs --output radius_spectrum.png
DynSysLabViewerCli render --source $E --view control_fourier --start 40 --stop 100 --coordinates dimensional --scale log --x-range 0:1e6 @runs --output drive_spectrum.png
# phase portrait of one run, colored by time, with target and estimate
DynSysLabViewerCli render --source $E --view "Phase space" --output phase.png
DynSysLabViewerCli render --source $E --view "Phase space" --start 40 --stop 100 --output phase_steady.png
# what the estimator measured and predicted
DynSysLabViewerCli render --source $E --view estimator --start 85 --stop 100 --output estimator.png
```

One eigenfunction along training and test trajectories on one phase plane (`observable_phase`,
one `coordinate_transform` eigenfunction series per trajectory). Overlaid runs keep their domain
colouring and grey paths. Brightness uses one modulus scale for all of them, so values compare
across runs. Each run is named at its initial state, and the legend lists one entry per kind of line:

```powershell
DynSysLabViewerCli render --source $PsiTrain --view observable_phase --label training `
  --compare-source $PsiTest1 --label "test (2.1, 0)" --compare-source $PsiTest2 --label "test (1.7, 0)" `
  --output psi_phase_plane.png
```

`inspect` takes the same selection and returns the payload as JSON (`--samples` adds
every plotted array).

## Arrays you saved yourself: `time-series`

One `.npz` with a strictly increasing clock (default key `time`) and one real array per
curve, each N or N × C on that clock. Each component becomes a panel.

```powershell
DynSysLabViewerCli time-series --source run.npz --output comparison.png `
  --curve "actual=RPE solution" --color "#D8242A" --line-width 2.7 `
  --curve "estimate=KLQE estimate" --color "#DE9B00" --line-style time-dash `
  --axis "Radius r" --axis "Radial velocity" --start 0.3 --stop 7.2
```

| Option | Meaning |
|---|---|
| `--curve key[=label]`, `--curve xkey,ykey[=label]` | a time curve, or a phase curve (y against x) |
| `--color`, `--line-width`, `--line-style` | per curve, after its `--curve` (`solid`, `dash`, `time-dash`, `sparse-dash`, `dot`, `dash-dot`, `stem`, `markers`, `open-markers`). A curve label of `_` omits it from the legend |
| `--integer-x-ticks` | ticks and grid at integer abscissas, for an index axis (singular-value index, mode number) |
| `--stem-baseline v` | baseline of `stem` curves: a vertical stem from it to every sample, with a tip marker. Default 0; with `--scale log`, one decade below the smallest positive stem value, shared by the figure |
| `--axis`, `--value-label`, `--x-label`, `--columns` | panel labels and tiling |
| `--time-key`, `--time-label`, `--start`, `--stop` | clock array, abscissa label, window |
| `--right-axis-label`, `--right-axis-scale s`, `--right-axis-offset c` | right ordinate for a curve stored as `c + s * value`; its ticks show `(plotted - c) / s`. In a `viewer_figure` job: `"right_axis":{"label":...,"scale":s,"offset":c}` |
| `--curve-axis left\|right` | per curve, after its `--curve`: a right-axis curve gets its own ordinate, fitted to span the left curves' range and labelled by `--right-axis-label` (do not give a scale or offset). In a job: `"axis":"right"` on the curve |
| `--stack-axes` | with right-axis curves: the left curves in an upper band and the right curves in a band of the same height below, each band with its own ticks and title (for example Re ψ above, Im ψ below). In a job: `"stack_axes":true` |
| `--legend auto\|top` | `auto` (default): inside the plot, in the first corner (top right, top left, bottom right, bottom left) with no data under it; a time or spectrum plot with data in every corner gets the legend as one row between the title and the plot instead. `top`: always that row. Entries keep the curve order, so list curves in the groups they belong to. Also for `render`. In a job: `"legend":"top"` |
| `--spectrum [--frequency-scale k]` | replace each time curve by its one-sided amplitude spectrum over the window (mean removed, Hann; the archive convention). Frequency is cycles per clock unit × k, without the zero bin. Use `--time-label` for the frequency label. With `--x-range`, only the bins in that band are kept, before display thinning, so a zoom keeps full resolution. In a job: `"spectrum":true`, `"frequency_scale":k` |
| `--scale`, `--x-range`, `--y-range`, `--title`, `--font-scale`, `--width`, `--height`, `--page` | as for archives. Plot text defaults to `--font-scale 2.3`; curves without `--line-width` use the Viewer default pens (3.6 pt curves, 3.0 pt targets) |

Arrays on different clocks of one archive can be plotted directly by naming the clock,
for example `--source RUN.npz --time-key estimator_time --curve applied_input=drive`.

## C++ API

`dynsyslab_viewer/ViewerSession.hpp` (target `dynsyslab_viewer_services`):

- `ViewerSession::workspace(ViewRequest)`: route or label, panel, window
  (`TimeWindow`), coordinates, comparisons (`compareIds`, `compareSources`), labels,
  colors, KROM options. Returns a `Payload`.
- `renderPlotPng(payload, output, PlotRenderOptions)`: PNG export, with `xRange` and
  `yRange`.
- `dynsyslab_viewer/RecordWindow.hpp`: `windowRecord(record, window)` slices every
  clock-aligned array of an archive record.
- `dynsyslab_viewer/ComparisonPayload.hpp`: `comparePayloads(entries)` (N runs, colors).
- `dynsyslab_viewer/TimeSeriesOverlay.hpp`: `timeSeriesOverlay` (with `spectrum`).

## Tests

`dynsyslab_native_payload_tests` (record windows on two clocks, controller-clock
control, estimate overlays, estimator panels, colored N-run comparisons with target
merging) and `dynsyslab_time_series_overlay_tests` (spectrum peak frequency and scale).
