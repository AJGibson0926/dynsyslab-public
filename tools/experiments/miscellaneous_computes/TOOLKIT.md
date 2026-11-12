# Using miscellaneous computes

Use this standalone side-car to fit dynamics in observed coordinates, replay a saved law on a new initial condition, compare delayed sensor readouts, and investigate evidence for closure. The [module inventory](TOOLKIT_REVIEW.md) traces these operations to the experiments. The [proposed pipeline split](PIPELINE_SPLIT.md) describes future ownership.

## Run a complete example

Start in the DynSysLab source root. Build the side-car with companions as described in the [README](README.md#build-and-run). The study client requires NumPy and SciPy. The plotting command uses the installed native Viewer and the trajectory-resolution library.

```powershell
$side = "tools/experiments/miscellaneous_computes"
$build = "build/miscellaneous-computes"
$demo = "$build/my-study"
python "$side/examples/make_demo.py" $demo
python "$side/toolkit.py" study "$demo/study.json" "$demo/results" --library "$build/misc_kernels.dll"
python "$side/toolkit.py" reconstruct "$demo/sensors.json" "$demo/reconstruction"
python "$side/toolkit.py" plot "$demo/plot.json" "$demo/viewer" --library "$build/resolution/trajectory_resolution.dll" --viewer "$env:LOCALAPPDATA/Programs/DynSysLab/bin/DynSysLabViewerCli.exe"
python "$side/toolkit.py" equations "$demo/results/three-taps/model.npz" "$demo/equations.json"
```

Use a new output directory for each run. The fixture supplies three independent records of a signal containing three harmonics, its rate, and state labels. It compares four and six signed coordinates on the same clock. Look at `results/results.json` for selection, coverage and errors, and `viewer/page-1.png` and `viewer/page-2.png` for train and held-out forecasts. Replace the record paths and candidate specifications in the generated JSON to use another system. [Example generator](examples/make_demo.py)

## Build a chart and fit its dynamics

One NPZ per trajectory contains a strictly increasing `time` and samples-by-channels `values` and `rates` of matching shape. For a joint echo and rate study, the two rows at a given instant are $e,e'$ in `values` and $e',e''$ in `rates`. Record the time variable and units, and how the derivative labels were obtained. For example RPE simulation can supply those labels analytically, while measured data require a rate estimation procedure.

An optional `observation_provenance` mapping records how each record's values and rate labels were obtained. For example, give each ID an entry such as `{"rates":"Lie derivatives of the observable evaluated along the RPE solve"}` or describe the differentiator, window and sampling interval used on measured data. The frozen model carries the training and validation entries. All declared entries stay with the study manifest.

Use these fields in a `dynsyslab.miscellaneous-computes/modal-study/v1` manifest.

- Assign separate whole-record IDs and files to `training`, `validation` and `test`. `basis_training` optionally selects a subset of training IDs. Compute the chart, scales and feature parameters on training data. Record-specific `windows` can restrict the scored time interval. All candidate geometries use the same maximum required history for their field training and evaluation.
- Set `geometries` to entries such as `{"id":"six","lags":[0,0.338,0.676],"rank":6}`. These numbers assume that the record clock is in microseconds. Convert lags if the clock uses nondimensional time. Lags begin at zero and increase. Tap order is newest first, with channels adjacent within each tap. For two channels, three taps provide at most six independent coordinates. Omit rank for the full dimension or choose a smaller supported rank. A rank-deficient training block cannot furnish all those coordinates.
- `interpolation: "linear"` supports general channels and irregular clocks. `"quintic_pair"` uses a uniformly sampled value, rate and acceleration to interpolate a derivative-consistent pair. `basis_scales` overrides training RMS channel scaling. `basis_stride`, `fit_stride`, and `speed_refinements` control training sampling. Reference curves retain the requested evaluation samples.
- Set `experiment: "direct"` to learn every component of the vector field on the retained modes. Candidates accept the global feature families `polynomial`, `bounded`, `rational`, `localized`, `gated`, `paired-gated`, and `harmonic`. Other direct candidates are `local`, `compact`, and `oscillator`. For global fields, specify ridge, threshold, and optionally `time_weighted`, `speed_scale` or `weak_durations`. Physical-duration weak fits regress integrated library columns against state increments within each record. Rational terms use finite reciprocal widths, which keep them defined through zero.
- Set `experiment: "structured"` to impose a value/rate pair at each tap on an invertible chart. `second-order` fits acceleration rows, `pair-network` shares maps by the physical separation of each tap pair, and `damped-interior` evolves interior taps from the endpoints with optional learned feedback. Passing `inherited_model` to the last family explicitly reuses an endpoint law. `"mixed"` permits a declared comparison of both classes. A fixed-support `calibration` on a direct sparse candidate optionally tunes its coefficients against short training rollouts. [Field definitions](structured.py) · [Training calibration](calibration.py)
- `modal_weight` sets the share of the validation score given to pooled relative modal error (default 0.5). The other share is the mean of relative current-value and rate errors. `validation_weights` weights trajectories by ID. Only candidates completing all validation windows with defined scores enter selection. All trials, including rejected candidates, are in the summary. The selected law is frozen before test scoring.

For a window matrix with normalized channels and windows in columns, write $H=U\Sigma V^T$. A new column $h$ has signed coordinates $\xi=\Sigma_k^{-1}U_k^Th$. On the training windows these are the rows of $V_k^T$. They retain their signs. The fitter uses $b=\sqrt{N}\xi$ for $N$ basis-training windows. Applying the same linear transform to window-rate labels gives $b'$. A direct SINDy fit regresses $b'$ against a fixed state-only library and converts the fitted field back to the signed convention. Orthogonality on the training windows is a statement about the sample covariance. Closure, state recovery and accuracy of a free forecast are separate tests. [Chart and model](modal.py)

The earlier `run_study.py` entry point truncates or extends one dense, fixed-spacing history basis. Use it for an old echo-only basis with locked signs and normalization. Use `toolkit.py study` to compare different physical delay layouts or full-tap models. [Fixed-basis manifest](README.md#rank-study)

## Replay, inspect and score

Each selected geometry has a self-contained `model.npz` and a `curves/<record-id>.npz` per training or test case. The `modal-model/v1` bundle includes the chart, lags, scaling, field type and coefficients, and any embedded endpoint model. Fitting provenance includes the chosen candidate and record IDs. The loader also recognizes the two historical RPE schema versions described in the [review](TOOLKIT_REVIEW.md).

A replay request contains `time` and exactly one of `initial` (a signed modal vector) and `window` (a taps-by-channels measurement array in original units and saved order). An adaptive solve accepts every supported field. The independent compiled RK4 path accepts a direct polynomial and requires a maximum step. The exact affine path uses a matrix exponential and takes no step or tolerance options.

```powershell
python "$side/toolkit.py" replay model.npz query.npz free.npz --method adaptive
python "$side/toolkit.py" replay model.npz query.npz fixed.npz --method rk4 --step 0.00005 --library "$build/misc_kernels.dll"
python "$side/toolkit.py" score reference.npz free.npz score.json --reference-columns 0 --prediction-columns 0 --time-weighted
```

Replay writes `time`, signed `xi`, decoded `windows` (samples × taps × channels), current channels `values`, completion status, solver message and evaluation count. The Python `replay` function uses `result["values"]` for the entire decoded window. For a scalar-only chart, obtain the rate of a decoded scalar by applying the decoder to the modal vector field. For a joint chart compare this chain-rule rate against the decoded rate channel. [Replay](replay.py) · [Kinematic checks](diagnostics.py)

`equations` exports signed-coordinate affine laws as $A,c$ with $\xi'=A\xi+c$, nonlinear polynomials as per-equation sparse monomials with one exponent per coordinate, and composite/local laws as their evaluator tree, numeric arrays and explicit coordinate transformation. Model time determines the derivative unit. [Equation export](equations.py)

Study curves store `time, xi, predicted, actual, chain_rate` for the available rollout prefix. Use `reference_time, reference` for the entire requested reference, especially for failed or partial held-out runs. `score` accepts separate time keys, value keys and column lists. It interpolates the forecast only within its own time span. Read the physical RMSE, per-channel relative RMS, sample and time coverage, and (with `--prominence`) peak counts, time shifts and relative heights. Phase alignment alone can conceal frequency error.

## Sensor and delay reconstruction

For fixed-support coefficient tuning, [waveform objective helpers](waveform_objective.py) supply quadrature-weighted residuals for physical output channels with global and local RMS normalization, and an SVD-based reduction of training-resolved coefficient perturbations. Give inactive coefficients zero perturbations to preserve their support. Evaluate the fitted coefficients with one unbroken training rollout and a separate held-out rollout.

When a fitted map returns derivatives of the **present observations** from their delay modes, a free forecast needs an initial history over the entire lag span. At each RK stage, assemble an observation window from the evolving present and the simulator's own earlier outputs, apply the frozen chart, evaluate the fitted derivative map, and advance just the present observations. [Delay rollout](delay_rollout.py) accepts a Python or native vector field and the previous measured samples. Its [compiled kernel](delay_rollout_kernel.hpp) supports fractional delays, several observation channels, magnitude guards and partial coverage. The uniform step must be no larger than the shortest positive delay. Refine the step to check interpolation and integration errors. Scoring the same map on fresh measured windows tests instantaneous derivative accuracy.

In a `reconstruction/v1` manifest supply `training` and `test` files with `time, values, state`. Each design is an ordered list of channel/physical-delay/optional-scale entries. For example, compare two physical receivers `[{"channel":0,"delay":0},{"channel":1,"delay":0}]` with one receiver and its delayed copy `[{"channel":0,"delay":0},{"channel":0,"delay":0.676}]` on a microsecond clock. Add another channel for corrected pressure or additional entries for longer histories. All scenarios begin after the longest requested history. [Sensor comparison](reconstruction.py)

The affine and inverse-distance neighbor readouts each estimate the state from fresh measurements at every held-out time. The default neighbor count is eight. Setting `steps_ahead: true` also fits a discrete delay-state advance from within-record snapshot pairs and decodes its free rollout with the affine readout. The forecast has its own clock. To run truncated Hankel-DMD through the Python API, first project each record with one frozen chart, then pass those projected records to `DiscreteDelayModel.fit`. Compare advances at the same physical step. The receiver companion generates retarded physical channels and reports ringwise errors. [Discrete model](delay_dmd.py) · [Receiver companion](../receiver_reconstruction/README.md)

## Crossings and closure

Search the full-resolution numerical path. For example:

```powershell
python "$side/toolkit.py" crossings sampled.npz visits.json --axes 0 1 --library "$build/misc_kernels.dll" --refine-with fine.npz
```

The native search returns isolated intersections of nonadjacent, nonparallel polygon segments and both visit brackets. The command can refine them on a spline through a separate reference run. The Python API also accepts a continuous observation, its velocity, coordinate scales and a higher-dimensional lift. Compare independent runs with `match_visits` to report one-to-one shifts and unmatched witnesses. Tangencies, collinear retracing and off-trajectory collisions require other tests. [Crossing API](crossings.py) · [Witness controls](witness_checks.py)

For $x'=F(x)$ and differentiable $q=H(x)$, a single-valued velocity on the entire image exists precisely when $H(x)=H(y)$ implies $DH(x)F(x)=DH(y)F(y)$ everywhere in the stated domain. The condition gives a set-theoretic law on the image. Smooth extension and uniqueness of its differential equation require additional hypotheses. Exact incompatible derivatives at equal coordinates disprove closure. A planar projection of a larger autonomous state can have crossings. For analytic jet questions, pass equations, observables and their regular domain to the existing native `dynsyslab-cli observable-closure` v2 service. Its report distinguishes order, global certificate, counterexample and unresolved search. [Global closure methods and RPE example](../../../docs/architecture/COORDINATE_DYNAMICS_SINDY.md)

The side-car `closure` command publishes the native report with its original status and certificate fields. Use the matching suite executable. The older installed CLI on the development machine lacks this subcommand.

```powershell
python "$side/toolkit.py" closure examples/native/duffing-position-global-jet.json duffing-proof.json --cli build/native-suite-fast/bin/dynsyslab-cli.exe
python "$side/toolkit.py" closure examples/native/rpe-echo-global-jet.json rpe-proof.json --cli build/native-suite-fast/bin/dynsyslab-cli.exe
```

The Duffing position query certifies minimum jet order one. For the specified autonomous viscous RPE echo on the positive-radius domain, the native query certifies a lower bound of three and a sufficient order of 2,672,672. The latter yields an implicit recurrence on realizable jets. The smallest working order remains unresolved by this query. These derivative-jet orders have a different meaning from the number of delay-SVD coordinates in a fitted model.

## Viewer and workflows


`comparison_pages` composes phase and time plots from independently clocked curves. A plot manifest supplies its reference and fit names, an optional baseline, axes, physical labels and styles. Use `axes: [null,1]` for the time series of channel 1. With `expanded: true` it includes a page on the fit's own scale. The full-resolution curves feed scoring and the native Viewer displays a geometrically bounded vertex selection. [Comparison API](comparison.py) · [Viewer adapter](viewer.py)

The bundled skills provide workflows for [delay identification](.agents/skills/dynsys-delay-identification/SKILL.md), [closure witnesses](.agents/skills/dynsys-closure-witness/SKILL.md) and [sensor reconstruction](.agents/skills/dynsys-sensor-reconstruction/SKILL.md). They live under this side-car's `.agents/skills`. Each can also be opened directly from these links.
