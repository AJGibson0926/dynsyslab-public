# Miscellaneous computes

A standalone computational side-car for testing coordinate changes on trajectories. Supply any modeled or measured real observables, their clock and their rates. The rank study fits an autonomous law for every retained coordinate, runs it from a single initial window, and compares the decoded observable and its derivative against independent reference data.

## Build and run

The [toolkit guide](TOOLKIT.md) includes a runnable, archive-independent example and the study, replay, diagnostic and Viewer commands. The [inventory](TOOLKIT_REVIEW.md) maps the RPE experiments to reusable APIs and records their verification. The [pipeline proposal](PIPELINE_SPLIT.md) recommends ownership and promotion order. [Bundled workflows](.agents/skills)

From the DynSysLab source root on Windows, with a matching native suite build, Eigen and the MSVC developer environment available:

```powershell
$src = "tools/experiments/miscellaneous_computes"
$build = "build/miscellaneous-computes"
$suite = (Resolve-Path build/native-suite-fast).Path
$eigen = (Resolve-Path build/vcpkg_installed/x64-windows/share/eigen3).Path
cmake -S $src -B $build -G Ninja -DCMAKE_BUILD_TYPE=Release "-DDYNSYSLAB_SUITE_BUILD=$suite" "-DEigen3_DIR=$eigen" -DMISC_BUILD_COMPANIONS=ON
cmake --build $build
$env:PATH = "$suite/bin;$env:PATH"
ctest --test-dir $build --output-on-failure
python "$src/run_study.py" study.json results --library "$build/misc_kernels.dll"
```

The compiled commands run independently of Python. The rank-study client uses NumPy and SciPy. On another platform, pass the path to the compiled shared library for that platform. The native kernels are documented in [hankel.hpp](hankel.hpp) and [local_rhs.hpp](local_rhs.hpp).

| Native command | Inputs | Result |
|---|---|---|
| `misc_computes gram signal.npz TAPS out.npz` | One sample-by-channel float matrix named `values` | Uncentered Gram matrix of complete causal windows, in delay-major then channel order |
| `misc_computes fit coords.npz DEGREE THRESHOLD model.npz` | `time`, coordinates `b`, matching rates `db` | Serialized native polynomial SINDy field |
| `misc_computes rollout model.npz request.npz forecast.npz` | Uniform output `time` and one `initial` state vector | Free state forecast on the requested clock |
| `misc_computes local-eval local-field.npz query.npz rates.npz` | Portable local model and sample-by-coordinate `b` | Predicted `db` for each state |
| `misc_computes rk4 model.npz request.npz MAX_STEP out.npz` | Native fitted model, increasing output clock, initial state and optional limits or step budget | Free forecast on recorded output times with a completion status and valid prefix |
| `misc_computes crossings plane.npz out.npz` | Increasing `time` and two columns in `values` | Crossing segment indices, visit times and positions on the supplied polygonal path |
| `misc_computes integral signal.npz DURATION out.npz` | `time` and scalar `values` in physical or model clock units | Current value and its integral over the complete trailing window, timestamped at the window end |

The native fit and integration reuse the existing identification library. The native standalone command creates a new result file. `MISC_BUILD_COMPANIONS=ON` also builds the existing [receiver reconstruction](../receiver_reconstruction/README.md) and [trajectory resolution](../trajectory_resolution/README.md) tools. The former handles retarded receiver signals, neighbor estimates and per-ring errors. The latter checks clocks and bounds path simplification for Viewer.

## Rank study

Use one NPZ per trajectory with `time`, `values` and `rates`. Arrays have shape (samples, observable channels). A one-dimensional array is also accepted for scalar data. Supply rates from measurements, analytic derivatives, or an independently chosen differentiator. Echo-only studies put only the echo in `values`. Joint echo and echo-rate studies put both in `values` and supply the corresponding rates in `rates`. Each trajectory must have a uniform clock. The sampling interval used for projection can divide the training tap interval by any positive integer.

Paths in a study file resolve relative to that file. The following example fits a basis on one training record, fits the field to both training records, selects candidates on a separate validation trajectory and scores one held-out trajectory. A new output directory receives the frozen basis, each selected field, scores and aligned phase curves for Viewer.

```json
{
  "schema": "dynsyslab.miscellaneous-computes/study/v1",
  "training": {"A": "train-a.npz", "B": "train-b.npz"},
  "basis_training": ["A"],
  "validation": {"C": "validation.npz"},
  "test": {"D": "held-out.npz"},
  "basis": {"taps": 64, "gap": 1},
  "channel": 0,
  "ranks": [2, 4, 8, 16],
  "candidates": [
    {"family": "polynomial", "degree": 2, "threshold": 0.004},
    {"family": "localized", "degree": 1, "centers": 20},
    {"family": "paired-gated", "weak_window": 8},
    {"family": "local", "stride": 24, "samples": 160, "neighbors": 12, "decay": 20}
  ],
  "solver": {"rtol": 1e-8, "atol": 1e-10, "max_abs": 400}
}
```

There must be enough independent training history for the requested rank and enough centers for the requested local neighborhood. `basis.locked_basis` accepts a previously saved side-car basis and checks its eigenpairs against the supplied training covariance. Pass matching scales via `basis.scales` if they were set explicitly. In Python, `fit_basis(..., locked=(U, sigma))` accepts an existing decomposition directly. `windows` can map a record ID to a forecast time interval, after preparing its full causal history. `validation_weights` follows the listed validation order. Within each rank the score averages half the relative modal-state error and one quarter each of the relative scalar observable and rate errors. A failed validation integration is excluded and recorded. Missing or zero reference energy has an undefined relative error. Specify a different score in a custom study driver if another scientific objective is desired.

## Python computations

`sparse_stlsq.sparse_stlsq` fits all derivative columns of a sparse feature dictionary from its weighted Gram matrix and feature-target cross products. Batched conjugate gradients and per-output thresholding bound memory without a dense or fill-heavy direct factorization. A failure to reach the requested tolerance is reported.

An optional `controls` object can run the comparisons from the same manifest. For example, `{"neighbor_count":8,"tangent_neighbors":[8,16,32],"restart_time":{"D":1.275},"compare_from":{"D":1.52},"visits":{"A":[2.145,6.182]}}` compares two readouts on measured and freely predicted modes, checks a fresh initialization and tighter integration, and measures the separation at two specified visits. The smooth readout chooses its support size on measured validation windows. Its additional phase curves are saved alongside the main forecast arrays.

Import `miscellaneous_computes` with `tools/experiments` on the module search path.

| API | Use |
|---|---|
| `history.fit_basis`, `HistoryBasis.project` | Compute multi-record, multi-channel Hankel covariance without materializing H. Extend a frozen basis while preserving its prefix exactly. Project values and optionally rates by overlap-add convolution on each admissible sampling phase. |
| `HistoryBasis.project_windows`, `window_encoder` | Evaluate a frozen rank prefix on individual or batched physical history windows. The encoder has axes (oldest-first tap, channel, mode) and returns unit-training-RMS coordinates. For a newest-first delay solver, reverse its tap axis before projection. |
| `HistoryBasis.endpoint`, `raw_scalar_xi` | Decode the newest sample in original units, differentiate that decoder, and convert unit-RMS coordinates to the original scalar-echo convention. |
| `history.separation` | Compare two visits at specified times, component by component and cumulatively. Provide unit-RMS coordinates for comparable component distances. |
| `tapped_history.sample_delay_jet`, `tapped_history.full_tap_svd` | Sample arbitrary physical delays on a uniform source grid using a quintic interpolant with consistent first and second derivatives. Fit and invert all signed H-SVD coordinates for disjoint blocks of multichannel windows. |
| `tapped_history.delay_consistency` | Check whether each delayed pair produced by an autonomous forecast tracks its own present output shifted by the prescribed lag, after enough simulated history exists. |
| `history.scalar_delay_pairs` | Form two scalar taps around each chosen lag. Return the raw history, the equivalent averages and secant slopes, and centered finite-difference curvature labels for training. Each pair has an explicit half-spacing latency. |
| `prehistory.backward_history` | Optionally generate missing history from a known reversible flow, then independently integrate forward to check the original state. Real measured history can instead be prepended directly. |
| `models.fit_field` | Fit normalized, ridge-regularized sequentially thresholded least squares. Families include polynomial, bounded polynomial, localized rational functions, global and paired energy gates, and regularized Fourier harmonics. A positive `weak_window` fits finite-window integrals separately within each record. |
| `paired.PairedField.fit` | Fit amplitude and angular velocity of each modal pair as functions of the full state. Set the winding numbers explicitly or estimate them on a training record. Options include phase locking and mixed instantaneous/integral regression. |
| `local_field.LocalField.fit` | Fit sparse local quadratic vector fields in two-dimensional PCA tangent charts, blend compactly weighted neighbors and optionally damp transverse departures. Tangent dimension two is an assumption of this candidate. Global libraries accept arbitrary state dimensions. |
| `models.free_rollout`, `study.rank_sweep` | Integrate from one state, select fields on validation trajectories and score free forecasts, instantaneous fields and decoding separately. Preserve coverage of partial integrations. |
| `readouts.NeighborReadout` | Predict arbitrary labels independently from a current measured or predicted state. The default uses eight inverse-square-weighted neighbors in scaled asinh coordinates. |
| `readouts.TangentReadout.phase` | Train a scalar output map on observed values and rates. Evaluate the value and its directional derivative using the supplied modal velocity. The local Taylor corrections saturate in the units of the scalar output. |
| `audit.replay_controls` | Repeat a forecast at tighter tolerances and after one later initialization. Compare earlier and restarted forecasts on identical physical times. |
| `viewer.dashboard`, `viewer.render` | Build multi-panel, multi-page plots with independent clocks, trace styles, explicit units and geometry-preserving display vertices. Render with the installed native Viewer CLI. |

For a retained basis of the normalized Hankel matrix H = UΣVᵀ with N training windows, the working state is b = √N Σ⁻¹Uᵀh, where h is the channel-normalized history. On the basis-training windows bᵀb/N = I when windows are rows. The covariance and the orthogonality check use uncentered data. For a scalar channel with training scale s, the old unscaled-window coordinates are ξ = s b/√N. Let c be the endpoint vector returned by the basis. The free phase pair is (cᵀb, cᵀf(b)). A nonlinear readout gives (g(b), Dg(b) f(b)). An independent-label readout of the rate has a different interpretation.

A sampled crossing audit and a successful free trajectory are empirical tests. For analytic equations and observables, see the existing [coordinate closure and global jet tools](../../../docs/architecture/COORDINATE_DYNAMICS_SINDY.md). These distinguish state reconstruction from the existence of a single-valued derivative on each coordinate fiber.

The [Colbrook spectral-method inventory](references/colbrook/README.md) records method summaries, pipeline integration recommendations, and a pinned acquisition manifest for external reference implementations. Its fetch helper prints a plan by default. No external algorithms are installed into the main pipeline by that reference preparation.
