# Identification options (native port, phase D)

The identification steps of the job pipeline and the settings phase D added to them. Every fitting
step publishes a `krom` ([DATA_OBJECTS.md](DATA_OBJECTS.md)); a step that adds to a krom publishes
its child.

## Pooled trajectories

- **H-DMD** (`hdmd`, `"additional_sources": [archives...]`, `workflow/native/src/HdmdPooled.cpp`):
  one Hankel DMD over the source and the listed trajectories of the same system and settings. Each
  trajectory contributes its own window pairs; the pooled matrices put their columns side by side,
  so no window and no shifted pair crosses a boundary (Arbabi & Mezić, *SIAM J. Appl. Dyn. Syst.* 16,
  2017; ensemble DMD: Tu et al., *J. Comput. Dyn.* 1, 2014). The krom's samples are every window of
  every trajectory (`training_sample_indices` index the input records concatenated in origin order),
  its interpolation is on the physical states there, and `trajectory_columns` counts each
  trajectory's windows. Modes keep DMD order unless `mode_sort` asks for the terminal-amplitude sort.
- **SINDy** (`sindy`, `"additional_sources"`, `workflow/native/src/SindyTraining.cpp`): the rows of
  all trajectories, each with the same `sample_stride` / `time_window`; finite differences stay
  within a trajectory (`KoopmanSINDyOptions::segments`). Eigenvalues from data need one trajectory.

## SINDy on a state-space domain

`"domain": {"lower": [...], "upper": [...], "density": n, "time": t}` (a regular grid, first
coordinate slowest, at most 250,000 states) or `"domain_points": [[...], ...]` trains on states with
the system's rates there; there is no time series, so `gamma_method` must be `system_rhs`. The
generator equation needs only the rate at each state (Kaiser, Kutz & Brunton, *Mach. Learn.: Sci.
Technol.* 2, 035023, 2021).

## SINDy fit options (`core/SINDy.hpp` `KoopmanSINDyOptions`)

| Setting | Effect |
|---|---|
| `rescale` | positive multipliers after normalization, q' = (q − c)/s · m (one value or one per component); archived as the expansion's scale |
| `sparsification_policy` | `iterative-relative` (default: sequential thresholding relative to the largest entry of the column-normalized generator) or `single-absolute` (one SVD of the unnormalized generator, entries below `threshold` zeroed, no refit) |
| `directional_derivative` | `analytic` (default, the dictionary gradient) or `finite-difference`: (Θ(q + h q̇) − Θ(q − h q̇))/2h with h = `directional_step` · max(1, maxᵢ\|qᵢ\|) |
| `dictionary` | fit in another dictionary: a list of terms (as in `eigenfunction_expansion`), an archived dictionary, or `{"kind": "gaussian", "centers": [[...]], "width": w, "include_constant": true}` (`workflow/native/src/DictionaryTerms.hpp`) |
| `readout_transfer` | `true` or `{"tolerance": 0.2}`: given a krom, keep its readouts C when every calibrated mode agrees with its reference within the tolerance; otherwise refit them. Recorded as `metadata.readout_transfer` |

## Vector-field SINDy (`sindy_dynamics`, `core/SindyDynamics.hpp`)

| Setting | Effect |
|---|---|
| `ridge` | Tikhonov ridge on the normalized problem (Hoerl & Kennard, *Technometrics* 12, 1970) |
| `sample_weights` | `uniform`, `trapezoid` (each sample's share of time) or a list, one per fit sample |
| `speed_scale` | divides each weight by 1 + \|q̇\|/s |
| `derivatives` | `finite_difference` (fourth order), `central_high_order` (eighth order, Fornberg, *Math. Comp.* 51, 1988) or `stored` |
| `monomials` | library of monomials with integer, possibly negative (Laurent) exponents, one row per term |
| `weak_durations`, `weak_stride` | weak (integral) form: rows q(t + D) − q(t) = ∫Θ dt ξ over windows of physical length D starting every stride samples, Θ piecewise linear (Schaeffer & McCalla, *Phys. Rev. E* 96, 2017; Messenger & Bortz, *Multiscale Model. Simul.* 19, 2021) |

## Delay coordinates at physical-time lags

`coordinate_transform` kind `delay_taps` (`core/DelaySampling.hpp`): coordinates
y_c(t − τ_l) at physical lags `lags` (or fractional `lag_samples`), linear between samples, causal
(no tap reads after its endpoint), for every sample with a complete window. With `"jet": [value,
rate, acceleration columns]` each tap is (y, y′, y″) from the C² quintic Hermite interpolant. Real-
valued delays: Sauer, Yorke & Casdagli, *J. Stat. Phys.* 65, 1991.

## Delay diagnostics (`core/DelayDiagnostics.hpp`, `core/TangentReadout.hpp`)

- **`delay_pairs`**: per lag, the taps p[e − l − 2h], p[e − l] available at each end e, their midpoint,
  secant slope and centred second difference, timed at e − h with latency h dt (`channel`, `lags` in
  samples, `half_width`, `stride`). Delay coordinates: Takens, *Lect. Notes Math.* 898, 1981.
- **`delay_consistency`**: a forecast of taps, (value, rate) per tap from lag 0, against the cubic
  Hermite shift of its current tap (Fritsch & Carlson, *SIAM J. Numer. Anal.* 17, 1980); errors are
  relative to `reference` (default the source).
- **`extendLockedBasis`**: extends a frozen prefix (U, s) of a window Gram matrix, kept bit-identical,
  with the leading eigenvectors of the complement-projected Gram (Broomhead & King, *Physica D* 20,
  1986), after checking that the prefix is orthonormal and still an eigenpair.
- **`historyEndpointDecoder`**: the exact linear decoder of the newest sample from history
  coordinates, d_i = U[end, i] s_i / gain; it maps coordinate rates to the endpoint rate too.
- **`TangentReadout`**: a derivative-consistent scalar readout: Wendland-weighted (*Adv. Comput.
  Math.* 4, 1995) bounded local approximants around k nearest neighbours in asinh coordinates, with
  the exact rate along the query velocity.

## Other options

- **HAVOK** `derivative_order` 1, 2 (default) or 4 (Brunton et al., *Nat. Commun.* 8, 2017).
- **Initial-condition tables**: `initial_conditions.mode: "file"` with `path`, optional `sha256`
  (refused on mismatch), `layout` `rows` or `columns`, `first`, and `norm` / `norm_tolerance` (every
  state on that norm). `examples/native/duffing-hdmd-sindy-sop.json` uses the retained Korda–Mezić
  table (`data/duffing_korda_mezic_initial_conditions.csv`; Korda & Mezić, *Automatica* 93, 2018).
- **Spectra of several observables**: one `hdmd` stage per scalar observable with identical
  settings, then `observable_spectrum_comparison` of their kroms (the spectra and their symmetric
  nearest-neighbour RMS distances; `examples/native/duffing-observable-spectra.json`).
- **Eigenvalue summaries**: every spectrum with continuous rates also stores `growth_rate`,
  `angular_frequency`, `cyclic_frequency` and `quality_factor` (Q = \|ω\|/2\|σ\|; `core/Eigenvalue.hpp`).
- **Empirical spectral weight**: `eigenfunction_errors` with `"spectral_weight": true` writes
  \|ψ_j(x_k)\|² / Σ_k \|ψ_j(x_k)\|² per mode, an exploratory localization, not an identified measure.
- **SVD driver**: any post-processing step's `svd_driver` (`gesdd`/`divide_and_conquer`, default, or
  `gesvd`/`jacobi`; Gu & Eisenstat 1995; Demmel & Veselić 1992). The dense kernels are Eigen, single
  threaded; `execution.max_threads` caps worker threads.
- **POD amplitudes as functions of the state**: `svd`/`pca`/`pod` `amplitude_queries` interpolate
  the modal amplitudes over the sample states (inverse distance, 8 neighbours; Shepard 1968) and
  write phi(x) and the reconstruction at each query (Holmes, Lumley & Berkooz 1996).
- **Archive coordinates**: `"analysis_coordinates": "archive"` analyses the stored values
  (dimensional for an RPE archive) instead of the integration coordinates; the output records it.
  A model fitted this way is not in the system's integration coordinates.
- **Stage timing and memory preflight**: done in phase C (derived archives record stage timing; plans
  report Hankel and SINDy geometry and refuse jobs over their memory budget).
