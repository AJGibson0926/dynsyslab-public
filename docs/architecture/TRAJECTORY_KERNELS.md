# Trajectory kernels, orbit geometry and model libraries

Small native tools that studies used to run from side-car scripts (`tools/experiments/`), now in
the C++ core with pipeline methods. Each is used from experiments like any other method.

| Tool | Core | Method |
|---|---|---|
| Cross-orbit injectivity and closure | `core/OrbitGeometry.hpp` | `orbit_geometry` (pooled) |
| Affine equations and exact propagation | `core/AffineDynamics.hpp` | `affine_dynamics` |
| Backward prehistory | `core/TrajectoryAudit.hpp` `backwardHistory` | `prehistory` |
| Clock/interpolation audit | `core/TrajectoryAudit.hpp` `compareClocks` | `clock_audit` |
| High-order central differences | `core/NumericalDifferentiation.hpp` `differenceWeights`, `highOrderCentralDifference` | derivative method `central_high_order` |
| Model library families | `core/ModelLibraries.hpp` | `sindy_dynamics` `library_family` |
| Conserved quantity from run identity | `core/RunInvariant.hpp` `fitRunInvariant` | `run_invariant` (pooled) |
| Crossing refinement on a fine record, visit matching | `core/OrbitGeometry.hpp` `refineVisits`, `matchVisits` | (core API) |
| Per-record fitting windows from metadata | — | `dmdc`/`edmdc` `time_window: {"from_metadata": ["drive_start", "drive_stop"]}` |
| Pooled and held-out DMDc/EDMDc | — | `dmdc`/`edmdc` `pool_sources` (one fit over several records; pairs never cross records) and `validation_sources` (frozen one-step score under `validation`) |
| Validation scores of KLQE subset benchmarks | — | `linear_analysis` `klqe-assessment` `role: validation` |
| Physical states and measurement along projected histories | — | `history_projection` keeps `source_states` and `measurement_values` |

Tests: `dynsyslab_orbit_geometry_tests`, `dynsyslab_trajectory_kernels_tests`,
`dynsyslab_model_libraries_tests`, `dynsyslab_experiment_sweep_tests`, `dynsyslab_dmd_control_postprocess_tests`
(windows from metadata).

## Orbit geometry (`orbit_geometry`)

A coordinate map (a delay embedding, an observable chart) is injective on an attractor only if
distinct states never share a point (Takens 1981; Sauer, Yorke and Casdagli 1991). Two orbits whose
projections cross at different physical states witness that it is not. The method pools every
source of its job (each one orbit in coordinates, for example a `history_projection` of a
validation record through a frozen basis, with its states in `source_states`).

```json
{"name":"orbit_geometry","settings":{"role":"validation","axes":[0,1],"dimensions":[2,4],
  "peaks":{"array":"measurement_values","prominence":0.15,"min_spacing_time":0.5},
  "cycle_samples":400,"cycle":0,"groups":["inside","inside","outside"]}}
```

- **Peaks** (`findPeaks`): local maxima of the peak signal (a column of an archive array, default
  the first plotted coordinate) whose topographic prominence is at least `prominence` times the
  signal's range, at least `min_spacing_time` apart (lower peaks nearby are dropped).
- **Closure** (`peakClosure`): the median distance between the orbit at consecutive peaks over its
  RMS radius. A periodic orbit in an injective chart closes; drift or a non-injective chart does
  not.
- **Representative cycle** (`representativeCycle`): the last complete period between peaks,
  resampled to `cycle_samples` points; its enclosed area (`polygonArea`, shoelace formula).
- **Cross-orbit crossings** (`findCurveIntersections`): proper crossings of every pair of
  representative cycles in the `axes` plane, by the parametric segment test (O'Rourke 1998,
  section 7.2), coordinates scaled by the range of all cycles. Each is lifted to the first
  `dimensions` coordinates and to the physical states (scaled by their ranges): a crossing whose
  states differ is a state ambiguity of the chart. `nearestSeparation` is the closest approach of
  two cycles.

The summary (under the `role` key, so a `select` stage can score it) counts crossings by group
pair, gives closure medians by group and the state-ambiguity median and maximum. Arrays `cycles`
and `crossings` (record a, record b, x, y, state separation, lifted separation per dimension).
This replaces `tools/sweep_rpe_echo_delay_injectivity.py` and the measurement code it imported
(`audit_rpe_echo_latent_injectivity.py` `proper_intersections`, `polygon_area`,
`pairwise_minimum_separation`; `project_rpe_inviscid_echo_basis_out_of_sample.py`
`closure_metric`); the delay sweep itself is an experiment (EXPERIMENTS.md, "Sweeps and
selection").

## Conserved quantities (`run_invariant`)

`{"degree": 3, "pair_stride": 1}` pools unforced runs (one per source) and finds the polynomial of
the coordinates whose mean differs most between runs while it varies least within each: the
leading generalized eigenvector of the between-run and within-run scatter of the monomials
(Fisher 1936), solved through a Cholesky factor. It is a nonconstant Koopman eigenfunction at
eigenvalue 0 (multiplier 1; Mezić 2005), found from run identity alone. The report keeps it as a
polynomial model (`invariant_model_json`), the separation and per run its mean, spread and sampled
invariance defect.

## Crossing refinement and visit matching

`refineVisits` locates the two visit times of a crossing witness on an independent, finer record by
Newton's method on its interpolant inside the witness's brackets; `matchVisits` pairs the visits of
two records one to one within time tolerances (maximum cardinality, then least cost; Hungarian
method, Kuhn 1955; Jonker and Volgenant 1987).

## Affine dynamics (`affine_dynamics`)

`{"model": "<dynamical_model.npz>"}` reads a fitted vector field (`dynamics_model_json`), extracts
A = Df(0) and c = f(0) and verifies the field is affine at probe points (an active nonlinear term
is an error: such a model has no A, c form). It publishes A, c and the exact propagation from the
source's first sample, exp(tM)[x0; 1] with M = [[A, c], [0, 0]] (Van Loan 1978), by repeated
doubling of one step's exponential on a uniform clock, with the relative error against the source.

## Prehistory (`prehistory`)

`{"duration": T, "step": dt, "maximum_step": h, "roundtrip_tolerance": 1e-8}` precedes a simulated
record with T of its own unforced system's past (RK4 backward, steps ≤ h), verified by integrating
forward again to the first sample. Causal history coordinates then have complete windows from the
record's start. It is refused for driven systems and recorded as simulated, never as measured.

## Clock audit (`clock_audit`)

`{"reference": "<dense.npz>", "window": [a, b]}` compares the source's samples, interpolated
linearly without extrapolation, with a dense reference at the reference's samples in their common
window: max |error|, RMS and the worst time per component.

## High-order central differences

The derivative method `central_high_order` uses the central stencil of `points` samples (odd,
default 9: eighth order) and one-sided stencils of the same width at the ends, with weights from
Fornberg's recursion (Fornberg 1988, Math. Comp. 51:699-706). `differenceWeights(nodes, m)` gives
the weights for any distinct nodes and derivative order.

## Model library families

`sindy_dynamics` accepts `library_family` `bounded`, `rational`, `localized`, `gated`,
`paired_gated` and `harmonic` besides `broad` and `polynomial` (`polynomial_degree`,
`library_centers`, `library_harmonics`, `rational_widths`). Each family is frozen on the training
coordinates and is a set of analytic expressions, so fits, derivatives and rollouts use the shared
expression engine. With z = x/√(1 + (x/3)²):

| Family | Terms |
|---|---|
| bounded | 1, x, monomials of z |
| rational | bounded, w²/(w² + x²) and w x/(w² + x²) per width (Kaheman, Kutz and Brunton 2020) |
| localized | bounded, a normalized partition of unity (1 + |z − c|²/s)⁻² / Σ over farthest-point centers (Gonzalez 1985; Shepard 1968), times [1, z] |
| gated | g(E) x with E the mean square of x and g ∈ {1, E/(E + s)} |
| paired_gated | gated plus p/(1 + p) per consecutive coordinate pair |
| harmonic | x and gated a(x) Re uᵏ, a(x) Im uᵏ, u = (x₁ + i x₂)/√(|·|² + 0.015) |

Ridge, sample weights and weak-form rows of the regression are SINDy options of phase D
(`SindyDynamicsOptions`).

## References

- F. Takens, Detecting strange attractors in turbulence, Lecture Notes in Math. 898:366–381, 1981.
- T. Sauer, J. A. Yorke and M. Casdagli, Embedology, J. Stat. Phys. 65:579–616, 1991.
- J. O'Rourke, Computational Geometry in C, 2nd ed., Cambridge, 1998.
- C. F. Van Loan, Computing integrals involving the matrix exponential, IEEE Trans. Autom.
  Control 23:395–404, 1978.
- B. Fornberg, Generation of finite difference formulas on arbitrarily spaced grids, Math. Comp.
  51:699–706, 1988.
- E. Hairer, S. P. Nørsett and G. Wanner, Solving Ordinary Differential Equations I, Springer, 1993.
- S. L. Brunton, J. L. Proctor and J. N. Kutz, Discovering governing equations from data by sparse
  identification of nonlinear dynamical systems, PNAS 113:3932–3937, 2016.
- K. Kaheman, J. N. Kutz and S. L. Brunton, SINDy-PI, Proc. R. Soc. A 476:20200279, 2020.
- R. A. Fisher, The use of multiple measurements in taxonomic problems, Ann. Eugenics 7:179–188, 1936.
- I. Mezić, Spectral properties of dynamical systems, model reduction and decompositions, Nonlinear
  Dyn. 41:309–325, 2005.
- H. W. Kuhn, The Hungarian method for the assignment problem, Naval Res. Logist. Q. 2:83–97, 1955;
  R. Jonker and A. Volgenant, Computing 38:325–340, 1987.
- T. F. Gonzalez, Clustering to minimize the maximum intercluster distance, Theor. Comput. Sci.
  38:293–306, 1985.
- D. Shepard, A two-dimensional interpolation function for irregularly-spaced data, Proc. 23rd ACM
  National Conference, 1968.
