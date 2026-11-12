# KLQE and Koopman observability: first-pass Python implementation

Date: 2026-08-30  
Status: implemented and covered by targeted regression tests  
Scope: Python DynSysLab only; the MATLAB suite remained read-only

## Outcome

DynSysLab now has an isolated first implementation of observable-driven
Koopman linear-quadratic estimation. It does not copy the MATLAB class layout.
Instead it preserves the useful mathematics behind clean ROM, estimation,
analysis, storage, CLI, and Viewer boundaries.

The implementation includes:

- multiple named observable reconstructions on one KROM;
- corrected H-DMD window/state alignment for scattered eigenfunctions;
- raw classical observability matrices and stable Gramian alternatives;
- right observability directions stored and plotted as columns;
- continuous-time steady-state KLQE design with explicit measurement and state
  reconstruction maps;
- covariance, shape, PBH detectability, CARE-residual, and pole-stability gates;
- exact augmented-exponential propagation for zero-order-held measurements;
- a standalone offline measurement runner that does not require a controller;
- separate estimate, measurement, predicted-measurement, innovation, truth, and
  state-error time-series objects;
- immutable numerical artifacts for observability, KLQE designs, and KLQE runs;
- report-only-by-default CLI commands with explicit output publication;
- Viewer-native observability and KLQE dashboards; and
- focused regression tests for all new numerical boundaries.

## Code boundaries

### ROM observable projections

`rom.observation.ObservableReconstruction` remains the authoritative
mathematical representation of

```text
h(x) = c0 + C Psi(x).
```

`KROM` now retains a named collection in `observable_reconstructions` while
preserving the earlier `observable_reconstruction` field as the selected
default. Existing single-map callers remain valid.

Stable calls are:

```python
krom.observable_map("far_field_acoustic_echo")
krom.observable_map("state")
krom.C("far_field_acoustic_echo", output="split")
krom.C("state", output="complex")
```

`with_observable_reconstruction(...)` returns an independent KROM derivative.
It can add a map without changing the default. KROM extraction restricts every
named C map to the same retained eigenfunctions, including the existing
conjugation policy. Named maps are stored as self-contained coefficient and
intercept arrays when KROM storage publishes a derived archive.

### H-DMD alignment gate

In scalar window-snapshot H-DMD, row `i` of a projected mode corresponds to the
physical state at index

```text
i + delays - 1,
```

the right/end edge of the delay window. The previous implementation attached
those values to `state[i]`. This has been corrected for individual and pooled
ensembles, and the anchor is now explicit in metadata as
`interpolation_state_anchor="window-end"`.

Delay-state H-DMD retains its existing window-start alignment and now records
that convention explicitly where appropriate.

### Observability analysis

`rom.observability` provides:

```python
classical_observability_matrix(A, C, horizon=...)
discrete_observability_gramian(A, C, horizon=...)
continuous_observability_gramian(A, C, horizon=...)
continuous_infinite_observability_gramian(A, C)
analyze_observability(A, C, ...)
analyze_krom_observability(krom, "far_field_acoustic_echo", ...)
```

The classical matrix is built by recurrence and has an overflow/safety gate.
Finite continuous Gramians use a deterministic Van Loan block exponential.
The infinite-horizon form is admitted only for a strictly stable generator.

`ObservabilityAnalysis` retains A, C, the Gramian, optional classical matrix,
singular strengths, right directions, numerical rank, tolerance, horizon,
grouped nonredundant-mode participation, and metadata. Directions are columns,
not rows. The artifact is saved and loaded without executable callables.

### KLQE design

`KLQEObservationModel` binds three explicit objects:

```text
KROM
measurement ObservableReconstruction (Cy and measurement intercept)
state ObservableReconstruction or StateObservationMap (Cx and state intercept)
```

`design_klqe(...)` and `design_klqe_from_krom(...)` use the continuous split
generator returned by `KROM.A(output="split")`. They validate:

- Q is symmetric positive semidefinite;
- R is symmetric positive definite;
- all shapes and entries are finite;
- every non-decaying A mode passes the PBH detectability test;
- the continuous algebraic Riccati residual is below a declared tolerance; and
- every pole of `A - L Cy` is decaying when the stability gate is enabled.

The result is an immutable `KLQEDesign` with P, L, observer poles, residual,
detectability report, observation metadata, and the exact Q/R used.

### Propagation and standalone filtering

For

```text
M = A - L Cy,
```

`KLQEEstimator` computes the exact held-measurement transition using

```text
exp([[M, L], [0, 0]] dt).
```

The measurement at sample `k` is held causally over `[t_k, t_{k+1})`. Matrix
exponentials are cached by time step. This replaces the earlier explicit Euler
update and remains valid when M is singular.

`run_klqe(...)` accepts measurement samples directly. Truth is optional and is
used only for benchmark scoring. It returns `KLQEResult`, containing separate
`TimeSeries` instances for:

- lifted estimate;
- reconstructed state estimate;
- measured observable;
- predicted observable;
- innovation;
- optional true state; and
- optional state error.

This makes a free echo-estimation experiment independent of controller
plumbing. A controller can still consume `KLQEEstimator.estimate(...)` through
the established estimator interface.

### Immutable artifacts

The public storage calls are:

```python
save_observability_analysis(...)
load_observability_analysis(...)
save_klqe_design(...)
load_klqe_design(...)
save_klqe_result(...)
load_klqe_result(...)
```

A saved KLQE design includes the exact KROM path and a streaming SHA-256
fingerprint. Loading refuses a changed KROM and recomputes the design from the
named maps before comparing A, Cy, Cx, P, L, and observer poles to the stored
arrays.

## Lightweight CLI

The module and installed entry point are:

```text
python -m estimation ...
dynsyslab-estimation ...
```

Available first-pass commands are:

```text
python -m estimation observability KROM.npz \
  --measurement far_field_acoustic_echo \
  --method continuous-gramian --horizon 10

python -m estimation design KROM.npz \
  --measurement far_field_acoustic_echo --state-output state

python -m estimation benchmark KROM.npz FREE_SOLUTION.npz \
  --measurement far_field_acoustic_echo --state-output state
```

All three commands print a JSON report and make no files by default. Use an
explicit `--output` to publish an immutable NPZ artifact. Observability and
benchmark commands accept `--viewer OUTPUT.html` for an explicit Viewer
report. Scalar process/measurement variance flags are first-pass conveniences;
the Python API accepts full Q and R matrices.

## Viewer behavior

`ObservabilitySource` and `PlotRequest("observability", ...)` render:

1. logarithmic observability strengths;
2. the signed split-coordinate coefficients of each right direction; and
3. grouped nonredundant Koopman-mode participation when available.

`KLQESource` and `PlotRequest("klqe", ...)` render:

1. each estimated state component with benchmark truth when present;
2. each measured observable with its reconstruction;
3. each innovation component;
4. a two-state phase portrait when available;
5. the state-error norm for benchmark runs; and
6. observer poles.

Viewer reads the stored analysis/result. It does not recompute an observer or
mutate the KROM. Existing automatic time-detail insets also apply to KLQE time
panels.

## Verification

Focused tests cover:

- single and ensemble H-DMD end-window alignment;
- classical matrix construction and right-direction orientation;
- equality of a discrete finite Gramian to `O.T @ O`;
- an analytic continuous-Gramian benchmark;
- observability artifact round trips;
- named C maps and their behavior under KROM extraction;
- named C-map archive restoration;
- observable-driven CARE design and stability diagnostics;
- exact held-measurement propagation against a closed form;
- offline histories and result round trips;
- Viewer observability/KLQE panels and routing; and
- estimation CLI discovery.

Final repository-wide verification: **368 passed, 1 skipped**. The only warning
is the pre-existing pytest collection warning for `tester.core.Tester`, which
has an explicit constructor and is exercised through its dedicated tests.

## Intentional first-pass limits

These boundaries are deliberately not hidden:

1. The implemented estimator design is continuous-time only. A genuinely
   discrete KROM needs a separate DARE design and discrete recursion.
2. The offline runner currently models free dynamics. A known lifted forcing
   input can later be added to the same augmented exponential without changing
   the observer law.
3. No MATLAB file or MATLAB source was modified, and no MATLAB behavior was
   copied wholesale.
4. No RPE echo benchmark is claimed yet. The echo and state projections must
   first pass held-out static reconstruction checks on the chosen frozen KROM.
5. H-DMD-native terminal modal amplitudes remain marked as native/unvalidated,
   not silently promoted as held-out fits.
6. Training-covariance coordinate whitening and measurement-covariance-weighted
   observability views remain a later refinement. The current artifact records
   the raw declared split coordinates and measurement units.
7. Noise generation, innovation whiteness testing, and Monte Carlo covariance
   calibration remain benchmark-layer work.
8. Joint plant/observer integration is not required for free estimation and was
   not coupled into the solver in this pass. The offline path gives the same
   one-way plant-to-observer experiment without controller dependencies.

## Recommended next refinement

Freeze one sufficiently rich RPE KROM, then fit and independently validate both
`far_field_acoustic_echo` and `state` projections on held-out free trajectories.
Run the noiseless offline benchmark with a deliberately displaced initial
estimate, inspect the observability artifact, and only then compare the same A,
Cy, Cx, Q, R, and initial lifted state against a coherent MATLAB export. This
keeps model/eigenfunction error separate from observer-integration error.
