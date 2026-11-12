# Joint Hankel DMD and echo-measured KLQE: experiment protocol

Work date: 2026-11-05 (America/New_York). Implementation and the first bounded experiment are complete; accurate state estimation was not achieved. See `docs/reports/JOINT_HDMD_KLQE_SELECTION_2026-11-05.html`.

## Data and controls on inference

Use the same archived viscous RPE ring-down described in `docs/reports/RPE_VISCOUS_HIRES_RINGDOWN_2026-11-04.md`: r(0)=2.333386598601735, v(0)=0, dt=0.0001, tau 0..80. Its resolved physical parameters, including vapor pressure 2330 Pa, are authoritative. The archive remains local. No new physical model or control law is introduced.

Joint inputs are radius and the native far-field acoustic echo. The training windows retain every echo sample. Window starts may be sparsely sampled across the decay; each DMD pair is still separated by dt=0.0001. Reference Psi is the fitted modal projection of an observed joint window, not a claim of an exactly known autonomous eigenfunction. Alignment is the window endpoint. State reconstruction can additionally fit velocity against the same modal coordinates with the existing least-squares reconstruction API.

Deployment receives only echo samples. Physical state and joint-window reference coordinates are reserved for scoring. A state interpolant may be saved in the normal KROM format but must not be used to supply the observer with truth or to initialize it at the correct state.

## Selection criteria

A mode's C_x column norm alone is not invariant: multiplying psi_j by a divides its reconstruction column by a. Normalize each mode by its training RMS, preserving the inverse rescaling of every reconstruction map. Scale state outputs by explicitly stated accuracy units, and whiten the measurement by its assumed sensor noise.

For a subset S, keep its original C_x and C_y columns. Do not refit them during the primary selection comparison. Report separately:

1. Actual normalized state reconstruction MSE from C_x,S Psi_S versus reference x. This measures truncation and the decoder together, including cancellation between modes.
2. Actual echo reconstruction MSE from C_y,S Psi_S versus y. Omitted echo dynamics can drive an observer innovation even with perfect estimates of its retained coordinates.
3. Autonomous propagation defect on time-aligned reference Psi. A good snapshot decoder does not establish useful dynamics.
4. A finite measurement-window observability Gramian in the normalized coordinates. Analyze singular values and directions, not only its diagonal. Directional information and correlations can make individually strong modes jointly hard to distinguish.
5. A state-weighted uncertainty benchmark. With an explicitly declared unit Gaussian prior in normalized modal coordinates and independent measurement noise, posterior initial-state covariance is P=(I+W_o)^-1. Propagate it to the observation-window endpoint and decode it with C_x. The trace, divided by the number of state outputs, measures expected normalized state error from measurement uncertainty. This is a linear-model benchmark, not an assertion that real model mismatch is Gaussian or independent.

A transparent candidate score adds state truncation MSE, state-weighted noise MSE, and separately weighted echo mismatch and propagation penalties. Retain the components and compare a Pareto frontier; do not hide tradeoffs behind one unexplained score. Compare rankings by reconstruction alone, observability alone, and the combined criterion. Test multiple subset sizes and horizons. Use held-out time samples and a separate same-physics trajectory where available; do not tune the ranking on final validation errors.

Conjugate pairs must be treated consistently. SVD directions are mixtures of eigenfunctions, not subsets. Preserve eigenfunction identities during subset selection. A real observer can represent each conjugate pair with one independent complex coordinate (two real coordinates), with the exact corresponding transformation of both reconstruction maps. Do not double-count the pair as four independent real unknowns.

## Required evidence

- Singular values of the joint Hankel fit, retained rank, and numerical threshold; never label unresolved numerical directions as reliable eigenfunctions just to exceed 3,000.
- Full/reference-Psi reconstruction of radius, velocity, and echo, before any observer.
- Retained-subset reconstruction using reference Psi.
- KLQE estimates from echo alone, with initial estimate visibly different from truth.
- Echo, radius, and velocity actual/reference reconstruction/observer overlays; Psi traces and errors at aligned comparison times.
- Observability spectra/directions, subset tradeoff plots, and errors before/after convergence, including early high-amplitude dynamics separately from late decay.
- Noise-free and stated-noise runs; held-out trajectory checks; estimation failure and unstable/poor-propagation candidates reported without concealment.

Echo agreement by itself does not guarantee Psi agreement or state agreement. Detectability and reconstruction relevance of weak directions, finite information windows, model defect, initial uncertainty, and omitted-mode forcing all matter.

## Scope of the executed first study

The requested large fit retained 3,188 modes. Fifteen subsets were assessed at T=2, and three were run through the standard KLQE. The full complex-bank Gramian uses T=0.001 because of unstable spurious rates. Separate 600- and 80-mode fits diagnose rank sensitivity; the full 80-mode observer was also evaluated with noise and on a held-out IC. Subset scores currently use fitted endpoints, not independent validation windows. A broad horizon sweep, Pareto optimization, held-out modal projection, and a validated successful subset remain future experiments; the report does not claim they were achieved.

Build this directory as a standalone CMake project with `DYNSYSLAB_SOURCE_DIR` and `DYNSYSLAB_SUITE_BUILD` pointing to the source and compiled suite. Put the suite's `bin` directory on PATH. `rpe_joint_klqe` supports `inspect`, `modal`, `screen`, `full-real`, `assess`, and `run`; all output directories should be new. `rpe_joint_klqe_plots` assembles native Viewer panels. `report.py` only assembles the HTML from those existing results. The exact local job specifications and archive pointers are listed in the report. No scientific data is bundled with this source.
