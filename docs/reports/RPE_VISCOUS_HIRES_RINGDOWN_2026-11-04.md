# High-resolution viscous RPE ring-down and echo-history coordinates

Date: 2026-11-04. This is a visualization and coordinate-extraction experiment,
not an eigenfunction fit or a control experiment. The physical trajectory was
already archived; no dynamics were rerun and no core pipeline code was changed.

## Source and physics

Source: `archives/solutions/rpe/rpe-viscous-echo-long-training-five-ics-trajectory-0__20261002T001350Z.npz`.
The record is unforced, starts at `(r,v)=(2.333386598601735,0)`, and contains
800,001 samples from tau 0 to 80 at spacing 0.0001. Its first 25 time units
contain 250,001 samples. This is ten times finer output sampling than the
original August spectral-training record, with the same physical parameters
and nondimensional groups.

Parameters: equilibrium radius 10 micrometers; density 998 kg/m3; dynamic
viscosity 0.001 Pa s; surface tension 0.0725 N/m; polytropic index 4/3;
ambient pressure 100000 Pa; vapor pressure 2330 Pa. One time unit is
3.0123833037972497 microseconds. Adaptive BDF tolerances are 1e-9.

## Native computation

The standard `history_basis` postprocessor evaluates the registered scalar
`far_field_acoustic_echo` observable and fits a Hankel basis using only that
scalar record. No radius, velocity, or echo-dot channel enters the decomposition.
The physical simulator supplies the synthetic echo; this is not laboratory data.

Settings: full interval 0–80, decimation 1, 12,751 delays, 1.275-tau oldest-to-newest
span, eight computed modes, two selected coordinates, oversampling 16,
32 subspace iterations, seed 20261104. The existing matrix-free native
Hankel/Rayleigh–Ritz implementation avoids materializing the enormous matrix.
Its standard convention uses RMS scaling without mean subtraction. The
normalization is recorded in the frozen transform.

The `history_projection` stage applies the frozen basis to causal echo windows:

\[
h_k=(e_{k-12750},\ldots,e_k)^T,\qquad
\xi_j(k)=U_j^T h_k/\sigma_j.
\]

The singular values correspond to the native RMS-scaled fitting matrix; the
runtime intentionally does not divide the raw echo by the fitting RMS a second
time. This is the existing native coordinate convention. Projection starts at
tau 1.275, with 787,251 valid coordinate samples. No padding or future samples
are used. The 1.275-tau span preserves the earlier single-trajectory visualization
window; it is distinct from the later 4.3215-tau pooled-control window.

## Results and checks

- First two modes capture 8.2508174325% of total Hankel squared singular-value
  energy; the first eight capture 17.8449466170%.
- Relative residuals of the leading two Gram eigenpairs: 4.35e-15 and 1.92e-15.
- Basis orthogonality error: 7.41e-15.
- Native projected coordinates versus direct history dot products at 19 times:
  maximum absolute difference 4.42e-17; relative difference 2.99e-15.
- The plotted analytic echo chart was checked against the native observable
  query at six states including the collapse/echo maximum: scaled error 2.22e-16.

Echo-dot is the analytic autonomous Lie derivative. Because this trajectory is
unforced, it equals the total time derivative of the same echo observable.
It is a comparison channel, not an input to the Hankel fit.

The first two coordinates form a contracting, distorted spiral. Their low
reconstruction-energy fraction means that this display is not evidence of a
two-coordinate lossless waveform representation, injectivity, or Markov closure.
Those properties require separate tests.

## Viewer output and local artifacts

All saved points are included in the native Viewer payloads, with no display
decimation. Greek labels use the native math renderer. IC/EC denote the physical
initial condition/equilibrium; `start` on the latent portrait denotes the first
complete history at tau 1.275.

![Full ring-down](../../viewer_outputs/rpe-viscous-hires-ringdown-full-20261104.png)

![Early 25 time units](../../viewer_outputs/rpe-viscous-hires-ringdown-early-20261104.png)

Frozen basis: `archives/derived/history-basis/rpe-viscous-echo-long-training-five-ics-trajectory-0__20261002T001350Z__history-basis__20261105T003259Z.npz`.

Projection: `archives/derived/history-projection/rpe-viscous-echo-long-training-five-ics-trajectory-0__20261002T001350Z__history-projection__20261105T003453Z.npz`.

Artifact names use UTC; the work date above uses America/New_York.
The two job JSON files, rendering helper, machine-readable checks, and Viewer
payloads are retained locally under `viewer_compute_archive/viscous-ringdown-20261104/`.
Reproduce the computation with the compiled `dynsyslab-cli run` workflow using
the settings above and the local source archive. The projection job's `transform`
must point to the newly published basis. Numerical archives and PNGs remain local.
