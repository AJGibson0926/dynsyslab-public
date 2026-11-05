# Poincare sections and bifurcation sidecar

DynSysLab's optional native phase-space sidecar exposes state-defined first
returns, fixed-forcing-phase maps, phase-conditioned periodic shooting, and
finite parameter continuation through `dynsyslab.phase-space-query/v1`. The
implementation is in the independently compiled C++ core; the workflow C ABI,
CLI, cache, and Viewer adapter consume the same versioned result.

This is numerical evidence for a declared model, section, branch, domain, and
tolerance set. A finite search is not a proof that all equilibria, cycles, or
bifurcations were found.

## Autonomous Poincare sections

The compact axis form uses `section_axis` and `section_value`. An affine
hyperplane may instead be supplied with full-state `section_normal` and
`section_reference` arrays:

`S(x) = section_normal . (x - section_reference) = 0`.

`section_direction` is `+1`, `-1`, or `0` for either orientation.
`section_minimum_time` excludes an explicitly declared initial interval. A seed
already on the section is not accepted as its own return: the event locator
requires departure and then selects the first oriented positive return.

For every `section_seeds` entry, the core integrates the state and variational
equation, brackets the return, refines its time, and repeats with a halved
maximum step and tighter tolerances. The result reports the section and return
states, return time, fixed-point residual, event orientation, endpoint and
derivative refinement errors, and the transversality denominator and angle.
Grazing returns are retained as `grazing` diagnostics without a false inverse.

For monodromy matrix `M`, final vector field `f`, and section covector `n`, the
reported derivative is

`DP = (I - f n^T / (n^T f)) M`,

with return-time derivative `D tau = -n^T M / (n^T f)`. The core constructs an
orthonormal tangent basis `Q` and emits `B = Q^T DP Q`. Matrices are row-major;
`tangent_basis` is an `n x (n-1)` matrix. Stability uses the eigenvalues of `B`,
not whichever full monodromy eigenvalue happens to be numerically closest to
one.

Autonomous periodic shooting uses the same configured section as its phase
condition. It emits both full Floquet multipliers and `section_multipliers`,
along with closure, phase, refinement, transversality, and tangent-identity
diagnostics. The neutral identity is checked directly as
`M f(x) = f(P(x))`; it is not inferred by deleting a near-unit eigenvalue.

## Fixed forcing phase

For `dynamics: periodic-forcing`, the section is the exact absolute forcing
time plus an integer number of forcing periods. The full physical-state
monodromy is the derivative of this stroboscopic map, so no autonomous neutral
multiplier is removed. See [periodically forced phase
context](PERIODIC_PHASE_CONTEXT.md) for clocks, subharmonics, waveform
reconstruction, and section-manifold semantics.

## Pseudo-arclength branches

The optional `bifurcation` object supports equilibrium and periodic-orbit
branches:

```json
{
  "enabled": true,
  "branch": "equilibrium",
  "parameter": "mu",
  "minimum": -0.1,
  "maximum": 0.8,
  "initial": 0.5,
  "step": -0.04,
  "minimum_step": 0.0001,
  "maximum_step": 0.08,
  "maximum_points": 55,
  "maximum_corrector_iterations": 20,
  "seed_state": [0.7071067811865476],
  "seed_period": 0.0
}
```

The signed first step selects the initial branch direction. Predictor and
bordered corrector operate in domain-scaled coordinates; failed correctors
halve the arclength step. Equilibrium branches solve `f_mu(x)=0`. Autonomous
periodic branches solve closure plus the configured section phase condition in
`(x,T,mu)`. Periodically forced branches solve `P_mu^m(x)-x=0` at the fixed
forcing phase. Parameter derivatives are centered finite differences of the
actual configured family, including the return endpoint.

Each point records parameter, state, period when applicable, scaled residual,
step/correction diagnostics, stability, and the relevant spectrum. Adjacent
spectral sign changes produce bracketed candidates:

- equilibrium real eigenvalue zero;
- equilibrium complex-pair real part zero (`hopf_candidate`);
- periodic multiplier `+1`;
- periodic multiplier `-1` (`period_doubling_candidate`);
- complex periodic multiplier modulus one (`neimark_sacker_candidate`).

These labels are candidate detectors, not nonlinear classification proofs.
For scalar zero-eigenvalue events only, the sidecar also estimates
`f_mu`, `f_xx`, `f_xmu`, and `f_xxx`, checks a preserved central branch and
reflection oddness, and may return `fold`, `transcritical`, or `pitchfork` as a
`structurally-classified-candidate`. Near-zero tests remain bounded by the
recorded finite-difference and continuation bracket errors.

## Exact model oracles and commands

The native catalog includes three small source-backed regression models:

- `quadratic_saddle_node`: `x' = mu + eta*x - x^2`;
- `cubic_bifurcation`: `x' = epsilon + mu*x + quadratic*x^2 + cubic*x^3`;
- `radial_limit_cycle`: the exact radial cycle family with
  `r' = r*(mu-alpha*r^2)` and `theta' = omega`.

Run the public CLI examples with:

```powershell
dynsyslab-cli phase-space examples/native/radial-limit-cycle-poincare-query.json
dynsyslab-cli phase-space examples/native/quadratic-fold-bifurcation-query.json
```

The first checks a return time of `2*pi` and transverse multiplier
`exp(-4*pi)` on the unit-cycle preset. The second crosses the fold with
pseudo-arclength continuation instead of stopping at the singular natural
parameterization.

The retained RPE demonstration
`examples/native/rpe-periodically-forced-mask-trajectory.json` drives the
10-micrometre air/water model at 100 kHz and 20 kPa. Its immutable archive is
included in the active Viewer Mask and supplies time, phase, frequency, and
control views. It is a response trajectory, not a claim of branch
continuation.

## Verification and provenance

Native tests cover the exact radial return, projected multiplier, affine
section geometry, neutral tangent identity, scalar fold/transcritical/
pitchfork structure, autonomous periodic continuation, fixed-phase forced
continuation, forced RPE response, tighter-integration checks, and the RPE
Liouville determinant identity. The exact consulted Grothendieck Labs record
IDs and SHA-256 digests are preserved in
`docs/provenance/phase-space-class-sources.json`; the peer repository remains
unchanged.

Legacy campaign bifurcation plots remain independent parameter-response
sweeps. They are not silently relabeled as continuation branches.
