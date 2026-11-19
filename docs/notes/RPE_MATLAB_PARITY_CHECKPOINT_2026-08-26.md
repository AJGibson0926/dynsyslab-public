# RPE MATLAB-parity checkpoint — 2026-08-26

This checkpoint isolates the autonomous training trajectory and H-DMD stages.
No SINDy fit was run. The retained artifacts are immutable inputs for later
comparison with results exported from the MATLAB suite.

## Physical model

The `air_water_20c` RPE preset now uses the unrounded values required by the
reference formulation:

| Parameter | Value |
| --- | ---: |
| equilibrium radius, `R0` | `1.0e-5 m` |
| polytropic exponent, `K` | `4/3 = 1.3333333333333333` |
| ambient pressure, `p0` | `100000 Pa` |
| vapor pressure, `pV` | `2330 Pa` |
| liquid density, `rhoL` | `998 kg/m^3` |
| surface tension, `sigma` | `0.0725 N/m` |
| liquid viscosity, `muL` | `0.001 Pa s` |
| equilibrium gas pressure, `pG0` | `112170 Pa` |

`CL` and `CG` are not RPE parameters and are not represented.

The exact seeded training IC is stored as the dimensionless double
`(r(0), u(0)) = (2.3333865986017348, 0)`. The physical archive request uses
`(R(0), Rdot(0)) = (2.333386598601735e-5 m, 0 m/s)`, chosen by converting the
dimensionless value through the model rather than typing a rounded physical
decimal. Converting it back inside the solver recovers the intended
dimensionless double exactly. This autonomous training request carries the
seeded override explicitly. Free solves otherwise use their system default or
request-specific IC; only controlled RPE solves resolve a missing IC to physical
equilibrium `(R0, 0)`.

## Dimensionless formulation and BDF solve

The run integrates `r = R/R0`, `u = dr/dtau`, and `tau = f0 t`, with

`f0 = sqrt(3 K pG0 - 2 sigma/R0) / (2 pi R0 sqrt(rhoL))`

and `T0 = 1/f0`. The realized values are:

| Quantity | Value |
| --- | ---: |
| `f0` | `331963.06683132035 Hz` |
| `T0` | `3.0123833037972497e-6 s` |
| pressure scale `rhoL (R0/T0)^2` | `10997.907878457558 Pa` |
| dimensionless gas pressure | `10.199212544752802` |
| dimensionless vapor pressure | `0.21185847578919523` |
| dimensionless ambient pressure | `9.092638445888207` |
| capillarity group | `1.31843257465379` |
| viscosity group | `0.12073680576341679` |

The reference request uses SciPy BDF, `rtol = atol = 1e-9`, the dimensionless
RPE RHS, `0 <= tau <= 25`, and PCHIP output onto `dtau = 0.001`. BDF produced
12,820 raw adaptive samples and the retained grid contains 25,001 samples.

Trajectory archive:

`archives/solutions/rayleigh-plesset/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z.npz`

SHA-256: `e4529d55dd549e12c66dbb78c9521c6bc273aa8bb135637319568a2edfdc9eba`

The provisional MATLAB `ode15s` reference available during this checkpoint was
`<local-data>/matlab_rpe_training_reference_tau25_dtau0p001.mat`
(SHA-256 `7a97c41a5e4a42afc2e3957f378255e42e9294dc03c115d215e7f0ccd4de0efe`).
Against that file, the SciPy-BDF trajectory differences are:

| State | RMSE | maximum absolute error |
| --- | ---: | ---: |
| `r` | `1.6140220610e-6` | `3.0610312661e-5` |
| `u` | `1.5225666940e-5` | `2.7035549339e-4` |

The final-state delta `(Python - MATLAB)` is
`(-2.1370650805e-8, -1.7484177806e-6)`. These residuals are the bounded
integrator/interpolator difference after aligning the equation, constants,
initial state, time convention, horizon, grid, and tolerances.

## Full-spectrum H-DMD

H-DMD consumes integration coordinates, not the dimensional archive view. It
uses scalar radius, a `20000 x 2000` shifted Hankel pair, `gesvd`, projected
modes, an absolute singular-value threshold of `1e-10`, and no explicit rank
cap. All 2,000 singular directions were above the threshold and were retained.

After mode construction, modes are ordered by the established terminal
amplitude-energy policy:

`b = pinv(projected_modes) @ X[:, -1]`

with descending `abs(b)`. The ordering is a complete permutation of all 2,000
modes; it does not truncate the spectrum. The stored metadata includes the
complete original-index permutation and all sorted amplitude magnitudes.

H-DMD archive:

`archives/derived/koopman-hankel-dmd/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-hankel-dmd__20260826T221322.756123Z__604a65f1.npz`

SHA-256: `c2fabe83b37b94ace1abb5b7d4fdd11cad14dee97c476d2fbb00e5e63eb0bf4b`

The archive is 790,524,387 bytes and contains the 2,000 eigenvalues, the
`2000 x 2000` reduced operator, and the `20000 x 2000` projected-mode matrix.
Its eigenfunction-value storage is an alias of the mode matrix rather than a
duplicate payload.

Measured H-DMD timings:

| Stage | Seconds |
| --- | ---: |
| SVD | `46.8391` |
| total factorization section | `73.6325` |
| projected operator | `3.1575` |
| eigendecomposition | `13.9228` |
| projected-mode construction | `9.7130` |
| terminal-amplitude energy ordering | `78.9440` |
| algorithm | `154.0375` |
| total including immutable archive publication | `210.5420` |

## Reproduction and validation

The trajectory is reproduced by:

`python scripts/generate_rpe_reference_training.py`

The named post-processing policy now includes
`mode_sort="terminal-amplitude-energy"`. The full repository test run after
this checkpoint passed 289 tests with one optional test skipped.

When new MATLAB exports arrive, compare them first to the trajectory archive in
dimensionless integration coordinates, then compare the unsorted spectral set
and finally the stored energy order. Do not start a SINDy parity comparison
until those three boundaries agree.
