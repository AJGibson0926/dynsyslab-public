# KROM observation maps and training-domain expansion

DynSysLab treats three Koopman properties as separate contracts:

1. eigenfunctions provide lifted coordinates and autonomous evolution;
2. an observation map reconstructs physical/model state from those coordinates;
3. the archived training domain states where polynomial evaluation is supported.

Having eigenfunctions does not imply that state is linearly observable from a
chosen subset. Likewise, a high-quality decoder inside the training domain does
not justify extrapolation beyond it.

## Validated observation-map workflow

`rom.fit_state_observation_map(...)` fits the affine real-coordinate decoder

```text
x = intercept + [Re(psi_1), Im(psi_1), ...] C
```

Complex eigenfunctions are indivisible real/imaginary coordinate groups. The
fit uses either an independent validation trajectory or a deterministic held-
out subset. It records per-state R2, RMSE, normalized RMSE, maximum error, the
regularized condition number, split provenance, and all promotion thresholds.

The default quality gate requires:

```text
held-out R2 for every state             >= 0.99
held-out normalized RMSE for every state <= 0.10
decoder condition number                 <= 1e12
```

Failure raises and publishes nothing. `--allow-quality-failure` is an explicit
diagnostic override; such a map remains `validated=false`. KROM capabilities
distinguish `state_reconstruction` from `validated_state_reconstruction`.
Extracting a different mode subset invalidates the decoder. A decoder can be
attached only through `KROM.with_observation_map(...)`, which checks lifted and
state dimensions.

CLI examples:

```powershell
python -m rom fit-observation selected-krom.npz --report-only
python -m rom fit-observation selected-krom.npz --validation-solution independent-solution.npz
python -m rom fit-observation selected-krom.npz --allow-quality-failure --report-only
```

Successful publication creates a new immutable
`koopman-krom-observation-map/v1` derivative. It never rewrites the source KROM.

## Current RPE result

The five-mode observable-selected RPE KROM failed the default held-out gate:

| state | fit R2 | held-out R2 | held-out normalized RMSE |
|---|---:|---:|---:|
| radius | 0.345551 | 0.344244 | 0.169835 |
| radial velocity | 0.815555 | 0.814929 | 0.027583 |

The normalized decoder condition number was `1.3244e16`. No observed KROM was
published. This independently confirms that changing retained-mode count cannot
make the present polynomial eigenfunction family a high-quality linear state
observation basis.

Polynomial KROM validation is vectorized through its dictionary coefficients;
the 25,001-sample five-mode audit takes about 1.7 seconds on the current host,
down from about 19 seconds for row-by-row callable evaluation.

## Generic KROM training-domain planner

`dynsyslab-krom-training` accepts any compatible ODE KROM plus a required
state-domain box, unions it with the current KROM domain, adds margin only on
axes that require expansion, and emits new Cartesian initial conditions
outside the existing domain. State dimension and coordinate convention come
from the KROM; the same implementation supports two-state Duffing, three-state
Lorenz, dimensionless RPE, and higher-dimensional registered systems. It is a
dry run unless `--publish` is supplied. Publication writes campaign/job plans
but does not execute them.

```powershell
dynsyslab-krom-training model.npz `
  --required-min -2 -1 --required-max 2 1 --points-per-axis 9 11
dynsyslab-krom-training lorenz-krom.npz `
  --required-min -20 -25 0 --required-max 20 25 50 `
  --points-per-axis 7 7 9 --preset lorenz_1963 --screen
```

The Python API additionally accepts named `StateEnvelopeRequirement` objects,
arbitrary anchor states, and `StateProbe` samples. `target_state_requirement`
converts any normal DynSysLab target into that representation. The former RPE
function remains only as a convenience preset: it supplies canonical harmonic
envelopes, phase probes, `air_water_20c`, and radius-specific safety limits to
the generic planner and screener.

For all current canonical RPE targets, the default plan reports:

```text
current domain:   r=[0.232875, 2.3334], u=[-26.909884, 26.910801]
required target:  r=[0.35, 6.833333],   u=[-19.373155, 19.373155]
planned domain:   r=[0.232875, 7.493379], u=[-26.909884, 26.910801]
new IC count:     72 (9 x 11 grid, existing-domain points removed, anchors kept)
integration:      tau=0..25, dtau=0.001
```

The planner converts model-coordinate initial conditions and time settings to
the solver's physical API boundary when the KROM is dimensionless. Dimensional
KROMs use identity conversion. Exact system name, preset, and constructor
parameters are carried into the resulting Campaigner trajectory ensemble.

## Executable screening gate

`--screen` runs the canonical ODE kernel without archive or Viewer side effects.
The public `solver.solve_ode_in_memory(...)` boundary rejects checkpoints and
postprocessors, so disposable screening cannot masquerade as a production
solution. The generic policy integrates each state for 0.25 model-time units,
limits the adaptive step to 0.002 model-time units, and rejects non-finite
states plus caller-specified component-wise bounds. The RPE preset supplies:

```text
r < 0.05
r > 15
abs(u) > 250
```

Screening covers proposed Cartesian ICs and every supplied probe. A rectangular
domain is not considered safe merely because accepted IC extrema span its
bounds. Campaign publication additionally requires every anchor and probe to
pass. `--publish` without `--screen` is rejected, as is publication from a
non-ready screen.

The retained RPE all-target preset screen produced:

```text
candidate ICs:       48 accepted, 24 rejected
target-path probes: 163 accepted, 35 rejected
publication-ready:  false
```

All 24 candidate failures crossed the radius floor. Rejected candidate runs
reached radii as low as approximately `6.66e-5` and dimensionless rebound
velocity above `3e8`. Target-path failures occurred only on inward-moving
states (`u < 0`): 35 of 198 probes failed, with radius as low as `8.44e-4` and
velocity above `2.0e6`.

Only `rpe_multisine_equilibrium_small` passed all 33 phase probes. Every harder
offset target contained at least four unsafe free-response probes. Therefore
the proposed all-target free-response campaign was not published. A box formed
from safe ICs is insufficient for those targets; the next training workflow
must use safely constrained open-loop/controlled excitation, local KROMs, or a
controlled-data identification method.

The immutable screening record is:

```text
archives/campaigns/screening/
rpe-krom-screen__fa375084285016c1e3d4a853e11abc1f101cf1eca9734c10e87c2eb1cb90c08b.json
```

Its plan fingerprint is
`d8e8cb0ff90b45f1610bf88121dd52da36b3e5b50c5a39df5bc1a05b4a66e2ac`.

After a successful safe campaign, fit H-DMD/SINDy using combined archives, select
modes, validate the state observation map on independent trajectories, and only
then repeat KLQR tuning and long-horizon control evaluation.
