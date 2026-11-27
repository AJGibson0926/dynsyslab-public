# KROM observation maps and training-domain expansion

DynSysLab treats three Koopman properties as separate contracts:

1. eigenfunctions provide lifted coordinates and autonomous evolution;
2. an observation map reconstructs physical/model state from those coordinates;
3. the archived training domain states where polynomial evaluation is supported.

Having eigenfunctions does not imply that state is linearly observable from a
chosen subset. Likewise, a high-quality decoder inside the training domain does
not justify extrapolation beyond it.

## Observation-map workflow

The native `observable_reconstruction` post-processing method fits the
ridge-regularized affine real-coordinate decoder

```text
x = intercept + [Re(psi_1), Im(psi_1), ...] C
```

for a named observable from a KROM's eigenfunction values (settings `krom`,
`observable` or `observable_array`, `ridge`). Complex eigenfunctions are
indivisible real/imaginary coordinate groups. Loading the KROM with the source
physics checks system, parameter, and coordinate compatibility. The result is a
new immutable KROM archive; the source KROM is never rewritten. Extracting a
different mode subset (`krom_extract`) invalidates a decoder.

Judge a decoder on independent, held-out trajectories. The thresholds used for
the result below were:

```text
held-out R2 for every state             >= 0.99
held-out normalized RMSE for every state <= 0.10
decoder condition number                 <= 1e12
```

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

## Generic KROM training-domain planner

A solve job may carry a `preflight` object with `"kind": "krom_training_domain"`
and `"enabled": true`. It names the KROM's `current_domain` and a
`required_domain` (each `minimum`, `maximum`, `coordinate_system`), and
`points_per_axis`, `margin_fraction`, optional `anchors`, and
`maximum_initial_conditions`. The planner unions the required box with the
current domain, adds margin only on axes that require expansion, and emits new
Cartesian initial conditions outside the existing domain. State dimension and
coordinate convention come from the system, so the same implementation serves
two-state Duffing, three-state Lorenz, dimensionless RPE, and
higher-dimensional registered systems. `dynsyslab-cli plan` checks that the
job's explicit initial conditions match the planned set; `run` screens them
before any solve is archived.

For all canonical RPE targets, the RPE study plan reported:

```text
current domain:   r=[0.232875, 2.3334], u=[-26.909884, 26.910801]
required target:  r=[0.35, 6.833333],   u=[-19.373155, 19.373155]
planned domain:   r=[0.232875, 7.493379], u=[-26.909884, 26.910801]
new IC count:     72 (9 x 11 grid, existing-domain points removed, anchors kept)
integration:      tau=0..25, dtau=0.001
```

The preflight coordinates must match the system's coordinates; a mismatch is
an error rather than a silent conversion.

## Executable screening gate

The preflight `screening` object runs the canonical ODE kernel without archive
or Viewer side effects, so disposable screening cannot masquerade as a
production solution. By default it integrates each state for 0.25 model-time
units (`time_end`), limits the adaptive step to 0.002 (`maximum_step`), and
rejects non-finite states plus the component-wise bounds `state_minimum`,
`state_maximum`, and `maximum_absolute_state`. Every initial condition must
pass (`require_all`). The RPE study used:

```text
r < 0.05
r > 15
abs(u) > 250
```

A rectangular domain is not considered safe merely because accepted IC
extrema span its bounds.

The RPE all-target screen (which also screened target-path probe states)
produced:

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
