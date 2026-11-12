# RPE echo observability and long-horizon KLQR stability

Date: 2026-09-14. Native DynSysLab on the installed 0.4.0 runtime.

## Result

The earlier instantaneous-echo KLQR multi-sine case was rerun for four times
the original horizon, from `tau=0` to `tau=40`, with all other model, target,
KROM, KLQR, solver and sampling settings unchanged. The radius does not level
off. It reaches 26.1683 equilibrium radii at the final sample and is still
increasing with `v=0.4741`. After the initial transient, the last nonpositive
velocity sample is at `tau=0.897`; every later saved velocity is positive.

Echo tracking remains extremely close while the state escapes. Echo NRMSE is
0.02394% over the full run and 0.02815% over `30 <= tau <= 40`. The radius is
10.2785 at `tau=10`, 16.5381 at `tau=20`, 21.7210 at `tau=30`, and 26.1683 at
`tau=40`. Thus the longer experiment strengthens the original conclusion:
accurate output tracking is not state regulation.

![Long-horizon stability diagnostic](../../viewer_outputs/rpe-echo-tau40-stability.png)

## Why the decreasing slope is not leveling

For this parameter set, write `S=R0*Pscale`, use `r=R/R0` and `v=dr/dtau`, and
let `nu` denote the dimensionless viscous coefficient. The stored echo is

```text
e/S = (2*r - 3/2)*v^2 - nu*v/r
      + p_g*r^(-4) + p_v - p_inf - c/r.
```

At large radius, bounded echo tracking on the positive-velocity branch gives

```text
v approximately sqrt((e/S + p_inf - p_v)/(2*r)).
```

The tracked target keeps the expression under the square root positive: its
observed range is 3.7467 to 23.8411. Consequently `d(r^(3/2))/dtau` tends to a
positive bounded oscillatory quantity and `r` grows approximately as
`tau^(2/3)`, rather than approaching a finite limit. On `10 <= tau <= 40`, a
linear fit to `r^(3/2)` has slope 3.37439 and R-squared 0.999913. This explains
the decreasing raw radius slope without mistaking it for stabilization.

The finite run by itself does not prove behavior at infinite time. The
large-radius balance is a conditional asymptotic statement: it applies while
the controller continues to track the bounded echo and remains on the observed
positive-velocity branch.

## What kind of observability failed

Instantaneous echo is a smooth scalar map from the two-dimensional physical
state `(r,v)` to one number. Its Jacobian has rank at most one, so it has an
instantaneously invisible tangent direction. It is also globally non-injective.
At `r=1`, the states `v=0` and `v=2*nu=0.2414736115` both have zero echo, but
their autonomous echo derivatives are 0 and -1.0461073. Their input gains
`L_g h` are also opposite, +0.1207368 and -0.1207368. A controller represented
only as a single-valued function of instantaneous echo therefore cannot assign
the correct local input action to both states.

This is weaker than saying that an echo time history can never observe the
bubble. For a two-state autonomous model, the first nonlinear differential
observability test stacks the gradients of `h` and `L_f h`. Exact native
analytic derivatives show that this two-row map is full rank at every one of
201 sampled states on one period of the multi-sine target. In native
dimensionless coordinates its minimum row-angle sine is 0.2545 and its maximum
condition number is 15.27. Thus a short, sufficiently resolved echo history is
locally informative near the target.

The failed controlled trajectory moves into a much worse region. Treating the
applied input as known, 401 samples along the long run give a median condition
number of 185.94 and a maximum of 18,710; the minimum row-angle sine is 0.00142.
This is an observability-conditioning collapse, not an exact rank loss at the
sampled points. Noise and modeling error would be amplified severely there.

The distinction also explains why classical observability of the two-state
real split KROM would be insufficient evidence. Such a test asks whether the
lifted eigenfunction coordinates can be reconstructed through the KROM output
matrix. It does not test whether the physical-state-to-lift map is injective,
well conditioned, or dynamically closed. Here the echo-only eigenfunction fit
already has a normalized generator defect near 0.97.

## Useful diagnostics for the next model

The next control study should retain separate diagnostics for:

1. static fiber ambiguity: state separation among pairs with nearly equal
   instantaneous measurements;
2. nonlinear differential observability: rank, smallest singular value,
   condition number, and normalized row angle of the Lie-derivative map;
3. a noise-weighted finite-horizon empirical observability Gramian computed
   from perturbed initial states under the same known input;
4. state-decoder error and uncertainty from held-out measurement windows;
5. sign and spread of `L_g h` within narrow measurement bins; and
6. Koopman generator closure, held-out rollout, and training-domain coverage.

Hermann and Krener's nonlinear observability framework supplies the local
Lie-derivative rank test ([IEEE TAC 22(5), 1977](https://doi.org/10.1109/TAC.1977.1101601)).
Takens supplies the qualified generic route from scalar delays to an embedded
state on a compact invariant manifold
([LNM 898, 1981](https://doi.org/10.1007/BFb0091924)). Neither theorem promises
good finite-data conditioning, correctness of this particular echo observable,
or closure in the chosen polynomial Koopman library.

## A laboratory-accessible missing observable

The most direct additional observable is bubble radius `R(t)`, obtained as a
time-resolved wall position. For the RPE state, `h_R(x)=r` and
`L_f h_R(x)=v`, so the radius observation map `(r,L_f r)=(r,v)` is exactly the
state in the ideal model. In practice, a state observer should estimate velocity
from a filtered radius history rather than naively differentiating noisy frames.

For a transparent water experiment, laser/Mie scattering with a photodiode or
high-speed shadowgraphy are established routes. Gaitan et al. measured
experimental radius-time curves for a single acoustically levitated bubble by
laser scattering ([JASA 91, 3166-3183, 1992](https://doi.org/10.1121/1.402855)).
For settings without optical access, active ultrasonic pulse-echo wall tracking
is a distinct measurement from passive far-field emission: Aglyamov et al.
used the reflected pulse to locate the wall quantitatively
([Optics Letters 33, 1357-1359, 2008](https://doi.org/10.1364/OL.33.001357)),
and Ghasemian, Reuter and Ohl demonstrated high-speed plane-wave imaging up to
74 kHz with simultaneous optical validation
([Physics of Fluids 33, 092006, 2021](https://doi.org/10.1063/5.0062260)).

For this 10-micrometre, 100-kHz-scale setup, a calibrated laser-scattering
photodiode is the most natural first laboratory design when optical access is
available. A multi-angle or dual-wavelength calibration would reduce the risk
that the optical transfer curve introduces another size ambiguity. Active
pulse-echo or ultrafast plane-wave imaging is the more realistic extension in
opaque media.

## Reproduction and claim boundary

The portable job is
`examples/native/case-studies/rpe-echo-control-prior-spectrum-tau40.json`.
It produced 40,001 saved samples, 40,386 accepted and 57 rejected internal BDF
steps. The immutable solution is
`archives/solutions/rpe/rpe-echo-control-r-hdmd-eigenvalues-exact-tau40__20260914T163042Z.npz`
with SHA-256
`bb54727d3e70d7235c990ba79e1944671c684fc027839adb0943a85dbc097e34`.
The machine-readable metrics and observability audit are in
`viewer_outputs/rpe-echo-tau40-stability-analysis.json`.

