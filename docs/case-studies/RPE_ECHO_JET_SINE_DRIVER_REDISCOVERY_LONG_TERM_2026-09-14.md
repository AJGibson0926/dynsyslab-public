# RPE echo-jet KLQR sine-driver rediscovery: complete pipeline and long-term result

Date: 2026-09-14  
System: dimensionless Rayleigh-Plesset equation (RPE), air bubble in water  
Control objective: recover the unknown pressure input that generated a prescribed
100 kHz periodic trajectory  
Final validation duration: `300 tau` = `903.715 microseconds` = `90.37` forcing cycles

## Executive result

The repaired controller succeeds in the simulation criterion that matters here:
after acquisition, it reproduces the correct long-term state orbit and rediscovers
the sine-wave pressure that generated that orbit.

The controller was **not given the pressure waveform**. Its reference contains
only the archived target state trajectory, mapped into the two measured-coordinate
variables `[echo, autonomous echo-dot]`. The original target archive's
`control_values` array was withheld from the controller and used only afterward
to score the recovered input.

Over the settled interval `100 <= tau <= 300`:

| Quantity | Result |
|---|---:|
| Known pressure | `20.000 kPa` at `100 kHz` |
| Recovered fundamental amplitude | `20.215 kPa` |
| Amplitude ratio | `1.01076` |
| Phase error | `-0.00418 deg` |
| Direct driver NRMSE | `2.1158%` |
| Driver correlation | `0.999838` |
| One-tone fit R-squared | `0.999676` |
| Radius RMSE | `2.199e-5 R0` = `0.220 nm` |
| Radial-velocity RMSE | `3.714e-4 R0/T` = `1.233 mm/s` |
| Echo NRMSE | `0.0640%` |
| Autonomous echo-dot NRMSE | `0.5520%` |
| Controlled radius range | `0.952782` to `1.064144 R0` |
| Largest settled pressure magnitude | `20.468 kPa` |

The cycle-by-cycle driver NRMSE for all 59 complete cycles after `tau=100`
lies between `2.1050%` and `2.1309%`; its median is `2.1170%`. The final
complete cycle has `2.1239%` NRMSE and `0.999837` correlation. There is no
late loss of lock or secular radius growth.

![Long-term sine-driver recovery](../../viewer_outputs/rpe-echo-jet2-sine-long-term-control.png)

## Pipeline at a glance

```text
large free RPE trajectory (r, v)
            |
            +--> radius H-DMD spectrum (2,000 delay modes)
            |                |
            |                v
            +--> degree-4 state-coordinate SINDy eigenfunctions
                             |
                             v
                  compact one-pair state KROM
                             |
                             +--> state-KROM KLQR follows a moderate
                             |    multisine state reference
                             |                |
                             |                v
                             |    forced branch-covering training states
                             |                |
                             |                v
                             +--> degree-4 SINDy refit in
                                  H(x) = [echo, L_f0 echo]
                                                |
                                                v
                                  repaired echo-jet KROM
                                                |
held-out 100 kHz sine-driven target state ------+
                                                |
                                                v
                                  exact-map echo-jet KLQR
                                                |
                                                v
                                  recovered pressure waveform
```

There are two distinct training datasets. The first identifies a useful Koopman
spectrum and state-coordinate eigenfunction. The second samples the moderate
state-space branch on which the echo jet is locally informative and refits that
same spectral mode as a function of echo and echo-dot. The sine trajectory is
neither one of those training sets; it is the held-out control reference.

## 1. Plant, units, and actuation

The dynamical state is

\[
  x=(r,v), \qquad r=R/R_0, \qquad v=dr/d\tau,
\]

with equilibrium radius `R0=10 micrometers`. The water parameters are
`rho=998 kg/m^3`, `sigma=0.0725 N/m`, `mu=0.001 Pa s`, vapor pressure
`2330 Pa`, ambient pressure `p0=100000 Pa`, and polytropic exponent
`kappa=4/3`. One model-time unit is

\[
  T=3.0123833038\ \mathrm{microseconds},
\]

and one model-input unit is `p0=100 kPa`.

Using the archived dimensionless groups

```text
G = 10.1992125448       gas-pressure ratio
V = 0.2118584758        vapor-pressure ratio
P = 9.0926384459        ambient-pressure ratio
C = 1.3184325747        capillarity
D = 0.1207368058        viscosity
```

the controlled RPE used by the native solver is

\[
\begin{aligned}
  \dot r &= v,\\
  \dot v &= \frac{G r^{-4}+V-P-C/r-Dv/r-\tfrac32v^2-Pu}{r}.
\end{aligned}
\]

The negative input sign means a positive acoustic-pressure command compresses
the bubble. The state-dependent input vector is therefore

\[
  g(x)=\begin{bmatrix}0 & -P/r\end{bmatrix}^{\mathsf T}.
\]

## 2. First training dataset: a large free transient

The upstream identification starts with the immutable free-response archive

```text
archives/solutions/rayleigh-plesset/
rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z.npz
```

SHA-256: `e4529d55dd549e12c66dbb78c9521c6bc273aa8bb135637319568a2edfdc9eba`

It contains 25,001 dimensionless samples over `0 <= tau <= 25` at `dtau=0.001`
(`75.310 microseconds` physical time), starting from `[2.3333866, 0]`. It is
unforced. The sampled state ranges are:

```text
r:  0.232721 to 2.333387
v: -27.0932   to 26.9480
```

This violent transient is intentionally broad. It supplies rich spectral data,
but it later proves too broad for a single well-conditioned polynomial inverse
from the echo jet to physical state.

## 3. H-DMD supplies the spectral candidates

Hankel DMD is applied to the scalar radius observable (`state index 0`), not to
echo. Its recorded settings are:

| Setting | Value |
|---|---:|
| Analysis samples | `22,000` |
| Sampling window | `20,000` samples |
| Delay coordinates | `2,000` |
| Delay stride | `1` |
| Hankel snapshot pair | shifted by one `dtau=0.001` sample |
| Hankel matrix | `20,000 x 2,000` |
| Retained rank | `2,000` |
| SVD | LAPACK `gesvd`, tolerance `1e-10` |
| DMD form | projected (`exact=false`) |
| Mode ordering | terminal-amplitude energy, descending |

The complete spectrum archive is

```text
archives/derived/koopman-hankel-dmd/
rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__
koopman-hankel-dmd__20260826T221322.756123Z__604a65f1.npz
```

SHA-256: `c2fabe83b37b94ace1abb5b7d4fdd11cad14dee97c476d2fbb00e5e63eb0bf4b`

This 790.5 MB artifact is a complete-spectrum intermediate: it stores the
`2000 x 2000` operator and `20000 x 2000` modes. It is not the compact KROM
used online.

## 4. State-coordinate SINDy turns the spectrum into eigenfunctions

The next stage searches for polynomial eigenfunctions in the physical state
coordinates `(r,v)`. It uses all 14 nonconstant monomials through total degree
four:

\[
  \Theta(r,v)=\{r,v,r^2,rv,v^2,\ldots,r^4,r^3v,r^2v^2,rv^3,v^4\}.
\]

For each H-DMD eigenvalue `mu`, the code converts to the continuous generator
eigenvalue

\[
  \lambda=\log(\mu)/\Delta\tau
\]

and finds a sparse coefficient vector `c` such that

\[
  \Gamma c \simeq \lambda\Theta c,
  \qquad
  \Gamma=\nabla\Theta\, f_0.
\]

The fit uses analytic RPE right-hand-side derivatives (`gamma_method=system_rhs`),
degree four, no constant term, threshold `1e-6`, one thresholding iteration,
no coordinate normalization, and conjugate-fit reuse. The full coefficient-first
archive contains 2,000 coefficient rows without materializing a
`25001 x 2000` eigenfunction-value array:

```text
archives/derived/koopman-with-sindy/
rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__
koopman-with-sindy__20260827T183645.394109Z__c4dde209.npz
```

SHA-256: `6e3725beac327fa749eeb0b0250cefa816d158cbbbc2d6fdea0d1b6bfa979119`

## 5. Compact parent KROM extraction

The extraction requests energy-sorted source mode index `1` and includes its
conjugate dual. The selected discrete pair is

\[
  \mu=0.99977779858\pm0.00481760114i,
\]

which, at `dtau=0.001`, corresponds to

\[
  \lambda=-0.21061644578\pm4.81863456156i.
\]

The resulting state-coordinate KROM retains one independent complex
eigenfunction plus its conjugate, their degree-four coefficient vectors, the
25,001 training states and evaluated training eigenfunctions. It supports exact
analytic input coupling but has no state reconstruction map:

```text
archives/derived/koopman-krom-extract/
rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__
koopman-krom-extract__20260827T183646.410439Z__0dea0923.npz
```

SHA-256: `8c3f221f563a82271e2995c8cc7e627b37f28b54f9feb93c59ae28050ae4ef86`

This is the compact state-basis parent used to generate the second training set.

## 6. Second training dataset: controlled coverage of the useful branch

The instantaneous echo map is non-injective over the large free transient, and
the two-coordinate echo jet is nearly singular at parts of that domain. Fitting
one global degree-four echo-jet eigenfunction there gave a normalized generator
defect of approximately 100%.

The repair therefore samples the moderate branch on which control is actually
desired. The compact state KROM drives the RPE from `[1,0]` along the canonical
multisine state target

\[
\begin{aligned}
  r_{ref}(\tau)&=1+0.09\sin(4\tau)-0.12\sin(3\tau),\\
  v_{ref}(\tau)&=0.36\cos(4\tau)-0.36\cos(3\tau).
\end{aligned}
\]

That state-KROM KLQR uses the exact lifted input map, `Q=5000 I`, `R=1e-4`,
continuous-RHS control, BDF integration, tolerances `1e-9`, maximum step
`0.001`, and produces 10,001 samples over `0 <= tau <= 10`:

```text
archives/solutions/rpe/
rpe-echo-study-state-basis-baseline__20260911T104742Z.npz
```

SHA-256: `d3e06fde49f7f3d704621fbb7443f83dd984ca6d7aa5a58eec3959eb73bfa40e`

The observed state range is `r=0.7975..1.2022` and
`v=-0.6444..0.7114`. This dataset is forced, but it is used only to select
states on the relevant branch. At each sampled state, the next fit evaluates
the **autonomous** RPE generator `f0`; it does not learn the preceding
controller's closed-loop vector field.

## 7. Refit the eigenfunction in physically accessible echo coordinates

The coordinate map is the two-component autonomous echo jet

\[
  H(x)=\begin{bmatrix}e(x)\\e_1(x)\end{bmatrix}
      =\begin{bmatrix}h(x)\\L_{f_0}h(x)\end{bmatrix}.
\]

The native far-field echo expression is

\[
  e=R_0P_s r\left(2v^2+\frac{Z}{r}\right),
\]

where

\[
  Z=G r^{-4}+V-P-C/r-Dv/r-\tfrac32v^2
\]

and `Ps=10997.9078785 Pa`. The second coordinate is obtained by exact symbolic
Lie differentiation through the autonomous RPE field. Under applied control,
the derivative measured from echo history is instead

\[
  de/d\tau=L_{f_0}h(x)+L_gh(x)u;
\]

recovering the autonomous coordinate in a laboratory therefore requires the
known applied input and an input-aware observer or derivative estimator.

The two coordinates are standardized with

```text
center = [ 0.1472323322, -0.0894546951 ]
scale  = [ 0.5530497264,  1.9641872085 ]
```

and the same 14-term, degree-four polynomial dictionary is built in these two
standardized variables. The fit reuses the parent spectral pair, uses analytic
coordinate Jacobians and system-RHS generator derivatives, disables input terms
(`include_input_terms=false`), and saves the coordinate velocities and
Jacobians used as evidence.

The repaired archive is

```text
archives/derived/koopman-with-sindy/
rpe-echo-study-state-basis-baseline__20260911T104742Z__
koopman-with-sindy__20260914T165124Z.npz
```

SHA-256: `e3e6301ed6ec6cb91e46540e20369bd47dc52cccfcb4502d78d741128ca30806`

Its local diagnostics are:

| Diagnostic | Repaired value |
|---|---:|
| Standardized degree-four feature condition | `71.04` |
| Minimum standardized observation-Jacobian singular value | `0.3307` |
| Median observation-Jacobian condition | `1.451` |
| Maximum observation-Jacobian condition | `2.786` |
| Generator residual RMS | `0.86965` |
| Residual / generator RMS | `3.7789%` |
| Generator-target alignment | `0.999288` |

This is the central observability repair: the same physical echo and echo-dot
coordinates that failed globally are locally informative and admit a good
eigenfunction fit on the moderate branch.

## 8. Held-out sine target

The target trajectory was generated separately from `[1,0]` by the prescribed
physical pressure

\[
  p_c(t)=20\ \mathrm{kPa}\,\sin(2\pi\,100\ \mathrm{kHz}\,t+\pi/2),
\]

or, in model coordinates,

\[
  u_*(\tau)=0.2\cos(1.8927362514011978\tau).
\]

Its period is `3.319630668 tau` (`10 microseconds`). The immutable target is

```text
archives/solutions/rpe/
rpe-periodic-100khz-20kpa-phase-portrait__20260911T214426Z.npz
```

SHA-256: `5c6aa54092618d6da58df4fcd384527e6a257b1b2f1b2af2c570c465ba67ec39`

The controller reads only its `time` and `values` arrays. At run time those
state values are mapped through `H` and then through the learned polynomial
eigenfunction to construct the moving lifted reference. The archived known
pressure is not part of that path.

## 9. KLQR construction

Only mode index `0` is selected from the conjugate pair. One complex
eigenfunction already becomes two real controller coordinates; selecting both
members separately would duplicate the same real subspace. For
`lambda=alpha+i beta`, the autonomous real lifted dynamics are

\[
  \frac{d}{d\tau}
  \begin{bmatrix}\Re\psi\\\Im\psi\end{bmatrix}
  =
  \begin{bmatrix}\alpha&-\beta\\\beta&\alpha\end{bmatrix}
  \begin{bmatrix}\Re\psi\\\Im\psi\end{bmatrix}
  +B_\psi(x)u,
\]

with `alpha=-0.2106164458` and `beta=4.8186345616`.

The controller uses the exact state-dependent lifted input map

\[
  B_\psi(x)=J_\phi(H(x))\,J_H(x)\,g(x),
\]

where both Jacobians are analytic. At each continuous RHS evaluation the native
KLQR path forms the local continuous-time Riccati solution for
`Q=5000 I`, `R=1e-4` and applies

\[
  u=-K(x)\left(z(x)-z_{ref}(\tau)\right).
\]

There is no feedforward copy of `u_*`, no input cap, and no pressure data from
the target archive in the feedback law. The final job is

```text
examples/native/case-studies/
rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300.json
```

It runs native compiled BDF from `tau=0` to `300`, publishes every `0.02`, uses
`rtol=atol=1e-9`, and restricts the internal step to `0.001`.

## 10. Long-term validation

The final immutable output is

```text
archives/solutions/rpe/
rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300__
20260914T171839Z.npz
```

SHA-256: `9f1611901c4a678308bf4e43263045f2b2f3406d9c93b1f9a67924c3e5461cda`

The important trend is convergence to a time-invariant periodic error, not a
short transient score:

| Window | Driver NRMSE | Correlation | Fundamental amplitude | Phase error | Radius RMSE |
|---|---:|---:|---:|---:|---:|
| `20..300` | `2.4528%` | `0.999761` | `0.202120` | `+0.00067 deg` | `2.201e-5` |
| `30..300` | `2.2400%` | `0.999811` | `0.202143` | `-0.00260 deg` | `2.200e-5` |
| `60..300` | `2.1193%` | `0.999837` | `0.202150` | `-0.00351 deg` | `2.198e-5` |
| `100..300` | `2.1158%` | `0.999838` | `0.202152` | `-0.00418 deg` | `2.199e-5` |
| `200..300` | `2.1108%` | `0.999839` | `0.202153` | `-0.00454 deg` | `2.199e-5` |
| `290..300` | `2.1201%` | `0.999838` | `0.202159` | `-0.00171 deg` | `2.197e-5` |

The near-identical `60..300`, `100..300`, `200..300`, and terminal-window
scores show that the result has settled. The small remaining waveform error is
mostly a roughly 1.08% amplitude excess plus weak non-fundamental content; it is
not phase drift.

## 11. Acquisition caveat

The full run remains numerically bounded (`r=0.91597..1.098997`), but the
unconstrained high-gain controller requests a brief early peak of
`28.606` model units, or `2.861 MPa`. After `tau=100`, the largest magnitude is
only `0.20468` model units (`20.468 kPa`).

This means the result establishes local long-term control and driver
rediscovery in simulation. It does **not** establish actuator-feasible capture
from startup. A physical design needs saturation-aware acquisition, preferably
hard-constrained KMPC or a scheduled handoff into local KLQR.

## 12. What was tried and rejected

Three negative branches are scientifically important:

1. **Instantaneous echo only.** The scalar echo eigenfunction tracked echo while
   the radius escaped to `21.33 R0`; its pressure NRMSE was `744%`. That was a
   genuine representation failure, not a tuning failure.
2. **Global echo-plus-echo-dot fit.** Adding the derivative did not help when a
   single polynomial was fit over the violent free transient. The map was nearly
   singular and the generator defect stayed near 100%.
3. **Sine-trajectory-specific refit.** A separate refit directly on the sine
   orbit had a larger relative generator defect (`8.72%`) and strongly varying
   input authority. Weight sweeps either failed the Riccati solve or produced
   poor driver recovery. It was rejected. The final successful controller uses
   the broader multisine-trained repaired KROM instead, which is also the cleaner
   held-out validation design.

During the rejected sine-specific branch, sourcing an older compact spectrum
archive without an explicit `continuous_eigenvalues` array caused its discrete
eigenvalues to be reinterpreted with the sine archive's `dtau=0.02` rather than
the original `0.001`. The preserved probe job documents this archive-compatibility
hazard. The accepted echo-jet archive contains explicit continuous eigenvalues,
and the successful control run uses those values directly.

## 13. Scientific interpretation and remaining boundary

The result supports a precise, local conclusion:

- instantaneous echo alone is not state identifying;
- `[echo, autonomous echo-dot]` is locally full-rank and well-conditioned on the
  moderate branch;
- a polynomial eigenfunction in that repaired coordinate chart can support
  long-term KLQR tracking and recover the generating pressure on a held-out
  sinusoidal target;
- this does not create a globally injective echo representation.

It is also not yet a laboratory sensor-only controller. Although the learned
eigenfunction depends only on the two echo-jet variables, the `exact` control
map evaluates `B_psi(x)` from the simulated physical state. A laboratory
implementation must estimate state or input-corrected autonomous echo-dot from
echo history plus known actuation, validate the estimator on noisy held-out
data, and monitor local observability before applying the same feedback law.

## 14. Reproduction artifacts

- Final job: `examples/native/case-studies/rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300.json`
- Final output: `archives/solutions/rpe/rpe-echo-jet2-repaired-klqr-rediscover-100khz-20kpa-tau300__20260914T171839Z.npz`
- Machine-readable metrics: `viewer_outputs/rpe-echo-jet2-sine-long-term-control.json`
- Figure: `viewer_outputs/rpe-echo-jet2-sine-long-term-control.png`
- Deterministic analysis script: `tools/analyze_rpe_echo_jet_sine.py`
- Echo-jet fit job: `examples/native/case-studies/rpe-echo-jet2-sindy-target-local.json`
- Rejected sine-local fit: `examples/native/case-studies/rpe-echo-jet2-sindy-sine-target-local.json`
- Preserved source-`dt` probe: `examples/native/case-studies/rpe-echo-jet2-sindy-sine-target-local-source-dt-probe.json`
- Earlier observability analysis: `docs/case-studies/RPE_ECHO_JET_OBSERVABILITY_REPAIR_2026-09-14.md`

Run the analysis from the repository root with:

```powershell
python tools/analyze_rpe_echo_jet_sine.py
```

The script checks that target and controlled archives share the same time grid,
re-evaluates the native RPE echo and autonomous echo-dot expressions, verifies
the mapped target, computes window and per-cycle metrics, hashes every pipeline
artifact, and regenerates the JSON evidence and figure.
