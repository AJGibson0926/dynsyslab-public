# Bubble models, model-owned observables, and Viewer support

Date: 2026-08-30  
Status: implemented and tested  
MATLAB suite role: read-only evidence

## Scope and evidence

This work traces the retained MATLAB cleanup suite without modifying it or
copying its architecture. The controlling evidence is:

- `scripts/dynsyslib/Bubble/Bubble.m`: model selection, parameters,
  nondimensionalization, RPE/Marmottant/Reddy equations, control maps, and
  observable definitions.
- `scripts/dynsyslib/Duffing/Duffing.m`: Duffing equation, parameters, and
  observables.
- `generate_training.m`: Marmottant and Reddy--Szeri training configurations.
- `playground.m`: common RPE and Duffing free/random-forced data recipes.
- `control_dynamics.m`: lower-amplitude Reddy--Szeri control initial state.

MATLAB's Code Analyzer was run read-only on the principal files. Its relevant
findings support treating the suite as evidence rather than a design template:
many unused intermediates, retained `keyboard` breakpoints, and unused Reddy
observable helpers were reported. No MATLAB source was changed.

## Shared bubble nondimensionalization

Let `R0` be equilibrium radius, `p0` ambient pressure, `pV` vapor pressure,
`rhoL` liquid density, `sigma` clean-interface surface tension, `muL` liquid
viscosity, and `K` the polytropic exponent. The MATLAB Bubble class first sets

\[
p_{G0}=p_0-p_V+\frac{2\sigma}{R_0}.
\]

For RPE and Reddy--Szeri,

\[
f_0=\frac{1}{2\pi R_0\sqrt{\rho_L}}
\sqrt{3Kp_{G0}-\frac{2\sigma}{R_0}},\qquad T=\frac1{f_0}.
\]

Marmottant deliberately uses half this frequency and therefore twice the
corresponding period:

\[
f_0^{M}=\frac12\frac{1}{2\pi R_0\sqrt{\rho_L}}
\sqrt{3Kp_{G0}-\frac{2\sigma}{R_0}}.
\]

The shared groups are

\[
M_L=\frac{R_0}{c_L T},\quad M_G=\frac{R_0}{3c_G T},\quad
Eu=\frac{p_0T^2}{\rho_LR_0^2},\quad p_V^*=\frac{p_V}{p_0},
\]

\[
We_{rc}=\frac{2\sigma}{p_0R_0},\quad
Re_{rc}=\frac{4\mu_L}{p_0T},\quad p_{G0}^*=\frac{p_{G0}}{p_0}.
\]

The code treats a zero sound speed as a request to remove the corresponding
compressibility terms, so its reciprocal is set to zero rather than infinity.

## Rayleigh--Plesset model

### Preset and scaling

The retained default is:

| parameter | value |
|---|---:|
| `R0` | `10e-6 m` |
| `K` | `4/3` |
| `p0` | `100e3 Pa` |
| `pV` | `2.33e3 Pa` |
| `rhoL` | `998 kg/m^3` |
| `sigma` | `0.0725 N/m` |
| `muL` | `0.001 Pa s` |
| `cL`, `cG` | disabled (`0`) |

MATLAB and Python agree on the derived values:

| quantity | value |
|---|---:|
| `f0` | `331963.06683132 Hz` |
| `T` | `3.01238330379725e-6 s` |
| `Eu` | `9.09263844588821` |
| `pG0*` | `1.1217` |
| `pV*` | `0.0233` |
| `We_rc` | `0.145` |
| `Re_rc` | `0.0132785226732528` |

### Equation

For `r=R/R0`, `v=dr/dtau`, and dimensionless external pressure `u`, define

\[
T_1=\frac32v^2,\quad T_2=Eu\,p_{G0}^*r^{-3K},\quad
LP=Eu\,We_{rc}/r,\quad VC=Eu\,Re_{rc}v/r,
\]

\[
Z=T_2-T_1+Eu\,p_V^*-(Eu+LP+VC).
\]

Then

\[
r'=v,\qquad v'=Z/r-Eu\,u/r.
\]

The pressure-input sign is the established external-pressure convention:
positive ambient/acoustic pressure compresses the bubble.

### Typical retained training recipes

- `playground.m` free trajectory: `tau in [0,10]`, `dtau=0.001`.
- `playground.m` random-forced ensemble: seed `2141442`, 10 trajectories,
  `r(0)=1+0.4 randn`, `v(0)=randn`, `tau in [0,5]`, `dtau=0.01`, and the
  suite's random controller.
- `generate_training.m` also retains the explicit free IC
  `(0.374736182055959, 0)`.
- The Lauterborn presets use `(0.9,0)`, `tau in [0,150]`, `dtau=0.001`, with
  the stored sinusoidal pressure settings.

### Observable catalog

The catalog keys are `radius`, `radial_velocity`, `state`,
`inertial_kinetic_energy`, `internal_gas_pressure`, `laplace_pressure`,
`viscous_damping`, `momentum_balance`, and `far_field_acoustic_echo`.
`r`, `R`, `rdot`, and `Rdot` are aliases.

The registered RPE far-field echo is the intrinsic, distance-scaled spherical
monopole: rho times (2 R Rdot squared + R squared Rddot), with zero-input
acceleration. In native coordinates it is S(2 r v^2 + r^2 a), S=R0 P*.
Echo derivatives are with respect to model time tau. Older assets with the same
name may lack one radius factor; compare the serialized formula and basis
identity before reuse. See the correction and provenance inventory (private development record `RPE_CORRECTED_MONOPOLE_RINGDOWN_2026-11-05.html`).

## Marmottant model

### State, preset, and hysteresis

The model remains exactly two-dimensional:

\[
x=(r,v)=(R/R_0,dr/d\tau).
\]

Shell rupture/reformation is a discrete hysteresis latch, as in `Bubble.m`;
it is not an added continuous state. DynSysLab exposes `reset_hysteresis()` so
separate runs start reproducibly while preserving the original hybrid model.

The retained default is:

| parameter | value |
|---|---:|
| `R0`, `Rbuck` | `0.975e-6 m` |
| `K` | `1.095` |
| `p0`, `pV` | `101e3 Pa`, `0 Pa` |
| `rhoL` | `1000 kg/m^3` |
| `sigma`, `muL` | `0.0725 N/m`, `0.001 Pa s` |
| `cL`, `cG` | `1480 m/s`, disabled |
| shell viscosity `Ks` | `15e-9 kg/s` |
| shell elasticity `xi` | `1 N/m` |
| critical/break tension | `1 N/m` |

After choosing `T` with the half-frequency rule, the shell equilibrium values
are recomputed:

\[
\sigma_0=\xi(R_0^2/R_{buck}^2-1),\quad
p_{G0}=p_0+2\sigma_0/R_0,
\]

\[
Res_{rc}=\frac{4K_s}{p_0R_0T},\quad
r_{rupt}=\frac{R_{buck}}{R_0}\sqrt{1+\sigma/\xi},\quad
r_{break}=\frac{R_{buck}}{R_0}\sqrt{1+\sigma_c/\xi}.
\]

MATLAB and Python agree on:

| quantity | value |
|---|---:|
| `f0` | `2115157.18279639 Hz` |
| `T` | `4.72778102797036e-7 s` |
| `ML` | `0.00139343125218005` |
| `Eu` | `23.7479895678326` |
| `pG0*` | `1` |
| `We_rc` | `1.47245493780147` |
| `Re_rc` | `0.0837686012988668` |
| `Res_rc` | `1.28874771229026` |
| `r_rupt`, `r_break` | `1.0356157588604`, `1.4142135623731` |

For the default linear shell law, surface tension is zero below buckling,
`xi(r^2/r_buck^2-1)` in the elastic regime, and the clean-interface value
after rupture. Between rupture and break radii it follows the hysteresis latch.
The exact switching inequalities from MATLAB are retained, including its
single equality branch at `r=r_break`.

### Equation

Define

\[
T_1=-\frac32v^2,\quad
T_2=(1-3KM_Lv)Eu\,p_{G0}^*r^{-3K},
\]

\[
LP=Eu\frac{2\sigma(r)}{p_0R_0}\frac1r,\quad
VC=Eu\,Re_{rc}\frac vr,\quad
VC_s=Eu\,Res_{rc}\frac{v}{r^2},
\]

\[
Z=T_1+T_2-(Eu+LP+VC+VC_s).
\]

Then

\[
r'=v,\qquad v'=Z/r-Eu\,u/r.
\]

No equation or continuous dimension was added or changed.

### Training and slow-manifold regime

`generate_training.m` fixes seed `2141442` and creates two 100-trajectory
ensembles:

- `r(0)` uniform on `[0,1]`, `v(0)` uniform on `[-10,10]`;
- `r(0)` uniform on `[0,1.5]`, reusing the same velocity range.

Both use `tau in [0,2]` and `dtau=0.001`. This broad, fast-velocity launch
followed by relaxation is the retained training setup that exposes the model's
slow attracting manifold. A separate commented/configuration experiment uses
`R0=0.82 micrometers`, `Ks=7.2e-9`, critical tension `0.13`, IC `(1,0)`, and
`tau in [0,600]`; its `200 kPa, 2 MHz` forcing call is commented out.

### Observable catalog

The catalog keys are `radius`, `radial_velocity`, `state`, `surface_tension`,
`shell_viscosity`, `inertial_kinetic_energy`, `internal_gas_pressure`,
`laplace_pressure`, `viscous_damping`, and `momentum_balance`.

## Reddy--Szeri model

### State and preset

The retained state is 12-dimensional:

\[
(r,r',x,x',a_2,a_2',a_3,a_3',a_4,a_4',a_5,a_5').
\]

`r` is the spherical radius mode, `x` is translation, and `a_2` through `a_5`
are shape modes. The default cleanup setup uses the RPE fluid preset with
`R0=2e-6 m`, translation enabled, acoustic amplitude zero, and
`fd=850.75 kHz`. Since the retained default has both sound speeds disabled,
`ML=MG=0`; compressibility-loss and traveling-wave pressure terms therefore
vanish unless the user supplies nonzero sound speeds.

Derived defaults are `f0=1.96444978751307 MHz`,
`T=5.09048389201115e-7 s`, `Eu=6.49123904178984`,
`pG0*=1.7017`, `We_rc=0.725`, `Re_rc=0.0785779915005227`, and
`omega=2.72107738974086`.

The typical training IC is

```
(0.524736182055959, 0, 0, 0.1, 0.1, 0,
 0.1, 0, 0.1, 0, 0.1, 0)
```

with `tau in [0,10]` and `dtau=0.001`. `control_dynamics.m` retains the lower
amplitude alternative `(1,0,0,0.01,0.01,0,...,0.01,0)`.

### Spherical and translation equations

For `R=r`, `V=r'`, the EOM switches the effective `K` to one when `R>1` and
defines

\[
t_1=-\frac32V^2,
\]

\[
t_2=Eu\,p_{G0}^*(1-3KM_LV)R^{-3K},
\]

\[
t_3=3K,Eu\,p_{G0}^*M_GV
\left(\frac32(K+1)M_LV-1\right)R^{-\frac32(K+1)},
\]

\[
lp=Eu\,We_{rc}R^{-1}(M_LV-1),\qquad
vc=Eu\,Re_{rc}(V/R)(1-M_LV),
\]

\[
D=R+EuM_LRe_{rc}+3K,Eu\,p_{G0}^*M_GM_LR^{-(3K+1)/2},
\]

\[
R''=(t_1+t_2+t_3-Eu-vc+lp)/D-Eu\,u/b(R).
\]

When translation is enabled,

\[
v_t=x'+P_a^*\sin(\omega\tau),
\]

\[
x''=-\left(3R'/R+\frac92EuRe_{rc}/R^2\right)v_t
      +3\dot P_a^*\cos(\omega\tau).
\]

### Shape-mode equations

With `c=3R'/R`, the retained coefficient functions are

\[
k_n=\frac12EuWe_{rc}\frac{(n-1)(n+1)(n+2)}{R^3}
 -(n-1)\frac{R''}{R}
 -\frac92\frac{n(n+1)^2}{(2n-1)(2n+3)}(v_t/R)^2,
\]

\[
d_n=\frac{9n(n-2)(n-1)(n+1)}{4(2n-3)(2n-1)}(v_t/R)^2,
\]

\[
f_n=\frac{3(n+1)[2(n-1)v_tR'+nR v_t']}{2(2n-1)R^2},
\quad g_n=\frac32(n+1)v_t/R,
\]

\[
h_n=\frac{3(n+1)^2[-2v_tR'+Rv_t']}{2(2n+3)R^2},
\quad
j_n=\frac{9(n+1)^2(n+2)(n+3)}{4(2n+3)(2n+5)}(v_t/R)^2.
\]

The four second-order shape equations are

\[
a_2'=\dot a_2,\quad
\ddot a_2=-k_2a_2-c\dot a_2+h_2a_3+j_2a_4-\frac94v_t^2/R,
\]

\[
a_3'=\dot a_3,\quad
\ddot a_3=f_3a_2+g_3\dot a_2-k_3a_3-c\dot a_3+h_3a_4+j_3a_5,
\]

\[
a_4'=\dot a_4,\quad
\ddot a_4=d_4a_2+f_4a_3+g_4\dot a_3-k_4a_4-c\dot a_4+h_4a_5,
\]

\[
a_5'=\dot a_5,\quad
\ddot a_5=d_5a_3+f_5a_4+g_5\dot a_4-k_5a_5-c\dot a_5.
\]

The state-dependent input influence on shape acceleration arises through the
controlled `R''` inside `k_n`; this is why the full analytical input map has
nonzero rows for the shape accelerations even though the EOM adds the explicit
input only to the spherical acceleration before computing `k_n`.

### Observable catalog

The guaranteed entries are all 12 state components plus `state`. The retained
Reddy observables add `compressibility_loss`, `radial_acceleration`, `c`,
`c_time_derivative`, `c_radius_derivative`, `c_velocity_derivative`, `k2`,
`k3`, `k4`, `k5`, `k5_time_derivative`, `k5_radius_derivative`,
`k5_velocity_derivative`, and `g1` (`r rdot`).

The MATLAB observable formulas intentionally use the base `K` and spherical
unforced acceleration even though the EOM switches `K` for `R>1`; that
provenance is retained rather than conflated with the EOM. Unused helper
expressions in `Bubble.m` that reference undefined local `v`/`vprime` are not
published as catalog observables.

## Duffing observables

For DynSysLab's equation

\[
x'=v,\qquad v'=-\delta v+\alpha x-\beta x^3+u,
\]

the MATLAB parameter mapping is `beta_MATLAB=-alpha` and
`alpha_MATLAB=beta`. The model-owned catalog contains `position`, `velocity`,
`state`,

\[
KE=\tfrac12v^2,\qquad
PE=-\tfrac12\alpha x^2+\tfrac14\beta x^4,\qquad TE=KE+PE,
\]

and the retained complex random observation

\[
h(x,v)=e^{ix}/\cosh(v).
\]

Aliases `x` and `xdot` select the first two state observables.

## Python API and Viewer

Every `ODESystem` now guarantees state-component observables and `state`.
Systems add domain observables through a small model-owned catalog:

```python
system = load_system("marmottant", preset="matlab_cleanup_linear")
catalog = system.observable_catalog("dimensionless")
momentum = catalog["momentum_balance"]
values = momentum.evaluate(states)       # shape: (samples, 1)
```

Viewer evaluates scalar observables read-only from an existing trajectory
archive:

```python
Viewer().add(PlotRequest(
    "observable_time", ArchiveSource(solution),
    observable="momentum_balance",
))

Viewer().add(PlotRequest(
    "observable_heatmap", ArchiveSource(solution),
    observable="momentum_balance",
    domain=((0.7, 1.5), (-5, 5)),
    domain_density=151,
))
```

Time plots evaluate the observable at every supplied trajectory sample. Heat
maps are limited to two-state systems and evaluate the same observable on the
requested phase-plane grid. If no domain is supplied, Viewer uses the pooled
trajectory bounds plus a five-percent margin. Complex scalar observables
support `real`, `imaginary`, `magnitude`, and `phase`; `auto` displays
magnitude. Both Plotly and dependency-free SVG renderers support the heat map.

## Code boundaries

- `dynsyslib/observables.py`: reusable `Observable` and `ObservableCatalog`.
- `dynsyslib/interfaces.py`: generic state/full-state observable guarantee.
- `dynsyslib/lib/rpe.py`, `duffing.py`: model-specific catalogs.
- `dynsyslib/lib/marmottant.py`: two-state model plus discrete hysteresis latch.
- `dynsyslib/lib/reddy_szeri.py`: 12-state Reddy--Szeri model.
- `viewer/core.py`: read-only observable display transforms.
- `viewer/sources.py`: model-aware axis labels.

The MATLAB suite remains a read-only source. These additions do not import,
call, or depend on MATLAB at runtime.
