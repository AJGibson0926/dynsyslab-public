# Native observable expressions, echo-only SINDy and control targets

The compiled core represents an observable as y=h(x) and a fitted eigenfunction
as psi(x)=phi(h(x)). H-DMD's scalar measurement and SINDy's polynomial coordinate
selection are independent. For an echo-only fit, select far_field_acoustic_echo
for both. H-DMD uses delays to estimate a spectrum; the fitted phi is still a
function of instantaneous echo alone. No delayed samples, radius or velocity
are additional polynomial variables.

## Native derivative API

AnalyticExpression is an immutable real expression graph supporting constants,
state components, arithmetic, constant real powers, exp, log, sin, cos and tanh.
Observable::fromExpressions creates both the value evaluator and exact Jacobian.
Observable::partialDerivative(axis) returns another expression observable;
Observable::flowDerivative(field, order) repeatedly applies L_f h = J_h f.
Differentiating that returned observable gives mixed spatial/flow derivatives.
These are symbolic chain-rule derivatives evaluated in floating point, without
finite-difference steps. An arbitrary opaque callback cannot be automatically
converted into an expression; it retains its supplied Jacobian or explicit
finite-difference fallback.

ODESystem::driftExpressions() describes autonomous zero-input drift in the
system's native coordinates. RPE now supplies this field and analytic expressions
for its pressure, inertia, viscous, momentum and echo observables. Systems without
an expression field reject analytic flow-derivative requests. This is not a total
derivative through an arbitrary time-dependent input or feedback law; such a
field and any extra time/input variables must be explicitly represented.

### Duffing force and pendulum torque

The existing `duffing` and `nonlinear_pendulum` systems expose force/torque through
the same `Observable::fromExpressions` catalog, exact Jacobian, expression JSON,
SINDy coordinates, reconstruction and lifted input-map machinery. They also supply
analytic drift/input fields for standard repeated flow-derivative queries. No new
observable class, solver or controller is involved.

- `duffing`: `force` (alias `force_dimensionless`) is
  `alpha*x - beta*x^3 - delta*v`, in dimensionless units. It equals the second
  zero-input RHS component in this normalized unit-mass model. The double-well
  case (`alpha > 0`, `beta > 0`) represents the cubic approximation to a magnetic
  pendulum attracted toward two magnets: the restoring contribution is
  `alpha*x - beta*x^3` and the dissipative contribution is `-delta*v`.
  It is not an exact magnet-geometry force law or an isolated magnetic force.
  See [Ma and Zhang, *The periodic solution of a pendulum under magnetic action*](https://doi.org/10.1177/14613484221122093)
  for this Duffing-Holmes reduction. Existing presets and equations are unchanged.
- `nonlinear_pendulum`: `torque` is
  `-mass_kg*gravity_m_s2*length_m*sin(theta) - damping_n_m_s*omega`, in **N m**.
  It equals `I*omega_dot` at zero input, with `I = mass_kg*length_m^2`.
  The acceleration RHS alone is not physical torque.
- `nonlinear_pendulum`: `torque_dimensionless` is `torque/(m*g*L)`.
  Its output is dimensionless, but its input remains the native dimensional
  state `(theta, omega)` with `omega = dtheta/dt` in rad/s. Output normalization
  does not change the solver clock or silently rescale the incoming state.
  The resolved configuration records `torque_scale_n_m = m*g*L` and
  `characteristic_time_s = T = sqrt(L/g)` in the standard scaling metadata.
  With `s=t/T`, `w=dtheta/ds=T*omega`, this gives
  `torque_dimensionless = -sin(theta) - (damping*T/I)*w = d2theta/ds2`
  on free dynamics. Convert dimensionless velocity back with `omega=w/T`
  before evaluating the native catalog. The normalized observable is available
  for positive mass, length and gravity with a finite positive torque scale;
  zero gravity retains physical `torque` but has no `m*g*L` normalization.

These are **intrinsic zero-input state observables**, including damping. An
applied driver remains separate: the full Duffing force is `force + input_gain*u`
in model units, and the full pendulum torque is `torque + u` in N m. Thus they do
not silently change definition in a controlled trajectory. Setting
`include_input_terms=true` for a flow derivative differentiates this intrinsic
observable along the driven dynamics; it does not add the driver to its value.
Use the standard velocity observable with `flow_order=1` and supplied input to
obtain the full acceleration, multiplying by inertia for physical pendulum torque.

The Duffing model does not prescribe a physical mass or length/time calibration.
For a chosen physical displacement `X=Lscale*x`, time `t=Tscale*s` and mass `m`,
physical force is `F=(m*Lscale/Tscale^2)*force`. Do not label the raw dimensionless
value as newtons, or use the input scale as a force scale without accounting for
`input_gain`. For a physical law `F=a*X-b*X^3-c*dX/dt`, the matching coefficients
are `alpha=a*Tscale^2/m`, `beta=b*Lscale^2*Tscale^2/m`, `delta=c*Tscale/m`.

Run `dynsyslab-cli observable examples/native/duffing-force-query.json` or
`dynsyslab-cli observable examples/native/pendulum-torque-query.json`.
For example, the standard observable query accepts:

```json
{"schema":"dynsyslab.observable-query/v1",
 "system":{"name":"nonlinear_pendulum","parameters":{"mass_kg":2,"length_m":0.5}},
 "state":[0.4,0.2],"observable":"torque","flow_order":0}
```

Select `force`, `torque` or `torque_dimensionless` by name in the existing
postprocessing observable settings. Native regression tests check non-unit
physical scalings, sign conventions, driver separation, symbolic derivatives,
and expression serialization.

The echo definition remains the existing state observable, with units Pa m:

    e = R0 * Pscale * r * (2*v^2 + Z/r)
    Z = pg*r^(-3*kappa) - 1.5*v^2 + pv - p0 - cap/r - visc*v/r

Here v=dr/dtau and Z/r is the zero-input acceleration. Thus this state observable
does not silently substitute the controlled acceleration for Z/r.

Use the headless CLI:

    dynsyslab-cli observable examples/native/rpe-echo-derivative-query.json

The public C++ workflow entry is evaluateObservableQuery. The C ABI entry
is dsl_workflow_observable_json, using the existing result handle functions.
The optional dependency-free client exposes NativeWorkflow.observable(query).
All three execute the same compiled implementation. The result reports values,
Jacobian, derivative order, backend and the retained expression definition.
Query partial_axes are applied in order, followed by flow_order (0–8).

Expression contracts are analytic-state-expression/v1 and
dynsyslab.observable-expression/v1. KROM archives retain coordinate_observable_json
alongside names, original training_states, transformed training_coordinates,
normalization and one exponent column per coordinate. These definitions survive
KROM extraction, save/reload and exact lifted input-map evaluation.

## Pipeline and control jobs

In the app's standard H-DMD → SINDy settings, choose far_field_acoustic_echo
as Scalar observable and Coordinate observable. Reconstruction observables may
also be set to far_field_acoustic_echo. The compiled postprocess settings are:

    hdmd:  {"observable":"far_field_acoustic_echo", ...}
    sindy: {"coordinate_observable":"far_field_acoustic_echo",
            "polynomial_degree":4,"include_constant":false,
            "gamma_method":"system_rhs", ...}

An explicit coordinate_expression can replace a named SINDy observable. The
system RHS and observable Jacobian produce ydot=J_h f; supplied eigenvalues are
converted to continuous-time generator eigenvalues. A standalone SINDy stage
may use spectrum_archive to reuse a saved spectrum. The archive records that
spectrum path and the actual coordinate gradient backend.

KLQR and KMPC accept target.space="observable". The callback then supplies the
coordinate dimension, rather than the dynamical-state dimension. For the echo
corresponding to a canonical radius/velocity target, use:

    "target": {
      "space":"observable", "kind":"from_state",
      "state_target": {
        "kind":"rpe_canonical",
        "name":"rpe_multisine_equilibrium_moderate_identification"
      }
    }

The target is evaluated through the exact coordinate observable retained by the
selected KROM. The app provides observable fixed-point, trajectory and waveform
choices for direct coordinate targets, and observable-from-state choices for
mapped references. Its observable canonical preset also maps through h.
State-space controllers reject observable targets. Existing state targets remain
compatible. Q is a cost in the KROM's real split coordinates; changing the fitted
basis changes the meaning of the same numerical Q.

Controlled archives store control_target_space, control_observable_target_values
and control_observable_target_names. A from_state target additionally retains
control_target_values as the state reference. Viewer draws the supplied echo
target on its observable time series; unrelated observables are not labeled as
control objectives. Direct scalar targets do not invent a radius/velocity target.

The exact map uses B_psi(x)=J_phi(h(x))*J_h(x)*B(x). Hence echo-only eigenfunctions
do not imply an echo-only exact-map controller: it still uses the plant state to
evaluate the map. Constant/bilinear maps can be trained in this new basis; a map
from a previous state-coordinate KROM must not be reused as if its coefficients
belonged to the echo basis.

## Legacy sources and scientific limits

Legacy Python-reference bubble trajectories store physical time/radius/velocity
even when integration_form says dimensionless. Native postprocessing now uses
archive coordinate tags and stored scales to convert these inputs explicitly,
records the conversion, and leaves source bytes untouched. Nested preset numeric
parameters are preserved. Legacy controlled inputs need declared units.

A scalar echo generally loses state information. A delayed echo signal can have
an informative H-DMD spectrum even when no nontrivial eigenfunction for that
spectrum can be represented accurately as phi(e) alone. Fit residuals and closed
loop behavior must therefore be assessed independently; success is not guaranteed
by a successful numerical fit. No radius, velocity or delay coordinates are added
to make an echo-only experiment succeed.

## Verification

Native tests cover analytic and repeated flow derivatives, chain-rule input maps,
expression archive round-trip, legacy physical-coordinate import, echo-only fit
shape, mapped and direct scalar control target archives, C ABI queries and Viewer
target identity matching. The MATLAB suite is read-only reference material.


## Echo and autonomous echo-dot coordinates

The RPE catalog provides echo, echo_dot, echo_ddot, echo_dot_dr and echo_dot_dv
aliases. Their canonical names begin far_field_acoustic_echo. The echo_jet2
coordinate preset contains exactly [echo, echo_dot]. With x=(r,v), the defining
expressions are e=h0(x), e1=L_f0 h0, e2=L_f0^2 h0, and the two components of
J_e1. These are autonomous, zero-input observables in the native clock tau;
physical-time derivatives are e1/T and e2/T^2. The spatial derivatives use r
and v=dr/dtau as their independent variables.

The SINDy library is theta(e,e1); its analytic generator library is
Gamma=theta_e*e1+theta_e1*e2 on free dynamics. All mixed monomials up to the
selected total degree are included. Neither physical state component becomes
a polynomial coordinate. Existing H-DMD spectra can be reused through
spectrum_archive; a saved compact KROM can supply its retained eigenvalue pair
without recomputing H-DMD.

Set archive_coordinate_derivatives=true to retain training_coordinate_velocities
(samples by coordinates) and training_coordinate_jacobians (samples by
coordinates by dynamical states). This requires system_rhs and an analytic
coordinate Jacobian. These arrays survive compact-KROM extraction. They record
the actual push-forward derivatives used by the configured source RHS.
Observable::pushForwardRealBatch is the shared native chain-rule API.

### Optional RPE receiver pressures

The default echo is the intrinsic leading monopole source in Pa m. A separate
receiver catalog exposes individual pressure terms in Pa for a stationary source
and an external point receiver. Select one with the normal `coordinate_transform`
workflow on a physical RPE trajectory, for example

```json
{"kind":"observable","observable":"spherical_near_field_pressure",
 "receiver":{"position_m":[0.005,0,0],"sound_speed_m_s":1480},
 "names":["spherical near-field pressure (Pa)"]}
```

Here D is the length of the receiver position vector in metres, n is its unit
direction from the center, c is the sound speed in m/s, R is the physical
radius and dots below use seconds. Each selection returns its own source-time
series. Add the recorded travel time D/c to physical emission timestamps to obtain
receiver timestamps. For native model time τ use τ + D/(c T*). The method
records the receiver and selected term in `acoustic_receiver` metadata.

| Selection | Pressure in Pa |
| --- | --- |
| `monopole_pressure` | Existing intrinsic `far_field_acoustic_echo` divided by D |
| `spherical_near_field_pressure` | −(ρ/2)(R/D)⁴ Ṙ², the spherical Bernoulli correction |
| `dipole_induction_pressure` | n·F/(4πD²) |
| `dipole_radiation_pressure` | n·Ḟ/(4πcD) |

F is the compact source’s net force **on the fluid** in N. Both dipole terms
vanish identically for fixed-center spherical RPE pulsation. For another model,
`compactForceDipolePressures` accepts separate Cartesian expressions for F
and its physical-time derivative in N/s. `rpeReceiverPressures` constructs the
standalone RPE catalog. The spherical correction describes the near-field
liquid pressure, a distinct effect from angular multipole radiation.
[Spherical pressure, Eq. 24](https://eprints.soton.ac.uk/45698/1/Pub9182.pdf)
and [point-force dipole, §1.7](https://assets.cambridge.org/97805218/12818/sample/9780521812818ws.pdf).

### Optional input-dependent terms

The derivative-query switch include_input_terms defaults to false. When true,
input supplies a scalar or vector of native input values and
input_time_derivative optionally supplies their rates in native time. The
system must provide analytic inputExpressions, matching its inputMatrix.
Observable::inputFlowDerivative computes the same result through C++.

For scalar affine actuation, with b=L_g h0, the included terms are

    d(e)/dtau = e1 + b*u
    d2(e)/dtau2 = e2 + (L_g e1 + L_f0 b)*u + (L_g b)*u*u + b*udot

The returned spatial Jacobian differentiates these expressions with input and
input rate held fixed. A missing input rate explicitly means held-constant
input; derivatives above order two reject nonzero supplied rates because higher
input derivatives would also be needed. With the switch off, supplied input
values/rates are ignored. Start a corrected echo derivative query from
observable="echo" and flow_order=1 or 2; the named catalog derivatives themselves
always retain their autonomous definitions.

For SINDy, include_input_terms controls whether known archived input enters the
source RHS used for Gamma. Its default true preserves previous source-derivative
behavior; set it explicitly to false for an autonomous fit. The standard
pipeline settings dialog exposes this as “Include known input in generator
derivatives.” This setting does not redefine the stored coordinate observable:
echo_jet2 remains the autonomous [e,e1] preset.

In the autonomous control experiment the fitted coordinates and mapped targets
use [e,e1], with input correction disabled. The plant still receives control,
and the exact lifted input map remains `J_phi(e,e1) * J_h(x) * g(x)`, where
`h(x) = (e(x), e1(x))`. Disabling
the observable correction does not disable actuation. A from_state target maps
both observables on the supplied state reference: its first component is the
same scalar echo target as before. Its second component is e1(x_ref), which is
distinct from differentiating the arbitrary reference waveform in time.
