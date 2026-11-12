# Observable-coordinate SINDy and KLQR

For the shipped C++ engine, CLI, echo pipeline and direct observable targets, see
[NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md](NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md).
The Python API examples below describe the historical reference implementation.

DynSysLab distinguishes the dynamical state `x` from the coordinates in which
an eigenfunction is represented. An eigenfunction may be fitted as

`psi(x) = phi(h(x))`,

where `h(x)` is an ordered collection of model observables. The standard state
case is retained as `h(x)=x`.

## Mathematical contract

For autonomous training dynamics `xdot=f(x)`, SINDy receives

- modeled coordinates `y=h(x)`;
- modeled coordinate derivatives `ydot=J_h(x) f(x)`;
- continuous-time eigenvalues for the generator equation;
- the selected dictionary evaluated in `y`, not in `x`.

For control-affine dynamics `xdot=f(x)+B(x)u`, the exact lifted input map is

`B_psi(x) = nabla_y phi(h(x)) J_h(x) B(x)`.

The identity observable makes `J_h=I`, so the standard state-coordinate path
is exactly the same implementation rather than a separate special algorithm.
The constant and bilinear learned control maps remain alternative KROM control
representations and do not use this exact chain-rule evaluation.

## API

Systems expose scalar observables through `observable(name)` and ordered maps
through `observable_coordinates(selection)`. SINDy accepts an `Observable`, a
named system coordinate preset, or an ordered sequence of observable names.

```python
system = load_system("reddy_szeri", preset="matlab_cleanup_reference")

radial = system.observable_coordinates("radial_state")
c_k2 = system.observable_coordinates("c_k2")
c_k5 = system.observable_coordinates("c_k5")

result = koopman_with_sindy(
    spectrum=eigenvalue_spectrum,
    trajectories=training_trajectories,
    system=system,
    coordinate_observable="c_k5",
    coordinate_system="dimensionless",
)
```

The standard H-DMD to SINDy CLI exposes the same boundary:

```text
dynsyslab-postprocess run TRAINING_ARCHIVE --execute \
  --sindy-coordinates c k5 --coordinate-system dimensionless
```

The Reddy--Szeri presets are:

- `radial_state`: `(radius, radial_velocity)`;
- `r_rdot`: alias of `radial_state`;
- `c_k2`: `(c, k2)`;
- `c_k5`: `(c, k5)`.

`KoopmanLQRController` needs no coordinate-specific branch. Its state target is
lifted through `KROM.psi`, while `KROM.B_psi` evaluates the composed Jacobian.
This permits the same controller API in radial or transformed coordinates.

## Analytic observable expressions

`AnalyticExpression` is a small serializable expression tree. It supports state
variables, constants, arithmetic, constant real powers, and common elementary
functions. `Observable.from_expressions` differentiates the expression exactly
with respect to every dynamical-state component and stashes each scalar partial
as an `Observable`. The main observable's Jacobian evaluates those exact
partials. Finite differences remain an explicit validation/fallback backend,
not the default for analytic observables.

The retained Reddy `c`, `k2`, `k3`, `k4`, `k5`, radial acceleration,
compressibility loss, and `g1` definitions use this representation. The `c`,
`k2`, and `k5` material derivatives and their radial/radial-velocity partials
are available in the system observable catalog.

## Eigenfunction and archive ownership

Each SINDy eigenfunction retains its coordinate observable, and a KROM verifies
that all of its eigenfunctions use the same map. Both objects expose the map in
metadata. The executable methods are:

- `eigenfunction.psi_from_state(x)`;
- `eigenfunction.nabla_psi_from_state(x)`;
- `krom.psi(x)`;
- `krom.nabla_psi(x)`;
- `krom.B_psi(x)`.

SINDy results and archives retain two separate arrays:

- `training_states`: the original dynamical states, used for provenance,
  state-domain gates, reconstruction, and control validation;
- `training_coordinates`: `h(training_states)`, used to evaluate stored SINDy
  coefficients and deferred eigenfunction values.

This distinction is required whenever `h` is not the identity. Analytic
expression observables and stacks are archive reconstructible, including their
exact derivative definitions.

## MATLAB parity boundary

The read-only MATLAB reference constructs SINDy inputs from each observable's
value and stored time derivative, then constructs its KLQR coordinate Jacobian
from stored spatial derivative observables. DynSysLab implements the same
mathematics with a single general chain rule. It does not reproduce the MATLAB
controller's fixed two-observable assumptions or its duplicated derivative
bookkeeping.
