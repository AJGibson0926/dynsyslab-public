# Observable-coordinate SINDy and KLQR

For the C++ engine, CLI, echo pipeline and direct observable targets, see
[NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md](NATIVE_OBSERVABLE_EXPRESSIONS_AND_TARGETS.md).

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

## Job settings

Systems own named scalar observables and ordered coordinate presets. The native
`sindy` stage selects the coordinate map with `coordinate_observable` (a named
observable or coordinate preset), `coordinate_observables` (an ordered list of
names), or an explicit `coordinate_expression`:

```json
{"name": "sindy",
 "settings": {"coordinate_observable": "c_k5", "polynomial_degree": 4}}
```

The Reddy--Szeri presets are:

- `radial_state`: `(radius, radial_velocity)`;
- `r_rdot`: alias of `radial_state`;
- `c_k2`: `(c, k2)`;
- `c_k5`: `(c, k5)`.

KLQR needs no coordinate-specific branch. Its state target is lifted through
the KROM's eigenfunctions, and `B_psi` evaluates the composed Jacobian. This
permits the same controller in radial or transformed coordinates.

## Analytic observable expressions

`AnalyticExpression` (`native/include/dynsyslab/core/AnalyticExpression.hpp`)
is a small serializable expression tree. It supports state variables,
constants, arithmetic, constant real powers, and common elementary functions.
An expression observable is differentiated exactly with respect to every
dynamical-state component, and its Jacobian evaluates those exact partials.
Finite differences remain an explicit validation backend, not the default for
analytic observables.

The retained Reddy `c`, `k2`, `k3`, `k4`, `k5`, radial acceleration,
compressibility loss, and `g1` definitions use this representation. The `c`,
`k2`, and `k5` material derivatives and their radial/radial-velocity partials
are available in the system observable catalog.

## Eigenfunction and archive ownership

Each SINDy eigenfunction retains its coordinate observable, and a KROM verifies
that all of its eigenfunctions use the same map. Both expose the map in
metadata, and evaluate `psi(x)`, its state gradient, and `B_psi(x)` from the
dynamical state.

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
