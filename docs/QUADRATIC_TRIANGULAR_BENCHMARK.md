# Quadratic triangular benchmark

Native system ID: quadratic_triangular. The dimensionless state is \(\bar{x}=(x_1,x_2)^T\):
\[
\dot{\bar{x}}=(a x_1,\ b x_2+c x_1^2)^T+(g_1,g_2)^T u.
\]

| Preset | a | b | c |
|---|---:|---:|---:|
| gradient_closed_reference (default) | -1 | -3 | 1 |
| brunton_2016 | -0.05 | -1 | 1 |
| kaiser_2021 | -0.1 | -1 | 1 |
| infinite_closure_saddle | 1 | -3 | 1 |

Explicit parameters override presets. The scalar-input direction uses input_gain_1 and input_gain_2 (defaults 1 and 0). The default initial state is (1,0); controlled runs default to the origin. Unknown presets and nonfinite parameters are rejected.

## Native properties

The system provides RHS, analytic Jacobian, input field, analytic expressions, exact unforced flow, configured dynamical facts and source references. The flow handles b=2a by the resonant limit and refuses nonzero prescribed drivers.

Observables include x_1, x_2, one, x_1_squared, x_1_dot and x_2_dot. For b!=2a, eta is \(x_2+c x_1^2/(b-2a)\), with eta_dot linked. When c=0, eta=x_2 even at resonance. All named rate links explicitly represent zero-input derivatives. The eigenfunctions coordinate preset selects [one,x_1,x_1_squared,eta].

The saddle preset additionally registers the analytic observables mixed_invariant \(I=x_1^3\eta\) and invariant_exponential \(e^I\). These are not misrepresented as polynomial modes.

## C++ API

~~~cpp
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <stdexcept>
namespace c = dynsyslab::core;
c::SystemConfiguration config;
config.identifier = "quadratic_triangular";
auto system = c::createSystem(config);
auto known = c::knownPolynomialEigenfunctions(*system);
if (!known.eigenfunctions) throw std::runtime_error(known.reason);
const double x[]{0.4,-0.2};
auto z = known.eigenfunctions->evaluate(x);
auto gradient = known.eigenfunctions->gradient(3,x);
auto flow = system->analyticSolutions().at(0).instantiate({0.4,-0.2},0);
auto stateAtOne = flow.evaluate(1);
~~~

The bank is the canonical `EigenfunctionFitResult` (a `FunctionExpansion` and its `Spectrum`) in physical coordinates, with four modes and rates 0,a,2a,b. It has no invented fit residuals or training data. derivativeOperators has row-major layout coordinate × output mode × predictor mode and implements \(\partial_jz=D_jz\). stateReadout is state × mode. Higher derivatives use products of these operators.

At b=2a,c!=0, the API reports why this diagonal family is unavailable. The finite invariant Jordan lift and exact flow still exist. Near resonance the eigenbasis can be poorly conditioned.

Systems, exact solutions and model facts belong to dynsyslab_foundation. The polynomial bank and derivative operators belong to dynsyslab_eigenfunctions. The workflow exposes these through the existing CLI; there is no reverse system-to-Koopman dependency.

## CLI and query contract

~~~text
dynsyslab-cli analytic resources/examples/triangular-analytic.json
~~~

~~~json
{
  "schema": "dynsyslab.analytic-query/v1",
  "system": {"name": "quadratic_triangular", "preset": "gradient_closed_reference"},
  "include_eigenfunctions": true,
  "initial_state": [1,0],
  "initial_time": 0,
  "time": [0,0.5,1,2]
}
~~~

The source-tree request is examples/native/triangular-analytic.json; installation
places it under resources/examples. Omit IC/time for metadata only. The C++ query
service is workflow::evaluateAnalyticQuery. Python clients can invoke the same
compiled CLI without duplicating numerical code.

The known_eigenfunctions response has schema dynsyslab.known-eigenfunctions/v1. It contains explicit dimensions, exponents, normalization, complex coefficients as [real,imag], rates, derivative operators, readout, conditions and sources. Ordinary solve jobs use the same system ID.

## Exact mathematical guarantees

At the default parameters the slow graph is \(x_2=x_1^2\); its residual decays as \(e^{-3t}\). The additional invariant graphs \(x_2=x_1^2+C x_1^3\) mean tangency alone does not identify a unique slow manifold.

The four eigenfunctions have Jacobian shape 4×2, Hessian shape 4×2×2, and zero spatial derivatives above order two. First-partial closure proves all-order closure. Three modes suffice with nonlinear state readout; four are minimal if both states must be linear readouts in a generator- and derivative-invariant space.

Adding \(x_1^4\) requires \(x_1^3\). Every finite collection of globally smooth eigenfunctions of the stable reference admits finite repair: they are polynomial sums of \(x_1^p\eta^q\) with fixed weight p+3q=n. Smoothness at the origin is essential. In the saddle preset, adding \(e^I\) forces the infinite independent sequence \(x_1^{3n}e^I\), despite global smoothness.

Sources: [Brunton et al. (2016)](https://doi.org/10.1371/journal.pone.0150171) and [Kaiser, Kutz and Brunton, Section 6](https://arxiv.org/pdf/1707.01146). Their quadratic form matches the default through \(X_1=x_1,X_2=3x_2,\mu=-1,\lambda=-3\). Validation lives in native/tests/triangular_system_tests.cpp. The comprehensive private derivations and functional-space report are indexed under koopman-triangular-gradient-closed.
