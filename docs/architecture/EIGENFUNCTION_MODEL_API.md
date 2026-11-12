# Eigenfunction and KROM API

The native runtime representation is now independent of the algorithm that
fitted it. This is a compatibility-preserving refactor, not a new eigenfunction
extraction method or controller. Numerical SINDy, DMDc, EDMDc and KLQR objectives
are unchanged. The MATLAB suite is unchanged.

## Ownership and vocabulary

| Type | Responsibility |
| --- | --- |
| `PolynomialDictionary` | Coordinate dimension, monomial powers, center and scale |
| `PolynomialEigenfunctionSet` | Immutable shared dictionary, coefficients, spectrum and coordinate contract |
| `PolynomialEigenfunctionView` | One mode; shares the bank's storage and owns its lifetime |
| `EigenfunctionCoordinates` | What the polynomial variables mean; optional physical dimension and frozen-provider pointer |
| `EigenfunctionSpectrum` | Primary continuous rates or discrete multipliers, sample interval and logarithm branches |
| `EigenfunctionFitResult` | A polynomial bank plus separate fitting diagnostics |
| `EigenfunctionValidationResult` | Errors measured along a particular supplied trajectory |
| `KoopmanModel` | The bank, optional physical-system adapter, input maps, reconstructions and dataset pointers |

No fit result, ODE, sampled training values or live history buffer is needed to
construct or evaluate a polynomial. A manually supplied expression is still a
candidate: construction validates its data contract, not its dynamical accuracy.

The representation remains

`psi_j(q) = sum_f coefficients[j,f] product_i ((q_i-center_i)/scale_i)^powers[f,i]`.

Coefficients are complex, mode-major; exponents are feature-major. An explicit
all-zero exponent row supplies the constant term. The dictionary is shared
across modes. Unique monomials, positive finite scales, finite coefficients and
nonzero polynomials are required. The evaluator and analytic gradient use the
same scalar kernels as the legacy API, preserving operation order.

## Coordinate contracts

`physicalState` means the polynomial variables are the state itself.
`analyticObservable` identifies a concrete expression or registered observable
names. `measured` and `causalHistory` require externally supplied coordinate
values. A history contract can retain `providerIdentity` and `providerArtifact`
pointing to its frozen projector; it does not own or update a running buffer.

The three counts remain distinct: `physicalStateDimension()` is optional;
`coordinateDimension()` counts independent polynomial variables;
`eigenfunctionCount()` counts complex eigenfunctions. `realLiftedDimension()`
counts the existing controller's split-real representation. A history length is
provider metadata, not one of these dimensions.

`evaluateFromCoordinates(q)` evaluates the polynomial only. It neither queries
an ODE nor infers hidden state. `evaluateFromState(x)` is available for physical
state coordinates and resolved instantaneous analytic observables. It rejects
history/measurement coordinates. Named analytic observables without a concrete
expression need a system catalog to resolve state evaluation; direct coordinate
evaluation remains available. An analytic expression carries its input state
dimension even when no ODE is attached.

## Manual construction, without a physical system

```cpp
#include <dynsyslab/core/KoopmanModel.hpp>
namespace c = dynsyslab::core;

c::EigenfunctionCoordinates coordinates;
coordinates.kind = c::EigenfunctionCoordinateKind::causalHistory;
coordinates.dimension = 2;
coordinates.names = {"xi_1", "xi_2"};
coordinates.providerArtifact = "frozen-echo-projector.npz";

// psi(q) = 2 + q1 + i*q2; no normalization, one sampled eigenvalue.
c::PolynomialEigenfunctionSet bank(
    {2, {0,0, 1,0, 0,1}, {}, {}},
    {2.0, 1.0, {0.0,1.0}},
    c::EigenfunctionSpectrum::fromDiscreteMultipliers({{0.99,0.06}}, 0.01),
    coordinates);
c::KoopmanModel model(bank, 1); // one input; no invented physical state
const double q[] = {0.3, -0.2};
const auto psi = model.evaluateFromCoordinates(q);
const auto gradient = bank.eigenfunction(0).gradient(q);
```

For a physical-system adapter, use the overload taking `SystemConfiguration`.
The existing constructor taking `KoopmanSINDyResult` still works and converts
into the same canonical bank. Scattered/interpolated models remain supported;
they are not falsely cast as polynomials.

## Fitting paths

`fitSindyEigenfunctions(...)` calls the existing generator-PDE fitter and
returns `EigenfunctionFitResult`. Its diagnostics explicitly identify
`sindy-generator` and the residual `Dpsi(q) qdot - lambda psi(q)`.
Supplied coordinate derivatives may be measured or model-derived; omitting
them retains the existing finite-difference behavior. Analytic state-to-observable
push-forward fitting remains available through the existing coordinate fitter
and `polynomialEigenfunctionsFromLegacy(...)` adapter.

`DmdControlResult::eigenfunctionSet()` and
`EdmdControlResult::eigenfunctionSet()` export the same canonical bank. Their
variables default to measured coordinates; callers must explicitly assert
physical-state or analytic/history semantics. Exporting a bank does not relabel
controlled-predictor residuals as autonomous eigenfunction errors. Original DMDc
and EDMDc results retain their own regression evidence and operators.

Legacy DTOs do not always identify the fitting algorithm or residual definition.
The compatibility adapter leaves those fields unspecified instead of inferring
"SINDy" from the C++ type name. Missing residuals mean unavailable, not perfect.
Optional large training-value arrays remain legacy evidence, not runtime-bank
storage; extraction does not duplicate them. Subsetting a bank does not count
as a new fitting run.

## Spectrum and input-map time domains

Use `fromContinuousRates(lambda, dt)` or
`fromDiscreteMultipliers(mu, dt, branches)`. The primary values retain their
declared domain. Discrete multipliers are not overwritten by a logarithm.
For nonzero multipliers, the selected continuous rate is
`lambda = (log(mu) + 2*pi*i*branch)/dt`. The interval is the observation sampling
interval in the model's time units, not an adaptive solver step.

Zero discrete multipliers are valid. Their `continuousRate(j)` is absent;
`continuousRates()` and continuous-generator consumers reject the request.
They admit only nonnegative integer-sample transitions. For nonzero sampled
modes, fractional-time transitions use the recorded logarithm branch: this is
an assumed interpolation, not evidence of a unique underlying generator.
Archive round-trips retain both primary values and validated cached rates.

`stateJacobianControlMap(x)` still computes the exact continuous input map from
the physical-system adapter and coordinate Jacobian. It is unavailable for
measurement/history coordinates or models with no physical input model.
Existing constant/bilinear `setControlMap` entries remain **continuous** maps.
The old optional `ConstantControlMap::discreteMatrix` field is retained for
compatibility, but new code should use the separately typed `DiscreteInputMap`:

`psi[k+1] = diag(mu) psi[k] + D u[k]`.

`setDiscreteInputMap` checks dimensions, finiteness and sample interval. Its map
is not discoverable as a continuous control map, and there is no implicit
division by dt or conversion to continuous B-psi. This refactor does not make a
discrete EDMDc predictor automatically ready for the existing continuous KLQR.

## Validation and provenance

New code calls `validateEigenfunctions` and explicitly chooses physical-state
or polynomial-coordinate inputs. The overload accepting a bank alone evaluates
coordinates directly. Metrics remain one-step, anchored, multiple horizons,
finite-difference generator defect, change-relative error and optional domain
coverage. Zero-multiplier modes have meaningful discrete defects but explicitly
unavailable generator statistics. Consult availability flags, not their empty
numeric storage. No errors are attached permanently to the expression itself.
Viewer omits unavailable generator samples and displays a warning rather than
plotting them as zero-error measurements.

The workflow continues to verify same-parameter **unforced** validation data
for autonomous defects, reuse existing data first, and offer—not automatically
run—new simulations. A coordinate-only model returns
`coordinate-provider-required`; it does not invent bubble parameters. The
generic low-level validator does not establish whether samples are unforced.
Its caller must retain dataset identity, time units, forcing classification,
sampling and coordinate-provider provenance with the result.

KROM training/validation/test pointers remain independent of forced/unforced
classification. They survive save/load and extraction. A dataset pointer is not
a certificate of accuracy. Completed source archives are not migrated in place.

## Naming migration and persistence

| Existing API | Preferred API |
| --- | --- |
| Runtime `KoopmanSINDyResult` | `PolynomialEigenfunctionSet` |
| `model.eigenfunctions()` legacy DTO | `model.polynomialEigenfunctions()` |
| `liftedDimension()` | `eigenfunctionCount()` |
| `splitDimension()` | `realLiftedDimension()` |
| Ambiguous `stateDimension()` | `physicalStateDimension()` or `coordinateDimension()` |
| `evaluate(x)` / `evaluateCoordinates(q)` | `evaluateFromState(x)` / `evaluateFromCoordinates(q)` |
| `generatorSplit()` | `continuousGeneratorReal()` |
| `exactControlMap(x)` | `stateJacobianControlMap(x)` |
| `validateKoopmanSemigroup(...)` | `validateEigenfunctions(...)` with explicit input kind |

Old names remain adapters for supported cases. Requesting a legacy
continuous-only DTO for zero-multiplier modes fails explicitly. This is source
compatibility, not binary ABI compatibility: rebuild C++ clients with the new
headers and matching libraries. Rebuilt Solver and Viewer keep their existing
workflow/UI behavior; no alternate controller or plot implementation is added.
Continuous-only downstream tools still require a finite continuous spectrum.

New KROM archives use `dynsyslab.koopman-model/v2`, with an embedded versioned
polynomial/coordinate contract. Established coefficient, normalization and
continuous-spectrum array names remain where valid. Old native and legacy
archives remain readable. Optional physical adapters, provider pointers, fit
evidence, primary spectrum, named continuous/discrete maps and dataset refs
round-trip. The raw `dynsyslab.edmdc/v1` fit archive remains a separate fitted
predictor contract; export it to the shared polynomial bank when needed.

Inspect an archive without fitting or mutating it:

```powershell
dynsyslab-cli krom-info C:/data/krom.npz
```

The JSON reports counts, coordinate semantics, primary spectrum, optional
continuous rates, input-map domains, fit-evidence availability and dataset refs.
The same service is `core/KoopmanInspection.hpp`.

## Runtime composition

Polynomial evaluation and legacy DTO adapters live in `dynsyslab_eigenfunctions`;
SINDy fitting lives in `dynsyslab_identification`. Runtime consumers can link the
former without the fitter. Controllers accept `valueFromCoordinates` or
`valueFromLifted`; estimators accept measured outputs and an initial lifted value.
The physical-state calls remain adapters.

Frozen history providers have a content identity, measurement clock, units and
coordinate selection. Per-stream buffers supply their current coordinates; they
are not stored in the eigenfunction bank. New input maps carry the identity of
the model whose coefficients, ordering and normalization they use. Archive
binding verifies physical compatibility and records any explicit legacy or
transfer acknowledgment.

See [pipeline composition and verification](PIPELINE_BOUNDARY_REFACTOR_2026-10-31.md)
for the native targets, job settings, legacy migration and deployment limits.

## Regression evidence

The case-study client in
`examples/native/case-studies/eigenfunction-api-refactor-20261028/` freezes the
pre-refactor results and compares 15,630 quantities: 15 modes' coefficients,
spectra, grid values, gradients, exact input maps and validation errors on five
unforced RPE trajectories. The local baseline is
`build/eigenfunction-api-parity/baseline.json`; large scientific assets remain
local and are only referenced. The comparison must not regenerate its baseline
after a code change. Existing KLQR case-study parity tests separately exercise
controller behavior. New `eigenfunction_contract` tests cover the public contracts,
zero modes, no-ODE deployment, provenance and archive round-trips.
