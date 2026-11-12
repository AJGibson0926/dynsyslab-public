# Coordinate datasets and dynamics SINDy

The identification layer supports two regression objectives:

- Koopman SINDy fits (D Theta(q) qdot - lambda Theta(q)) a approximately 0,
  with a nonzero normalization.
- Original dynamics SINDy fits qdot approximately Theta(q) C by sequential
  thresholded least squares. Each column predicts one coordinate derivative.

These share a framework. A vector field is not an eigenfunction, and no
eigenvalue is invented for a vector-field fit.

## Reused components

| Responsibility | Existing implementation |
| --- | --- |
| Analytic observations and chain-rule derivatives | Observable / AnalyticExpression |
| Frozen causal history projection and online ring buffer | FrozenHistoryTransform / projectHistory / CausalHistoryProjector |
| Polynomial monomial order and normalization | PolynomialDictionary / polynomialExponentTable |
| Least squares and rank tolerance | numerics::solveLeastSquares |
| Sample matrices and clock | TimeSeries |
| Adaptive integration | integrate / ODESystem |
| Immutable archives, compatibility, saved jobs | Existing workflow post-processing |
| Plots | Viewer Figure, Trace, PlotCanvas and CLI renderer |

PolynomialDictionary::asObservable exposes the same dictionary as a real
candidate library. The history_projection job now calls the common
coordinate-series adapter; its frozen scaling and identity are unchanged.
The unfinished standalone RPE fitter was removed.

## Public native APIs

CoordinateTransform.hpp belongs to the observation module:

- transformCoordinates evaluates an Observable and optionally pushes supplied
  source rates through its Jacobian.
- transformHistoryCoordinates projects scalar measurements and optional
  measurement rates through the frozen U/sigma. Timestamps mark window ends.
  It does not recompute a decomposition.
- differentiateTimeSeries uses a centered five-point fourth-order derivative
  on a uniform clock, dropping two samples at each end.
- trimTimeSeries and requireAlignedTimeSeries enforce explicit alignment.
- `windowIntegralCoordinates(signal, W)` returns the current scalar measurement
  and its signed integral over the preceding duration W. The trapezoidal rule
  integrates the piecewise-linear signal, including fractional first intervals
  on irregular clocks. Output starts with the first complete window in each record.

SindyDynamics.hpp belongs to the existing identification module:

- polynomialDynamicsLibrary uses the shared dictionary (15 terms including
  the constant for total degree four in two coordinates).
- fitSindyDynamics accepts aligned coordinates/rates, a candidate Observable
  library and sparsity/rank options.
- SindyDynamicsModel stores the library and feature-by-output coefficients.
  evaluate predicts rates; vectorField returns a standard analytic Observable.
  definition/fromDefinition use the portable dynsyslab.sindy-dynamics/v1 schema.
- dynamicsError reports componentwise absolute RMS, relative RMS and maximum
  errors. Relative errors for zero-reference components are JSON null.
- integrateSindyDynamics adapts the model to the existing ODE solver.

Thresholding follows training feature/output RMS normalization. Saved
coefficients are in the original candidate-library units. Coordinate
normalization is embedded in the expressions and frozen for validation.
Zero derivative outputs and empty supports are valid zero vector-field
components. They would not be valid nontrivial Koopman eigenfunctions.

## Standard jobs

Use dynsyslab-cli plan JOB.json and dynsyslab-cli run JOB.json --archive-root ARCHIVES.

coordinate_transform accepts kind identity, observable, history or window_integral. Analytic
maps use a registered observable/preset or an expression definition. A history
map takes transform_archive and a scalar measured array, with an optional
registered observable as a simulator adapter. system_derivatives defaults to
false; true explicitly requests an offline EOM reference from confirmed
unforced physical-state data.

For example, `{"kind":"window_integral","observable":"far_field_acoustic_echo","window_duration":1.275,"names":["echo","echo integral"]}`
uses the RPE simulator observable on a physical-state trajectory. For recorded
measurements, use `"measurement_array":"echo_values"` in place of `observable`
(the default array is `values`). The window uses the same clock as the scalar
signal. The second column has units of signal times clock unit. The transform
computes values only and records its window and input definition for compatibility checks.
Viewer reads the axis labels for transformed trajectories from `coordinate_names`
and plots the saved values on their original scale. Physical equilibrium markers
and phase overlays require the source state trajectory.

sindy_dynamics accepts derivatives finite_difference (default) or stored
(coordinate_derivatives). use_previous_output consumes a preceding coordinate
transform. The existing sindy stage still fits Koopman eigenfunctions.

Dynamics SINDy defaults to library_family broad. Sparse thresholding selects terms
across families; it does not presume that one family is correct. The broad library
contains total-degree polynomials, sin/cos at several frequencies, reciprocal
powers through degree five, log(q) where positive, log(1+z^2), exponentials at
multiple rates, tanh, sqrt(1+z^2), and mixed terms with other coordinates. Here z
is the frozen normalized coordinate. Singular raw reciprocals and invalid logs
are excluded using the training range; very large exponential arguments are
excluded before fitting. Exclusions and selected terms are recorded. These checks
do not guarantee validity outside the training range. library_family polynomial
keeps the smaller baseline, and library_families can enable/disable individual
families. Explicit expression libraries remain supported. Koopman SINDy and its
existing options are unchanged.

Settings:

- polynomial_degree 4, threshold 0.001, normalize true, max_iterations 20.
- library: optional standard observable-expression candidate library.
- fit_stride: thin fitting rows after full-resolution differentiation.
- validation_sources: separate unforced records at the same physics and
  with the same frozen coordinate transform.
- evaluation_intervals: ranges for regional derivative errors.
- rollout_starts / rollout_duration: free simulations from observed coordinates.
- diagnostic_stride: thin saved display arrays only, never error calculations.

Archives retain the portable model, fitting settings, source pointers,
full-resolution error summaries and sampled display/rollout arrays. They do not
copy the original physical training archive. Viewer has sindy_dynamics and
sindy_dynamics_transient views with separate training/validation pages.
The RPE comparison uses saved jobs in
examples/native/case-studies/rpe-vector-field-sindy-20261105. Scientific archives
remain local. There is no experiment-specific fitter or integrator.

## Interpretation

For q=h(x), the chain rule gives qdot=Dh(x)f(x). Expressing that as f_q(q)
requires an inverse or consistent derivatives on all fibers of h. Exact
trajectory derivatives do not establish global closure. The same caveat
applies to two-dimensional projections of history buffers.

Small derivative errors do not guarantee good free rollouts or generalization.
Polynomial library inadequacy can also cause failure in perfectly observed
physical-state coordinates: RPE contains inverse powers of radius.

## Observable-jet closure side-car

`dynsyslab-cli observable-closure QUERY.json` builds the analytic Lie-derivative
jet

`H_k(x) = (h(x), L_f h(x), ..., L_f^k h(x))`

for any registered autonomous system and analytic observable. The query schema
is `dynsyslab.observable-closure-query/v1`; see
`examples/native/rpe-echo-closure-query.json`. It requires a finite state box
and reports, for every requested derivative order:

- sampled Jacobian rank and the least sampled singular value;
- a nonlocal coordinate collision, if refinement finds distinct states with
  matching jets;
- a closure obstruction, if such a collision also has different
  `L_f^(k+1) h` values.

A coordinate collision disproves use of the jet as an injective coordinate
map. It obstructs autonomous jet dynamics only when the successor derivative
also differs. The result therefore keeps `coordinate_collision_witness` and
`closure_obstruction_witness` separate.

This side-car is deliberately labeled `bounded_domain_numerical` with
`search_is_proof: false`, `global_guarantee: false`, and a null
`proof_certificate`. Matching within floating-point tolerances is numerical
counterexample evidence; failure to find a collision is not a global theorem.
The RPE example includes reproducible seeds for the known order-one and
order-two echo-jet obstructions.

The theorem-grade propositions are different and should not be inferred from
the search:

- global injectivity asks whether any off-diagonal pair satisfies
  `H_k(x) = H_k(y)`;
- global closure asks whether any such pair also satisfies
  `L_f^(k+1)h(x) != L_f^(k+1)h(y)`.

For rational systems on semialgebraic state sets these are real-algebraic
decision problems after clearing denominators. Exact quantifier elimination,
a real-radical/ideal certificate, or a valid global-univalence argument can
prove them. Interval Newton/Krawczyk isolation can certify a counterexample in
a bounded box. Sampled rank, unconstrained nonlinear least squares, and a
finite box sweep are not proof certificates. Use the v2 query below for global analysis.

### Global queries and equations

Run the same CLI with `"schema":"dynsyslab.observable-closure-query/v2"`.
The global query uses the regular real domain of the supplied expressions.
`positive_axes` restricts specified zero-based axes to the entire positive
half-line. Simulator stopping conditions and invariant submanifolds require
separate interpretation. Inputs are autonomous in native model time.

For a specified time law, include a clock coordinate with `t'=1` in the
supplied equations. A known sinusoidal drive can instead be added as two
states with `a'=omega*b, b'=-omega*a`. The solver tests the full declared
state space, including all driver phases and amplitudes unless encoded in
the coordinate expressions. An arbitrary external input can prevent any
autonomous jet closure. For `x'=u(t), h=x`, admissible histories
`x=0` and `x=t^(k+1)` agree through order k at zero and differ at the
next derivative, for every k.

| Input | Choices |
| --- | --- |
| Equations | Registered `system.name` with preset/parameters, or `system.equations`, an array of native analytic state-expression nodes in state order |
| Measured coordinates | A registered `observable` or preset, an ordered list of registered `coordinates`, zero-based `coordinate_axes`, or a scalar/vector `expression` in the existing observable-expression schema |
| Explicit search | `highest_derivative_order` (default 4, maximum 12), `polynomial_degree` (default 2, maximum 6), optional `pair_seeds` |
| Exact real feasibility | `emit_smtlib: true` exports each two-state polynomial problem. `solver_timeout_ms` runs unresolved tests with Z3 if the process environment names an absolute trusted Z3 shared library in `DYNSYSLAB_Z3_LIBRARY`. Zero disables the solver. Each nonzero timeout applies per decision. |
| Evaluation | `evaluate_at_jet` evaluates an explicit certified RHS at a supplied realizable jet. The C++ gates are `globalJetClosure` and `evaluateJetRhs` in `GlobalJetClosure.hpp`. |
| Higher orders | `search_all_orders: true` continues after the first closed law up to the requested order or algebra budget. The report records the first tested injective jet, smooth RHS and smooth state-recovering jet it certifies. |

Components are ordered by derivative, then by observable component. For
`h=(a,b)` the order-one coordinates are `q=(a,b,Da,Db)`, where `D=L_f`.
Each row of `rhs_in_jet` lists exact rational coefficients and signed integer
powers of q. Empty rows are zero. `rhs_observable` is the corresponding
native expression for floating-point evaluation. A real-cube-root law stores
the exact radicands and root degrees instead. `state_in_jet`, when present,
gives a left inverse on realizable jets. Examples:

- `examples/native/duffing-position-global-jet.json` proves the minimum order
  is one and gives `q0'=q1, q1'=q0-q1/2-4*q0^3`.
- `examples/native/cubic-observation-global-jet.json` supplies both the ODE
  and sensor directly. For `x'=1, q=x^3`, it finds `q'=cbrt(27*q^2)` and
  evaluates it at -8. This is a single-valued derivative. At zero this ODE
  admits waiting solutions. The same query continues through order 2 and
  finds the smooth law `q'''=6` and the inverse `x=q''/6`.
  This distinguishes single-valued derivatives from uniqueness of a rollout.

The native algebra differentiates sparse Laurent polynomials with exact
rational arithmetic starting from the stored binary64 constants. It compares
coefficients to construct polynomial laws, real cube-root laws and inverse
maps. A literal state chart also permits inverse powers. Direct state-axis
substitution supports analytic drift expressions outside this algebra.

The same analysis works after adjoining variables for reciprocals,
sines/cosines, exponentials, logs, hyperbolic tangents, and powers on a
certified positive base. The chain rule supplies a polynomial ODE in the
enlarged state. For example, with `x'=1` and `q=exp(x)*sin(x)`, the code
constructs `q''=2*q'-2*q`. An arbitrary rational denominator p becomes
an auxiliary z satisfying `p*z=1` and `z'=-z^2*p'`. These graph equations
are imposed for both states in the solver. Polynomial division by known graph
identities can simplify the resulting laws. Auxiliary variables stay internal
to the analysis. The explicit RHS uses only the original observable jet.
[Hemery, Fages and Soliman on polynomialization](https://arxiv.org/abs/2106.15884).

The report distinguishes exact algebraic graphs from supersets of the
physical graph. For example `sin(x)^2+cos(x)^2=1` alone permits phases
inconsistent with a particular x. A closure identity proved on the superset
holds for every physical state. A solver collision on the superset leaves
the physical question open. Rational powers with small dyadic exponents
use a polynomial graph and a positive-root condition. A positive fractional
power at a possible zero of its base requires separate branch analysis.
The symbolic budget is 32 state/observable components and 32 lifted state
variables. Reports distinguish budget limits, solver outcomes, and proofs.

At a generic constant-rank point a rank increase from `H` to `(H,Dnext)`
proves that the successor varies along a fiber. At nonlocal intersections,
seed validation fixes enough state components to give a square system and
uses rational interval bounds to prove contraction of the preconditioned
Newton map into its box. A strictly signed successor difference over the box
certifies an actual counterexample. The two-state SMT problems independently
test for an unequal successor or unequal physical state at equal jets.
`unsat` for the first test certifies the single-valued implicit graph
`{(q,s) : exists x, H(x)=q and successor(x)=s}` on the declared domain.
The earlier blocks of q shift to the next block. `unknown` keeps that order
open. `minimality_proved` requires a refutation at every earlier order.

State recovery has its own tests. A dependency walk follows the observable's
state axes through the equations for those axes. A state axis outside this
closure is unobservable at every order on a nonempty open domain.
A field-equivariant sign reflection that
preserves the measurement and domain produces indistinguishable trajectories
at every derivative order. A smooth closed law whose jet has generic rank
below the state dimension also has positive-dimensional indistinguishable
fibers on an open set. These permit exact reduced dynamics even when the
original state cannot be recovered.

### Why the global search terminates in principle for rational models

In the ring of polynomials in two copies of the state, invert the finitely
many denominators of f and h on both copies. Every Lie derivative is in this
same ring. Form the ascending chain of ideals generated by componentwise
differences `D^j h(x)-D^j h(y)` through order k. The ring is Noetherian, so
the chain eventually stabilizes. At that order, equality of the jets implies
equality of every higher derivative. This proves a single-valued successor
on the real jet image. The ideal argument alone supplies neither state
injectivity nor smoothness of a formula on the whole image.
[Polynomial observability by ideals](https://research.rug.nl/en/publications/observability-conditions-by-polynomial-ideals/)
and [localization of Noetherian rings](https://stacks.math.columbia.edu/tag/00FM).

The side-car also computes a sufficient derivative order from rational input.
In the doubled state space of dimension m, clear the denominators of the
paired field with one common, nonvanishing time multiplier. Clear the output
difference denominators as well. Both operations preserve finite contact
order along the paired trajectory. If the resulting field and output have
degree bounds delta and d (padded to delta >= 1 and d >= m), every finite
contact order is at most `K=2^(m+1)*(d+(m-1)*(delta-1))^m`. A stationary
pair with equal output has constant equal output. Thus derivatives agreeing
through K agree forever. For Laurent input the code also tries reciprocal
coordinates on each excluded pole axis to lower the degree bound. When the
explicit search is unresolved, `guaranteed_order_dynamics` records this K,
the input equations and the implicit recurrence on the realizable jet image.
For a minimum, test each order in increasing order by real polynomial feasibility.
In principle this procedure terminates by K. The side-car's explicit algebra
budgets and solver timeout can stop an individual run much earlier.
[Binyamini, Corollary 1](https://arxiv.org/pdf/1406.1858).

Analyticity alone permits a genuine global obstruction. For example, take
`x'=1` and the observable `h(x) = product_{n>=1} (1-exp(x-n))^n`.
The product converges locally uniformly since the tails are controlled by
`sum n exp(Re(x)-n)`. At each positive integer n it has a zero of exact
multiplicity n, coming from the nth factor alone. At x=k+1 and
x=k+2 every derivative through order k vanishes, but the derivative of order
k+1 is nonzero at the first point and zero at the second. This gives a
counterexample for *each* finite k for one fixed analytic system and observable.
An unbounded family of such pairs rules out every finite closure. A finite
list of pairs provides only a lower bound.

### Global RPE echo result

`examples/native/rpe-echo-global-jet.json` checks the autonomous viscous
model at kappa=4/3. It certifies unequal successors for orders 0, 1 and 2.
At order 2 the two boxes are near (0.7500373, 5.6275344) and
(1.1887084, 7.5753553). The third-derivative gap is about 14625.
The original chart gives m=4, delta=11 and d=9. The global change u=1/r
gives delta=d=5, reducing the bound to K=32(5+3(5-1))^4=2,672,672.
The global minimum lies between 3 and K. The supplied six-second solver
budget returned `unknown` for orders 3 and 4.

There is a stronger state-recovery proof on that entire half-plane. The code
recognizes the field by its coefficients and verifies the identities below.
Write `r'=v` and `v'=g/r^5+b/r-c/r^2-d*v/r^2-3*v^2/(2*r)`, with
g>0, b<0, c>=0 and d>0. The measured echo is a nonzero constant times
`D^2(r^3)`. The energy and its exact derivative are

    E = r^3*v^2/2 + g/r - b*r^3/3 + c*r^2/2
    DE = -d*r*v^2.

Energy sublevels are compact inside r>0. On v=0 the only invariant motion
is the unique positive equilibrium, since `g+b*r^4-c*r^3` decreases from
positive to negative. LaSalle's invariance principle gives forward
convergence to this same equilibrium from every positive-radius state.
Two equal infinite echo jets give identical echo histories by analyticity.
The difference of their cubed radii is then affine in time and tends to zero,
so it vanishes identically. Radius and velocity therefore coincide. By the
finite contact bound, H_K already recovers the state. For realizable q its
exact implicit dynamics shifts the blocks and uses `D^(K+1)h(x)` at the
unique positive-radius solution of `H_K(x)=q`. The example does not expand
this many derivatives. This proof uses the mathematical continuation on
r>0 and also distinguishes initial states above the simulator's radius cutoff.
[LaSalle's invariance principle](https://www.egr.msu.edu/~khalil/NonlinearSystems/Sample/Lect_10.pdf).

The dynamics_validation job evaluates a saved model against a coordinate archive,
without fitting. model_archive is required; compare_model_archive optionally
compares two existing models through native evaluation. Coordinate transforms and
physical system compatibility are checked. Fitted archives record
sparse_regression_steps: support size and relative residual before every pruning
step, separately for each output, so sparsification losses are visible.
