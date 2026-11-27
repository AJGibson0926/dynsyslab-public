# Lagrangian coherent structures and spectral submanifolds

Two families of invariant geometry are computed natively in the C++ core and run as post-processing
methods: finite-time Lyapunov exponent (FTLE) fields with their ridges and variational Lagrangian
coherent structures (LCS), and spectral submanifolds (SSMs) with their reduced dynamics and backbone
curves, from a known vector field or from trajectory data. No external toolbox (LCStool, SSMtool,
SSMLearn) is used; the methods below are implemented from the papers cited.

| Piece | Where |
|---|---|
| Flow map on a grid, Cauchy–Green tensor, FTLE, ridges, strainlines | `core/LagrangianCoherence.hpp` (`src/FiniteTimeLyapunov.cpp`, `src/LcsExtraction.cpp`) |
| Parameterization method, backbone curves, data-driven SSMs | `core/SpectralSubmanifold.hpp` (`src/SsmParameterization.cpp`, `src/SsmLearning.cpp`, internal `src/PowerSeries.*`) |
| Library | `dynsyslab_invariant_geometry` (foundation; eigenfunctions for the polynomial table; numerics) |
| Methods `ftle`, `ssm` | `workflow/native/src/LcsPostprocess.cpp`, `SsmPostprocess.cpp` |
| Tests | `native/tests/lagrangian_coherence_tests.cpp`, `native/tests/spectral_submanifold_tests.cpp`, `workflow/native/tests/geometry_postprocess_tests.cpp` |

## FTLE and LCS

For x′ = f(t, x) the flow map Φ = Φ_{t0}^{t0+T} is integrated from every node of a `DomainGrid`
of initial conditions (classical RK4, steps ≤ `maximum_step`; T > 0 forward, T < 0 backward;
time-dependent and forced systems use their own clock). Its gradient comes from central differences
between neighbouring nodes (Shadden, Lekien and Marsden 2005) or, with `auxiliary_spacing` δ > 0,
from an auxiliary pair of orbits ±δ about each node along each gridded axis (Farazmand and Haller
2012, §4.1). The Cauchy–Green tensor C = DΦᵀDΦ is restricted to the gridded coordinates — a grid
that anchors the other coordinates is a slice — and its eigenvalues λ₁ ≤ … ≤ λ_k come from cyclic
Jacobi rotations. The FTLE is

$$\sigma = \frac{\log \lambda_k}{2|T|}$$

(Haller 2001; Shadden et al. 2005). On two-axis grids:

- **Ridges** (`ridgeNodes`): nodes above `ridge_threshold` × max σ where σ is a maximum along the
  Hessian eigenvector of most negative curvature — the second-derivative ridges of Shadden et al.
  2005, §4. Forward ridges approximate repelling LCS, backward ridges attracting LCS.
- **Variational LCS** (`strainlines`): repelling LCS are strainlines r′ = ξ₁(r) of the forward
  tensor on which λ₁ ≠ λ₂ > 1 and ⟨ξ₂, ∇²λ₂ ξ₂⟩ < 0 (Haller 2011). Seeds are λ₂ ridge nodes; each
  strainline is integrated both ways (RK4 with the bilinear eigenvector field sign-aligned to the
  direction of travel) until it leaves the grid or the conditions fail over the failure length
  l_f, and the failing ends are trimmed (Farazmand and Haller 2012). Candidates shorter than
  `minimum_length` are dropped, and of curves within `separation` of each other only the one with
  the largest mean λ₂ is kept. Attracting LCS are the repelling LCS of the backward flow (Haller
  and Yuan 2000; Haller 2015).

### The `ftle` method

```json
{"name": "ftle", "settings": {
   "domain": {"lower": [0, 0], "upper": [2, 1], "points": [201, 101], "axes": [0, 1]},
   "initial_time": 0, "integration_time": 15, "directions": ["forward", "backward"],
   "maximum_step": 0.01, "auxiliary_spacing": 0, "ridge_threshold": 0.3,
   "lcs": {"enabled": true, "seed_stride": 3, "step": 0.25, "minimum_length": 0,
           "failure_length": 0, "separation": 0}}}
```

The domain defaults to the source trajectory's bounds over its first two coordinates and the
integration time to its duration. It publishes one `field` object (`dynsyslab.grid-field/v1`) with
outputs `ftle_forward`, `ridge_forward`, `lcs_repelling`, `ftle_backward`, `ridge_backward` and
`lcs_attracting` (the masks are 0/1 at nodes a strainline crosses), so Viewer draws every one as a
heat map through the domain-field payload. The strainlines themselves are the arrays `lcs_points`
(samples × 2), `lcs_offsets` (curve starts), `lcs_direction` (+1 repelling, −1 attracting),
`lcs_length` and `lcs_mean_stretching`; `cauchy_green_eigenvalues_<direction>` keeps λ per node.

Tests: the linear saddle (σ = 1 exactly, both gradient schemes and both directions), the steady
double gyre (its ridges and strongest strainline lie on the separatrix x = 1), the time-periodic
double gyre of Shadden et al. 2005 (FTLE of the published magnitude; the backward field equals the
forward one reflected through (1, ½), the flow's symmetry under (x, y, t) ↦ (2 − x, 1 − y, −t); the
strongest LCS crosses the domain beside the midline), and through the job, the double-well Duffing oscillator, whose forward and
backward ridges lie on the homoclinic loop v²/2 − x²/2 + x⁴/4 = 0.

## Spectral submanifolds

An SSM W(E) of an equilibrium x₀ of x′ = f(x) is the smoothest invariant manifold tangent to a
spectral subspace E of A = Df(x₀); it exists and is unique among C^{σ(E)+1} manifolds, where the
spectral quotient σ(E) = Int(min_{λ∉E} Re λ / max_{λ∈E} Re λ), when E is non-resonant with the rest
of the spectrum (Haller and Ponsioen 2016). The parameterization method writes it as
x = x₀ + W(p), p ∈ ℂ^d, with reduced dynamics p′ = R(p), and solves the invariance equation
f(x₀ + W(p)) = DW(p) R(p) (Cabré, Fontich and de la Llave 2003; Haro et al. 2016). At order k, for
each monomial p^m with μ = m·λ_E,

$$(A - \mu I)\, W_m = V_E R_m + [DW\,R]_m - [f(x_0 + W)]_m ,$$

where the brackets are the order-k Taylor coefficients of what is already known. The composition
f(x₀ + W(p)) is evaluated exactly by truncated power-series arithmetic over the expression tree of
f (automatic differentiation; Griewank and Walther 2008), so any analytic drift built from
polynomials, powers, exp, log, sin, cos, tanh and sech is supported. The equation is solved in the
eigenbasis of A. In normal-form style R_m keeps exactly the inner resonances λ_j ≈ μ; for lightly
damped modes the comparison is of imaginary parts (near-resonance; Ponsioen, Pedergnana and Haller
2018; Jain and Haller 2022). In graph style R_m is the whole master projection. A near outer
resonance (a non-master λ_l ≈ μ) is reported as a warning: the SSM does not exist at that order.

For a two-dimensional SSM of a complex pair the normal form is
p′ = λp + Σ_j γ_j p^{j+1} p̄^j, and with p = ρ e^{iθ}

$$\dot\rho = \operatorname{Re}\lambda\,\rho + \sum_j \operatorname{Re}\gamma_j\,\rho^{2j+1}, \qquad
\dot\theta = \operatorname{Im}\lambda + \sum_j \operatorname{Im}\gamma_j\,\rho^{2j},$$

so the backbone curve is the frequency θ′(ρ) and damping ratio −ρ′/(ρ|θ′|) against the amplitude
max_θ |x_i(ρ e^{iθ}) − x₀,i| of chosen coordinates (Szalai, Ehrhardt and Haller 2017; Breunung and
Haller 2018). The real chart q (p = q_a + i q_b for a pair) gives W and R as real polynomials over
the canonical monomial table.

**Data-driven SSMs** (Cenedese, Axås, Bäuerlein, Avila and Haller 2022, SSMLearn): from
trajectories decaying to an equilibrium y₀ (each source one trajectory, its components the
observables; `skip_time` drops the off-manifold transient), the tangent space V is spanned by the
leading right singular vectors of the data, the reduced coordinates are η = Vᵀ(y − y₀), the
manifold is the polynomial graph y = y₀ + Vη + w(η) fitted by least squares (degrees 2…order,
optional ridge), and the reduced dynamics η′ = r(η) are regressed on numerical derivatives
(`derivative_method`, any method of [numerical differentiation](NUMERICAL_DIFFERENTIATION.md)).
The normal form of r — the parameterization method above with every mode a master mode — gives its
backbone. Each source is then predicted from its first post-transient sample and the normalized
trajectory error recorded.

### The `ssm` method

```json
{"name": "ssm", "settings": {
   "mode": "parameterization",
   "vector_field": {"expressions": ["v1", "-x1 - 0.5*x1^3", "v2", "-7.29*x2 - x1^2"],
                    "variables": ["x1", "v1", "x2", "v2"]},
   "equilibrium": [0, 0, 0, 0], "master_modes": [0], "order": 7, "style": "normal_form",
   "backbone": {"radius": 0.2, "points": 101, "coordinates": ["x1"]}}}
```

Without `vector_field` the source system's analytic drift is used; without `equilibrium`, Newton's
method from the first source's last sample. Master modes index the eigenvalues ordered by
decreasing real part (slowest first) and a complex one brings its conjugate. The data mode takes
`"mode": "data"`, `dimension`, `order` (manifold), `dynamics_order`, `normal_form_order`, `ridge`,
`skip_time`, `equilibrium` and `derivative_method`, and pools all sources of the job.

It publishes a `dynamical_model` object (schema `dynsyslab.spectral-submanifold/v1`): the reduced
vector field in real chart coordinates as `dynamics_model_json` and `coefficients` (the
dictionary-and-coefficients form shared by every polynomial dynamical model), `ssm_parameterization`
(x as a polynomial in the chart, features × state over `polynomialExponentTable(d, order, true)`),
the eigenvalues, master modes, spectral quotient, warnings, invariance error (parameterization) or
fit errors, tangent space, singular values and predictions (data). When the SSM is two-dimensional
and in normal form, its backbone is a second, one-axis `field` object over ρ with outputs
`frequency`, `damping_ratio` and `amplitude_<coordinate>`, which Viewer draws as line plots.

Tests: the Duffing oscillator x″ + 2ζx′ + x + κx³ = 0, whose backbone matches the
exact conservative frequency ω(a) = π√(1 + κa²) / (2K(m)), m = κa²/(2(1 + κa²)) (Kovacic and
Brennan 2011), to 3·10⁻⁶ at amplitudes up to 0.3; a two-mode
oscillator whose slow mode is that Duffing oscillator driving a faster mode, so its slow SSM
(spectral quotient ≥ 100) has the same backbone, in both styles with small invariance error; the
real chart against the complex parameterization; and SSMLearn-style fitting of that two-mode
system's trajectories (identified spectrum, backbone and held-out prediction). Through the job, the
hardening Duffing oscillator from its own drift and from three decaying trajectories.

## References

- G. Haller, Distinguished material surfaces and coherent structures in three-dimensional fluid
  flows, *Physica D* 149:248–277, 2001.
- S. C. Shadden, F. Lekien and J. E. Marsden, Definition and properties of Lagrangian coherent
  structures from finite-time Lyapunov exponents in two-dimensional aperiodic flows, *Physica D*
  212:271–304, 2005.
- G. Haller and G. Yuan, Lagrangian coherent structures and mixing in two-dimensional turbulence,
  *Physica D* 147:352–370, 2000.
- G. Haller, A variational theory of hyperbolic Lagrangian coherent structures, *Physica D*
  240:574–598, 2011.
- M. Farazmand and G. Haller, Computing Lagrangian coherent structures from their variational
  theory, *Chaos* 22:013128, 2012.
- G. Haller, Lagrangian coherent structures, *Annu. Rev. Fluid Mech.* 47:137–162, 2015.
- G. Haller and S. Ponsioen, Nonlinear normal modes and spectral submanifolds: existence,
  uniqueness and use in model reduction, *Nonlinear Dyn.* 86:1493–1534, 2016.
- X. Cabré, E. Fontich and R. de la Llave, The parameterization method for invariant manifolds I,
  *Indiana Univ. Math. J.* 52:283–328, 2003; A. Haro et al., *The Parameterization Method for
  Invariant Manifolds*, Springer, 2016.
- S. Ponsioen, T. Pedergnana and G. Haller, Automated computation of autonomous spectral
  submanifolds for nonlinear modal analysis, *J. Sound Vib.* 420:269–295, 2018.
- S. Jain and G. Haller, How to compute invariant manifolds and their reduced dynamics in
  high-dimensional finite element models, *Nonlinear Dyn.* 107:1417–1450, 2022.
- R. Szalai, D. Ehrhardt and G. Haller, Nonlinear model identification and spectral submanifolds
  for multi-degree-of-freedom mechanical vibrations, *Proc. R. Soc. A* 473:20160759, 2017;
  T. Breunung and G. Haller, Explicit backbone curves from spectral submanifolds of forced-damped
  nonlinear mechanical systems, *Proc. R. Soc. A* 474:20180083, 2018.
- M. Cenedese, J. Axås, B. Bäuerlein, K. Avila and G. Haller, Data-driven modeling and prediction
  of non-linearizable dynamics via spectral submanifolds, *Nat. Commun.* 13:872, 2022.
- A. Griewank and A. Walther, *Evaluating Derivatives*, 2nd ed., SIAM, 2008 (power-series
  arithmetic); G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., 2013, §8.5 (Jacobi).
- I. Kovacic and M. J. Brennan (eds.), *The Duffing Equation*, Wiley, 2011, §2.2 (exact frequency).
