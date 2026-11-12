# Analytic trajectories and phase-space context

Analytic solution families are properties of the configured native system. The solver evaluates an applicable family alongside its numerical trajectory by default. Set solver.analytic_solution to off to disable this, auto to select an applicable family, or a registered family identifier to request one. Applicability failures leave the numerical solution intact and retain a reason. Archived analytic_values and analytic_valid contain raw samples in the same coordinates as values. The sampled maximum absolute difference is a comparison diagnostic, not an error certificate.

Bibliographies, paper catalogues, and source literature belong to dynamical systems. Solution archives contain computational metadata and raw data; they do not copy paper references. See [System research metadata](SYSTEM_RESEARCH_METADATA.md) for the additive catalogue API.

## Available families

- Harmonic oscillator: free and sinusoidally forced linear motion, including damping, critical damping, resonance, and the zero-stiffness limit.
- Duffing: conservative quartic-energy trajectories using the Weierstrass representation; bounded wells and separatrix limits are checked against initial energy.
- Pendulum: undamped libration, rotation, separatrix, and zero-gravity limits using Jacobi elliptic functions.
- RPE: autonomous inviscid bounded gas-filled motion for kappa=4/3 with capillarity, or kappa=3/2 with zero capillarity. The radius is an elliptic-function parametrization; physical time uses adaptive quadrature and monotone inversion. This is not a general viscous RPE closed form.
- Gravity-free double pendulum: conservative real hyperelliptic reduction, including relative libration, rotation, separatrix, rest and rigid rotation. See [native equations, Lauricella support and parameter conditions](GRAVITY_FREE_DOUBLE_PENDULUM.md).
- Additional ODE catalogue cases: the zero-nonlinearity Van der Pol oscillator, Lorenz invariant-axis and zero-sigma limits, Rossler invariant-plane limit, affine Chua and Hindmarsh–Rose limits, harmonic Henon–Heiles and its nonlinear elliptic symmetry-axis orbits, linear FPUT normal modes, uncoupled Kuramoto, and uncoupled spring pendulums.
- Two-oscillator Kuramoto: exact Adler reduction with locked, critical and slipping branches; projective exponentials retain unwrapped phase winding and negative-time evaluation.
- Hodgkin–Huxley: exact fixed-voltage gate relaxation only when channel-current coefficients vanish and the bias balances leak, making constant voltage an invariant set of the native free-voltage equations. General firing trajectories do not use a voltage-clamp formula.
- Spring pendulums: equal-length synchronous motion with an unstretched spring and zero damping/input follows the full nonlinear pendulum, including unequal masses and nonzero spring stiffness.
- Marmottant: exact stationary states after checking the native pressure residual and reset shell history. Published linearized oscillations remain approximate and are not selected as exact references.
- Reddy–Szeri: stationary-radius shape oscillators for modes 2–5 when radial pressure balance, zero translation and effective forcing conditions hold. This is exact within the native perturbative shape closure; it is not an exact finite-amplitude hydrodynamic solution.

- Heat equation: continuum Fourier evolution of the supplied samples' trigonometric interpolant, with periodic, homogeneous Dirichlet, or Neumann boundaries and optional linear decay. A separately selectable family evaluates the exact semigroup of the actual finite-difference stencil. Continuum comparison errors include spatial and temporal discretization; the discrete family isolates temporal error. The edge/zero-flux stencil reflects at half-grid boundaries and is explicitly distinguished from the continuum interval.
- Mackey–Glass: constant equilibria, feedback-free decay, the exponent=0 linear-delay limit using a finite delayed-exponential sum, and the nonlinear constant-history formula through the first requested delay interval. After that interval, invalid reference samples are NaN with analytic_valid=0; no continuation is invented. Numerical delay rounding remains recorded separately.

The native special-functions API supplies real Jacobi sn/cn/dn, elliptic integrals, Jacobi amplitude, theta functions, and Weierstrass p and its derivative. The delay module also exposes a finite delayed-exponential evaluator for a scalar linear delay equation with constant history. Families state their conditions, formula, evaluation method, version, and system-owned references. Further catalogue research and broader special-function domains remain extensible through system/extension implementations.

## PDE and delay API conventions

PDE and DDE descriptors own an analyticSolutions factory bound to their resolved configuration. C++ clients can also call heatEquationAnalyticSolutions or mackeyGlassAnalyticSolutions from the corresponding native family headers. The existing analytic-query/v1 CLI/C ABI supports these families without an ODE adapter. For PDEs, initial_state is the complete mesh field; for Mackey–Glass it contains the one constant prehistory value.

Heat family identifiers are heat-continuum-fourier (automatic default) and heat-discrete-semigroup. The continuum initial function is explicitly the finite periodic Fourier, cosine, or sine interpolant; samples alone do not establish an unknown continuous initial profile. Zero-value boundaries require zero endpoint data. Viewer overlays the final spatial profile with the same muted purple dashed Analytic curve used for ODE trajectories.

Mackey–Glass family identifiers are mackey-glass-linear-delay, mackey-glass-no-feedback, mackey-glass-constant-equilibrium, and mackey-glass-first-delay. Automatic selection prefers global families before the interval-limited formula. The first interval is determined by the requested continuous delay, rather than the Euler history grid's rounded lag. Viewer reports partial validity and terminates the reference at the last valid sample. Setting solver.analytic_solution to off is honored for both PDE and DDE solves.

The special-function-query/v1 operation linear_delay accepts beta, gamma, delay, history, and elapsed, and evaluates x'=-gamma*x+beta*x(t-delay). Its beta is the actual delayed coefficient; Mackey–Glass with exponent=0 uses beta/2. Evaluation uses 50-digit working arithmetic, a finite sum, and explicit finite-domain/precision checks. It does not step a numerical DDE solver.

## Explicit phase-space analysis

The solver sidecar is called only by an explicit compute request. Opening a Viewer phase plot performs a cache lookup. A Phase context action above the workspace computes missing context on a worker and offers cancellation. Viewer also checks cheap system-provided equilibrium candidates against the native zero-input equations before marking them.

The native query schema is dynsyslab.phase-space-query/v1 with mode lookup or compute, a system configuration, and options containing full-state lower and upper bounds. Compute also checks the cache before running a finite multistart analysis. Results include equilibrium residuals and eigenvalues, one-dimensional stable/unstable saddle branches with seed-refinement checks, shooting-based periodic orbit candidates with Floquet multipliers, invariant-guide diagnostics, and a sparse vector field. Higher-dimensional invariant surfaces are explicitly reported as unsupported. Center-family periodic orbits are distinguished from isolated limit cycles. Searches do not establish global completeness or basin boundaries.

Autonomous zero-input analysis remains the default API mode. The explicit periodic-forcing mode follows the configured nonautonomous RHS at its actual time, including RPE's state-dependent pressure gain and Reddy–Szeri's internal acoustic clock. A common forcing period is resolved in native model time; an explicit period must be an integer multiple of every active clock. Incommensurate clocks are rejected. See [periodic forcing conventions, numerical validation and limits](PERIODIC_PHASE_CONTEXT.md).

Viewer automatically requests periodic context for a stored sinusoidal driver or a continuous single-tone waveform controller. It translates the archived controller's explicit sine formula to the native cosine-driver convention, with physical/model time and input conversions. Windowed, sampled, offset and multiple-tone controller forms are explicitly rejected by this adapter. Other feedback or random-control trajectories may still use labeled autonomous zero-input context; that context does not represent their controlled dynamics.

Marmottant still requires the entire domain to remain strictly within one memory-independent shell regime; domains crossing shell history or switching boundaries are rejected. Reddy–Szeri domains may not cross its gas-exponent boundary unless kappa=1. For higher-dimensional state projections, vector arrows require an explicit full-state slice; otherwise Viewer omits them.

Caches live under system_metadata, keyed by resolved equations, parameters, coordinates, phase domain, method version, tolerances, options, forcing waveform/phase, and numerical guide metadata. Bibliography additions do not invalidate the numerical cache. Cache records use a distinct .dslphase extension and do not enter the trajectory data tree. Existing solution archives are never rewritten.

## Viewer and CLI

Viewer renders analytic state trajectories as muted purple dashed curves labeled Analytic. Phase context uses thin subdued manifold/periodic-orbit curves and sparse gray direction arrows behind the principal trajectories. Equilibria are closed red dots; initial states are open red dots. Context cannot change trajectory-based axis limits and does not appear in time series, frequency plots, observable heat maps, or Koopman eigenfunction heat maps.

DynSysLabViewerCli phase-space --source trajectory.npz --mode lookup returns a read-only cache query using the same selected data-tree object and persistent Mask as Viewer. Use --mode compute to compute explicitly; --options options.json supplies numerical/domain/slice overrides. The existing inspect and render commands display cached context automatically without opening windows. The solver CLI exposes phase-space query.json --archive-root PATH and the workflow C ABI exposes dsl_workflow_phase_space_json with progress/cancellation callbacks.

The Viewer uses the installed native workflow C ABI without starting another application or shell. Standalone Viewer installations lacking the solver module can still show trajectories and analytic arrays; an explicit analysis request reports the missing capability.


## Catalogue-limit CLI examples

The same native query used by solver/Viewer integrations can be inspected without a window:

    dynsyslab-cli.exe analytic examples/native/kuramoto-pair-analytic-query.json
    dynsyslab-cli.exe analytic examples/native/hh-constant-voltage-analytic-query.json
    dynsyslab-cli.exe analytic examples/native/spring-pendulum-synchronous-analytic-query.json
    dynsyslab-cli.exe analytic examples/native/marmottant-stationary-analytic-query.json
    dynsyslab-cli.exe analytic examples/native/reddy-szeri-shapes-analytic-query.json

Each example uses explicit initial data and returns complete raw sample values. Catalogue-only queries omit the initial state and return system-owned families, applicability metadata and references. Solution archives retain the sampled arrays and computational identity without those references.

For variable-length state arrays, native presentation labels are expanded by systemStateLabels from SystemLabels.hpp. Kuramoto phases retain one-based mathematical subscripts at any configured oscillator count. This changes display labels only; state names and indices remain stable.

Two-dimensional autonomous saddle manifolds use validated adaptive surface meshes; see [surface continuation and limits](MANIFOLD_SURFACES.md). The remaining nonlinear Hénon–Heiles family is documented [here](HENON_HEILES_ANALYTIC.md).
