# Gravity-free double-pendulum reference

The native double_pendulum system owns the double-pendulum-gravity-free analytic family. It applies when gravity_m_s2, damping1_n_m_s, damping2_n_m_s and the external driver amplitude are all zero. Masses and lengths must be finite and positive, with a numerically resolvable inertia determinant. Ordinary nonzero-gravity double-pendulum trajectories do not receive this reference.

The native state is [theta1, omega1, theta2, omega2], with both angles absolute. In the Bocci–Mingari Scarpello paper, the second angle and velocity are relative: convert them by theta2=theta+psi and omega2=theta_dot+psi_dot. The 70°, 88° and 160° examples are included in native numerical-comparison tests with this conversion.

## Reduction and evaluation

Let delta=theta1-theta2, A=(m1+m2)l1², B=m2*l1*l2 and C=m2*l2². The two conserved quantities are

- H=(A*omega1²+2B*cos(delta)*omega1*omega2+C*omega2²)/2;
- L=(A+B*cos(delta))*omega1+(C+B*cos(delta))*omega2.

With D=A+C+2B*cos(delta), E=A+B*cos(delta) and Q=AC-B²*cos²(delta), the exact reduced equations are delta_dot²=(2H*D-L²)/Q and omega2=L/D-(E/D)*delta_dot. They determine a real hyperelliptic time integral and a second angle integral. Native evaluation uses adaptive quadrature and monotone inversion of these analytic relations; it does not call an ODE integrator to generate the reference.

For positive H, epsilon=(L²-2H(A+C))/(4HB) classifies relative motion. Libration has -1<epsilon<1, rotation has epsilon<-1, and epsilon=-1 is the relative separatrix. Rest and collinear rigid rotation are evaluated directly. The separatrix approaches folded rigid rotation and must not be labeled as a manifold of an isolated equilibrium in the original four-state coordinates.

Turning points are regularized using delta=2*asin(k*sin(u)), k²=(1-epsilon)/2. The separatrix uses a hyperbolic coordinate and an unbounded clock; rotation uses a winding coordinate. Angles remain unwrapped. Periodic primitives use 512 quadrature intervals and checked local inversion, rather than interpolating the numerical trajectory. Evaluation rejects unresolved finite-domain or convergence failures; invalid analytic samples retain the common analytic_valid mask. Quadrature tolerance is 2e-12, and periodic winding beyond 1e12 cycles is rejected because its phase cannot be resolved reliably in double precision.

## Native special function and interfaces

special::lauricellaFD(a,b,c,x) in dynsyslab/core/SpecialFunctions.hpp evaluates the real Euler-integral branch of Lauricella FD. It requires c>a>0, matching finite b/x vectors and x_j<=1. If arguments equal one, their combined endpoint exponent must remain integrable: c-a-sum(b_j for x_j=1)>0. Negative arguments are supported. Complex continuation through x_j>1 and the general genus-two Kleinian sigma function are not exposed by this API.

The versioned special-function-query/v1 command uses function=lauricella_fd and fields a, b, c, x. For example, a=1, b=[1], c=2, x=[0.5] returns 2*log(2). The existing analytic-query/v1 and compute-job/v1 interfaces discover the family automatically; solver.analytic_solution=off disables archived references.

- [Analytic query](../examples/native/double-pendulum-analytic-query.json)
- [Solve job](../examples/native/double-pendulum-gravity-free.json)
- [Special-function query](../examples/native/lauricella-query.json)

Viewer uses the usual muted purple dashed Analytic curves on state time series and phase plots. Paper catalogues remain on the system; solve archives contain raw numerical and analytic samples and computational metadata only.

## Research and checked guides

[Enolskii, Pronine and Richter (2003)](https://doi.org/10.1007/s00332-002-0514-0) describe the genus-two theta-divisor construction. [Bocci and Mingari Scarpello (2022, arXiv v2)](https://arxiv.org/abs/2104.13211) express the real integrals through Lauricella functions. Their publisher/primary records are retained in the extensible system bibliography. The native implementation uses the real integral reduction rather than claiming a general complex Kleinian implementation.

System metadata also supplies nonzero-gravity vertical equilibrium candidates, the zero-gravity equilibrium continuum, mechanical energy (including gravity) and angular momentum under the appropriate conditions. The sidecar checks all candidate states against the configured native field and checks invariants along computed branches. Tests independently check the gradients of the invariants against the native equations and the published Lauricella period formula against the analytic trajectory.
