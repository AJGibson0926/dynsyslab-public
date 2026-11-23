# Analytic tori and perturbation tools

The qualification criterion is a known invariant-set parametrization with a
valid harmonic expansion in **native model time**. Elliptic functions are not
required. Elementary, elliptic and other analytic representations can implement
the same interface. A closed-form expression alone does not certify a Fourier
expansion, integrability, a global Koopman eigenfunction or a KAM torus.

## Native ownership

- `AnalyticTrajectory` in foundation exposes a time evaluator, an optional
  physical period, optional angular frequencies, and a torus map. Torus angles
  are radians relative to the supplied initial state. A multifrequency torus
  does not receive an invented common period.
- `sampleInvariantTorus` checks map values, coordinate wrapping and the sampled
  invariance equation \(DK(\theta)\Omega=f(K(\theta))\). The derivative check
  uses a fourth-order directional difference. It is a numerical diagnostic,
  not a rigorous enclosure or an integrability theorem.
- `TorusFourier.hpp` owns multidimensional FFT, Fourier evaluation, Parseval,
  angular/time derivatives and a finite cohomological solve. It uses the
  existing native arbitrary-length FFT; it owns no systems or archives.
- `Perturbation.hpp` owns finite small-divisor diagnostics, a homogeneous
  normal-form step, a sampled planar Melnikov integral and a constant-form
  conformal symplectic defect. Workflow binds systems; Python only calls CLI.

## Registered analytic examples

These are ordinary DynSys library systems, with state, energy, energy-rate,
input-power and time-derivative observables. Their preset is `regular_bounded`.
All use zero damping and forcing for the analytic invariant-set evaluator.
Nonzero configured damping remains available for numerical experiments.

| System key | Hamiltonian / state | Analytic scope |
|---|---|---|
| `morse_oscillator` | \(H=p^2/(2m)+D(1-e^{-aq})^2\) | Elementary bounded orbit, \(0<H<D\). |
| `isotonic_oscillator` | \(H=(p^2+\omega^2q^2+g/q^2)/2\) | Positive component \(q>0\), regular orbit, frequency \(2\omega\). |
| `poschl_teller` | \(H=p^2/2-D\operatorname{sech}^2(aq)\) | Elementary bounded orbit, \(-D<H<0\). |
| `cosh_oscillator` | \(H=p^2/2+2a\cosh q\) | Quartic reduction, evaluated in physical time. |
| `radial_twist` | \(I=(q^2+p^2)/2,\ H=\omega I+\alpha I^2/2\) | Elementary circles; frequency depends on action. |
| `kepler_planar` | \(H=|p|^2/2-\mu/|q|\), state \((q_1,q_2,p_1,p_2)\) | Bound noncollision ellipse; safeguarded numerical inverse of Kepler's analytic clock equation. |
| `euler_top` | \(H=\tfrac12\sum M_j^2/I_j\), state \((M_1,M_2,M_3)\) | Euler Lie–Poisson flow, \(I_1<I_2<I_3\), both regular Jacobi branches; excludes separatrix and principal-axis equilibria. |
| `coupled_duffing_hamiltonian` | \(H=\sum_j[p_j^2/2+\omega_j^2q_j^2/2+\beta_jq_j^4/4]+\epsilon q_1^2q_2^2\) | Product torus at \(\epsilon=0\), positive component energies; nonzero coupling is a Hamiltonian perturbation, not certified integrable. |

The radial-twist observable catalog includes `action`, `hamiltonian_frequency`,
`phase_cosine` and `phase_sine`, with analytic derivative links. Its canonical
phase convention is \(q=\sqrt{2I}\cos\theta\),
\(p=-\sqrt{2I}\sin\theta\). The normalized phase observables are undefined at
the equilibrium and must be used on a punctured domain. `hamiltonian_frequency`
is \(\omega+\alpha I\); with nonzero damping it is not the entire instantaneous
phase derivative. The generated derivative observables use the full drift.

Kepler also exposes angular momentum and eccentricity-vector components;
Euler exposes its Casimir; the Duffing pair exposes component and coupling
energies. The one-pair systems expose `position_squared` and
`position_momentum`. Equilibria are not assigned arbitrary phase circles.

The parameter \(m=k^2\), not modulus \(k\), is used by all native Jacobi APIs.
Even confining quartics crossing zero use a stable `cn` representation. Other
quartics use the existing Weierstrass representation. Tests include exact
quarter periods, nonzero initial velocity and both directions of time.

## CLI and optional client

`dynsyslab-cli analytic request.json`, schema `dynsyslab.analytic-query/v1`:

```json
{"schema":"dynsyslab.analytic-query/v1","operation":"torus",
 "system":{"name":"coupled_duffing_hamiltonian"},
 "initial_state":[0.7,0.4,0,0],"shape":[32,32],
 "maximum_harmonics":[12,12],"observables":["q_1","q_2"]}
```

The product mesh has no repeated endpoints; the last axis varies fastest.
Outputs include states, complex measurements, multiindex Fourier coefficients,
physical angular frequencies and invariance/wrap residuals. Defaults are 256
samples in one angle and 32 per angle in higher dimensions. The mesh is capped
at one million points; each retained harmonic is below its axis Nyquist value.
Refine the mesh to assess aliasing. `sampled_checks_pass` is a diagnostic;
`integrability_certified` and `kam_certified` remain false. Lack of a registered
evaluator returns `not_applicable`, not a nonintegrability proof.

```python
from dynsyslab_client.periodic import analyze_torus
result = analyze_torus({"system":{"name":"radial_twist"},
                        "initial_state":[1,0],"shape":[128]})
```

## Mathematics operations

Use `dynsyslab-cli mathematics request.json`, schema `dynsyslab.math-query/v1`.
Torus series contain `angular_frequencies` and `coefficients`, each record
`{"index":[k1,k2],"coefficient":[real,imag]}`.

| Operation | Input and behavior |
|---|---|
| `torus_fourier`, `from_samples` | Complex `samples`, `shape`, `angular_frequencies`, `maximum_harmonics`; product FFT with \(c_k e^{ik\cdot\theta}\) convention. |
| `torus_fourier`, `evaluate` | `series`, `phases`; evaluate one point. |
| `torus_fourier`, `time_derivative` | `series`; multiply by \(ik\cdot\Omega\). |
| `torus_fourier`, `partial_derivative` | `series`, zero-based `axis`; multiply by \(ik_j\). |
| `torus_fourier`, `mean_square` | `series`; sum squared coefficient magnitudes. |
| `torus_fourier`, `cohomological` | `series`, `divisor_tolerance`; solve \(\Omega\cdot\partial_\theta\chi=g-g_{retained}\). Mean and small-divisor terms remain in `retained`. |
| `small_divisors` | `angular_frequencies`, `maximum_order`, `exponent`, optional `retain`; search finite \(0<|k|_1\le N\), one representative of \(\pm k\). |
| `normal_form` | `eigenvalues`, homogeneous `terms` with `component`, `powers`, `coefficient`, `divisor_tolerance`; diagonal semisimple linear part only. |
| `melnikov` | Increasing `time`, flat two-component `drift` and `perturbation`; trapezoid integral of their determinant. No infinite-tail bound. |
| `symplectic_defect` | Flat `jacobian`, constant skew `form`, `dimension`, optional `damping`; Frobenius norm of \(A^TJ+JA+\gamma J\). Caller must supply a nondegenerate symplectic form. |
| `hamiltonian_first_order` | Torus-series `perturbation`, list `action_gradient`, flat symmetric `hessian`, `divisor_tolerance`; canonical first-order generator and action/angle corrections at one action. |
| `hamiltonian_twist` | `angular_frequencies`, flat symmetric `hessian`, optional `relative_tolerance`; singular values of the Hessian and isoenergetic bordered matrix. Requires the optional native linear-algebra backend. |

Normal-form convention is \(x=y+h(y)\). The divisor is
\(m\cdot\lambda-\lambda_j\); the change coefficient is the nonlinear
coefficient divided by this divisor. This operation treats one homogeneous
degree after any lower-degree changes have been applied. It does not perform
a complete multidegree normalization or certify convergence.

KAM proofs require more than small sampled residuals and a finite resonance
search: appropriate nondegeneracy, infinite Diophantine bounds, function-space
norms and quantitative smallness estimates must also hold. The current tools
are building blocks and diagnostics, not an automatic KAM theorem prover.

## First-order Hamiltonian torus corrections

For \(H=H_0(I)+\varepsilon H_1(I,\theta)\), supply the Fourier coefficients of
\(H_1\), of its action derivatives, and \(D^2H_0\) at the selected action.
The API solves \(\Omega\cdot\partial_\theta\chi=H_1-H_{1,\mathrm{retained}}\).
Its convention is
\[
 I_{\mathrm{old}}=I-\varepsilon\partial_\theta\chi,\qquad
 \theta_{\mathrm{old}}=\theta+\varepsilon\partial_I\chi.
\]
The returned corrections are coefficients of \(\varepsilon\). Action
differentiation includes the derivative of \(k\cdot\Omega(I)\) in the
denominator. A coefficient that vanishes at the chosen action may still have
a nonzero action derivative, which the implementation retains.

Mean frequency shift is \(\partial_I\langle H_1\rangle\). Other retained
resonant terms can still produce angle-dependent dynamics, so this mean
alone is not a complete resonant normal form. The routine does not compose
higher Lie-transform orders or bound their remainder. Tests compare the
corrected embedding directly with Hamilton equations and check second-order
residual scaling.

The twist query returns singular values of \(D^2H_0\) and
\[
 \begin{pmatrix}D^2H_0&\Omega\\\Omega^T&0\end{pmatrix}.
\]
These are different nondegeneracy conditions. Relative numerical thresholds
depend on the supplied action units and are not interval certificates.

## Algorithm references

- Arnold, *Mathematical Methods of Classical Mechanics*, action–angle theory:
  https://doi.org/10.1007/978-1-4757-2063-1
- Arnold, *Geometrical Methods in the Theory of Ordinary Differential Equations*:
  https://doi.org/10.1007/978-1-4612-1037-5
- Krajňák and Wiggins, Morse trajectories and action–angle variables:
  https://arxiv.org/abs/1905.04059
- Pinney, isotonic/Ermakov equation: https://doi.org/10.1090/S0002-9939-1950-0038586-4
- NIST Jacobi identities and series: https://dlmf.nist.gov/22.6 and https://dlmf.nist.gov/22.11
- Teschl, planar Melnikov theorem, §13.3:
  https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf
- Figueras, Haro and Luque, rigorous a-posteriori KAM verification for symplectic
  maps: https://arxiv.org/abs/1601.00084

