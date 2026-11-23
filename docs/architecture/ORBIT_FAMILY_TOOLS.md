# Orbit-family tools

Reusable native tools for one-degree-of-freedom period annuli, their weakly damped
continuation, separatrix integrals and orbit-resolved Gramians. Every routine takes
caller-supplied callables or sampled tables; nothing depends on a registered system.
Results are numerical diagnostics, not enclosures or theorems.

## Native ownership

Namespace `dynsyslab::math`; headers in `native/include/dynsyslab/math`.

| Module | Responsibility |
|---|---|
| `OrbitIntegrator` | Adaptive Gragg-Bulirsch-Stoer step (order 16 extrapolation) and local table quadrature. |
| `OrbitFamily` | \(T(E)\), \(\Omega(E)\), action, equal-phase samples, Fourier data, \(\Omega'(E)\), monotone inverse, shell ladder, logarithmic separatrix fit. |
| `DampedOrbit` | Averaged \(D(E)\), \(\chi\), \(S(E)\), WKB and Laplace-limit principal eigenfunctions, isochrons, isostables, annular eigenfunctions, shell weight. |
| `MelnikovTools` | Connection integrals with exponential tails, damping and forcing integrals, threshold, lobe flux, capture probabilities. `sampledMelnikov` is unchanged. |
| `OrbitGramian` | Modal Cauchy Gramians, oscillator Hankel law, undamped and damped orbit Gramians, Jordan K-ROM response. |
| `OrbitQuery` | JSON operations for `dynsyslab-cli mathematics`. |

The core `integrate` API binds `ODESystem` and output grids, while section events need
local steps, so `OrbitIntegrator` is a small separate stepper.

## Conventions

- **Phase.** Orbits start on the ray `centre + r * section_direction`, which is phase
  zero. Use a reversibility axis to measure phase from a turning point. The first
  bracketed energy crossing on the ray is used.
- **Action.** \(I=\frac1{2\pi}|\oint p\,dq|\) for `canonical_pair`, or the area of
  \(\rho\,dx\wedge dy\) through the Poincaré-lemma primitive around the centre when a
  `density` is declared. The density must be smooth on the star-shaped region.
- **Separatrix fit.** \(T\sim\lambda^{-1}\log(C/\delta)\). The fit \(\lambda\) is the
  saddle exponent divided by the number of saddle passages per period.
- **Damped phase.** With \(\dot E=-dQ\), \(D=\langle Q\rangle\) and
  \(\chi=\Omega^{-1}\int(Q-D)\,d\phi\),
  \[
  S(E)=\int_{E_{\min}}^{E}\frac{\Omega-\Omega_0}{D}\,dE,\qquad
  \psi\approx\exp\bigl(a(E)+i[\phi+S/d+S'\chi]\bigr),\qquad a'=-\frac{\operatorname{Re}\lambda}{dD}.
  \]
  The measure is \(dE\): then \(\Omega-S'Q+\Omega S'\partial_\phi\chi=\Omega_0\), and the
  shell weight \(\sqrt{2\pi dD/|n\Omega'|}\) follows from \(S''=\Omega'/D\) at resonance.
  An integral in \(dJ\) with an energy-rate \(D\) would be off by a factor \(\Omega\).
  The amplitude is normalized to zero at the last table energy.
- **Annular eigenfunction.** \(W=\psi^n|\psi|^s\), \(s=(iw-n\lambda)/\operatorname{Re}\lambda\),
  has eigenvalue \(iw\) and \(|W|=1\) for real \(w\). It is defined only where \(\psi\ne0\).
- **Melnikov.** \(M(t_0)=-dD+F|A|\cos(wt_0+\arg A)\), \(D=\int Q\,dt\),
  \(A=\int f_0\wedge Y\,e^{iwt}dt\). The threshold is \(F/d=D/|A|\).
- **Gramians.** Entry order is `[EE, Ephi, phiphi]`. Modal matrices are row-major.

## JSON operations

Schema `dynsyslab.math-query/v1`. Expressions use the suite encoding
`analytic-state-expression/v1` with `state_dimension` 2; below, `<H>` and `<Q>` stand
for such objects. Complex values are `[real, imag]`.

| Operation | Required input | Optional input |
|---|---|---|
| `orbit_family` | `hamiltonian`, `centre`, `maximum_radius`, `energies`; or `table` with `energies`, `frequencies` | `section_direction`, `canonical_pair`, `density`, `phase_samples`, `maximum_harmonic`, `include_samples`, `include_fourier`, `target_frequency`, `drive_frequency`, `maximum_shell`, `separatrix_energy`, `separatrix_points` |
| `damped_averages` | `damping_scale`, `equilibrium_frequency`, `real_eigenvalue`, and either `table` (`energies`, `frequencies`, `damping`) or a declared system with `dissipation` | `q_samples` with `frequency`, `include_chi` |
| `melnikov_homoclinic` | `hamiltonian`, `midpoint`, `forward_saddle`, `backward_saddle`, `forward_rate`, `backward_rate`, plus `dissipation` and/or `forcing` with `forcing_frequency` | `damping_scale`, `forcing_scale`, `tail_radius`, `maximum_time` |
| `capture_probabilities` | `lobe_integrals` (nonnegative) | |
| `modal_gramian` | `eigenvalues`, `inputs`, `outputs` (rows are channels) | `horizon` |
| `orbit_gramian` | `stiffness`, `damping`; or `collocated`, `input_variance`, `action_twist`, `horizon`; or `input_variance`, `frequency_slope`, `horizon` | `output_variance`, `input_cross`, `input_phase`, `output_cross`, `output_energy` |
| `damped_orbit_gramian` | `energies`, `input_variance`, `output_variance`, `damping`, `frequencies`, `damping_scale`, `initial_energy`, `target_energy`, `final_energy` | |
| `jordan_response` | `modes` with `harmonic`, `amplitude`, `slope` | `real_observable`, `times` with `frequency` and `action_twist`, `lift_harmonic` |

```json
{"schema":"dynsyslab.math-query/v1","operation":"orbit_family",
 "hamiltonian":"<H>","centre":[0,0],"maximum_radius":3.14159,
 "energies":[-0.5,0,0.5],"drive_frequency":1.6,"maximum_shell":3}
```

```json
{"schema":"dynsyslab.math-query/v1","operation":"melnikov_homoclinic",
 "hamiltonian":"<H>","midpoint":[0,2],
 "forward_saddle":[3.141592653589793,0],"backward_saddle":[-3.141592653589793,0],
 "forward_rate":1,"backward_rate":1,"dissipation":"<Q>"}
```

```json
{"schema":"dynsyslab.math-query/v1","operation":"modal_gramian",
 "eigenvalues":[[-0.1,1],[-0.1,-1]],"inputs":[[1,1]],"outputs":[[0.5,0.5]],"horizon":20}
```

## Scope and limits

- The energy list must lie in one period annulus. Crossing a separatrix or another
  centre is not detected except through bracketing or closure failure.
- Averaging and WKB forms are first order in \(d\) and degrade near the separatrix and
  at the equilibrium. No remainder is computed.
- The Laplace limit stops at `settleRadius`, leaving a relative error of that order.
- Melnikov tails are leading-order exponential estimates. The midpoint must lie on the
  connection. A simple zero is a first-order prediction for small perturbations.
- Capture probabilities are ensemble ratios under the adiabatic phase assumption.
- Sampled-table inversion is linear. Declared Hamiltonians refine roots on the flow.
- Modal Hankel values use dense Jacobi, limited to 128 modes. Very small values of
  ill-conditioned Cauchy matrices are not resolved to high relative accuracy.
- Jordan response is the linear response at one action, for a real observable.

## Algorithm references

- Gragg, SIAM J. Numer. Anal. 2 (1965): https://doi.org/10.1137/0702030
- Bulirsch and Stoer, Numer. Math. 8 (1966): https://doi.org/10.1007/BF02165234
- Hairer, Nørsett and Wanner, *Solving Ordinary Differential Equations I*, §II.9:
  https://doi.org/10.1007/978-3-540-78862-1
- Arnold, *Mathematical Methods of Classical Mechanics*, §50:
  https://doi.org/10.1007/978-1-4757-2063-1
- Arnold, Kozlov and Neishtadt, *Mathematical Aspects of Classical and Celestial
  Mechanics*, chapter 5: https://doi.org/10.1007/978-3-540-48926-9
- Mauroy and Mezić, IEEE Trans. Automat. Control 61 (2016):
  https://arxiv.org/abs/1408.1379
- Melnikov, Trans. Moscow Math. Soc. 12 (1963).
- Guckenheimer and Holmes, *Nonlinear Oscillations, Dynamical Systems, and Bifurcations
  of Vector Fields* (1983), §4.5: https://doi.org/10.1007/978-1-4612-1140-2
- Neishtadt, J. Appl. Math. Mech. 39 (1975):
  https://doi.org/10.1016/0021-8928(75)90060-X
- Henrard, Celestial Mechanics 27 (1982): https://doi.org/10.1007/BF01228946
- Rom-Kedar and Wiggins, Arch. Rational Mech. Anal. 109 (1990):
  https://doi.org/10.1007/BF00375090
- MacKay, Meiss and Percival, Physica D 13 (1984):
  https://doi.org/10.1016/0167-2789(84)90270-7
- Moore, IEEE Trans. Automat. Control 26 (1981):
  https://doi.org/10.1109/TAC.1981.1102568
- Antoulas, *Approximation of Large-Scale Dynamical Systems* (2005):
  https://doi.org/10.1137/1.9780898718713
- Golub and Van Loan, *Matrix Computations*, cyclic Jacobi, §8.5.
