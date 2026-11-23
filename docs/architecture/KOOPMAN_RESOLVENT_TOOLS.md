# Koopman resolvent tools

Native, stateless tools in `dynsyslab::math` for the Laplace-domain view of Koopman dynamics. For
an orbit `y(t)=f(S^t x0)`, the one-sided Laplace transform is the resolvent acting on `f`:
`Y(s)=[(sI-L)^{-1} f] at x0` (Susuki, Mauroy & Mezic 2021, eq. 4). The tools read numbers and return
reports. They own no archive and fit no model.

| Module | Responsibility |
| --- | --- |
| `ResolventLinearAlgebra` | Eigen-free complex LU, SVD and eigenvalues |
| `LaplaceTransform` | Sampled transform, tail bound, region-of-convergence estimate, z-transform, Talbot and Bromwich inversion |
| `KoopmanResolvent` | Pole/residue models with Jordan blocks, Cauchy residues, pole lattices, resolvent norm, pseudospectrum |
| `GeneralizedLaplace` | GLA projections, harmonic estimates, stationary/transient split, pole-order probe |
| `KoopmanFrequencyResponse` | Harmonic response `H_n`, symmetric Volterra kernels, real-input harmonic balance, Floquet check |
| `ShellResolvent` | Orbit-table boundary values, Jordan orbit transfer function, damped limit, rational-kernel projector |
| `FrequencyGramian` | Frequency-domain Gramians, H2 by Parseval, H-infinity peak |
| `ResolventQuery` | JSON operations for `dynsyslab-cli mathematics` |

Classical channel transfer functions remain in `core/TransferFunction.hpp`. The H-infinity peak
reuses `core::hInfinityNorm` through a real embedding, so it needs the optional numerics library;
without it the JSON field `hinfinity` is `null`.

## Conventions

- **Transform.** `Y(s)=int_0^inf e^{-st} y(t) dt`. Complex values in JSON are `[real, imag]`; a bare
  number is real. Matrices are arrays of rows.
- **Pole model.** `{"eigenvalue":z,"laurent":[a1,a2,...]}` means `sum_q a_q/(s-z)^q`, which is
  `sum_q a_q t^(q-1)/(q-1)! e^(z t)` in time.
- **Region of convergence.** The abscissa of a given orbit and observable is the largest `Re(lambda)`
  with a nonzero residue. It can lie left of the principal exponent when a residue vanishes, as on a
  strong stable manifold. A value fitted from samples is an estimate.
- **Lattices.** Equilibrium `k.lambda`; limit cycle `k.nu + i m Omega`; torus `k.nu + i m.Omega`, with
  `k` in `N0`. A complex `nu` is defined modulo `i Omega`, and the lattice set does not change.
- **Harmonic response.** Under the complex input `u=u0 e^{iwt}` the order-`n` output is
  `H_n(w) u0^n e^{inwt}`. For an analytic plant `H_n(w)` is the diagonal symmetric Volterra kernel
  `H_n^sym(iw,...,iw)`, so it does not depend on `u0`.
- **Real input.** Under `u=a cos(wt)` harmonics mix: the second harmonic is
  `(a^2/4) H_2(iw,iw)` and the DC term is `(a^2/2) H_2(iw,-iw)`, plus higher-order corrections.
  Use `real_input` for this case.
- **Shell pairing.** `Y(s)=sum_n int rho(E)|g_n(E)|^2/(s-i n Omega(E)) dE`, the autocorrelation of the
  observable over an orbit family with density `rho`.

## JSON operations

Every request has `"schema":"dynsyslab.math-query/v1"` and an `"operation"`.

| Operation | Required | Optional |
| --- | --- | --- |
| `laplace_transform` | `times`, `values`, `s` | `derivatives`, `interpolant` (`linear`, `cubic`), `growth_rate` with `growth_constant`, `z` with `sample_interval` |
| `inverse_laplace` | `poles`, `times` | `method` (`talbot`, `bromwich`), `nodes`, `shift`, `damping_parameter` |
| `koopman_resolvent` | `poles` with `s`, or `lattice` | both |
| `resolvent_residue` | `poles`, `center`, `radius` | `maximum_order` |
| `generalized_laplace` | `times`, `values`, `eigenvalues` | `refine`, `amplification_limit`, `decomposition` |
| `pseudospectrum` | `matrix`, `real_axis`, `imag_axis` | |
| `harmonic_response` | `A`, `B`, `C`, `frequency` | `N`, `maximum_order`, `direction`, `amplitude`, `real_input`, `maximum_harmonic` |
| `volterra_kernels` | `A`, `B`, `C`, `arguments` | `N`, `channels` |
| `shell_resolvent` | `mode` and its fields | see below |
| `frequency_gramian` | `A`, `B`, `C` | `tolerance` |

`shell_resolvent` modes:

- `orbit`: `series` (`angular_frequency`, `coefficients` rows `[n,re,im]`), optional `phase`, `s`, `times`.
- `jordan`: `modes` (`harmonic`, `amplitude`, `slope`), `frequency`, `action_twist`, optional `s`, `times`.
- `boundary`: table (`energies`, `frequencies`, `densities`, `harmonics`), `frequency`, optional
  `slope_tolerance`, `damping` with `damping_scale`.
- `damped`: table, `s`, `damping_rate`, `damping_scale`.
- `projector`: `bandwidth`, `angle`, optional `order`, `poles` with `sample_interval`.

Laplace transform of a sampled signal, with a declared growth bound:

```json
{"schema":"dynsyslab.math-query/v1","operation":"laplace_transform",
 "times":[0,0.5,1,1.5,2],"values":[1,0.6065306597,0.3678794412,0.2231301601,0.1353352832],
 "s":[[1,0],[0.5,2]],"interpolant":"cubic","growth_rate":-1,"growth_constant":1,
 "z":[[1.2,0]],"sample_interval":0.5}
```

Inversion of a model with a double pole:

```json
{"schema":"dynsyslab.math-query/v1","operation":"inverse_laplace",
 "poles":[{"eigenvalue":[-0.5,2],"laurent":[1,0.5]},{"eigenvalue":[-0.5,-2],"laurent":[1,0.5]}],
 "times":[0.5,1,2],"method":"talbot","nodes":48}
```

Resolvent values and a limit-cycle pole lattice:

```json
{"schema":"dynsyslab.math-query/v1","operation":"koopman_resolvent",
 "poles":[{"eigenvalue":[0,1],"laurent":[1,2]}],"s":[[1,1]],
 "lattice":{"kind":"limit_cycle","exponents":[[-0.3,0]],"frequency":1,"maximum_order":2,"maximum_harmonic":2}}
```

Residue extraction on a circle:

```json
{"schema":"dynsyslab.math-query/v1","operation":"resolvent_residue",
 "poles":[{"eigenvalue":[0,1],"laurent":[2,3]},{"eigenvalue":[-1,0],"laurent":[1]}],
 "center":[0,1],"radius":0.5,"maximum_order":2}
```

Generalized Laplace Analysis:

```json
{"schema":"dynsyslab.math-query/v1","operation":"generalized_laplace",
 "times":[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
 "values":[2,1.60653,1.36788,1.22313,1.13534,1.08208,1.04979,1.0302,1.01832,1.01111,1.00674,1.00409,1.00248,1.00151,1.00091,1.00055,1.00034,1.0002,1.00012,1.00007,1.00005,1.00003,1.00002,1.00001,1.00001],
 "eigenvalues":[[0,0],[-0.5,0]]}
```

Pseudospectrum of a Jordan block:

```json
{"schema":"dynsyslab.math-query/v1","operation":"pseudospectrum",
 "matrix":[[-1,1],[0,-1]],"real_axis":[-1.5,-1,-0.5],"imag_axis":[-0.5,0,0.5]}
```

The worked example of arXiv:2603.05771, `x1'=a1 x1+x2^2`, `x2'=a2 x2+u`, with lift
`z=(x1,x2,x2^2)`, `a1=-1.3`, `a2=-0.7`:

```json
{"schema":"dynsyslab.math-query/v1","operation":"harmonic_response",
 "A":[[-1.3,0,1],[0,-0.7,0],[0,0,-1.4]],"B":[[0],[1],[0]],"C":[[1,0,0]],
 "N":[[[0,0,0],[0,0,0],[0,2,0]]],"frequency":1.9,"maximum_order":4,
 "amplitude":0.5,"real_input":true,"maximum_harmonic":6}
```

Symmetric Volterra kernels of the same lift:

```json
{"schema":"dynsyslab.math-query/v1","operation":"volterra_kernels",
 "A":[[-1.3,0,1],[0,-0.7,0],[0,0,-1.4]],"B":[[0],[1],[0]],"C":[[1,0,0]],
 "N":[[[0,0,0],[0,0,0],[0,2,0]]],"arguments":[[[0,1.9],[0,1.9]],[[0,1.9],[0,-1.9]]]}
```

Shell boundary value at `s=iw+0`:

```json
{"schema":"dynsyslab.math-query/v1","operation":"shell_resolvent","mode":"boundary",
 "energies":[1,1.25,1.5,1.75,2],"frequencies":[1,1.25,1.5,1.75,2],"densities":[1,1,1,1,1],
 "harmonics":[{"harmonic":1,"values":[1,1,1,1,1]}],"frequency":1.4}
```

Jordan-block orbit transfer function:

```json
{"schema":"dynsyslab.math-query/v1","operation":"shell_resolvent","mode":"jordan",
 "modes":[{"harmonic":1,"amplitude":[0.6,0.1],"slope":[0.2,-0.3]},
          {"harmonic":-1,"amplitude":[0.6,-0.1],"slope":[0.2,0.3]}],
 "frequency":1.1,"action_twist":-0.4,"s":[[0.2,1]],"times":[0,1]}
```

Damped-limit pairing:

```json
{"schema":"dynsyslab.math-query/v1","operation":"shell_resolvent","mode":"damped",
 "energies":[1,1.25,1.5,1.75,2],"frequencies":[1,1.25,1.5,1.75,2],"densities":[1,1,1,1,1],
 "harmonics":[{"harmonic":1,"values":[1,1,1,1,1]}],
 "s":[[0,1.4]],"damping_rate":[1,1,1,1,1],"damping_scale":0.01}
```

Rational-kernel projector as a resolvent combination:

```json
{"schema":"dynsyslab.math-query/v1","operation":"shell_resolvent","mode":"projector",
 "bandwidth":0.2,"order":2,"angle":0.5,"sample_interval":0.5,
 "poles":[{"eigenvalue":[0,1],"laurent":[1,0.2]}]}
```

Frequency-domain Gramians:

```json
{"schema":"dynsyslab.math-query/v1","operation":"frequency_gramian",
 "A":[[-0.3,2],[-2,-0.3]],"B":[[0],[1]],"C":[[1,0]]}
```

## Scope and limits

- **Sampling.** Exact panel integration removes the oscillation error, not the interpolation error of
  `y`. A finite window never proves a region of convergence. The tail bound is only as good as the
  declared or fitted growth rate.
- **Talbot** needs singularities inside the contour, roughly `|Im p| t < 0.83 N`, and reports
  `scope_verified`. **Bromwich** needs `f` continuous at `t`; its error is about `e^{-A}`, and it loses
  accuracy for strongly growing signals.
- **Residues.** A contour must enclose only the intended pole. On a torus with rationally independent
  frequencies, poles are dense on the imaginary axis, so a small contour there is not meaningful.
- **GLA** needs the eigenvalues as input. Non-stationary coefficients are ill-conditioned, because
  `e^{-lambda t}` amplifies noise and errors in slower modes. The window is limited by
  `amplification_limit`. `sequential_coefficient` is the classical ordered pass; `coefficient` also
  solves the coupled window equations unless `refinement_condition` exceeds `1e12`.
- **`H_n`** exists only where a steady periodic response exists. It is invalid at external resonance
  (`i n w` in the spectrum of `A`), past a bifurcation, and beyond the series radius in `u0`. The
  Floquet check covers the bilinear model only, and a lifted model covers only its closure error.
- **Volterra kernels** describe the weakly nonlinear regime around the expansion point.
- **Harmonic balance** truncates at `|k|<=K` and reports the dropped edge term; it does not bound the
  error or detect coexisting responses.
- **Shell values** assume one period annulus, a smooth table and simple roots. A root with vanishing
  twist or at a table edge makes the result `valid:false`. The Hermite interpolant assumes `Omega` is
  monotone within each panel.
- **Damped pairing** is a Lorentzian closure with a declared rate `gamma(E)`. It is not derived from a
  specific damped flow. Panels are exact only for piecewise-linear data.
- **Gramians** need a Hurwitz `A`. Cost is one dense factorization per quadrature node, so they suit
  small models.

## References

- Y. Susuki, A. Mauroy, I. Mezic, Koopman resolvent: a Laplace-domain analysis of nonlinear
  autonomous dynamical systems, SIAM J. Appl. Dyn. Syst. 20(4) (2021) 2013-2036, arXiv:2009.11544.
- Y. Susuki, N. Katayama, A. Mauroy, I. Mezic, On Koopman resolvents and frequency response of
  nonlinear systems, arXiv:2603.05771 (2026).
- I. Mezic, Spectral properties of dynamical systems, model reduction and decompositions, Nonlinear
  Dynamics 41 (2005) 309-325.
- M. Budisic, R. Mohr, I. Mezic, Applied Koopmanism, Chaos 22 (2012) 047510.
- S. Das, J. A. Yorke, Super convergence of ergodic averages for quasiperiodic orbits, Nonlinearity
  31 (2018) 491-501.
- W. J. Rugh, Nonlinear System Theory: The Volterra/Wiener Approach, Johns Hopkins (1981).
- C. Hayashi, Nonlinear Oscillations in Physical Systems, McGraw-Hill (1964).
- T. Kato, Perturbation Theory for Linear Operators, Springer (1995).
- L. N. Trefethen, M. Embree, Spectra and Pseudospectra, Princeton (2005).
- L. N. Trefethen, J. A. C. Weideman, The exponentially convergent trapezoidal rule, SIAM Rev. 56
  (2014) 385-458.
- A. Talbot, The accurate numerical inversion of Laplace transforms, IMA J. Appl. Math. 23 (1979)
  97-120.
- L. N. Trefethen, J. A. C. Weideman, T. Schmelzer, Talbot quadratures and rational approximations,
  BIT 46 (2006) 653-670.
- J. Abate, W. Whitt, Numerical inversion of Laplace transforms of probability distributions, ORSA
  J. Comput. 7 (1995) 36-43.
- L. N. G. Filon, On a quadrature formula for trigonometric integrals, Proc. R. Soc. Edinburgh 49
  (1930) 38-47.
- N. I. Muskhelishvili, Singular Integral Equations, Noordhoff (1953).
- M. J. Colbrook, C. Drysdale, A. Horning, Rigged Dynamic Mode Decomposition, arXiv:2405.00782.
- K. Zhou, J. C. Doyle, K. Glover, Robust and Optimal Control, Prentice Hall (1996).
- A. C. Antoulas, Approximation of Large-Scale Dynamical Systems, SIAM (2005).
- S. Boyd, V. Balakrishnan, A regularity result for the singular values of a transfer matrix and a
  quadratically convergent algorithm for computing its L-infinity norm, Systems Control Lett. 15
  (1990) 1-7.
- N. A. Bruinsma, M. Steinbuch, A fast algorithm to compute the H-infinity norm of a transfer
  function matrix, Systems Control Lett. 14 (1990) 287-293.
- G. H. Golub, C. F. Van Loan, Matrix Computations, 4th ed., Johns Hopkins (2013).
- M. R. Hestenes, Inversion of matrices by biorthogonalization, J. SIAM 6 (1958) 51-90.

## Data-driven resolvent tools

These modules work from samples or a declared vector field. Their operations are served by
`evaluateResolventDataQuery` in `ResolventDataQuery`.

| Module | Responsibility |
| --- | --- |
| `ResolventLeastSquares` | Thin complex SVD and minimum-norm least squares |
| `DiscreteResolvent` | z-domain pole/residue model, exact sampling of Jordan models, alias report |
| `PronyIdentification` | Matrix-pencil poles and residues, order scan, cluster and near-continuous warnings |
| `ForcedResponseEstimation` | `H_n` by harmonic projection and skew-product EDMD; subharmonics on an n-fold cover |
| `IsochronDecomposition` | Stationary/transient split by asymptotic phase, decay rate |
| `HarmonicContinuation` | Pseudo-arclength harmonic balance, folds, Floquet stability, branch loss |
| `ContractionCheck` | Logarithmic norms on a sampled box, or a bound from a caller-supplied interval |

| Operation | Required | Optional |
| --- | --- | --- |
| `discrete_resolvent` | `poles` with `sample_interval`, or `discrete_poles` | `z`, `indices`, `branches` |
| `prony_identification` | `values` | `sample_interval`, `order`, `pencil_parameter`, `rank_tolerance` |
| `forced_harmonic_estimate` | `frequency`, and `values` with `times`, `amplitude`, and/or `dmd` | `maximum_harmonic`, `complex_input`, `transient_fraction` |
| `subharmonic_response` | `times`, `values`, `frequency`, `amplitude`, `cover` | `maximum_index`, `transient_fraction`, `detection_tolerance` |
| `isochron_split` | `times`, `values`, `cycle` | `asymptotic_phase`, or `phase_map` with `initial_state`; `transient_fraction`, `phase_tolerance`, `include_series` |
| `harmonic_continuation` | `field` or `expressions`, `frequency`, `amplitude`, `variable`, `lower`, `upper` | `parameter`, `harmonics`, `step`, `minimum_step`, `maximum_step`, `maximum_steps`, `direction`, `output_state`, `tolerance`, `compute_stability`, `initial_coefficients` |
| `contraction_check` | `field` or `expressions` with `lower`, `upper`; or `interval_lower`, `interval_upper` | `weights`, `samples_per_axis`, `parameter` |

- **Fields.** `field` has one list of terms per state. Each term is
  `{"coefficient":c,"powers":[x...,u,p]}`. `expressions` are `analytic-state-expression/v1` objects of
  dimension `states+2`, with arguments `(x...,u,p)`.
- **Discrete poles.** `{"multiplier":mu,"laurent":[a1,...]}` means `sum_q a_q z/(z-mu)^q`, or
  `a_q C(k,q-1) mu^(k-q+1)` in the sequence. `X(z)=[z(zI-U)^{-1}f] at x0`, `U=exp(Lh)`, `mu=exp(lambda h)`.
- **Aliasing.** Sampling fixes `lambda` only modulo `2 pi i/h`. `aliasing.entries` marks poles above the
  Nyquist rate `pi/h`, and `collisions` lists distinct eigenvalues that sampling merges.
- **Subharmonics.** `H_{1/n}` is reported through `|H|`, `arg(H)` modulo `2 pi/n`, `H^n`, and all `n`
  branches. A shift by one forcing period rotates it by `e^{2 pi i/n}`.
- **Continuation coefficients** are for real forcing `a cos(wt)`. They equal `H_n (a/2)^n` only in the
  small-amplitude limit. `termination` is `parameter_bound`, `step_limit`, `corrector_failed`,
  `non_finite_field`, `tangent_jump` or `singular_tangent`; `branch_lost` marks the failures.
- **Contraction labels.** `sampled` is evidence at grid points. `caller_interval` holds for every
  Jacobian inside an enclosure the caller supplies and the tool trusts.

z-domain model of a sampled double pole, with its alias report:

```json
{"schema":"dynsyslab.math-query/v1","operation":"discrete_resolvent",
 "poles":[{"eigenvalue":[-0.5,1],"laurent":[1,0.5]},{"eigenvalue":[-0.5,13.566370614359172],
 "laurent":[0.2]}],"sample_interval":0.5,"z":[[1.2,0]],"indices":[0,3]}
```

Matrix-pencil identification:

```json
{"schema":"dynsyslab.math-query/v1","operation":"prony_identification",
 "values":[1.3,0.96026849,0.53190896,0.18367253,-0.00203647,-0.02440074,0.06492016,0.19718485,
 0.31467989,0.38518704,0.4030975,0.38162436,0.34180558,0.30296623,0.27716144,0.26794245,
 0.2722835,0.28391617,0.29649926,0.3056708,0.30972897,0.30921816,0.30594637,0.30195279],
 "sample_interval":0.5,"rank_tolerance":1e-06}
```

`H_n` by harmonic projection, and by EDMD on the skew lift of `x'=-x+u`:

```json
{"schema":"dynsyslab.math-query/v1","operation":"forced_harmonic_estimate",
 "times":[0.0,0.3926991,0.7853982,1.1780972,1.5707963,1.9634954,2.3561945,2.7488936,3.1415927,
 3.5342917,3.9269908,4.3196899,4.712389,5.1050881,5.4977871,5.8904862,6.2831853,6.6758844,
 7.0685835,7.4612826,7.8539816,8.2466807,8.6393798,9.0320789,9.424778,9.817477,10.2101761,
 10.6028752,10.9955743,11.3882734,11.7809725,12.1736715,12.5663706,12.9590697,13.3517688,
 13.7444679,14.1371669,14.529866,14.9225651,15.3152642,15.7079633],"values":[[0.7,0.0],
 [0.4949747,0.4949748],[-1e-07,0.7],[-0.4949747,0.4949748],[-0.7,0.0],[-0.4949748,-0.4949747],
 [0.0,-0.7],[0.4949748,-0.4949747],[0.7,1e-07],[0.4949748,0.4949747],[0.0,0.7],[-0.4949747,
 0.4949747],[-0.7,-0.0],[-0.4949747,-0.4949748],[-1e-07,-0.7],[0.4949747,-0.4949748],[0.7,-0.0],
 [0.4949747,0.4949748],[-0.0,0.7],[-0.4949748,0.4949747],[-0.7,0.0],[-0.4949748,-0.4949747],[0.0,
 -0.7],[0.4949748,-0.4949747],[0.7,1e-07],[0.4949748,0.4949747],[0.0,0.7],[-0.4949747,0.4949748],
 [-0.7,-0.0],[-0.4949747,-0.4949748],[1e-07,-0.7],[0.4949747,-0.4949748],[0.7,-0.0],[0.4949747,
 0.4949748],[-0.0,0.7],[-0.4949748,0.4949747],[-0.7,1e-07],[-0.4949748,-0.4949747],[-0.0,-0.7],
 [0.4949748,-0.4949747],[0.7,0.0]],"frequency":2,"amplitude":0.5,"maximum_harmonic":2}
```

```json
{"schema":"dynsyslab.math-query/v1","operation":"forced_harmonic_estimate",
 "frequency":1.5,"dmd":{"sample_interval":0.2,"states":[[[0.8,0.0]],[[0.7441978403,
 0.0139402466]],[[0.6904066951,0.0510952492]],[[0.6310168007,0.1037122171]],[[0.5611689281,
 0.1634609078]],[[0.478780513,0.2220320532]],[[0.3843976987,0.2717603217]],[[0.2808731068,
 0.3062113799]],[[0.1728871432,0.3206766971]],[[0.0663466167,0.3125286809]],[[-0.0322929684,
 0.2814017966]],[[-0.1167269302,0.2291813382]],[[-0.1814651811,0.159799094]],[[-0.2224050022,
 0.0788527176]]],"inputs":[[0.5,0.0],[0.4776682446,0.1477601033],[0.4126678075,0.2823212367],
 [0.3108049841,0.3916634548],[0.1811788772,0.466019543],[0.0353686008,0.4987474933],
 [-0.1136010473,0.4869238154],[-0.2524230523,0.4316046833],[-0.3686968578,0.3377315903],
 [-0.452036071,0.2136899401],[-0.4949962483,0.070560004],[-0.493739885,-0.0788728471],
 [-0.4483792082,-0.2212602216],[-0.3629661521,-0.3438830796]],"maximum_harmonic":1,
 "dictionary":[[1,0],[0,1]]}}
```

Subharmonic on the two-fold cover:

```json
{"schema":"dynsyslab.math-query/v1","operation":"subharmonic_response",
 "times":[0.0,0.3926991,0.7853982,1.1780972,1.5707963,1.9634954,2.3561945,2.7488936,3.1415927,
 3.5342917,3.9269908,4.3196899,4.712389,5.1050881,5.4977871,5.8904862,6.2831853,6.6758844,
 7.0685835,7.4612826,7.8539816,8.2466807,8.6393798,9.0320789,9.424778,9.817477,10.2101761,
 10.6028752,10.9955743,11.3882734,11.7809725,12.1736715,12.5663706,12.9590697,13.3517688,
 13.7444679,14.1371669,14.529866,14.9225651,15.3152642,15.7079633],"values":[[0.7,0.0],
 [0.4949747,0.4949748],[-1e-07,0.7],[-0.4949747,0.4949748],[-0.7,0.0],[-0.4949748,-0.4949747],
 [0.0,-0.7],[0.4949748,-0.4949747],[0.7,1e-07],[0.4949748,0.4949747],[0.0,0.7],[-0.4949747,
 0.4949747],[-0.7,-0.0],[-0.4949747,-0.4949748],[-1e-07,-0.7],[0.4949747,-0.4949748],[0.7,-0.0],
 [0.4949747,0.4949748],[-0.0,0.7],[-0.4949748,0.4949747],[-0.7,0.0],[-0.4949748,-0.4949747],[0.0,
 -0.7],[0.4949748,-0.4949747],[0.7,1e-07],[0.4949748,0.4949747],[0.0,0.7],[-0.4949747,0.4949748],
 [-0.7,-0.0],[-0.4949747,-0.4949748],[1e-07,-0.7],[0.4949747,-0.4949748],[0.7,-0.0],[0.4949747,
 0.4949748],[-0.0,0.7],[-0.4949748,0.4949747],[-0.7,1e-07],[-0.4949748,-0.4949747],[-0.0,-0.7],
 [0.4949748,-0.4949747],[0.7,0.0]],"frequency":4,"amplitude":1,"cover":2,"maximum_index":4}
```

Isochron split with a declared asymptotic phase:

```json
{"schema":"dynsyslab.math-query/v1","operation":"isochron_split",
 "times":[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,
 10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,16.5,17.0,17.5,18.0,18.5,19.0,19.5,
 20.0,20.5,21.0,21.5,22.0,22.5,23.0,23.5,24.0,24.5,25.0,25.5,26.0,26.5,27.0,27.5,28.0,28.5,29.0,
 29.5,30.0,30.5,31.0,31.5,32.0,32.5,33.0,33.5,34.0,34.5,35.0,35.5,36.0,36.5,37.0,37.5,38.0,38.5,
 39.0,39.5,40.0,40.5,41.0,41.5,42.0,42.5,43.0,43.5,44.0,44.5,45.0,45.5,46.0,46.5,47.0,47.5,48.0,
 48.5,49.0,49.5,50.0,50.5,51.0,51.5,52.0,52.5,53.0,53.5,54.0,54.5,55.0,55.5,56.0,56.5,57.0,57.5,
 58.0,58.5,59.0,59.5,60.0],"values":[[1.4945759,0.0],[0.2824336,0.2772417],[-1.0030843,
 0.3650196],[-1.74494,0.2752033],[-1.6155113,0.0919227],[-0.712908,-0.0828491],[0.5051018,
 -0.1771781],[1.4494651,-0.1718986],[1.6793382,-0.0950669],[1.1029639,0.0021794],[0.0105058,
 0.0732968],[-1.0712244,0.0948888],[-1.6304794,0.0706321],[-1.4099326,0.0226997],[-0.5230312,
 -0.0224389],[0.6045248,-0.0463596],[1.4383379,-0.0444145],[1.5872209,-0.0241326],[0.9854568,
 0.0011298],[-0.0789848,0.0193654],[-1.1021516,0.0246597],[-1.6021168,0.018121],[-1.3452632,
 0.00559],[-0.4549299,-0.0060646],[0.6479347,-0.0121258],[1.4436046,-0.0114723],[1.5581438,
 -0.0061215],[0.938815,0.0004392],[-0.1218179,0.0051132],[-1.1240715,0.0064067],[-1.5964005,
 0.0046472],[-1.3170683,0.0013723],[-0.4181111,-0.0016359],[0.6771079,-0.0031705],[1.4532302,
 -0.0029624],[1.5453133,-0.0015516],[0.910347,0.0001518],[-0.1527015,0.0013493],[-1.1436462,
 0.001664],[-1.5963905,0.0011913],[-1.2981101,0.0003358],[-0.3892638,-0.0004405],[0.7025572,
 -0.0008287],[1.463787,-0.0007647],[1.53643,-0.000393],[0.8863993,4.92e-05],[-0.1804993,
 0.0003558],[-1.162431,0.0004321],[-1.5975688,0.0003053],[-1.2812895,8.18e-05],[-0.3623892,
 -0.0001184],[0.7269211,-0.0002165],[1.4743054,-0.0001974],[1.5282635,-9.94e-05],[0.8634386,
 1.53e-05],[-0.2074694,9.38e-05],[-1.1807815,0.0001122],[-1.5987317,7.82e-05],[-1.2647593,
 1.99e-05],[-0.3359483,-3.18e-05],[0.7508571,-5.66e-05],[1.4845113,-5.09e-05],[1.519967,
 -2.51e-05],[0.8405543,4.6e-06],[-0.2341827,2.47e-05],[-1.1987747,2.91e-05],[-1.5995586,2e-05],
 [-1.2480414,4.8e-06],[-0.3095503,-8.5e-06],[0.7745252,-1.48e-05],[1.4943264,-1.31e-05],
 [1.5113201,-6.3e-06],[0.8175153,1.4e-06],[-0.2607794,6.5e-06],[-1.216424,7.6e-06],[-1.599964,
 5.1e-06],[-1.2310149,1.2e-06],[-0.2831002,-2.3e-06],[0.7979604,-3.9e-06],[1.503727,-3.4e-06],
 [1.5022667,-1.6e-06],[0.7942666,4e-07],[-0.2872894,1.7e-06],[-1.2337283,2e-06],[-1.5999252,
 1.3e-06],[-1.213652,3e-07],[-0.2565793,-6e-07],[0.8211665,-1e-06],[1.5127047,-9e-07],[1.492794,
 -4e-07],[0.7707989,1e-07],[-0.3137149,5e-07],[-1.2506836,5e-07],[-1.5994362,3e-07],[-1.1959489,
 1e-07],[-0.2299881,-2e-07],[0.8441396,-3e-07],[1.5212553,-2e-07],[1.4829007,-1e-07],[0.7471148,
 0.0],[-0.3400509,1e-07],[-1.2672853,1e-07],[-1.5984956,1e-07],[-1.1779084,0.0],[-0.2033325,
 -0.0],[0.8668738,-1e-07],[1.5293759,-1e-07],[1.4725885,-0.0],[0.7232198,0.0],[-0.3662905,0.0],
 [-1.2835287,0.0],[-1.5971032,0.0],[-1.1595352,0.0],[-0.1766196,-0.0],[0.8893629,-0.0],
 [1.5370642,-0.0],[1.4618601,-0.0],[0.6991204,0.0],[-0.3924265,0.0],[-1.2994092,0.0],[-1.5952594,
 0.0]],"cycle":{"angular_frequency":1.4,"coefficients":[[1,0.8,0],[-1,0.8,0]]},
 "asymptotic_phase":0.9}
```

Duffing frequency sweep through both folds:

```json
{"schema":"dynsyslab.math-query/v1","operation":"harmonic_continuation",
 "field":[[{"coefficient":1,"powers":[0,1,0,0]}],[{"coefficient":-0.1,"powers":[0,1,0,0]},
 {"coefficient":-1,"powers":[1,0,0,0]},{"coefficient":-0.4,"powers":[3,0,0,0]},{"coefficient":1,
 "powers":[0,0,1,0]}]],"frequency":0.6,"amplitude":0.3,"variable":"frequency","lower":0.6,
 "upper":2.2,"harmonics":1,"step":0.02,"maximum_step":0.04}
```

Sampled contraction check, and a bound from a supplied interval:

```json
{"schema":"dynsyslab.math-query/v1","operation":"contraction_check",
 "field":[[{"coefficient":-1,"powers":[1,0,0]},{"coefficient":-1,"powers":[3,0,0]},
 {"coefficient":1,"powers":[0,1,0]}]],"lower":[-2,-1],"upper":[2,1],"samples_per_axis":9}
```

```json
{"schema":"dynsyslab.math-query/v1","operation":"contraction_check",
 "interval_lower":[[-3,-1],[0,-4]],"interval_upper":[[-2,1],[0.5,-3]]}
```

### Scope and limits of the data tools

- **Prony** assumes simple poles and low noise. Set `rank_tolerance` above the data precision: the
  default `1e-9` overestimates the order on data rounded to 8 digits. Poles closer than `2 pi/N`
  are unreliable. A warning is a heuristic; its absence does not prove a discrete spectrum.
- **Harmonic projection** needs a settled, resolved record of whole periods. The bump window leaks
  slowly for short records.
- **Skew DMD** needs transient data that excite every dictionary function. Steady-state data are reported
  as infeasible. The result carries any closure error of the dictionary.
- **Subharmonics.** Detection is an energy threshold at the declared cover. It does not find the cover,
  and PhaseCoverFourier labels exist only for prime-power degrees.
- **Isochron split** trusts the supplied phase map, and checks it only against the tail of the same
  record. With no map, the phase is fitted from that record and is not an independent test.
- **Continuation** uses finite-difference Jacobians, a truncated harmonic basis and a fold location
  accurate to about the squared step. It detects folds only. Branch points, period doubling and
  Neimark-Sacker points are not classified, though Floquet radii are reported. Autonomous limit cycles
  are out of scope.
- **Contraction.** A sampled pass can miss a narrow expanding region. Only diagonal metrics are
  supported. A negative log norm gives a unique steady response only on a forward-invariant set;
  invariance is tested at samples and is not proved. Interval bounds use no outward rounding.

### References for the data tools

- G. de Prony, Essai experimental et analytique, J. Ecole Polytechnique 1(22) (1795) 24-76.
- Y. Hua, T. K. Sarkar, Matrix pencil method for estimating parameters of exponentially damped/undamped
  sinusoids in noise, IEEE Trans. Acoust. Speech Signal Process. 38(5) (1990) 814-824.
- T. K. Sarkar, O. Pereira, Using the matrix pencil method to estimate the parameters of a sum of
  complex exponentials, IEEE Antennas Propag. Mag. 37(1) (1995) 48-55.
- Y. Susuki, A. Mauroy, Z. Drmac, Koopman resolvents of nonlinear discrete-time systems: formulation and identification, Proc. ECC
  (2024) 627-632.
- E. I. Jury, Theory and Application of the z-Transform Method, Wiley (1964).
- M. O. Williams, I. G. Kevrekidis, C. W. Rowley, A data-driven approximation of the Koopman operator,
  J. Nonlinear Sci. 25 (2015) 1307-1346.
- A. T. Winfree, Biological rhythms and the behavior of populations of coupled oscillators, J. Theor.
  Biol. 16 (1967) 15-42.
- J. Guckenheimer, Isochrons and phaseless sets, J. Math. Biol. 1 (1975) 259-273.
- A. Mauroy, I. Mezic, J. Moehlis, Isostables, isochrons, and Koopman spectrum for the action-angle
  representation of stable fixed point dynamics, Physica D 261 (2013) 19-30.
- H. B. Keller, Numerical solution of bifurcation and nonlinear eigenvalue problems, in Applications
  of Bifurcation Theory, Academic Press (1977).
- M. S. Nakhla, J. Vlach, A piecewise harmonic balance technique for determination of periodic response
  of nonlinear systems, IEEE Trans. Circuits Syst. 23 (1976) 85-91.
- T. M. Cameron, J. H. Griffin, An alternating frequency/time domain method for calculating the
  steady-state response of nonlinear dynamic systems, J. Appl. Mech. 56 (1989) 149-154.
- W. Lohmiller, J.-J. E. Slotine, On contraction analysis for non-linear systems, Automatica 34 (1998)
  683-696.
- A. Pavlov, N. van de Wouw, H. Nijmeijer, Uniform Output Regulation of Nonlinear Systems: A Convergent
  Dynamics Approach, Birkhauser (2006).
- A. Pavlov, N. van de Wouw, H. Nijmeijer, Frequency response functions for nonlinear convergent
  systems, IEEE Trans. Automat. Control 52 (2007) 1159-1165.
- G. Soderlind, The logarithmic norm: history and modern theory, BIT 46 (2006) 631-652.
