# Native mathematics toolkit

`dynsyslab_math` is an independently linkable C++ library with public headers
under `native/include/dynsyslab/math/`. It uses foundation's existing
arbitrary-length FFT and JSON, and privately uses Boost integer arithmetic.
It owns no trajectories, archives, fitting or UI. Unpublished case-study
derivations and data remain outside this module.

## Components

| Header | Operations | Limits of the result |
|---|---|---|
| Arithmetic.hpp | Exact rationals, prime valuations, p-adic distance/digits, compatible phase residues, continued fractions, Stern–Brocot runs, Ford circles, mediants, multiplicative word distance | No inference of rationality from measurements. |
| Topology.hpp | Components, spanning forest, first Betti number, degrees, Dijkstra distances, lower-star H0 persistence | Caller supplies an undirected multigraph; no certified Reeb/separatrix extraction. |
| HarmonicAnalysis.hpp | FFT coefficients, evaluation, derivatives, periodic primitive, convolution, cyclic sectors, Parseval norm, omitted-term bound | Finite coefficient bank; the bound excludes the unknown remainder of a sampled function. |
| Resonance.hpp | Detuning, phase drift, beat time, record-small primitive candidates | Floating-point near-resonance diagnostics, not exact rationality certification. |
| Distributions.hpp | Delta jets, test pairing, smooth multiplication, orbit/coarea quadrature, regular spectral density, Gaussian regularization | No singular distribution products; critical values and atoms need separate treatment. |
| SpectralFilters.hpp | Shared periodic rational kernels, Fourier packets, weighted orbit spectra, finite-frequency prediction bounds | Supplied finite coefficients and quadrature; no automatic continuous-spectrum certification. |
| PeriodicAtlas.hpp | Frequency and Fourier-coefficient interpolation across an orbit parameter | Consistent phase origins required; extrapolation rejected. |
| AnalyticFourier.hpp | Closed-form Jacobi and logarithmic-cosine coefficients and infinite-tail majorants | Bounds omit floating-point error; Jacobi parameter is m=k². |
| PeriodicSnapshots.hpp | Signed-delay pairs from Fourier banks and explicit orbit masses | Finite supplied series and quadrature; no inference of an invariant measure. |
| TorusFourier.hpp / Perturbation.hpp | Product Fourier analysis, cohomological equations, finite small divisors, homogeneous normal forms, Melnikov and symplectic defects | See [analytic tori and perturbations](ANALYTIC_TORI_AND_PERTURBATIONS.md); no automatic KAM certificate. |
| InverseScattering.hpp | Cell transfer matrices, finite-box eigenvalues, KdV scattering evolution, reflectionless inverse | See [inverse scattering](INVERSE_SCATTERING.md); boundary conditions and Lax/Koopman distinction are explicit. |
| CyclicCocycle.hpp | Periodic matrix characters, quotient projections, ordered products, transverse rollouts, finite residue-cell geometry | No nonlinear closure or infinite-limit certificate; computed eigenvalues use the optional native numerics backend. |
| Query.hpp | Versioned JSON adapter | Stateless, with no file or archive access. |

Rationals use integer or numerator/denominator strings, at most 4096 characters
per supplied integer. Denominators normalize positive. Prime arguments are
checked and capped at one million; digit precision and phase depth at 4096.
Infinite zero valuation is an empty optional, or JSON null accompanied by
`valuation_is_infinite:true`. For x=p^v*u, digits describe u in increasing
powers of p. Distance is p^(-v_p(a-b)), returned as an exact rational.
Phase compatibility checks upper mod p^depth=lower with valid residue ranges.
It tests one finite level, not a whole infinite solenoid. The logarithm of
`multiplicativeDistanceProduct` is a prime-weighted word metric, not a p-adic
metric.

Graphs support loops, parallel edges, isolated vertices and zero lengths.
b1=edges-vertices+components. Unreachable distances and essential H0 deaths
are null. Lower-star edges enter at max(endpoint values); H0 uses the elder
rule and vertex-index ties, retaining zero-length bars. Limits: one million
vertices, ten million edges in C++; each JSON array is capped at one million.

Fourier convention: c_n exp(i*n*phase), phase=Omega*t+phase0, Omega>0.
Periodic primitives reject nonzero mean. Samples cover one physical period
without repeating the endpoint. The requirement 2*maximum_harmonic<N removes
the Nyquist sign ambiguity, not aliasing from unresolved harmonics.
Harmonic magnitude is capped at one million, derivative order at 64, and
convolution at ten million term products. The omitted-term bound is the sum
of absolute values of omitted *known* coefficients.

Delta jets mean sum c_j*delta^(j)(x-a). Test jets contain actual derivatives,
not factorial-scaled Taylor terms. Pairing is sum (-1)^j*c_j*phi^(j)(a).
Smooth multiplication follows the distributional Leibniz rule; maximum
order is 64. Gaussian derivatives allow order 32, use log scaling, and reject
overflow. A Gaussian regularization is not a pointwise value of delta.

Orbit pairing is bilinear with exp(+i*n*phase), no test conjugation:
normalized integration against dt/T, or against dt otherwise. Coarea sums
phi*surface_weight/|grad H| and requires positive gradient norms.
Spectral density sums rho*|c_n|^2/|n*Omega'| over caller-supplied regular roots.
It does not locate roots, handle critical roots or add spectral atoms.

Resonance orders m,n mean detuning n*B-m*A for angular frequencies A,B.
Phase drift is |detuning|*horizon and beat time 2*pi/|detuning|.
Candidate search retains record-small detunings from nearest numerators
as denominators increase, with both orders bounded; it is not exhaustive.

## CLI and Python

```text
dynsyslab-cli mathematics query.json
python -m dynsyslab_client mathematics query.json
```

```json
{"schema":"dynsyslab.math-query/v1","operation":"padic",
 "a":"4/9","b":"1/9","prime":3,"precision":8}
```

Responses contain schema `dynsyslab.mathematics/v1`, operation and result.
This example returns distance "3", valuation -2 and the unit digits of 4.
The existing CLI error/exit contract applies.

```python
from dynsyslab_client.mathematics import query
answer = query({"operation": "continued_fraction", "value": "355/113"})
print(answer["result"]["convergents"])  # ["3", "22/7", "355/113"]
```

The optional Python client writes a temporary request and invokes the native
CLI. Pass `executable="build/native-suite-fast/bin/dynsyslab-cli.exe"` for a
source build. Python contains no duplicate scientific algorithms.

## Request fields

Every request needs schema and operation. Complex values are a real number
or [real,imag]. Numbers must be finite and integer fields integral.

| operation | Fields |
|---|---|
| rational | a,b: exact strings; b nonzero because ratio is also returned |
| padic | a,b: exact strings; prime,precision |
| phase_tower | prime; digits: least-significant-first branch choices (at most 512) |
| phase_cover_flow | prime,digits; base_phase_turns and elapsed_base_periods: exact rational strings |
| phase_cover_character | prime,digits; base_phase_turns and frequency_ratio: exact rational strings |
| phase_compatibility | lower,upper: integer strings; prime,lower_depth |
| continued_fraction | value: positive rational string |
| ford | a,b: rational strings |
| multiplicative_distance | a,b: nonzero rational strings with positive ratio |
| resonance | frequency_a,frequency_b,maximum_order; optional horizon=0 |
| graph | vertices,edges as [from,to] or [from,to,length]; optional source,vertex_values |
| fourier | action and fields below |
| analytic_fourier | kind: sn/cn/dn/log_cosine/kepler_cos_eccentric/kepler_sin_eccentric; parameter, maximum_harmonic; optional angular_frequency=1 |
| distribution | action and fields below |
| orbit_pairing | test_samples,harmonic,period; optional normalized=true |
| coarea | test_values,surface_weights,gradient_norms of equal length |
| spectral_density | roots: source_density,frequency_derivative,harmonic,coefficient |
| gaussian_delta | x,epsilon; optional order=0 |

Fourier series objects contain angular_frequency and coefficients as
[harmonic,real,imag] rows; duplicates are rejected. Actions:
evaluate (series,phase); differentiate (series,optional order=1);
primitive (series); multiply (series,other); cyclic_sector (series,degree,residue);
norms (series,maximum_harmonic); from_samples (samples,angular_frequency,maximum_harmonic).

The two Kepler kinds expand cosine or sine of the eccentric anomaly in the
**mean anomaly**, which advances uniformly in physical time. Here `parameter`
is eccentricity, `0 <= e < 1`; `angular_frequency` is the mean motion. The native
Fourier–Bessel coefficients use Boost.Math Bessel functions, with the constant
coefficient `-e/2` for cosine and zero for sine. These queries return
`uniform_tail_bound: null`: no certified truncation bound is implemented for
this family. A null bound must not be interpreted as a zero error. This is the
classical Kepler expansion; see [Borghi (2024)](https://doi.org/10.3390/math12010154).

Distribution objects contain location and complex coefficients. Actions:
pair (distribution,test_derivatives); differentiate (distribution,optional
order=1); multiply_smooth (distribution,smooth_derivatives).

## Sources and verification

### Sampled cycles and finite cyclic spectra

`math/DiscreteCycles.hpp` provides reusable diagnostics for uniformly sampled
return maps and arbitrary finite periods (no dyadic assumption). The existing
compiled mathematics CLI and `dynsyslab_client.mathematics.query` expose:

- `sampled_cycles`: sample-major `values`, `dimension`, `candidate_periods`,
  `tolerance`; optional positive coordinate `scales`, `coordinate_periods`
  (zero for a linear coordinate), and `minimum_repeats` (default 4).
  Reports maximum/RMS recurrence defects, the smallest accepted **supplied**
  period, a folded cycle, and its RMS defect. `certified=false`: finite
  recurrence does not prove periodicity, stationarity, or a limiting attractor.
- `cycle_spectrum`: sample-major `values` for one supplied cycle and `dimension`.
  Returns native DFT coefficients, angles, observable masses and reduced
  denominators of the cyclic eigenangles. Uses the uniform measure on the
  cycle; `total_mass` is its mean squared observable norm. Coefficients are
  harmonic-major, with each complex value encoded as `[real, imaginary]`.
  `conductor_masses` groups squared masses by the reduced denominator of each
  character: its smallest cyclic period. These groups are invariant under
  cyclic shifts of the phase origin.
- `cycle_projection`: `values`, `dimension`, `quotient_period` dividing the
  cycle length. Averages over equal residues to give uniform conditional
  expectation onto the quotient characters. Returns quotient and lifted values,
  retained DFT indices, and errors. Arithmetic coordinates are required; unwrap
  angles appropriately before averaging. Projected values need not follow the
  underlying nonlinear physical dynamics.
- `cycle_alignment`: equal-size `values` and `candidate`, with `dimension`.
  Returns the best forward cyclic shift, aligned values and RMS/maximum error.
  Uses bounded exhaustive least squares (at most 20 million coordinate terms),
  without reversal or geometric deformation. This is registration, not conjugacy.
- `sampled_hankel_pairs`: `values`, `dimension`, `delays`; optional
  `sample_stride`, `delay_stride`, `snapshot_stride`, `columns`.
  Calls the canonical native Hankel builder and returns sample-major `x`,`y`
  for the spectral API. The successor advances one sample-stride tick, even
  when snapshot columns have a larger spacing. The first block is the oldest
  sample. Observable selection and physical time units remain caller-owned.

The finite-cycle character basis and DFT masses are standard harmonic analysis
on a cyclic group; see Eisner et al., *Operator Theoretic Aspects of Ergodic
Theory*, chapters 14, 17–18
([author text](https://www.math.uni-leipzig.de/~eisner/book-EFHN.pdf)).
Tests cover non-dyadic cycles, corrupted recurrence, angular wrapping,
insufficient repeats, Parseval mass and independent Hankel successor spacing.
The native Viewer presents already computed spectra, sections, regularized
measures and free predictions (`dynsyslab.research-figure/v1` documents
`cycle_spectrum_sequence`, `return_sections`, `regularized_measures`,
`prediction_comparison`; [Viewer reports and export](VIEWER_REPORTS_AND_EXPORT.md));
it does not fit models or infer periods.

Notebook packaging (standalone fragments, iframes, executed notebooks) is retired;
the native Viewer exports any payload as self-contained interactive HTML or SVG.

### Exact finite phase-cover clock

`advancePhaseCover` in `math/Arithmetic.hpp`, exposed as `phase_cover_flow`,
advances a finite prefix of a classical circle-cover translation flow.
Base phase is a rational number of turns in `[0,1)`; elapsed time is a signed
rational number of base periods. For phase `s`, address prefix `a` modulo
`p^d`, and elapsed time `t`, the result has base phase `s+t-floor(s+t)` and
address `(a+floor(s+t)) mod p^d`. `prefix_carry` is the signed integer added
to the unspecified higher digits. `return_count` counts base-section crossings.
Each returned level includes its exact `phase_turns` and `frequency_ratio`.

```json
{"schema":"dynsyslab.math-query/v1","operation":"phase_cover_flow",
 "prime":2,"digits":[1,1],"base_phase_turns":"5/6",
 "elapsed_base_periods":"1/3"}
```

This gives phase `1/6`, digits `[0,0]`, return count `1`, carry `1`, and level-2
phase `1/24`. After `p^d` base periods the displayed prefix repeats while the
unspecified higher address changes. `finite_prefix:true` deliberately does not
assert return of an infinite state or existence of any physical cascade.
Arithmetic remains exact for up to 512 supplied digits; no floating-point
phase-to-rational inference occurs. Source: the translation and return-section
description in [Clark (1999)](https://arxiv.org/html/math/9906150).

`phaseCoverCharacter` / `phase_cover_character` evaluates a character whose
`frequency_ratio` belongs to `Z[1/p]`. It reports `minimum_depth`: the first
factor circle determining this character. If the supplied prefix is too short,
`resolved:false` and `phase_turns:null` explicitly preserve the ambiguity.
Otherwise the complex value is `exp(2*pi*i*phase_turns)`. The reduced denominator
must be a power of the supplied prime; other denominators are rejected.
For `p=2`, base phase `1/6` and frequency ratio `1/2`, zero digits do not
determine a value; digit `[0]` gives `1/12` turns and `[1]` gives `7/12`.
The ratio `1` instead gives `1/6` turns in both cases. This is exact factor
character arithmetic, not a claim that a measured physical state supplies the
missing phase choices. The native tests also verify character evolution against
the exact suspension clock, including negative time.

The native mathematics tests check gluing, negative time, compatibility and
inverse evolution. Python contracts also check the CLI schema and hidden carry.
`tests/test_phase_cover_viewer.py` checks the analytic `phase_cover` query and
the native Viewer's phase-cover pages through the compiled CLIs. The library
and CLI have no browser dependency.

Fourier action `evaluate_many` evaluates a supplied `series` at `phases`.
Action `evolve` takes `series`, `time`, optional `initial_phase` and
`initial_time`, and returns complex values after native character evolution
\(a_n e^{in(\phi_0+\Omega(t-t_0))}\). Both retain the supplied finite series;
neither estimates missing harmonics or fits a law. The C++ API is
`evolveFourierSeries`; the optional thin client is
`mathematics.evolve_fourier_channels`. Direct evaluation is bounded to 100 million
term evaluations. The native Viewer plots already computed decoders, character masses and errors
(`coordinate_reconstruction`, `character_mass`, `error_comparison` research figures)
without estimating them.

See [Periodic spectral analysis and estimation](PERIODIC_SPECTRAL_ANALYSIS.md)
for the additional native APIs, CLI operations, thin Python adapters and sources.

- Gouvêa, *p-adic Numbers*, chapters 1–2:
  [Springer](https://doi.org/10.1007/978-3-030-47295-5).
- Ford, “Fractions,” *Amer. Math. Monthly* 45 (1938), 586–601:
  [paper](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/ford.pdf).
- Short, “Ford circles, continued fractions, and best approximation”:
  [author paper](https://arxiv.org/abs/0912.1997).
- Hatcher, *Algebraic Topology*, §1.A:
  [author text](https://pi.math.cornell.edu/~hatcher/AT/ATpage.html).
- Dijkstra, "A note on two problems in connexion with graphs":
  [DOI](https://doi.org/10.1007/BF01386390), nonnegative shortest paths.
- Edelsbrunner, Letscher and Zomorodian, “Topological Persistence and
  Simplification”: [DOI](https://doi.org/10.1007/s00454-002-2885-2);
  Edelsbrunner and Harer, *Computational Topology* (AMS, 2010).
- Katznelson, *An Introduction to Harmonic Analysis*, chapter I. FFT
  provenance remains in [TimeSeries and Fourier](TIME_SERIES_AND_FOURIER.md).
- Hörmander, *The Analysis of Linear Partial Differential Operators I*,
  distribution differentiation and smooth multiplication; Federer,
  *Geometric Measure Theory*, coarea formula.

The `invariant_readout_floor` query takes complex `target_component_means`,
nonnegative `component_weights`, and optional `allow_constant` (default true).
Assuming all source observables have exactly zero conditional mean on each
invariant component, it computes the squared-norm floor for a linear readout
of their time translates. With a constant this is the weighted variance of
the target means; without one it is their weighted second moment.
Weights are normalized. Approximate component means or quadrature weights
give a numerical estimate, not a certified continuous-measure bound.
Nonlinear readouts are outside this restriction. The C++ API is
`math/InvariantProjection.hpp`. This is orthogonal-projection contractivity;
see [Eisner et al., Remark 8.9](https://www.math.uni-leipzig.de/~eisner/book-EFHN.pdf).

`dynsyslab_math_tests` checks exact identities, ultrametric inequalities,
graph topology and persistence, Fourier calculus and FFT recovery, weak
distribution identities, coarea normalization, density, resonance and bad
inputs. `dynsyslab_multiwell_system_tests` separately checks
[native multiwell models](../MULTIWELL_SYSTEMS.md). No general theorem prover
or automatic global certification is claimed.

## Finite periodic matrix cocycles and residue partitions

`math/CyclicCocycle.hpp` exposes `cyclicCocycle` and `cyclicPartitionGeometry`;
both are available through `dynsyslab-cli mathematics` and the existing thin
`dynsyslab_client.mathematics.query` adapter.

For `operation: cycle_cocycle`, provide `dimension: d` and phase-major,
row-major `values` containing N matrices of shape d by d. Optional fields are
`quotient_period` (default N, must divide N), `initial` (d entries, default zero),
`steps` (default zero), and `phase` (initial cycle index, default zero). The
matrix sequence is projected by uniform conditional expectation onto residues
modulo the quotient period. Its normalized character coefficients are computed
by the native DFT and synthesized back by the native inverse DFT. Returned
synthesis errors check this convention. No state-space derivative is estimated
by this operation.

The update is v[k+1] = A[(phase+k) mod N] v[k]. `monodromy` is the ordered
product of N matrices with the last update multiplying on the left, even when
the retained quotient period is smaller. Outputs include `multipliers` when
the numerics backend is enabled, `trace`, `determinant_2d`, flattened
`deviations` including the initial vector, `projection`, and
`character_spectrum` including masses grouped by reduced denominator.
Matrix averaging does not commute with products and need not preserve stability.
Very small multipliers or a determinant obtained by cancellation may be below
floating-point resolution. Bounds limit dimension to 16, allocations to one
million input values and matrix/rollout work to 20 million scalar loop steps.

For `operation: cycle_partition_geometry`, provide sample-major physical
`values`, `dimension`, and optional integer `base` (default 2). N must be a
power of the base. At depth n, samples are grouped by their chronological index
modulo base^n. Each level returns cell diameters, maximum diameter, minimum
intercell sample distance, and a `singleton_terminal_level` flag. Optional
`images` with the same shape contain independently supplied one-step return
images; `sample_transition` reports their cyclic point-defect RMS and maximum.
Coordinate units/scaling and any angular unwrapping are the caller's choice;
the distances are Euclidean. Work is bounded to 4096 points and 30 million
pairwise component comparisons.

Alternatively supply `partition_periods`, an explicit nonempty list of at most
32 strictly increasing periods, each dividing the next and each dividing N.
For example, `[1,3,6,12]` diagnoses an odd-base doubling hierarchy; `[1,3,9]`
diagnoses two tripling levels. Do not supply `base` with this list. The list
may stop before N to retain non-singleton cells. The corresponding C++ entry
point is `cyclicTowerPartitionGeometry`. Outputs include `partition_periods`
and the integer `branching_factor` of each level (relative to one at the first
level). No prime-adic limit or continuation beyond supplied levels is inferred.

This operation always reports `certified_infinite_limit: false`. Index residues
enforce a finite combinatorial hierarchy by construction. Point images are not
set-image enclosures, and a terminal singleton diameter of zero is automatic.
An infinite odometer conclusion requires compatible partitions for a single
invariant set and map at every depth, together with a uniform shrinking proof.
See de Carvalho, Lyubich and Martens, *Renormalization in the Hénon Family I*,
§5 ([author manuscript](https://www.math.stonybrook.edu/preprints/ims05-07.pdf))
for a theorem obtaining such geometry under substantive renormalization
hypotheses. This finite diagnostic does not verify those hypotheses.

The exact adding-machine arithmetic is the native `phase_cover_flow` and
related mathematics queries; the native Viewer draws their results.

### Circle-cover characters with periodic forcing

`phase_cover_resonance` takes `prime` and an explicitly rational positive
`base_to_drive_ratio` string. It returns an exact generator
(`character_ratio`, `drive_harmonic`) of the resonant integer-character
relations, `minimum_depth`, `branch_classes`, and `fiber_blind`. Input
rationality is asserted by the caller, never inferred from floating samples.
The independently compilable C++ API is `phaseCoverResonance` in
`dynsyslab/math/Arithmetic.hpp`.

`phase_cover_forcing` takes `prime`, common `depth` (0..32), a finite positive
`base_angular_frequency`, up to seven `forcing_frequencies`, sparse
`coefficients` in the existing torus Fourier format, and nonnegative
`divisor_tolerance`. The first integer Fourier index k represents the
circle-cover character k/p^depth; subsequent indices refer to driving angles.
The C++ API `phaseCoverForcing` is declared in
`dynsyslab/math/PhaseCoverFourier.hpp` and composes the existing
`solveCohomological` implementation with exact character labels.

Outputs include `forcing`, `solution`, `retained`, `maximum_inverse_divisor`,
`retained_mean_square`, and per-mode reduced labels, minimal depths,
detunings, forcing magnitudes, and correction magnitudes. A retained mode has
null correction magnitude in the JSON output; its required division was not
performed. Retention is based on floating detuning and the supplied tolerance,
not a proof of exact resonance. `finite_truncation` is true and `kam_certified`
is false: no unknown Fourier tails, nonlinear remainders, physical invariant
embedding, or torus breakup are certified.

Both use the existing compiled `math-query` CLI and optional
`dynsyslab_client.mathematics.query` Python client. Standard sources are
[Burgos and Verjovsky, §2](https://arxiv.org/abs/1603.05676) for circle-cover
characters and [Pöschel](https://arxiv.org/abs/0908.2234) for the small-divisor
perturbation framework. The analytic query `phase_cover` samples a caller-chosen orbit tower and checks
its clock ratios (`|omega_j p^j / omega_0 - 1|` within a tolerance, else an error);
the native Viewer draws its `dynsyslab.phase-cover-geometry/v1` result, one panel per
level with the compatible phase points.


### Positive Fourier clock changes

`dynsyslab/math/TorusClock.hpp` supplies `analyzeTorusClock` and
`evolveTorusClock`; compiled operation `torus_clock` exposes both through the
existing mathematics client. Input `clock_rate` uses the sparse torus Fourier
contract. It represents $\rho=dt/ds$ for $d\theta/ds=\omega$, so physical
velocity is $d\theta/dt=\omega/\rho$. This orientation is explicit.

The analysis composes `solveCohomological` to obtain
$\omega\cdot\nabla h=\rho-\bar\rho$ and returns the straightened phases
$\eta=\theta+\omega h/\bar\rho$, with frequencies $\omega/\bar\rho$.
The algorithm is the finite Fourier specialization of the positive time-change
criterion in [Clark, Theorem 1.1](https://arxiv.org/abs/math/9906150).
It applies to torus clocks and finite circle-cover factors. It complements
`PeriodicClock.hpp`, which transforms observable coefficients using clock
quadrature/Taylor expansions, rather than computing time evolution and the
straightening map.

Required conditions: a positive real mean, exact Hermitian coefficient pairs,
and a positive lower bound $\bar\rho-\sum_{k\ne0}|\rho_k|$.
This positivity condition is sufficient, not necessary; failure rejects the
request without claiming that the clock is negative. Nonzero nonconstant terms
retained by `divisor_tolerance` also reject the request; they are never silently
removed. All claims concern the supplied finite polynomial. Floating bounds
are not interval certificates and do not control an unknown infinite tail.

For evolution, supply `initial_phases`, elapsed physical `time` values (signed
and unordered are accepted), and optional integer `characters`. The routine
inverts the monotone scalar clock integral using bounded bisection, reporting
the achieved physical-time residual. `relative_time_tolerance` defaults to
1e-12; the stopping bound is tolerance times max(1, absolute elapsed time).
Outputs include bare and corrected character values, original and straightened
phases, mean clock rate, sufficient lower/upper bounds, and the primitive's
Fourier l1 and mean-square norms. Samples share the same explicit time origin.

The corrected characters are orthonormal with the normalized clock-weighted
measure $\rho\,dm/\bar\rho$, where $m$ is Haar measure, not generally with
unweighted Haar measure. No infinite-limit, Hamiltonian preservation, physical
embedding, KAM persistence, or spectral-type certificate is asserted.

The native Viewer's `clock_character` research figure renders a native query
response against a supplied linear reference on the same mesh. The existing native
Fourier evolution operation can supply that reference; Viewer contains no
clock integration or regression.
