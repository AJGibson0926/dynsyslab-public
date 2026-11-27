# Spectrum and eigenfunctions of a K-ROM

A Koopman reduced-order model is a spectrum and a set of eigenfunctions
(`native/include/dynsyslab/core/KoopmanModel.hpp`, the one class `core::KoopmanModel`), independent
of the method that found them; input maps B and readouts C are optional parts
([DATA_OBJECTS.md](DATA_OBJECTS.md)). This page describes the eigenfunctions as data. H-DMD is the
first method that produces them.

## Input: time-series records

A time-series method takes one or more `TimeSeriesRecord`s (`core/TimeSeriesRecord.hpp`) and
nothing else. A record holds:

- the samples, as a `core::TimeSeries`;
- one `ChannelDescriptor` per column: the name, and the catalog observable that produced it, such
  as `state` for a solve or `echo` after a coordinate change. The observable is empty for a raw
  measurement;
- an `origin` (`DataOrigin`): the producing method, its settings and its input `DataReference`s
  (path plus SHA-256);
- the producing system's configuration, if one exists;
- `stateSource`: the solver record the channels were computed from. It is absent for measured
  data, and then no state domain exists.

`workflow/native/src/TimeSeriesSource.cpp` builds a record from any archived solve or coordinate
change. H-DMD reads extra records from `additional_time_series`, and
`joinTimeSeriesRecords` joins them. The records must share one clock, because nothing is
resampled.

## Domains

Every representation is evaluated on a domain: the coordinates q it takes as input. `FunctionDomain`
(`core/FunctionDomain.hpp`) is the single type for this. It is used by expansions, interpolations,
interpolation domains and K-ROMs. A domain is:

- **Identity:** the dimension and `names`, one per coordinate. These are the catalog observables or
  channel names that `recognizeDomain` matches.
- **State map:** one optional `Observable` q(x) of the physical state. The identity (q = x, built by
  `FunctionDomain::physicalState`), a formula of x (`fromStateMap`) and catalog observables
  (`fromCatalog`) are all just observables. When the map is absent, q cannot be computed from x and
  the caller supplies it (measured or history coordinates, `FunctionDomain::supplied`).
- **Catalog maps** are stored by name (`stateMapFromCatalog`, with `names` the catalog observables)
  and bound from the system's `ObservableCatalog` by `bound(catalog)` when the domain is loaded
  without one. `Observable::fromRepresentation` and `KoopmanModel` bind them automatically when a
  system is attached.
- **Measurement map:** an optional frozen provider that produces q online from measurements,
  recorded as `providerIdentity` (SHA-256) and `providerArtifact`.

The archive format (`dynsyslab.function-domain/v1`) is unchanged: `state_map` is written as
`identity`, `expression` (with the formula), `catalog` or `none`. The former four exclusive
coordinate kinds map as follows:

| Former kind | State map | Measurement map |
|---|---|---|
| physical-state | the identity | — |
| analytic-observable | a formula, or catalog observables | — |
| causal-history | none | yes |
| measured | none | — |

Old archives that use the former kinds still read correctly.

`workflow/native/src/MeasurementMap.cpp` is the single loader for measurement maps. It reads either
a history-basis archive or an H-SVD archive and checks the identity the same way for every
consumer: controllers, input maps, target lifts, certificates and trajectory audits. An end-aligned,
single-channel, uncentered H-SVD archive is itself a measurement map. Coordinate SINDy fitted on its
coordinates records that archive as its provider, so a controller projects live windows through the
same basis.

## Output: spectrum and eigenfunctions

- **`Spectrum`** (`core/Spectrum.hpp`). The one spectrum class, used by every producer and by
  `KoopmanModel`. Each mode is stored in the time domain it was produced in: a continuous rate
  λ, or a sampled multiplier μ = e^{λΔt} with its logarithm branch. Conversions are explicit and
  lossless, and a zero multiplier is allowed. `eigenvalue(j)` (or `spectrum[j]`) and
  `eigenvalues()` give `Eigenvalue` objects with conjugate partners paired. The spectrum also
  carries the origin that produced it, each eigenpair's residual, its continuous part and an
  optional pseudospectrum ([below](#eigenpairs-spectrum-and-pseudospectrum)).
- **`KoopmanModel`** (`core/KoopmanModel.hpp`): the eigenfunctions of one spectrum, stored
  together because they share their samples and each representation's dictionary, window layout or
  base points. For each training sample it holds the eigenfunction values (`samples()`,
  `EigenfunctionSamples`), the label time and the input-record row the sample belongs to. It also
  carries the origin (method, runtime settings, input records) and an ordered list of named
  representations, each one multi-output `FunctionRepresentation`; the first is the one the model
  evaluates (`useRepresentation` chooses another).
- **`Eigenfunction`** (`core/Eigenfunction.hpp`): one eigenfunction, an `Observable` with an
  eigenvalue. `KoopmanModel::eigenfunction(j)` gives mode j with every representation;
  `psi.as("interpolation")` makes another the default, and `psi.representationNames()` lists them.
  Any `Observable` can hold several named representations the same way
  (`Observable::withRepresentation`).

### Eigenpairs, spectrum and pseudospectrum

| Class | Header | What it is |
|---|---|---|
| `Eigenvalue` | `core/Eigenvalue.hpp` | One eigenvalue: λ, μ, Δt, the logarithm branch through the spectrum. |
| `Eigenfunction` | `core/Eigenfunction.hpp` | One function with its named representations, sharing dictionaries and interpolation data with its siblings. |
| `Eigenpair` | `core/Eigenpair.hpp` | An `Eigenvalue`, its `Eigenfunction` (absent on a bare spectrum) and an optional `EigenpairResidual {value, definition, scope}`. |
| `Spectrum` | `core/Spectrum.hpp` | The point part (the eigenvalues and their pairs' residuals), the continuous part (`SpectralMeasure` components, `core/SpectralMeasure.hpp`; they have no eigenfunctions) and an optional `Pseudospectrum`. |
| `Pseudospectrum` | `core/Pseudospectrum.hpp` | Grid points z with the minimum residual at each (and the minimizers); `levelSet(ε)` is the sampled ε-pseudospectrum. |
| `KoopmanModel` | `core/KoopmanModel.hpp` | The K-ROM: a spectrum whose point part are its eigenpairs (`eigenpair(j)` joins eigenvalue, eigenfunction and residual; `for (const auto& pair : krom.pairs())` visits them all in mode order, assembling each on access, so nothing is copied or stored twice), optional B and C, and a spectral fit's `spectralDiagnostics()`. |

The residual is a property of the pair: it measures how well this function and this value satisfy
the eigen-equation together, so it sits on neither the eigenvalue nor the eigenfunction. Each
method stores its own definition; definitions are never mixed or converted.

| Definition | Producer | Value |
|---|---|---|
| `one-step` | the `koopman_spectral` family (ResDMD; Colbrook & Townsend, *Comm. Pure Appl. Math.* 77, 2024, in the EDMD dictionary of Williams, Kevrekidis & Rowley, *J. Nonlinear Sci.* 25, 2015) | ‖ψ(y) − μψ(x)‖ / ‖ψ(x)‖ as a weighted quadrature over the training snapshot pairs (scope `empirical-training`) |
| `generator` | `sindy` (Kaiser, Kutz & Brunton, *Mach. Learn.: Sci. Technol.* 2, 2021) | RMS over the training rows of ∇ψ·f − λψ for the published (calibrated) eigenfunction, in the units of ψ |
| `one-step-successor-mean` | history EDMD (`echo_history_koopman`) | the one-step residual over ‖ψ(y)‖, averaged over the training runs (its mode score) |

H-DMD (Arbabi & Mezić, *SIAM J. Appl. Dyn. Syst.* 16, 2017), DMD, EDMD, DMDc/EDMDc, boundary
eigenfunctions, `eigenfunction_expansion` and polynomial fits compute none: their pairs' residuals
are empty. H-DMD's per-representation value error and Koopman residual stay its method properties.

The pseudospectrum is computed evidence about the spectrum, not part of it: σ ⊂ σ_ε for every ε and
σ = ∩_ε σ_ε (Trefethen & Embree, *Spectra and Pseudospectra*, Princeton 2005). ResDMD evaluates it
on the `residual_grid` you give, in the plane of μ, with the definition and scope of its eigenpair
residuals. Spectral measures of measure-preserving systems (mpEDMD, Colbrook, *SIAM J. Numer.
Anal.* 61, 2023; Rigged DMD and ResDMD smoothed measures, Colbrook, Drysdale & Horning,
arXiv:2405.00782) are the continuous components; the trial space's invariance, wave packets,
stochastic residuals and RKHS adjoints are the krom's spectral diagnostics.

Selecting modes (`extract`, `krom_extract`) keeps the selected pairs' residuals and the operator's
continuous part and pseudospectrum; it clears the diagnostics of the original trial space. The
Viewer's K-ROM spectrum shows each eigenvalue with its residual and, when the grid is a rectangular
lattice, the pseudospectrum as a heat map of log₁₀ of the minimum residual under the eigenvalues.

### Representations: two kinds

Every representation is a `FunctionRepresentation` that evaluates all eigenfunctions at a point of
its own domain. There are two kinds, and they are the same two that any `Observable` uses. Each has
a name, unique within the object (`representationName`): `interpolation`; `time-delay` for an
expansion on series windows; `polynomial` for monomials; else the dictionary's one term kind.

**Expansion: ψ_j(q) = Σ_k c_jk φ_k(q).** The dictionary (`core/ExpansionDictionary.hpp`) is one
ordered collection of basis functions, and it may mix any term kinds:

- monomials;
- orthogonal polynomials: Chebyshev, Laguerre, Hermite and Fourier–Hermite;
- radial kernels: Gaussian, Laplacian, Lorentzian and Matérn;
- sine, cosine and exponential terms;
- Voronoi partitions;
- any analytic expression.

Each term may carry its own centre and scale. A method only chooses terms and supplies
coefficients:

- SINDy, EDMD and EDMDc fit monomials;
- EDMD may also use Gaussian radial-basis terms;
- DMD and DMDc use degree-one monomials;
- spectral fits convert their dictionary term for term, so the values are identical;
- H-DMD's window map, `time-delay`, is one coordinate term per window entry on the window domain
  (`windowMapFromModes`, `core/HankelRepresentations.hpp`), stored as the matrix T;
- `eigenfunction_expansion` takes the terms and coefficients from the job.

The method is recorded in the origin, not in the type. Gradients are exact for monomials, Gaussian,
trigonometric, exponential and expression terms, and finite differences otherwise. A dictionary
made only of monomials sharing one normalization uses the established polynomial arithmetic, so
SINDy results are unchanged.

**Interpolation: anchored to sampled values.** There is one framework for every observable built
from samples, eigenfunctions included, and no method owns its own:

- `core::interpolateSamples(states, samples, domain, values, method, neighbors)`
  (`core/FunctionInterpolation.hpp`) takes the full states x_k at the samples. It maps them
  through the domain's state map to the base points q_k = q(x_k), which the interpolation stores.
  The default domain is the full state (identity map). A domain without a state map takes the
  given points as q_k (measured or history coordinates).
- In the workflow, `sampledDomain` (`workflow/native/src/SampledDomain.hpp`) finds x_k by following
  the sampled series' origin back to the solver record of the dynamical system, matched by time.
  The setting `interpolation_domain` picks the domain:
  - `null` or `"state"`: the full state (the default);
  - a list of catalog observables, e.g. `["far_field_acoustic_echo", ...]`;
  - an expression definition q(x);
  - `"input"`: the series' own channels, computed from x when every channel is a catalog
    observable, otherwise the sampled values themselves.

  Without a producing system, only `"input"` exists.
- `Observable::fromRepresentation` and `Eigenfunction::fromRepresentation` wrap the result. The
  eigenfunctions hold it directly, one domain for every mode, named `state` or after its
  coordinates.

The framework is shared, but each observable chooses its own interpolant. The method and its
neighbour count belong to the interpolation object, so two observables of the same samples may
use different ones. The methods (`interpolationMethods()`) are:

| Method | Interpolant over the k nearest samples |
|---|---|
| `inverse_distance` (default) | Shepard weights 1/dᵖ, p = 2 |
| `nearest` | the nearest sample's value |
| `gaussian_kernel` | Nadaraya–Watson weights exp(−d²/2h²) |
| `local_linear` | weighted affine least squares at the query (moving least squares) |
| `rbf_gaussian`, `rbf_multiquadric`, `rbf_thin_plate` | local radial-basis interpolation through the k samples; the thin-plate spline adds an affine term |

A method is `{"name", "neighbors", "power", "width"}`, where a `width` of 0 means the mean
neighbour distance. To sweep, call `setInterpolationMethod(set, method, neighbors)`, evaluate
`set.errors()`, and repeat. The job method `eigenfunction_errors` does the same on an archive.
`interpolation_methods` lists the candidates, and defaults to every method. With
`include_expansion` it also scores the expansion. The output has one row of value error and one row
of Koopman residual per method. The method is set by H-DMD's `interpolation_method`, stored in the
archive, and can be overridden with `krom_extract`'s `interpolation_method`.

`KoopmanModel` has exactly these two backends. A spectral K-ROM step publishes a krom with both
representations of its eigenfunctions: the frozen-dictionary expansion and an interpolation over
the training states. DMD, EDMD and DMDc kroms store their expansion, so `loadKoopmanModel`
evaluates them in their own dictionaries. A step that adds to a K-ROM publishes a child krom that
references it and its modes ([DATA_OBJECTS.md](DATA_OBJECTS.md)).

**Eigenfunctions given by hand.** `eigenfunction_expansion` is a delivery method, not a kind of
eigenfunction. You choose dictionary terms and coefficients:

```json
{"name": "eigenfunction_expansion", "settings": {
  "domain": "state", "variables": ["r", "v"], "eigenvalues": [[-0.1, 1.7]],
  "dictionary": [{"monomial": [1, 0]}, {"exponential": [0, -1]}, {"expression": "sin(r)*v"},
                 {"gaussian": {"node": [1, 0], "width": 0.3}}, {"cosine": [2, 0], "phase": 0.1}],
  "coefficients": [[2, [0, -1], 0.1, [0, 1], 0]]}}
```

There is one coefficient row per eigenvalue, and each entry is real or `[re, im]`. The shorthand
`"eigenfunctions": [{"real": "r", "imag": "-v/1.7"}]` decomposes each text into dictionary terms
(polynomial parts become monomials). `"library": true` or `"library": ["x_1", "eta"]` loads the
source system's own library entries instead (see below); `"reference_state": [...]` picks the orbit
whose eigenvalues shell entries take. `domain` may instead list catalog observables, such as
`["echo", "echo_dot"]`. Expressions support `+`, `-`, `*`, `/`, `^` with a constant exponent, `exp`,
`log`, `sin`, `cos`, `tanh`, `sech`, `sqrt`, `pi`, `atan2(y, x)` and the elliptic integrals in the
parameter m = k² (`ellipk(m)`, `ellipe(m)`, `ellipf(phi, m)`, `ellipe(phi, m)`; Carlson 1995,
*Numer. Algorithms* 10:13–26). Every function has an exact symbolic derivative.

- **One domain, chosen where the eigenfunctions are made.** H-DMD's `interpolation_domain` sets it
  (the former `interpolation_domains` list is read as its first entry). The K-ROM loader uses
  the domain it finds and refuses a request for a different one; it never re-domains.
- **Domain recognition.** `KoopmanModel::recognizeDomain(observables, representation)` names the
  representation's domain and checks it against the observables of the data handed in. For
  example, echo windows are the `time-delay` representation's `window` domain.
- **No copies of the input.** H-DMD stores neither its input channels nor the states. The
  interpolation stores its base points, the one copy the framework allows, because they define the
  function. Stages that need the training states, such as refit reconstructions, find them
  through the K-ROM's origin.

### Errors

`KoopmanModel::errors(representation, input)` evaluates one representation at the training
samples. With no arguments it uses the first representation on its domain. Per
mode, it reports:

- the value error ‖R(p_k) − ψ_k‖ / ‖ψ_k‖;
- the Koopman residual ‖R(p_{k+1}) − μ^Δt R(p_k)‖ / ‖R(p_k)‖.

The interpolation reproduces its own samples, so its value error is zero by construction; for
interpolation, the Koopman residual is the meaningful number. H-DMD writes the default value error as
`eigenfunction_value_error` and attaches the Koopman residual to each eigenpair as its `one-step`,
`empirical-training` residual (`core/Eigenpair.hpp`): a residual belongs to an eigenvalue together
with the eigenfunction approximation it was measured on, not to the eigenvalue alone.

## Observables and eigenfunctions

An eigenfunction is an observable with an eigenvalue. Along a solution it evolves as
ψ(x(t)) = e^{λt} ψ(x(0)). The classes say exactly that.

| Class | Header | What it is |
|---|---|---|
| `Observable` | `core/Observable.hpp` | A function of the state. It has a `representation()` property (`representationKind()`: `expansion`, `interpolation` or `evaluator`). A formula observable (`fromExpressions`) is an expansion with one expression term per output; `expressions()` reads the closed form back from the expansion. |
| `Eigenfunction : Observable` | `core/Eigenfunction.hpp` | One scalar observable plus its `Eigenvalue` (`eigenvalue()`, `propagate(value, t)`). |
| `Eigenvalue` | `core/Eigenvalue.hpp` | λ, the multiplier e^{λΔt}, Δt, and the conjugate partner. `transition(t)` = e^{λt}. |
| `FunctionRepresentation` | `core/FunctionRepresentation.hpp` | The interface both representations implement: domain, outputs, `evaluate`, `jacobian`, `exactJacobian`, `selectOutputs`. It holds no eigenvalues. |
| `FunctionExpansion` | `core/FunctionExpansion.hpp` | Dictionary terms and coefficients, g_r(q) = Σ_k A_rk φ_k(q). |
| `FunctionInterpolation` | `core/FunctionInterpolation.hpp` | Sampled values at scattered points, with a selectable interpolant. |
| `Spectrum` | `core/Spectrum.hpp` | The eigenvalues of a set of eigenfunctions, stored with their time domain and branches. `eigenvalue(j)` gives one `Eigenvalue`. |

How the pieces connect:

- A representation may hold several outputs, so that one dictionary evaluation or one neighbour
  search serves every eigenfunction. Each `Eigenfunction` points at the shared representation and
  selects its own output.
- `Observable::fromRepresentation(descriptor, representation, catalog)` gives g(q(x)) on the
  physical state through the domain's state map (the optional `catalog` binds a catalog map); with
  no map it is a function of q. When the representation and the map are exact, the Jacobian is the
  exact chain rule, so partial derivatives, Jacobians and the Koopman generator ∇ψ·f work the same
  way for observables and eigenfunctions.
- `Eigenfunction::fromRepresentation(descriptor, representation, output, eigenvalue, catalog)`
  builds one eigenfunction the same way.
- Producers (SINDy, EDMD, DMD, DMDc, EDMDc, spectral fits, H-DMD, the library, jobs) return
  `EigenfunctionFitResult {functions, spectrum, diagnostics}`. That is the shared `FunctionExpansion`
  and, separately, its `Spectrum`, with the pairs' residuals when the fit computed them.
- `KoopmanModel(functions, spectrum, …)`, `KoopmanModel(InterpolationData, spectrum, …)` or
  `KoopmanModel(spectrum, samples, origin, representations)` holds them. `model.eigenfunction(j)`
  hands out an `Eigenfunction` with every representation attached (the evaluated one first). Each
  representation is a `FunctionRepresentation` on its own domain: an interpolation, a dictionary
  expansion, or the window map as a `FunctionExpansion` in coordinate terms of the window.

**Expansion tools** (`FunctionExpansion`, `core/ExpansionProjection.hpp`).

- Expressions are not a separate representation. Every expression is an expansion: at worst one
  expression term with coefficient 1. The expansion classes therefore live in the foundation
  library beside `Observable`.
- `FunctionExpansion::fromExpressions` turns analytic expressions into dictionary terms.
  Polynomial parts become monomials, and other pieces become expression terms. A hand-written
  eigenfunction and a fitted bank therefore share terms.
- `canonical()` rewrites normalized monomials as raw monomials, so equal functions have equal
  coefficients.
- `realExpressions()` and `expressionObservable(expansion, descriptor)` give the closed form as an
  expression `Observable`, when every term has one (every kind except Matérn kernels and
  partitions). That observable has `partialDerivative`, `flowDerivative` and
  `inputFlowDerivative`.
- `stackExpansions` and `mergeDictionaries` put several expansions in one dictionary.

**Projection onto eigenfunctions.** `projectObservable(g, eigenfunctions, samplePoints, tolerance)`
finds B in g_r ≈ Σ_j B_rj ψ_j, for an observable and eigenfunctions whose representations are
expansions. `projectExpansion` does the same directly on expansions. Both are canonicalized and
written in the union of their dictionaries.

- With no sample points, the fit is least squares on the dictionary coefficients. When every g_r
  lies in span{ψ_j}, this is an identity of functions, `exact` is set, and B is the analytic
  K-ROM reconstruction.
- With sample points, the fit is least squares over those points.

The result also holds the residual g − Bψ as an expansion in the merged dictionary, plus a
relative residual per output. For example, on the triangular system 2x₁ + 3x₂ − 3x₁² projects onto
(x₁, η) exactly with B = (2, 3), and x₂ alone leaves the residual −k x₁².

**System eigenfunction library.** Each system can carry hand-derived eigenfunctions in its
`DynamicalKnowledge::eigenfunctions`. They are stored as data, so systems do not depend on the
Koopman libraries. Each `EigenfunctionEntry` has:

- an id, a title and a kind;
- real and imaginary expression text in the state and parameter names;
- the eigenvalue as rate text;
- conditions, a derivation note and references;
- an `applicable` flag computed from the configuration.

There are three kinds of entry:

- **exact**: one eigenvalue on the stated region;
- **principal**: the linearization wᵀ(x − x*) at an equilibrium, where w is a left eigenvector of
  Df(x*). It is accurate near x* only;
- **shell**: exact on the stated region, but the eigenvalue depends on the orbit through a conserved
  quantity, for example i Ω(H) for a phase e^{iφ}. Along every orbit ψ(x(t)) = e^{iΩ(H)t} ψ(x(0)).
  Loading one needs a reference state, and the eigenvalue is that orbit's. (Across an energy
  band the fixed-eigenvalue object is the distribution δ(H − E) e^{iφ}; Mezić 2020, *J. Nonlinear
  Sci.* 30, section 11.)

`libraryEigenfunctions(system, ids, reference)` substitutes the configured parameters and returns
`Eigenfunction` objects on the physical state, named by entry id. `libraryEigenfunctionBank` returns
the same entries as one shared expansion plus spectrum, for building a K-ROM or projecting.
`eigenfunctionGeneratorResidual` checks any eigenfunctions against the configured equations,
using |∇ψ·f − λψ|. The job method `eigenfunction_expansion` with `"library": true` (or a list of
ids) publishes them as a K-ROM. The triangular benchmark registers x₁, η = x₂ + c x₁²/(b − 2a),
and the principal pair x₁ and x₂ at the origin; x₂ is the linearization of η (Brunton et al.
2016, PLOS ONE 11(2):e0150171; Kaiser, Kutz and Brunton 2021, Mach. Learn.: Sci. Technol.
2:035023). The other catalog entries live in `native/src/SystemEigenfunctionLibrary.cpp`, and
`system_eigenfunction_library_tests` checks every one against the configured equations:

| System | Entries |
|---|---|
| `nonlinear_pendulum` | `energy`; `libration_phase` e^{iφ}, φ = πF(α\|k²)/(2K), α from sin(θ/2) = k sin α, Ω = π√a/(2K) (shell); `rotation_phase` (shell, H > 2a); `principal` and `principal_conjugate` damped modes at the bottom |
| `duffing` (δ = 0) | `energy`; `outer_phase` (H > 0, x = A cn) and `well_phase` (either well, x = ±A dn) (shell) |
| `harmonic_oscillator` | `mode`, `mode_conjugate`: (λ + c/m)x + v, exact for any underdamped damping |
| `isotonic_oscillator` | `energy`; `isochronous_mode` Z = q² − H/ω² − iqp/ω with λ = 2iω (exact on the whole half-plane) |
| `morse_oscillator`, `poschl_teller`, `cosh_oscillator` | `energy`; `phase` (shell, closed form; the cosh phase is elliptic) |
| `radial_twist` | `action`; `rotation` q − ip, λ = i(ω + αI) (shell) |
| `kepler_planar` | `energy`, `angular_momentum`; `mean_anomaly` e^{iM}, M = E − e sin E from the state, Ω = √(μ/a³) (shell) |
| `euler_top`, `coupled_duffing_hamiltonian`, `triple_well_duffing`, `three_leaf_hamiltonian_clover` | conserved quantities (eigenvalue 0) |

Conservative entries require zero damping and no driver; their `applicable` flag is false otherwise.
Sources: Brizard (arXiv:1108.4970) and Jaffé and Brumer (1984, *J. Phys. Chem.* 88) for the
pendulum; DLMF §22.19 for the elliptic solutions; DeMarcus (1978), Kuru and Negro (2008,
arXiv:0709.4649), Delisle-Doray et al. (2019, arXiv:1812.11582) and Henrici and Kappeler
(arXiv:0802.4252) for the Morse, isotonic, Pöschl–Teller and cosh oscillators; Murray and Dermott
(1999), chapter 2, for Kepler's equation; Mezić (2005, *Nonlinear Dyn.* 41) for linear modes.

## H-SVD: Hankel basis and coordinate record

`hankel_svd` (also `hankel_pod`) takes the same time-series records as H-DMD and produces two
objects:

- **`HankelBasis`** (`core/HankelBasis.hpp`, `dynsyslab.hankel-basis/v1`): the window layout, U,
  Σ, the optional mean and the origin. `coordinates(window, scaling)` gives either
  q = Uᵀ(h − mean) (`energy`) or ξ = Σ⁻¹q (`whitened`, the history-basis convention).
- **A coordinate time series.** The output archive stores `time` and `values` = coordinates, one
  row per window, labelled by the window's newest sample when `window_alignment` is `end`. Its
  channels are `hsvd[k]`, `value_space` is `history-coordinates` and `state_source_archive` names
  the original solve. It loads as an ordinary time series, so it can feed H-DMD, SINDy or a
  controller.

For the same uncentered window matrix and rank, H-DMD's eigenfunctions are ψ = W⁻¹q, so they are a
fixed linear change of basis of the H-SVD coordinates. States are matched to the samples by time,
so an H-DMD fitted on H-SVD coordinates still interpolates on the original solve's full state.

### Window settings

Use the `window` object, `{"delays", "columns", "sample_stride", "delay_stride",
"snapshot_stride"}`, together with `window_alignment` (default `"end"`: a window is labelled by its
newest sample, so values use only history), for both methods. The legacy keys disagree:
`matrix_size` means the number of delays in H-SVD but the number of columns in H-DMD. They remain
readable only for old jobs.

## Adding representations with SINDy

Representations have names, unique within an object: `interpolation`, `time-delay` (H-DMD's window
map) and `polynomial` (a monomial dictionary expansion; other dictionaries are named after their
term kind). Jobs, the Viewer and K-ROMs choose a representation by name. A kind such as `expansion`
still works when only one representation has it.

The `sindy` method publishes a krom. What it is given decides what it does:

| Given | SINDy does | Publishes |
|---|---|---|
| `"eigenvalues": [[re, im], …]` | fits each eigenvalue's eigenfunction and normalizes it | a new krom |
| a spectrum (any archive with eigenvalues that is not a krom) | the same, for all modes or `modes` | a new krom |
| a krom (`krom`, `eigenfunctions`, `spectrum_archive`, or the preceding method's object) | fits the selected modes, calibrates them to the eigenfunctions' own values, and adds the polynomial | the krom's child |
| nothing, or `"eigenvalues": "from_data"` | also finds the eigenvalues from the training data | a new krom |

In H-DMD → SINDy, SINDy receives H-DMD's krom. It keeps its eigenvalues and publishes a child with
three representations. The polynomial is evaluated; `make_default: false` keeps H-DMD's.

**Fit.** The default `fit: "nullspace"` solves the generator equation alone, as published (Kaiser,
Kutz & Brunton, *Mach. Learn.: Sci. Technol.* 2, 035023, 2021). Given eigenfunctions, it then takes
the least-squares complex scale c_j = p_jᴴψ_j / p_jᴴp_j, where ψ are the eigenfunctions' values at
their samples. Eigenfunctions are defined up to such a factor. `fit: "anchored"` (opt-in) solves, for
each mode, min ‖Θ(X_a)ξ − ψ‖² + w‖(Γ − λΘ)ξ‖². Here X_a are the states at the samples (through the
origin), Θ the polynomial library and Γ its rate along the dynamics. The weight w is `anchor_weight`
(default 1) times tr(Θ_aᵀΘ_a)/tr(G(λ)). The anchored fit trades generator residual for agreement, so
when the library cannot represent the references it no longer returns an eigenfunction. H-DMD's
values, for example, cover only its last window (near rest for a decaying trajectory).

**Normalization.** Without reference values, or for a mode whose polynomial disagrees with its
references beyond `agreement_tolerance`, each eigenfunction has unit RMS over the training rows and
is real and positive at the first row. A disagreeing reference fixes no meaningful factor.
`metadata.agreement.calibration` records `reference` or `normalized` for each mode.

**Agreement.** `agreement_mismatch` records ‖c_j p_j − ψ_j‖/‖ψ_j‖ for each mode, at the samples.
Modes above `agreement_tolerance` (default 0.1) are listed in `metadata.agreement.flagged_modes`, and
SINDy warns.

**Eigenvalues from data** follow Kaiser, Kutz & Brunton (2021), Algorithm 1:
1. The initial guesses are the eigenvalues of the library generator K = Θ⁺Γ, the generator EDMD
   matrix (Klus et al., *Physica D* 406, 132416, 2020).
2. For each guess, SINDy takes the sparse null vector of (Γ − λΘ). It uses sequential thresholding,
   not their alternating-directions search.
3. It updates λ ← ξᴴKξ/ξᴴξ until λ settles.
4. It validates each eigenpair along the training series by the error of φ(x(t)) against c e^{λt}
   (their Eq. 47), with c fitted by least squares. The trivial λ = 0 and eigenvalues above
   `validation_tolerance` are dropped; `metadata.eigenvalue_estimate` records the candidates.

A finite library may not be closed under the dynamics, so this can miss eigenvalues. Use one
trajectory.

**Children.** Archives are immutable, so SINDy publishes a new krom: a child of the one it was
given, `krom.parent = {path, sha256, modes}`. It holds its own coefficients, its readouts, the
calibration and the agreement. Readers rebuild the parent's representations, samples, readouts and
maps from the parent archive, refusing it if it is missing or changed. Nothing large is copied:
H-DMD's window map and samples stay in H-DMD's krom.

**Representation choice.** `krom_extract` (and `control_regression`) take
`eigenfunction_representation`; the default is the krom's first. The child krom lists it first in
its `order` (also written as `eigenfunction_representation_name`). A K-ROM's input map needs ∇ψ
with respect to the state, so a `time-delay` K-ROM (a function of past windows) refuses a state and
its exact map. To compare representations in control, change this setting in the K-ROM stage and
rerun the experiment.

## Archive: the krom (`dynsyslab.krom/v1`)

The arrays are:

- `continuous_eigenvalues`, `discrete_eigenvalues`, and the spectrum's evidence (`eigenpair_residual`,
  `pseudospectrum_*`, `continuous_spectrum_json`; `core/SpectrumArchive.hpp`);
- `eigenfunction_values`, `training_time`, `training_sample_indices`;
- `interpolation_points`, the interpolation's base points q_k (samples × domain dimension);
- `window_map` (modes × rows): row j is eigenfunction j's own linear map from a raw window, so each
  eigenfunction stands alone.

`loadKoopmanModel(path).eigenfunction(j).as(representation)` gives one eigenfunction on its own,
with every representation it has and `representation` as its default. Its window map reads raw
windows built by `HankelWindowLayout::window_at` from the channels and layout its description
names, so new windows of the right length evaluate without the producing method. The archive is
the `krom` data object ([DATA_OBJECTS.md](DATA_OBJECTS.md)), which also holds the readouts and
input maps.

`metadata_json.krom` holds the origin (the input series, by path and SHA-256), the
input channels, the representations in order and any parent; H-DMD's default errors are in
`metadata_json.eigenfunction_errors`. The interpolation's
description holds its domain, with its state map, and its interpolant. The window map's
description holds the window layout: the record columns, lags, scales, delays and strides. Together
with `training_sample_indices`, this says exactly which input values each eigenfunction value came
from, without copying them. `loadKoopmanModel` uses the interpolation as stored. Archives written
before the shared interpolation (with `domain__<name>` arrays) are refused and must be regenerated.

## Planned: SINDy model access

SINDy still reaches past its inputs. On a coordinate-change archive, it reloads the solve and
recomputes q(x). Its default `gamma_method: system_rhs` also evaluates a fixed model f(x, u).

The fix: SINDy must take the dynamical system only from what it is handed, through the origin
pointer of its time-series record, of its spectrum or of another input object. When a system is
found, it may use that system's equations for Γ. When none is found, it falls back to finite
differences or another numerical derivative of the record. It must never assume a model.

## Known limits

- **An interpolated scalar domain can be multivalued.** ψ(h) depends on the whole window, so a
  scalar channel such as echo does not determine it. On the RPE baseline, the best single-valued
  function of e leaves about 70% of ψ's variance unexplained.
- **Start-aligned windows are not causal.** Legacy H-DMD jobs align windows at their start, so
  ψ_k depends on later samples. The window map records `causal: false` for them. Use
  `window_alignment: "end"` for an online map.
- **Full-set window-map evaluation is expensive.** It costs one ΦᴴH product and a Cholesky solve
  of the rank × rank Gram matrix, which is factored once. Extraction materializes T rows only for
  the kept modes.
- **Interpolant gradients are finite differences** of the chosen interpolant.

References: Arbabi & Mezić, *SIAM J. Appl. Dyn. Syst.* 16 (2017); Brunton et al., *Nat. Commun.* 8
(2017); Schmid, *J. Fluid Mech.* 656 (2010).
