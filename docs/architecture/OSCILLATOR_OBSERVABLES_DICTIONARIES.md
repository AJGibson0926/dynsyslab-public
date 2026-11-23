# Oscillator observables and reusable spectral dictionaries

Implemented 2026-11-21. These are native pipeline features, not experiment-only
scripts. Scalar/vector observables belong to the system catalog; finite trial
spaces belong to the shared spectral dictionary. A chosen dictionary is not an
analytically known Koopman eigenbasis: the estimator computes approximate
eigenfunctions as linear combinations of its features.

## System catalog

Write the state as $(q,v)$. Duffing and Van der Pol already use normalized state
coordinates. For the dimensional pendulum the native state is
$(\theta,\omega=d\theta/dt)$; the reference formulas below use
$q=\theta$, $v=\sqrt{\ell/g}\,\omega=d\theta/d\tau$ with
$\tau=t\sqrt{g/\ell}$. This scaling is recorded in each real component's analytic
expression. The new pendulum entries require positive length and gravity.

| System | Catalog entries | Formula / meaning |
|---|---|---|
| Pendulum | Existing `angle`, `velocity`, `state` | Native angle, angular velocity, their pair |
| Pendulum | Existing `torque`, `torque_dimensionless` | $-mg\ell\sin\theta-b\omega$ and division by $mg\ell$; excludes applied input torque |
| Pendulum | New `phase_sech` | $\exp(iq)\operatorname{sech}v$ |
| Pendulum | New `phase5_sech` | $\exp(5iq)\operatorname{sech}v$ |
| Pendulum | New `phase_gaussian_velocity` | $\exp(iq)v\exp(-v^2/2)$ |
| Duffing | Existing `position`, `velocity`, `state` | $q,v,(q,v)$ |
| Duffing | Existing `force`, `force_dimensionless` | $\alpha q-\beta q^3-\delta v$, without applied drive |
| Duffing | Existing energy observables | `kinetic_energy`, `potential_energy`, `total_energy` |
| Duffing | Existing `random_observation`; new canonical `phase_sech` | Both $\exp(iq)\operatorname{sech}v$; deterministic, despite the old name |
| Duffing | New `gaussian_sum` | $\exp(-(q+v)^2)$ |
| Van der Pol | Existing `position`, `velocity`, `state` | $q,v,(q,v)$; use any generic dictionary below on these coordinates |

Every new complex name has three additional catalog forms: `_components`
returns `[Re, Im]`, `_real` returns Re, and `_imag` returns Im. They have analytic
Jacobians and serialized real expressions, including analytic partial and flow
derivatives. The complex form has an analytic complex Jacobian; use its real
components for higher expression derivatives or real-valued fitting. The new
`sech` expression primitive evaluates stably without overflowing `cosh`.
`gaussian_sum` also has analytic partial/flow derivatives.
The global jet analyzer also supports `sech` through a polynomial lift with
its `tanh` companion, positivity and $\operatorname{sech}^2+\tanh^2=1$.
Its transcendental source graph remains explicitly marked as a relaxation.

For example, with $a=\sqrt{\ell/g}$ and $g_k=\exp(ik\theta)\operatorname{sech}(a\omega)$,
$\partial_\theta g_k=ikg_k$ and
$\partial_\omega g_k=-a\tanh(a\omega)g_k$.
The catalog's `angle`/`position` to `velocity` time-derivative links remain
available. Analytic flow derivatives of the new real components use the supplied
system vector field, with the existing explicit zero-input / supplied-input
contract; their derivatives are not assumed to be other catalog scalars.

The downloaded Van der Pol examples supplied dictionary constructions, not an
additional system-specific scalar measurement. Laguerre and radial dictionaries
are therefore reusable options, not artificial Van der Pol-only catalog names.

## Dictionary families

For all families, $s_j=(x_j-c_j)/a_j$. `center` and `scale` explicitly specify
$c,a$; otherwise both default to identity. `normalize:true` learns mean/std from
training X only and cannot be combined with explicit center/scale. Frozen KROMs
store the resulting arrays, indices, centers, bandwidth and kernel parameters;
held-out evaluation does not refit them.

| `kind` | Feature | Settings |
|---|---|---|
| `linear` | $1,s_1,\ldots,s_d$ | Existing; `constant` |
| `polynomial` | $\prod_j s_j^{\alpha_j}$ | Existing; now also tensor-product indices |
| `chebyshev` | $\prod_j T_{\alpha_j}(s_j)$ | `degree`, `index_set` |
| `laguerre` | $\prod_j L_{\alpha_j}(s_j)$, ordinary Laguerre | `degree`, `index_set` |
| `hermite` | $\prod_j h_{\alpha_j}(s_j)$, normalized Hermite **functions** | `degree`, `index_set` |
| `fourier_hermite` | $F_k(s_1)h_n(s_2)$ | Dimension 2; `angular_order`, `hermite_order` |
| `gaussian` | $A\exp(-r^2/(2b^2))$ | Existing; `bandwidth`, `amplitude` |
| `laplacian` | $A\exp(-r/b)$ | `bandwidth`, `amplitude` |
| `lorentzian` | $A/(1+(r/b)^2)$ | `bandwidth`, `amplitude` |
| `matern` | $A\,2^{1-\nu}z^\nu K_\nu(z)/\Gamma(\nu)$, $z=r/b$ | `bandwidth`, `smoothness` $\nu\in[0.5,8]$, `amplitude`; value $A$ at $r=0$ |
| `partition` | Indicator of the nearest center's Voronoi cell | Existing; first center wins distance ties |

Here $r$ is Euclidean distance to a center in the normalized coordinates;
`nodes` are stored in original input coordinates. Amplitude defaults to 1 and
bandwidth to 1. Our Matérn bandwidth absorbs the common $\sqrt{2\nu}$ convention:
with the usual length parameter $\ell_K$, set $b=\ell_K/\sqrt{2\nu}$.

Product indices use `index_set:"total_degree"` (default, $|\alpha|\le p$) or
`"tensor_product"` ($0\le\alpha_j\le p$). Auto-generated degree is at most 32;
explicit frozen indices support orders through 512. `maximum_features` defaults
to 4096 and bounds generated dictionaries. `constant:false` omits the all-zero
multi-index: for Hermite functions this omits the Gaussian ground state, **not**
a constant function. It does not change radial or partition features.

Hermite normalization is
$h_n(s)=H_n(s)e^{-s^2/2}/\sqrt{2^n n!\sqrt\pi}$, using physicists' Hermite
polynomials. The real Fourier basis is
$F_0=1/\sqrt{2\pi}$,
$F_{2k-1}=\cos(k\theta)/\sqrt\pi$,
$F_{2k}=\sin(k\theta)/\sqrt\pi$.
It spans the same complexified space as normalized $e^{ik\theta}$ for
$-K\le k\le K$, and yields $(2K+1)(N+1)$ features before optional zero-index
omission. Frozen `exponents` store `(Fourier code, Hermite order)` in that order.
The angular scale must be 1 (radians); automatic normalization is rejected.
For raw pendulum states use `scale:[1,sqrt(g/length)]` to normalize only velocity.
These functions are orthonormal under their defining continuous measures, not
necessarily under an arbitrary trajectory's empirical measure. Chebyshev values
outside [-1,1] and Laguerre values on negative arguments are allowed; conditioning
and the appropriateness of the weighting measure remain scientific choices.

Radial/partition dictionaries take `centers` and `center_selection:"farthest"`
(default) or `"kmeans"`. K-means uses deterministic farthest-point initialization,
then Lloyd iterations (default 100, maximum 1000), with empty centers retained.
It uses training X only. Upstream scripts sometimes pool X and Y and use randomized
k-means; supply explicit frozen `nodes` to reproduce a particular upstream set.
Repeated input points can result in fewer centers than requested.

## CLI, Python and pipeline use

The same native endpoint builds a dictionary without fitting a Koopman operator:

```python
from dynsyslab_client.workflow import NativeWorkflow
api = NativeWorkflow()
result = api.koopman_spectral({
    "action": "dictionary",
    "data": {"dimension": 2, "x": [0, 0, 1, -1, 2, .5]},
    "settings": {"dictionary": {"kind": "chebyshev", "degree": 9,
        "index_set": "tensor_product", "scale": [2.5, 2.5]}}
})
# 3 samples x 100 features, flattened in row-major order.
frozen = result["dictionary"]
held_out = api.koopman_spectral({"action": "dictionary", "dictionary": frozen,
                                "points": [.4, -.2]})
observable = api.observable({"system": {"name": "nonlinear_pendulum"},
    "observable": "phase_sech_components", "state": [.7, -.3], "flow_order": 1})
```

Pass the same dictionary settings to `action:"fit"`, or to the existing
Post-Processor `resdmd` and other spectral stages. For EDMD select
`koopman_spectral` with `settings.method:"edmd"`; the older `edmd` stage
retains its separate existing contract. The Spectral Koopman dialog
lists all families. Its advanced settings accept index sets, Fourier/Hermite
orders, center selection, explicit normalization and Matérn parameters.
A frozen dictionary contains `nodes` and/or `exponents`; when these fields are
present the API uses the serialized basis instead of generating another one.

Runnable, self-contained examples:

```powershell
dynsyslab-cli koopman-spectral examples/native/observable-dictionaries/chebyshev.json
dynsyslab-cli koopman-spectral examples/native/observable-dictionaries/fourier-hermite.json
dynsyslab-cli koopman-spectral examples/native/observable-dictionaries/laguerre.json
dynsyslab-cli koopman-spectral examples/native/observable-dictionaries/matern.json
```

The `dictionary` result includes `dictionary`, `samples`, `features`, and flat
row-major `values`. Complex observable-query values and Jacobian entries use
`[real,imag]` pairs and `value_type:"complex"`; real responses retain flat numeric
arrays with `value_type:"real"`. Jacobian shape is outputs by state dimension.
Use `_real`/`_imag` as separate named channels in the
[observable/Hankel pipeline](OBSERVABLE_HANKEL.md), for example:

```json
{"observables":[{"observable":"phase_sech_real","taps":[0,2]},
                {"observable":"phase_sech_imag"}],
 "window":{"delays":16,"columns":100}}
```

The delay dictionary $g,g\circ T,\ldots$ in the reference pendulum examples is
already handled by that existing Hankel module. Do not invent another windowing
implementation for these observables.

## Sources and scope

For the pendulum, the exact derivative-closure check
rules out a nontrivial finite bank closed under both drift and physical
derivatives. The new dictionaries are useful approximation spaces; their finite
size or orthogonality does not bypass that obstruction.

The pinned source copies are in ignored `build/reference_sources/colbrook`,
with revisions in the [manifest](../../tools/experiments/miscellaneous_computes/references/colbrook/repositories.json).
The formulas above were independently implemented using standard recurrences,
not by embedding MATLAB or copying unlicensed source.

| Reference example | Observable/dictionary now available |
|---|---|
| [Rigged DMD pendulum](https://github.com/MColbrook/Rigged-Dynamic-Mode-Decomposition/blob/013f5a429868/PENDULUM_example.m) | Phase/sech scalar and its existing Hankel delays |
| [spectralcomputations pendulum](https://github.com/MColbrook/spectralcomputations/blob/8e7d5e469168/chapter_11/pendulum_example.m) | Fifth phase harmonic/sech scalar |
| [mpEDMD pendulum](https://github.com/MColbrook/Measure-preserving-Extended-Dynamic-Mode-Decomposition/blob/e473dee7b32d/mpEDMD_pendulum.m) | Phase times Gaussian-weighted velocity scalar |
| [ResDMD pendulum](https://github.com/MColbrook/Residual-Dynamic-Mode-Decomposition/blob/0f7f1bf67832/Examples_gallery_1/pendulum_example.m) | Real-equivalent Fourier-Hermite trial space |
| [ResDMD Duffing](https://github.com/MColbrook/Residual-Dynamic-Mode-Decomposition/blob/0f7f1bf67832/Duffing_example.m) | Laplacian RBF sections and k-means centers |
| [spectralcomputations Duffing](https://github.com/MColbrook/spectralcomputations/blob/8e7d5e469168/chapter_11/Duffing_oscillator_convergence.m) | Gaussian-sum scalar; exponential RBF with bandwidth $d/2$ |
| [ResDMD kernels](https://github.com/MColbrook/Residual-Dynamic-Mode-Decomposition/blob/0f7f1bf67832/main_routines/kernel_ResDMD.m) | Gaussian/Laplacian/Lorentzian sections; polynomial trial spaces already existed |
| [TrustKoopman Duffing](https://github.com/GustavConradie1/TrustKoopman/blob/ee5769f3b44d/examples/duffing_principal_obs.m) | $T_j(q/2.5)T_k(v/2.5)$, $j,k=0,\ldots,9$, 100 features |
| [Koopman-PageRank source](https://github.com/HyukpyoHong/Koopman-PageRank/blob/6408ef364d60/src/koopman_source.py) | Laguerre products through total degree 12, 91 features in 2D |
| [SpecRKHS Duffing](https://github.com/GustavConradie1/SpecRKHS/blob/8b888786c3e7/examples/duffingoscillator.m) | $z^2K_2(z)$ with $z=6r$: Matérn smoothness 2, bandwidth 1/6, amplitude 2 |
| [MultDMD basis](https://github.com/NBoulle/MultDMD/blob/85277b046a9a/src/construct_basis.m) | Nearest-center indicators, now with optional deterministic k-means |

Kernel-section EDMD/ResDMD with these features is supported. This does **not**
port the full upstream kernel-whitening algorithm or every theorem/workflow.
The separate `specrkhs` adjoint endpoint remains Gaussian-only; choosing a Matérn
dictionary for EDMD is not the Matérn SpecRKHS adjoint method. Upstream Gaussian
$e^{-r^2/d^2}$ corresponds to our bandwidth $d/\sqrt2$. Polynomial kernels have
the same finite polynomial span under appropriate nondegenerate conditions,
but the raw monomial coordinates do not reproduce their kernel metric.

## Ownership and validation

- Foundation: [OscillatorObservables.hpp](../../native/include/dynsyslab/core/OscillatorObservables.hpp),
  [BasisFunctions.hpp](../../native/include/dynsyslab/core/BasisFunctions.hpp), analytic expressions and system registrations.
- Koopman: [SpectralDictionary.cpp](../../native/src/SpectralDictionary.cpp), frozen evaluation and serialization.
- Spectral: [KoopmanDictionaryFit.cpp](../../native/src/KoopmanDictionaryFit.cpp), training-only construction and center fitting.
- Workflow: [KoopmanSpectralQuery.cpp](../../workflow/native/src/KoopmanSpectralQuery.cpp),
  [ObservableQuery.cpp](../../workflow/native/src/ObservableQuery.cpp), existing CLI/C ABI and Python marshalling.
- Workbench: [SpectralKoopmanDialog.cpp](../../compute_workbench/native/src/SpectralKoopmanDialog.cpp), settings only.

[Native tests](../../native/tests/observable_dictionary_tests.cpp) check formulas,
analytic gradients, flow derivatives, mixed partials, Hermite orthogonality,
periodicity, center selection, bounds and serialization.
[Python integration tests](../../tests/test_native_observable_dictionaries.py)
compare against independent SciPy formulas, fit/archive/evaluate every new
family, and exercise the CLI and Post-Processor. This validates implementation;
it does not assert spectral convergence for a particular oscillator dataset.

Validation on 2026-11-21: rebuilt the native CLI, workflow C API and workbench;
all eight focused native test targets and ten Python integration tests passed,
including the existing pinned MATLAB mpEDMD/Rigged parity fixture. All four
standalone dictionary examples ran through the compiled CLI. The strict
architecture audit passed with zero boundary violations and broken links
(pre-existing unindexed documentation remains). The built applications are in
`build/native-suite-fast/bin`; the installed application was not replaced.
