# Echo-only sampled Koopman eigenfunction results

**Date:** 2026-10-01
**Branch:** `echo-only`
**Status:** reproducible research result; no KROM promotion and no core-pipeline changes

## Question

Can causal echo history alone provide online coordinates from which useful
Koopman eigenfunctions can be evaluated, without the bubble equations, radius,
radial velocity, echo derivatives, or an online SVD?

The answer from the first controlled experiment is qualified:

- **yes on one invariant orbit:** an independently validated orbit-local
  complex eigenfunction exists and is accurately represented by a polynomial
  of the first two causal Hankel coordinates;
- **no as one fixed-frequency oscillatory eigenfunction over the tested
  amplitude family:** the inviscid nonlinear bubble is non-isochronous, so its
  phase speed changes with action;
- **possibly for a global invariant:** an approximate action-like
  (mu=1) function was extracted, but its held-out accuracy is not yet strong
  enough for promotion.

## Echo-only representation

Five equal-duration unforced inviscid echo records were used to fit one shared
Hankel basis. Every trajectory produced its own boundary-safe Hankel block;
no history crossed an initial-condition boundary. The blocks were concatenated
horizontally before the SVD.

The retained representation used:

| Quantity | Value |
|---|---:|
| Echo sample interval | (0.0001\tau) |
| Embedding decimation | 5 |
| Latent sample interval | (0.0005\tau) |
| Delay samples | 8,644 |
| Delay duration | (4.3215\tau) |
| Reference-period duration | approximately 3 periods |
| Initial retained dimension | 2 |

At evaluation time the only input is the current causal echo ring buffer
(d_k). The latent coordinate is

\[
\xi_i(k)=\frac{u_i^{\mathsf T}d_k}{\sigma_i}.
\]

The basis (u_i), singular-value scaling, echo mean, and echo scale are fixed
training artifacts. Radius and radial velocity were used only as simulation
audit channels and to choose independent phase points for one validation
experiment.

## Sampled polynomial eigenfunction fit

The latent coordinates were standardized using training-only statistics. A
total-degree-five dictionary in two variables was used:

\[
\Theta(z)=
\begin{bmatrix}
1&z_1&z_2&z_1^2&z_1z_2&z_2^2&\cdots
\end{bmatrix},
\]

containing 21 terms. Paired samples were separated by
(Delta\tau=0.01). Dense EDMD first fitted

\[
\Theta_+\approx\Theta_-K.
\]

For a right eigenvector (c),

\[
Kc=\mu c,
\qquad
\psi(z)=\Theta(z)c,
\]

so the sampled eigenfunction condition is

\[
\psi(z_{k+1})\approx\mu\psi(z_k).
\]

Dense fitting was intentional: representability was tested before introducing
sparsification error.

## One-orbit result

Fitting only the original unforced orbit produced

\[
\mu=0.99904758+0.04360206i,
\]

with implied period (1.440570\tau). The echo-measured period was
(1.440531\tau).

On the fitting record:

| Metric | Result |
|---|---:|
| (0.01\tau) sampled-step relative error | 0.2030% |
| (0.1\tau) semigroup error | 0.3797% |
| (1\tau) semigroup error | 0.3283% |
| (4\tau) semigroup error | 0.3579% |
| Eigenfunction-magnitude coefficient of variation | 0.0513% |
| Unwrapped-phase linearity (R^2) | effectively 1 |

### Independent same-orbit validation

Four fresh unforced solves were initialized at 20%, 40%, 60%, and 80% phase
points from the original orbit. These state values were used only to construct
simulation validation records. The retained basis and polynomial coefficients
were not refitted. All four records had the same echo-measured period at the
retained time resolution.

| Validation | Step error | (0.1\tau) | (1\tau) | (4\tau) | Fitted period |
|---|---:|---:|---:|---:|---:|
| Phase 20% | 0.2178% | 0.3931% | 0.3421% | 0.3835% | 1.440534 |
| Phase 40% | 0.2169% | 0.4055% | 0.3424% | 0.3858% | 1.440534 |
| Phase 60% | 0.2035% | 0.3837% | 0.3291% | 0.3672% | 1.440533 |
| Phase 80% | 0.2033% | 0.3819% | 0.3221% | 0.3414% | 1.440532 |

The median independent errors were 0.2102% for one sampled step and 0.3754%
at four (	au). This rules out simple memorization of the original time
record. The result is an orbit-local eigenfunction of the causal echo-history
coordinates.

## Sparsification

The successful dense multiplier was held fixed. Candidate supports were
selected from the dense coefficient magnitude, and coefficients on each
support were refitted by minimizing

\[
\left\|\Theta_+c-\mu\Theta_-c\right\|_2
\]

subject to a nonzero training-function normalization. Selection used training
errors only; the four phase-shifted trajectories remained independent tests.

The sparsest candidate satisfying the declared training tolerances retained
19 of 21 terms. Its independent median step error was 0.2052%, and its
four-(	au) error was 0.6707%. More aggressive supports of 18, 16, 14, or 11
terms substantially worsened the training semigroup error.

Thus the function can be sparsified slightly, but the evidence does not
support calling it strongly sparse. The 21-term dense expression remains the
preferred scientific reference; the 19-term version is the current compact
candidate.

## Why pooled oscillatory EDMD fails

The five unforced training periods ranged from approximately (1.020\tau) to
(1.532\tau). A pooled EDMD fit selected an implied period near
(1.433\tau), but the phase of the resulting polynomial followed each
trajectory's actual period rather than this fixed multiplier.

The mean training sampled-step residual was only about 2.53%. This apparently
small number was misleading because a (0.01\tau) sampled map is close to the
identity. Held-out four-(	au) errors ranged from approximately 0.79 to 1.78.

This failure persisted after the three-period top-two chart eliminated the
proper cross-orbit intersections detected in the finite injectivity audit.
Therefore two different obstructions have been separated:

1. insufficient history can fold distinct physical states onto the same
   latent point;
2. even an empirically unfolded chart cannot give a continuum of
   non-isochronous invariant orbits one common nonzero Koopman eigenvalue.

For action-angle coordinates ((I,\theta)), the dynamics have the form

\[
I_{k+1}=I_k,
\qquad
\theta_{k+1}=\theta_k+\omega(I_k)\Delta\tau.
\]

Then (e^{i\theta}) evolves with the state-dependent multiplier
(e^{i\omega(I)\Delta\tau}), not one constant (mu), unless the family is
isochronous or the fit is restricted to one action level.

## Global invariant and latent-rank audit

A separate (mu=1) fit maximized between-trajectory polynomial variation
while minimizing within-trajectory variation. It uses run identity and latent
coordinates, not hidden physical state. This produces an approximate
action-like invariant.

The controlled rank comparison retained the same training records, validation
records, three-period parent SVD, normalization, polynomial degree, and scoring
rules:

| Rank | Dictionary terms | Median held-out invariant std. | Maximum std. | Held-out (omega(I)) error |
|---:|---:|---:|---:|---:|
| 2 | 21 | 0.0610 | 0.1696 | 16.3% |
| 3 | 56 | **0.0548** | **0.1630** | **11.9%** |
| 4 | 126 | 0.0756 | 0.2171 | 15.5% |

Rank three improves the preliminary invariant. Rank four is worse because the
fixed fifth-degree dictionary expands to 126 terms and becomes more weakly
constrained by five training orbits. This is an observed conditioning and
overfitting cost, not evidence that the fourth coordinate is intrinsically
harmful.

The invariant and the quadratic frequency diagnostic are not yet accurate
enough for promotion. In particular, (omega(I)) is a parameter-varying phase
law, not a fixed-eigenvalue Koopman KROM.

## Viscous RPE comparison

The conservative result should not be interpreted as a special failure of
echo-history coordinates. A conservative nonlinear RPE oscillator has an
amplitude-dependent frequency, just as a conservative pendulum does. Its
nonzero-frequency Koopman eigenfunctions have orbit/frequency-selective
support; one fixed-frequency eigenfunction is not expected to extend smoothly
over the full family of energy levels. The same obstruction appears in
((r,\dot r)) coordinates.

To separate that spectral obstruction from the echo representation, the
experiment was repeated with the standard air-in-water viscosity
((\mu=0.001\,\mathrm{Pa\,s})). Five long unforced training trajectories and
five long unforced held-out trajectories were integrated to (80\tau) at
(0.0001\tau) output resolution. The echo-history basis was fitted from the
training records only. The held-out records were projected through this fixed
basis without a new SVD.

The equilibrium linearization has the continuous eigenvalue

\[
\lambda_{\mathrm{eq}}
=-0.0603684+6.2828953i.
\]

A single degree-five polynomial fit over the original (12\tau) records was
not successful. Those records cover less than one linear damping e-fold and
remain dominated by amplitude-dependent nonlinear frequency. Importantly, a
matched direct-state EDMD audit failed over the same interval, establishing
that the early failure was not caused specifically by the latent variables.

The long-record experiment therefore retained one fixed Hankel basis and
refitted the eigenfunction on suffixes beginning at progressively later times.
The data-discovered results were:

| Fit begins at | Coordinates | Discovered (lambda) | Median held-out four-(tau) error | Median held-out anchored error |
|---:|---|---:|---:|---:|
| (5\tau) | echo latent, rank 2 | (-0.02180+5.96971i) | 0.6700 | 1.1947 |
| (30\tau) | echo latent, rank 2 | (-0.05073+6.24850i) | 0.09082 | 0.2742 |
| (40\tau) | echo latent, rank 2 | (-0.05734+6.28184i) | 0.007733 | 0.02069 |
| (50\tau) | echo latent, rank 2 | (-0.060245+6.282949i) | 0.000397 | 0.000778 |
| (60\tau) | echo latent, rank 2 | (-0.060366+6.282896i) | (1.68\times10^{-5}) | (3.00\times10^{-5}) |
| (50\tau) | echo latent, rank 4 | (-0.060347+6.282875i) | 0.000145 | 0.000332 |
| (50\tau) | direct state audit | (-0.060327+6.282654i) | 0.000656 | 0.00156 |

Thus the viscous trajectories do converge to one common oscillatory Koopman
mode, and the first two causal echo-history coordinates recover it as
accurately as the direct-state audit in the asymptotic regime. The held-out
time histories from the rank-two (50\tau) fit are visually indistinguishable
from the predicted exponentially decaying sinusoids over the remaining
(30\tau).

This does **not** yet demonstrate an accurate global eigenfunction over the
entire basin. It demonstrates a common asymptotic eigenfunction and a clear
transition into its domain of accurate polynomial representation. The early
large-amplitude transient remains difficult for the same physical reason in
both coordinate systems: its instantaneous frequency is still strongly
amplitude dependent, and a low-degree global polynomial has to approximate a
highly nonlinear isochron map.

The fixed full-transient Hankel basis also makes an important distinction
between energy and dynamical sufficiency. Across the complete training Hankel
matrix, the first two modes retain only 6.14% of Frobenius energy and the first
four retain 9.03%, because short, violent nonlinear collapse features dominate
the energy. On suffixes beginning at (50\tau), the same fixed modes retain
about 20.3% and 59.1% of training-window energy, respectively (20.5% and 58.8%
on held-out data). Nevertheless, rank two is already sufficient for the
asymptotic eigenfunction. Singular-value energy is therefore not a valid proxy
for Koopman-coordinate sufficiency.

The conservative and viscous conclusions are consistent:

1. conservative RPE supports accurate orbit-local oscillatory eigenfunctions,
   but not one regular fixed-frequency eigenfunction across different actions;
2. viscous RPE retains conservative-like frequency dispersion during its
   large-amplitude transient;
3. dissipation eventually collapses all tested trajectories into a common
   stable-focus spectrum;
4. echo history alone recovers that common mode online from a fixed training
   basis, with no radius, radial velocity, echo derivative, model equation, or
   held-out decomposition required.

## Formal eigenfunction and native-evaluator reproduction

The rank-two viscous fit beginning at (50\tau) was recomputed from the retained
training-only Hankel basis and the archived unforced trajectories. It recovered
the previously reported multiplier and continuous eigenvalue exactly:

\[
\mu=0.9974257989803501+0.0627503477792618i,
\qquad
\lambda=-0.0602447488246558+6.28294918416634i.
\]

The 21-term degree-five polynomial was then wrapped in the canonical
`KoopmanEigenfunction` contract and evaluated with the canonical
`FreeTrajectoryDataset` error framework. On the five training suffixes, the
pooled one-step and four-(\tau) relative (L^2) errors are
(1.27\times10^{-6}) and (4.46\times10^{-4}). On the five held-out suffixes,
they are (4.30\times10^{-6}) and (1.72\times10^{-3}). The last number is a
pooled norm over all held-out samples. The earlier value
(3.97\times10^{-4}) is the median of five independently normalized
per-trajectory errors; both are retained because they answer different
aggregation questions.

The canonical generator residual, computed from numerical derivatives along
the latent trajectories, is (4.04\times10^{-4}) on training data and
(1.37\times10^{-3}) on held-out data. Only 0.908% of held-out late-time
samples lie outside the training bounding box, with a maximum normalized
excursion of 8.90% of a box width.

The polynomial coefficients were also passed to the native
`KoopmanModel::evaluateCoordinates` implementation without changing core
code. Four probes spanning the fitted coordinate domain agreed with the
canonical Python values to a maximum absolute difference of
(8.88\times10^{-16}) and maximum relative difference of
(2.43\times10^{-16}). Thus the standard analytic-polynomial representation
is already numerically compatible with a polynomial in ((\xi_1,\xi_2)).

This is not yet a complete KROM reproduction. The missing contract is the
stateful coordinate map

\[
(e_{k-d+1},\ldots,e_k)\longmapsto(\xi_1(k),\xi_2(k)),
\]

including its ring-buffer initialization policy. The current KROM
`coordinate_observable` is an instantaneous map of the current dynamical
state, and the standard Viewer phase-map path likewise assumes instantaneous
state coordinates. Direct coordinate evaluation is exact once (\xi) is
supplied, but the core cannot yet obtain (\xi) from echo history by itself.
For that reason the reproduced object remains an isolated formal
eigenfunction artifact rather than a promoted KROM.

The accuracy maps also delimit what was produced. The fit uses only
(50\tau\leq\tau\leq80\tau). Four-(\tau) error is order one in the earlier
large-amplitude transient and collapses only as the viscous trajectories enter
the late-time annulus. This is a common asymptotic eigenfunction, not a global
basin-wide eigenfunction.

## Accepted and rejected conclusions

Accepted as research evidence:

- the pooled three-period echo basis supports an accurately evaluable
  top-two orbit-local eigenfunction;
- that function generalizes across independent phase shifts on its invariant
  orbit;
- online evaluation requires only the echo ring buffer;
- long-horizon testing exposes failures that one-step residuals conceal.

Rejected:

- the pooled oscillatory EDMD pair as a global fixed-eigenvalue KROM;
- the proposition that higher latent rank automatically improves an
  eigenfunction fit;
- the proposition that high singular-value energy or a small one-step error
  establishes Koopman validity;
- promotion of the provisional action invariant or (omega(I)) law.

## Viewer evidence

- `viewer_outputs/rpe-inviscid-echo-sampled-koopman-summary.png`
- `viewer_outputs/rpe-inviscid-echo-sampled-koopman-evolution.png`
- `viewer_outputs/rpe-inviscid-echo-sampled-koopman-invariant.png`
- `viewer_outputs/rpe-inviscid-echo-latent-rank-koopman.png`
- `viewer_outputs/rpe-inviscid-echo-local-koopman-validation.png`
- `viewer_outputs/rpe-inviscid-echo-local-koopman-phase-replay.png`
- `viewer_outputs/rpe-viscous-echo-asymptotic-koopman-summary.png`
- `viewer_outputs/rpe-viscous-echo-asymptotic-koopman-evolution.png`
- `viewer_outputs/rpe-viscous-echo-asymptotic-eigenfunction-field.png`
- `viewer_outputs/rpe-viscous-echo-asymptotic-eigenfunction-trajectories.png`
- `viewer_outputs/rpe-viscous-echo-asymptotic-eigenfunction-accuracy.png`

Every listed PNG was rendered at 1900 by 1900 using Viewer PlotCanvas.

## Reproduction

The experiments are intentionally separate from the core pipeline:

```powershell
python tools/compare_rpe_echo_pooled_training_injectivity.py
python tools/analyze_rpe_echo_sampled_koopman.py
python tools/analyze_rpe_echo_latent_rank_koopman.py
python tools/validate_rpe_echo_local_koopman.py
python tools/analyze_rpe_viscous_echo_koopman.py
python tools/analyze_rpe_viscous_echo_asymptotic_koopman.py
python tools/plot_rpe_viscous_echo_eigenfunction_domain.py
```

Machine-readable results are retained in the corresponding
`viewer_outputs/*.metrics.json` and `viewer_outputs/*.npz` artifacts. The
phase-shift validation solve specification is
`viewer_outputs/experiments/rpe-inviscid-echo-local-eigenfunction-phase-validation.json`.

## Next control-facing step

If this eigenfunction is used for KLQR, it must be treated as local to its
target action level. Controlled data may fit a statistical (B_\psi) map in
the usual way, but the control study must monitor departure from the validated
orbit. A controller that moves substantially across action levels cannot rely
on the same fixed autonomous multiplier without a local-bank,
gain-scheduled/parameter-varying, or nonlinear phase-law extension.
