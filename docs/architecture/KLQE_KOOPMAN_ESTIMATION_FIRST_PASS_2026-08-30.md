# KLQE / Koopman observer architecture: MATLAB trace and first Python design pass

Date: 2026-08-30  
Status: investigation complete; the candidate architecture now has a tested first-pass implementation

Implementation record: [KLQE_FIRST_PASS_IMPLEMENTATION_2026-08-30.md](KLQE_FIRST_PASS_IMPLEMENTATION_2026-08-30.md)

The source trace below remains the design rationale and MATLAB comparison. The
linked implementation record identifies exactly which proposals have now been
implemented, their tests, and the boundaries intentionally deferred.

## Purpose

This report traces the existing MATLAB code for Koopman linear-quadratic
estimation (KLQE), Koopman-Luenberger/Kalman filtering, observable
reconstruction through Koopman modal amplitudes, and classical observability
analysis. It then compares those mechanisms with the current Python DynSysLab
abstractions and proposes a controlled first implementation and benchmark
sequence.

The goal is not to copy the MATLAB suite. The MATLAB suite is being used as a
read-only record of the successful experimental pathway and of the intended
mathematics. Several parts of that implementation are incomplete,
coordinate-inconsistent, or plainly experimental. Those defects are called out
explicitly below so they are not canonized in Python.

## Executive findings

1. The MATLAB suite does contain the requested RPE acoustic-echo estimator
   experiment. The primary script is `KLQE.m`; the reusable estimator is
   `scripts/comp_modules/ESTM.m`; joint plant/observer integration is in
   `scripts/comp_modules/compute_node.m`.

2. The intended lifted observer is

   ```text
   zhat_dot = A zhat + Bpsi(xhat) u + L (y - Cy zhat)
   xhat     = Cx zhat
   yhat     = Cy zhat
   ```

   For the free RPE experiment, `Bpsi u = 0`. The measured scalar `y` is the
   far-field acoustic echo. The lifted state `zhat` can have hundreds of real
   coordinates, while the two physical states are reconstructed with the first
   two observable rows of `C`.

3. MATLAB computes the observer gain through dual continuous-time LQR,

   ```text
   L = lqr(A', Cy', Qe, Re)'
   ```

   which is equivalent to a steady-state continuous LQE/Kalman-Bucy gain when
   `Qe` is the lifted process-noise covariance (or `G W G'`) and `Re` is the
   measurement-noise covariance. The active examples do not actually inject
   noise: both disturbance and measurement-noise signals are multiplied by
   zero. The current evidence therefore demonstrates a deterministic
   Luenberger observer designed by an LQE CARE, not yet a validated noisy
   Koopman-Kalman filter.

4. There are two meaningful MATLAB routes to a Koopman observable matrix `C`:

   - H-DMD supplies modal amplitudes `nu` for every observable included in its
     horizontally stacked Hankel construction.
   - `KROM.compute_nu` fits arbitrary observable rows from sampled data through
     the least-squares projection `Y = C Psi`.

   The method `observable.project` is not the working projection route: it
   references an undefined `x` and does not integrate, sum, or normalize the
   proposed inner product.

5. The MATLAB observability routine constructs the classical matrix

   ```text
   O = [Cy; Cy A; ...; Cy A^(n-1)]
   ```

   using `obsv(A, Cy)` and computes `O = U S V'`. The right singular-vector
   columns `V(:,i)` are lifted-state directions ordered by observability over
   that algebraic horizon. The experiment plots `V(i,:)` instead, which is a
   row and does not represent the corresponding singular direction.

6. The MATLAB reusable KLQE path has a serious time-domain defect. `ESTM.train`
   requests a continuous model, but `plant.model` overwrites that request with
   the KROM's stored mode. An H-DMD KROM is stored as discrete, so the class path
   can feed discrete eigenvalues into continuous `lqr` and continuous ODE
   integration. The manual block in `KLQE.m` explicitly requests continuous
   eigenvalues and is closer to the intended formulation.

7. The current Python suite already has strong pieces of the desired design:

   - formal `Observable` objects and system-owned observable catalogs;
   - `ObservableReconstruction` for a declared `h(x) ~= C Psi(x)`;
   - native H-DMD modal amplitudes archived as a reconstruction object;
   - validated least-squares observable reconstruction and H-DMD-to-SINDy
     normalization transfer;
   - a continuous-time `KLQEEstimator` using a CARE gain;
   - internally consistent nonredundant-complex and real/imaginary-split KROM
     conventions.

   What is missing is the complete product path: named multiple observable
   reconstructions on one KROM, classical/finite-horizon observability
   analysis, a numerically sound high-dimensional observer propagator, an
   estimator-only solver/offline workflow, archived estimate/innovation
   histories, CLI/API tooling, Viewer displays, and benchmark parity.

8. One current Python H-DMD detail needs correction before direct H-DMD
   eigenfunctions are trusted for estimation. In the window-snapshot layout,
   each projected-mode row corresponds to the state at the end of its delay
   window, but interpolation is currently attached to the state at the start of
   that window. MATLAB uses the end-of-window indices. This is a delay of
   `matrix_size - 1` samples in the Python interpolation-state association. The
   terminal modal-amplitude `C` remains internally associated with the DMD mode
   rows, but evaluating those modes as functions of physical state is shifted.
   This must be isolated by a test before an H-DMD KROM is used directly in a
   filter.

## MATLAB source map

| Concern | MATLAB location | Relevant behavior |
|---|---|---|
| RPE KLQE experiment | `KLQE.m:12-37`, `65-176`, `224-342` | Builds free RPE data, H-DMD KROM, observable-selected sub-KROM, echo observer, and plots |
| Reusable estimator | `scripts/comp_modules/ESTM.m:38-112`, `138-160` | Builds `Cy`, dual-LQR gain, innovation injection, and `Cx` reconstruction |
| Joint integration | `scripts/comp_modules/compute_node.m:73-105`, `208-259` | Packs plant and lifted observer state into one ODE and records estimates |
| H-DMD and native `nu` | `scripts/comp_modules/DMD.m:265-433` | Horizontal Hankel blocks, projected modes, modal amplitudes, energy sorting |
| Arbitrary `C` fit | `scripts/data_structs/Koopman_tree/KROM.m:655-756` | Fits `Y = C Psi` and attaches each `nu` to each eigenfunction |
| KROM `A`, `C`, and observability | `KROM.m:387-448`, `608-653`, `758-771` | Real/complex matrix assembly and `obsv`/SVD |
| Mode selection | `scripts/data_structs/Koopman_tree/bank.m:213-356` | Observable-energy selection and an unfinished `maxO` experiment |
| Observable class | `scripts/data_structs/Koopman_tree/observable.m` | Observable functions and an incomplete `project` method |
| RPE echo definition | `scripts/dynsyslib/Bubble/Bubble.m:180-211` | Defines radius, radial velocity, balance terms, and far-field echo |
| Standalone prototype | `KoopmanKalman/KKF.m` | Older CSV-driven LQE/Luenberger test |
| Standalone failsafe | `KoopmanKalman/working example - failsafe/KKF.m` | Dimensionally coherent 100-state prototype input set |

## Exact RPE experiment traced from `KLQE.m`

### Plant and trajectory

The active `sys` is the first `Bubble` system:

```matlab
rpe = dynsys(dynamics=Bubble,trajectories=[1]);
rpe.dyn.prms.CL = 0;
rpe.dyn.prms.CG = 0;
tres = 0.001;
T = 10;
bubbles.compute("soft");   % ode15s
sys = rpe;
```

The initial condition is not explicitly assigned in this block; it comes from
the suite's then-current default initialization machinery. This must be exported
rather than guessed for parity. The default Bubble preset is nondimensional
Rayleigh-Plesset with the now-familiar physical values
`R0=10e-6`, `K=4/3`, `p0=100e3`, `pV=2.33e3`, `rhoL=998`,
`sigma=0.0725`, and `muL=0.001`; the experiment explicitly zeros `CL` and `CG`.

The fitted observables are:

```matlab
r     = sys.dyn.obsv(1);
rdot  = sys.dyn.obsv(2);
echo  = sys.dyn.obsv(end);
```

For Rayleigh-Plesset, the echo is defined as

```text
Z  = T2 - T1 + Eu*pVstar - (Eu + LP + VC)
e  = (p0*R0/Eu) * r * (2*rdot^2 + Z/r)
```

where `r` and `rdot` are the nondimensional bubble state and the intermediate
terms use the suite's nondimensional groups. The state is nondimensional but
the echo includes a dimensional prefactor. Its exact measurement units and
normalization therefore need to be stored explicitly in a parity artifact; they
should not be inferred later from a plot label.

### H-DMD construction

The call is:

```matlab
Koop.process(
    data=sys,
    eigenvalue_method='HDMD',
    eigenfunction_method='HDMD',
    observable=[r,rdot,e],
    rows=10000,
    cols=2000,
    stack_method='horizontal',
    sortmode='Energy');
```

For each scalar observable, MATLAB builds a Hankel block whose rows are sliding
windows and whose 2,000 columns are delays. The blocks are concatenated
horizontally, so the requested full matrix is approximately

```text
10,000 rows x (3 * 2,000) columns.
```

In the ten-time-unit script, the source trajectory is too short to supply all
10,000 rows after reserving 2,000 delay samples, so the actual row count is
about 8,001. The later `playground.m` version uses a 50-time-unit trajectory and
can fill the requested row count.

The DMD kernel is:

```text
X = U S V'
Atilde = U' Y V S^-1
Atilde W = W Lambda_d
Phi = U W
```

`Phi` is called `ProjectedModes`. Because the rows are time-delay windows, these
mode values are treated as samples of Koopman eigenfunctions along the training
trajectory.

### Why H-DMD supplies `C`

For horizontal stacking, MATLAB selects the final delay column of each
observable's Hankel block and computes

```text
nu_h = pinv(Phi) h_window_end.
```

With three scalar observables, `nu` has shape

```text
3 observables x full complex DMD mode count.
```

Each row is the modal coefficient vector for that observable in the H-DMD
eigenfunction normalization. After sorting, the code attaches column `n` of
`nu` to eigenfunction `psi_n`. Thus the intended reconstruction is

```text
h_j(x) ~= sum_n nu[j,n] psi_n(x).
```

This explains the user's recollection accurately: H-DMD gives the `C` rows for
the observables used to form the horizontal Hankel matrix. It does not
automatically give a `C` row for an arbitrary observable that was absent from
that matrix.

The coefficient is tied to both the observable and the eigenfunction
normalization. Rescaling `psi_n -> alpha_n psi_n` requires
`nu_n -> nu_n/alpha_n`. A `C` matrix cannot be transferred from H-DMD
eigenfunctions to separately fitted SINDy eigenfunctions without a normalization
alignment or a fresh projection fit.

### Energy/observable selection

The initial H-DMD spectrum is reduced with:

```matlab
rom = Koop.ROM(1).extract(
    method="observable",
    threshold=[0.99,0.99,0.99],
    caps=[30,3,200]);
```

This takes the union of modes needed to represent 99 percent of the modal
amplitude energy for radius, radial velocity, and echo, subject to distinct
per-observable caps. It is not an observability-optimal selection; it is
observable-reconstruction-energy selection followed by a union.

The later call

```matlab
rom.extract(method="maxO", estm_O=1, measurable=3, threshold=0.9)
```

is an unfinished experiment, not a reliable algorithm. The code computes
candidate observability matrices but does not complete a sound scalar ranking
or replace the selected eigenfunction set. Specific defects include:

- candidate pairs from `nchoosek` are indexed as individual elements;
- zero entries are removed from `C`, destroying coordinate alignment;
- `V(:,1)' * diag(S)` produces a vector, not a scalar score;
- the selected candidate set is not cleanly promoted into `keigf` before the
  extracted KROM is constructed.

It should not be ported.

## Observable projection and reconstruction

### Working arbitrary-observable projection

`KROM.compute_nu` evaluates every requested observable on sampled states,
evaluates the full complex eigenfunction basis including conjugate duals, and
solves

```text
Y = C Psi,
C = Y pinv(Psi).
```

This is the MATLAB function call the user remembered for obtaining modal
amplitudes by projecting observables onto a supplied eigenfunction basis. It is
more accurately a global least-squares reconstruction fit than an orthogonal
projection, because the learned eigenfunctions are not assumed orthonormal.

The mathematical shapes in MATLAB's sample-major language are:

```text
Psi:  number_of_full_complex_modes x number_of_samples
Y:    number_of_observable_outputs x number_of_samples
C:    number_of_observable_outputs x number_of_full_complex_modes.
```

The first two rows can be used as `Cx` for radius and radial velocity; the echo
row can be used as `Cy`.

### Incomplete alternatives in MATLAB

`observable.project` is not executable as written. It uses an undefined `x`:

```matlab
vpt = obsv.g(x) * conj(obj.g(x));
```

and has no sampling measure, quadrature, sum, or normalization. Likewise, the
fallback in `ROM.add_observable` is disabled by `elseif ~1`. If enabled, it
would sum trajectory products and divide only by trajectory count, so its scale
would depend on the number of samples. Neither is the production path.

`compute_nu` itself also has fragile large-dimension code: it chunks the mode
row dimension at 50,000 and concatenates independent pseudoinverse solutions,
which is not equivalent to one global least-squares problem. This is harmless
for the observed KROM sizes because the chunk count is one, but it should not be
copied.

### Complex and real-coordinate convention warning

MATLAB's expanded eigenfunction coordinate is

```text
[Re(psi), -Im(psi)].
```

The corresponding dynamic block is the transpose of the usual
`[Re(psi), Im(psi)]` block, which explains the explicit transpose in the manual
observer:

```matlab
A = rom.get_A(form="expanded",mode="Continuous")';
```

For a real observable reconstructed from a nonredundant complex conjugate pair,
the exact contribution is

```text
nu psi + conj(nu psi) = 2 Re(nu psi).
```

`KROM.get.C` emits `[Re(nu), Im(nu)]` for the primary complex mode, with no
visible factor of two. Unless `nu` was doubled elsewhere—and the H-DMD path does
not visibly do that—the expanded `C` reconstructs only one real part rather
than the full conjugate-pair contribution. The full-complex path that explicitly
includes both duals does not have this ambiguity. This factor must be checked
against exported MATLAB `C Psi` values before using an expanded MATLAB matrix as
a benchmark.

The current Python convention is cleaner:

```text
z_split = [Re(psi), Im(psi)]
A_pair  = [[a, -b], [b, a]]
C_pair  = [2 Re(nu), -2 Im(nu)]
```

for a retained eigenvalue `a + i b`. Python's `ObservableReconstruction` makes
this conversion explicitly.

## KLQE equations actually implemented in MATLAB

Let

```text
z = Psi(x),        z_dot = A z
y = h_echo(x),     yhat = Cy zhat
xhat = Cx zhat.
```

The class-based free-dynamics estimator computes

```text
innovation = y - Cy zhat
zhat_dot   = A zhat + L innovation.
```

The forced form adds the lifted control input:

```text
zhat_dot = A zhat + Bpsi(xhat,t) u + L innovation.
```

In `ESTM.m`, `Cx` is not a declared object. State estimation is hard-coded as
the first two rows of the plant `C`:

```matlab
xhat = [Cs{1}*zhat; Cs{2}*zhat];
```

and `Cy` is selected through an integer `Cidx`. This only works while the
observable order remains exactly `[radius, radial velocity, ..., echo]`.

The gain is the dual CARE solution:

```text
A P + P A' - P Cy' Re^-1 Cy P + Qe = 0
L = P Cy' Re^-1.
```

MATLAB evaluates the same result as `lqr(A',Cy',Qe,Re)'`. The manual RPE case
uses highly ad hoc tuning:

```text
Qe = 10000 I
Qe[1,1] = 1000
Qe[2,2] = 0.01
Re = 1000000.
```

The comments call these disturbance and noise covariances, but they have not
been identified from data. Because the echo has a large dimensional scale while
the lifted coordinates are arbitrarily normalized, these values currently act
as observer tuning weights as much as physical covariance estimates.

The observer starts from `zhat(0)=0`. A zero lifted vector is generally not on
the Koopman image of any physical state and is incompatible with a constant
eigenfunction equal to one. It is acceptable as a convergence stress test, but
it must be labeled and compared with physically consistent initialization.

## Continuous/discrete-time inconsistency in the MATLAB paths

There are three distinct paths:

1. **Manual `KLQE.m` path:** explicitly calls `get_A(..., mode="Continuous")`
   and transposes for the `[Re,-Im]` convention. This is the closest to the
   intended continuous observer.

2. **Reusable `ESTM` path:** calls `plant.model(..., mode="Continuous")`, but
   `plant.model` immediately assigns `mode = sys.mode`. An H-DMD KROM's
   `sys.mode` is discrete. The resulting `A` can therefore contain discrete
   eigenvalues even though it is passed to continuous `lqr` and an ODE.

3. **Old `KoopmanKalman/KKF.m` path:** loads `A.csv` and uses it directly in
   continuous `ss` and `lqr`. The current `702 x 702` file has eigenvalue
   magnitudes between approximately `0.99693` and `1.00002`, and all eigenvalues
   have positive real parts between about `0.695` and `1.000`. It is clearly a
   discrete transition matrix, not a continuous generator. That stored
   prototype is not a valid continuous LQE benchmark.

The separate failsafe `cA.csv` is `100 x 100` and has oscillatory
continuous-looking blocks (frequencies up to about `353.3` in its stored time
units). Its `cC.csv` is one row with 100 columns and is dimensionally compatible
with the failsafe script. It is the better standalone seed, although its exact
provenance and expected output still need to be exported.

## Classical observability and its SVD

For an `n`-coordinate discrete lifted model and a selected measurement row
`Cy`, MATLAB calls:

```text
O_n = [Cy;
       Cy A;
       ...;
       Cy A^(n-1)]

O_n = U Sigma V'.
```

Interpretation:

- `rank(O_n) = n` is the exact-arithmetic algebraic observability condition.
- `sigma_i` measures the amplification of the lifted direction `V(:,i)` into
  the stacked measurements over the chosen algebraic horizon.
- `V(:,1)` is the most observable right-singular direction and `V(:,end)` the
  least observable direction under the selected Euclidean coordinates.
- `U(:,i)` describes the corresponding time/measurement-history signature.

The existing plot

```matlab
stem(V(n,:)); title(sigma_n)
```

uses row `n` of `V`, not column `n`; it therefore does not plot the direction
associated with `sigma_n`.

### Important limitations of raw `obsv` SVD

The SVD is coordinate- and scale-dependent. Koopman eigenfunctions can be
rescaled arbitrarily, and `C` rescales inversely. Raw singular values and vector
entries change under that otherwise equivalent representation. Echo units and
state units add another scale dependence.

For a diagonal complex KROM with distinct eigenvalues, the observability matrix
is a coefficient-weighted Vandermonde matrix. In exact arithmetic a scalar
measurement can observe every distinct mode for which `Cy_i != 0`. In floating
point, clustered eigenvalues, tiny echo coefficients, and a long horizon can
make those directions effectively unobservable even when the matrix has formal
full rank.

Powers of a continuous generator are especially unsafe. The saved
`RPE_nonvisc_good_Oy_correct_echo.csv` is `159 x 159` but currently contains
5,887 nonfinite entries, with finite magnitudes reaching approximately
`6.94e307`. The saved viscous matrix is `22 x 22`, has a raw 2-norm condition
number around `3.67e49`, and is numerically rank 4 under NumPy's default
tolerance. These artifacts strongly suggest overflow and conditioning problems
from unscaled algebraic powers; they are not evidence that the physical system
has only four observable directions.

### Recommended observability products

The future analysis should produce all of the following, with explicit time
mode and horizon metadata:

1. **Discrete finite-horizon matrix**

   ```text
   O_N = [Cy; Cy Ad; ...; Cy Ad^(N-1)].
   ```

2. **Discrete finite-horizon Gramian**

   ```text
   Wo(N) = sum(k=0..N-1) (Ad^k)' Cy' R^-1 Cy Ad^k
         = O_weighted' O_weighted.
   ```

3. **Continuous finite-horizon Gramian**

   ```text
   Wo(T) = integral(0..T) exp(A' t) Cy' R^-1 Cy exp(A t) dt.
   ```

4. **Infinite-horizon continuous Gramian**, only when `A` is asymptotically
   stable:

   ```text
   A' Wo + Wo A + Cy' R^-1 Cy = 0.
   ```

The finite-horizon definitions remain meaningful for marginal or weakly
unstable Koopman modes. The measurement-noise weighting gives the singular
directions a direct signal-to-noise interpretation.

For display, each complex conjugate pair should be treated as one indivisible
mode group. A direction's pair participation should be shown as
`sqrt(v_re^2 + v_im^2)`, which is invariant to rotations within the real/imag
pair. Useful native plots are:

- log singular-value scree with rank/threshold annotations;
- singular-direction-by-mode heat map;
- grouped eigenfunction participation bars for selected directions;
- spectrum colored by modal observability participation;
- optional cumulative observable energy and effective-rank plots.

## Current Python capability audit

### What already exists

#### Observable contract

`dynsyslib.observables.Observable` formalizes `y=h(x)`, dimensions, optional
Jacobian, value type, descriptor metadata, and coordinate-system ownership.
RPE now exposes `far_field_acoustic_echo` in its observable catalog along with
radius, radial velocity, and the balance-term observables.

#### Koopman `C`

`rom.observation.ObservableReconstruction` stores:

- the declared `Observable`;
- nonredundant complex modal coefficients;
- an intercept;
- validation quality;
- provenance.

It can return the correct real split matrix and reconstruct the observable from
lifted coordinates. `fit_observable_reconstruction` performs a held-out affine
least-squares fit. `transfer_observable_reconstruction` aligns H-DMD and SINDy
eigenfunction normalizations before transferring a `C` matrix.

The candidate full-complex fitter in
`rom/candidates/observable_projection.py` reproduces MATLAB's no-intercept
`Y = C Psi` calculation without promoting it into the stable KROM API.

#### H-DMD modal amplitudes

The window-snapshot H-DMD path computes

```text
nu = pinv(projected_modes) @ terminal_delay_column
```

and attaches it as an `ObservableReconstruction` for the scalar source
observable. The result is intentionally marked unvalidated until held-out
reconstruction is checked.

The Python H-DMD path currently accepts one scalar observable in this layout,
not MATLAB's horizontal multi-observable stack. Additional `C` rows must be fit
on the resulting KROM, or the H-DMD interface must later be extended to a
declared observable stack.

#### KLQE class

`estimation.koopman.KLQEEstimator` already:

- obtains the continuous generator with `rom.A(output="split")`;
- solves the continuous dual CARE;
- evolves a lifted estimate using the innovation;
- accepts an explicit physical measurement function;
- accepts an explicit lifted-to-state reconstruction function.

This is a useful prototype, but it is not yet a production RPE KLQE path.

### Critical gaps before a benchmark claim

1. **No native classical observability analysis.** The current
   `select_observable_modes` routine uses state-reconstruction regression and
   labels that concept “observability,” but it does not construct `obsv(A,C)` or
   a Gramian. The two analyses must remain separately named.

2. **Only one attached observable-reconstruction object.** A KLQE case needs at
   least a state reconstruction (`Cx`) and echo measurement (`Cy`) with stable
   named row selection. A stacked `[radius, radial_velocity, echo]` observable
   can represent this today, but the KROM API does not expose a convenient
   named registry or row-selection contract.

3. **The default KLQE measurement is unsuitable for the RPE case.** If no
   `measurement_map` is supplied, the estimator measures `rom.psi(x)`, not the
   acoustic echo. RPE KLQE must require the echo `Observable` explicitly.

4. **Explicit-Euler observer propagation.** The estimator advances with one
   Euler step per increasing sample-clock time. Explicit Euler is intrinsically
   unstable for undamped imaginary modes and can be poor for a large KROM with
   high frequencies or a high observer gain. It must be replaced by an exact
   or stable propagator before large-mode tests.

5. **Estimator execution is coupled to controller evaluation.** The ODE solver
   currently invokes the estimator only while obtaining a controller signal.
   A free-system state-estimation experiment with no controller has no native
   estimator-only execution and archive path.

6. **No estimate/innovation archive.** There is no canonical time series for
   measured echo, predicted echo, innovation, lifted estimate, reconstructed
   state, or estimation covariance/quality.

7. **No KLQE API/CLI or Viewer report.** The class is importable but there is no
   lightweight command for analysis/design/run, no reproducible spec artifact,
   and no standard Viewer dashboard.

8. **No KLQE-specific tests.** The generic one-dimensional LQE has a smoke test;
   `KLQEEstimator` itself has no direct unit or integration test.

9. **No high-dimensional CARE gates.** Shape, symmetry, definiteness,
   detectability, stabilizability, observer-pole stability, and conditioning
   are not checked or reported before accepting a gain.

10. **H-DMD interpolation-state delay misalignment.** For a delay count `d`,
    the DMD mode row for window `i` corresponds to the window-end state
    `x[i+d-1]`. Python currently associates it with `x[i]`. MATLAB associates
    the rows with its `idx2 = n:M`, i.e. the end-of-window states. This should be
    fixed and tested independently of KLQE.

## Candidate clean architecture

This is a design proposal only.

### 1. Observable reconstruction remains a ROM concern

A KROM should retain a named collection of observable projections, each with:

```text
observable identity and metadata
C in nonredundant complex coordinates
C in real split coordinates (derived, not separately authoritative)
intercept
fit/validation diagnostics
source method: HDMD-native, least-squares, transferred, or user-provided
training/test dataset provenance
coordinate system and units
eigenfunction normalization fingerprint
```

Suggested stable calls:

```python
krom.observable_map("far_field_acoustic_echo")
krom.observable_map("state")
krom.C("far_field_acoustic_echo", output="split")
krom.C("state", output="split")
```

The existing single `ObservableReconstruction` can be evolved into this
collection without changing its coefficient mathematics.

### 2. Observability analysis is a separate immutable analysis artifact

Suggested call:

```python
report = analyze_observability(
    krom,
    measurement="far_field_acoustic_echo",
    time_mode="discrete",
    horizon_steps=500,
    measurement_covariance=R,
    coordinate_scaling="training-covariance-whitened",
)
```

The report should retain `A`, `Cy`, the time convention, horizon, scaling,
singular values/eigenvalues, right directions, numerical rank under declared
tolerances, grouped mode participation, and source fingerprints. Viewer should
render the artifact and never recompute the analysis silently.

### 3. Estimator design is an estimation concern

A declarative specification should identify:

```text
KROM archive and fingerprint
measurement observable and Cy projection
state output observable and Cx projection
continuous or discrete estimator formulation
Q/GWGT and R with units and provenance
initial lifted estimate and covariance policy
sample time / interpolation policy
gain and observer-pole diagnostics
```

The design boundary should refuse a continuous CARE with a discrete `A`, or a
discrete DARE with a continuous `A`.

### 4. Propagation should match the measurement model

For continuous `A` and piecewise-constant measurements over `dt`, define

```text
M = A - L Cy.
```

Then propagate with an augmented matrix exponential, which remains valid when
`M` is singular:

```text
[F G; 0 I] = exp([[M, L], [0, 0]] dt)
zhat[k+1] = F zhat[k] + G y[k].
```

If forcing is later included, its held input matrix can be appended to the same
augmented exponential. For a genuinely discrete KROM, use a discrete Kalman or
Luenberger recursion and a DARE-derived gain. Do not convert back and forth
inside the time loop.

### 5. Solver and offline filtering are separate runners

- **Simulation runner:** integrate plant and observer together, but archive
  plant truth, measurement, estimate, and innovation as separate time-series
  objects.
- **Offline runner:** accept a measured echo `TimeSeries` and a KROM, with no
  access to true state except optional benchmark scoring. This is the path that
  most faithfully represents an acoustic sensor.
- **Controller composition:** a controller may consume the estimator's state,
  but free-estimation experiments must not require a dummy controller.

### 6. Viewer is a renderer

A standard KLQE Viewer report should show:

1. measured echo, reconstructed echo, and innovation;
2. true versus estimated radius;
3. true versus estimated radial velocity;
4. true versus estimated phase portrait;
5. selected lifted-coordinate truth/estimate panels when truth is available;
6. estimation error norms and convergence envelope;
7. observer poles;
8. observability singular-value scree;
9. grouped most/least observable KROM directions.

Truth-dependent panels should be marked benchmark-only and omitted for real
measurement-only runs.

### 7. Lightweight API/CLI boundary

Candidate commands, not yet implemented:

```text
python -m estimation klqe analyze-observability KROM.npz --measurement far_field_acoustic_echo
python -m estimation klqe design KROM.npz --measurement far_field_acoustic_echo --state-output state --q ... --r ...
python -m estimation klqe run DESIGN.json ECHO_SOLUTION.npz --report-only
python -m estimation klqe benchmark DESIGN.json TRUTH_SOLUTION.npz --viewer OUTPUT.html
```

Every command should be report-only unless an explicit publish/output path is
given, and should create immutable derived artifacts rather than mutate a KROM.

## Recommended isolation and benchmark sequence

No RPE filter result should be interpreted until the preceding gates pass.

### Gate 0: coordinate and matrix identity

For a frozen KROM, verify:

- `Ad = exp(Ac dt)` mode by mode;
- complex and split-real `A` propagate the same eigenfunctions;
- complex and split-real `C Psi` give the same observable;
- conjugate-pair factors and signs are exact;
- every `C` is tied to a KROM normalization fingerprint;
- state, time, echo, `Q`, and `R` units are explicit.

### Gate 1: H-DMD window alignment

Use a known analytic oscillator and a scalar observable. Confirm that H-DMD
mode row `i` is attached to the end-of-window physical state. Test both one
trajectory and an ensemble. Only then enable direct H-DMD interpolation in an
estimator.

### Gate 2: static observable reconstruction

Before running any observer, evaluate on both training and held-out free RPE
trajectories:

```text
r_actual       versus Cr Psi(x)
rdot_actual    versus Crdot Psi(x)
echo_actual    versus Cecho Psi(x).
```

Report RMSE, normalized RMSE, maximum error, `R2`, conditioning, and
out-of-domain status. A filter cannot recover information that its `Cy` cannot
represent, and it cannot reconstruct a state that `Cx` cannot represent.

### Gate 3: synthetic lifted observer

Generate `z(t)` directly from the frozen KROM `A`, generate exact
`y(t)=Cy z(t)`, and recover `z` from a deliberately incorrect initial estimate.
This removes all eigenfunction and nonlinear-plant errors. Compare exact
propagation with the previous Euler prototype.

### Gate 4: standalone MATLAB failsafe parity

Use the coherent `cA.csv`, `cC.csv`, and `cz0.csv` set only after a wrapper
exports:

```text
time
input/measurement time series
Q, R
L
observer poles
zhat(t)
Cy zhat(t)
solver/time convention.
```

Python must reproduce `L`, poles, and `zhat(t)` within declared tolerances. The
current 702-state `A.csv/C.csv` prototype is rejected as a benchmark until its
time mode and dimensional mismatch are repaired or regenerated.

### Gate 5: noiseless nonlinear RPE echo estimation

Freeze one RPE free trajectory, the exact KROM, `Cx`, `Cy`, `Q`, `R`, and
`zhat0`. Drive the observer with `echo(x(t))` only. Compare:

- MATLAB manual path;
- Python offline exact-propagation path;
- Python joint simulation path.

Measure echo innovation, state RMSE, phase error, convergence time, and terminal
error. Keep the plant trajectory fixed across implementations.

### Gate 6: known-noise benchmark

Inject seeded measurement noise and, separately, seeded process disturbance.
Use physically declared covariance and compare empirical innovation/error
statistics. Report whether the innovation is approximately zero mean and
temporally white. This is the first gate that justifies calling the result a
Koopman-Kalman filter rather than only a Koopman-Luenberger observer.

### Gate 7: robustness and extrapolation

Evaluate new free initial conditions inside and outside the eigenfunction
training domain. Separate:

- eigenfunction evolution error;
- `Cx`/`Cy` reconstruction error;
- numerical observer error;
- lack of observability;
- extrapolation beyond the KROM domain.

## Minimal MATLAB parity payload needed next

A noninteractive wrapper should export one self-contained MAT file containing:

```text
RPE physical parameters and nondimensional groups
initial condition, time vector, and integration method/tolerances
truth states r and rdot
measured acoustic echo
full complex and expanded-real eigenvalues/A matrices, with time mode
eigenfunction normalization identifiers or sampled Psi on the trajectory
Cr, Crdot, and Cecho in full complex and expanded-real conventions
Q, R, gain L, observer poles, and zhat0
zhat(t), xhat(t), yhat(t), and innovation(t)
observability horizon, O or its singular values, and V columns
mode selection indices and conjugate-pair map
```

The wrapper should call existing MATLAB code and save outputs only; it should
not modify the suite. Exporting both complex and expanded representations is
essential because the current MATLAB sign/factor conventions are otherwise
ambiguous.

## Decisions recommended for review

1. Treat the manual continuous-time block in `KLQE.m` as the intended observer
   equation, not the mode-overwriting reusable path or the stale 702-state CSV
   prototype.
2. Treat `KROM.compute_nu` and H-DMD's per-observable modal amplitudes as the two
   MATLAB `C` sources; do not use `observable.project` as a reference.
3. Preserve Python's explicit nonredundant-complex convention and conjugate-pair
   factor of two.
4. Add named observable reconstruction maps before adding a public KLQE command.
5. Implement finite-horizon, noise-weighted observability analysis as an
   immutable ROM analysis, while retaining raw classical `obsv` only for parity.
6. Replace explicit-Euler lifted observer updates with exact/stable propagation.
7. Decouple free estimation from controllers and archive estimator histories.
8. Fix and unit-test H-DMD end-of-window interpolation alignment before direct
   H-DMD KLQE.
9. Validate `Cx` and `Cy` separately before attributing any failure to the
   observer gain.
10. Do not promote a “Kalman” benchmark until seeded nonzero-noise tests pass.

## Bottom line

The MATLAB suite contains the full conceptual chain the user remembered:

```text
free RPE trajectory
  -> [r, rdot, acoustic echo] horizontal H-DMD
  -> eigenvalues/eigenfunction samples plus modal amplitudes C
  -> observable-energy sub-KROM
  -> echo row Cy and state rows Cx
  -> dual-LQR/LQE gain
  -> lifted Luenberger observer
  -> reconstructed r and rdot
  -> classical observability matrix and SVD.
```

It does not yet provide a single trustworthy reference implementation of every
step. The manual observer is the best statement of intent; the reusable class
has a discrete/continuous overwrite; the old CSV prototype mixes a discrete
matrix with continuous tools; observable-SVD plotting uses rows instead of
right-singular-vector columns; the `maxO` selector is unfinished; and raw
observability powers can overflow catastrophically.

Python already has better boundaries for observables, KROM coordinates, `C`
objects, and continuous eigenvalue handling. The next work should be a narrow,
reversible KLQE addition built around those boundaries and verified gate by
gate against a frozen MATLAB payload—not a port of the MATLAB estimator classes.
