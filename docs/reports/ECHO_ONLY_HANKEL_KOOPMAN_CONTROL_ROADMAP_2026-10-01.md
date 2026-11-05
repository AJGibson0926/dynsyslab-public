# Echo-only causal Hankel–Koopman control roadmap

**Status:** approved future-work plan and experimental policy

**Branch:** `echo-only`

**Date:** 2026-10-01

**Implementation state:** planning plus model-blind research evidence; this
document does not introduce a production numerical pipeline

## Executive decision

DynSysLab will investigate a model-free control pathway in which the bubble's
physical state, governing equations, analytic echo map, and echo derivatives
are unavailable to the learning and online-control system. The only plant
measurement required by the proposed controller is a sampled echo signal.
Quantities derived causally from the past echo, including Hankel/POD latent
coordinates, are permitted. Controlled identification additionally requires a
time-synchronized record of the applied actuator input; this is control data,
not a hidden-state measurement.

The primary formulation will be discrete time. A causal delay vector will be
formed from a rolling echo buffer, projected through a fixed basis learned from
training data, and evaluated by sparse sampled-Koopman eigenfunctions. The
control response will be learned statistically from paired input/output data.
The eventual controller will therefore require neither radius nor wall
velocity, neither an RPE model nor an analytic echo expression, and neither
echo-dot nor any online numerical differentiation.

The central hypothesis is:

> Over a deliberately bounded operating regime, the current echo together
> with a compact coordinate extracted from its causal history contains enough
> information to support a well-conditioned sampled Koopman representation and
> an identifiable control response.

This is a falsifiable hypothesis, not an architectural assumption. The first
work on this branch must measure whether a two-variable chart is sufficient.
The implementation must retain the ability to use more than one latent
coordinate when the data reject the two-variable model.

## Scope and non-goals

### In scope

- Scalar echo measurements sampled on a documented clock.
- Causal functions of echo history, including delay vectors, Hankel/POD
  coordinates, filters, and learned encoders.
- Sparse polynomial eigenfunctions of measured/history-derived coordinates.
- Sampled Koopman and input/output models learned from data pairs.
- Statistical constant, state-dependent, or bilinear control maps.
- Discrete KLQR, LQR, or a closely related constrained sampled-data controller.
- Simulation experiments in which the training algorithm is deliberately
  blinded to the available physical model.
- Later laboratory identification from echo and synchronized actuator data.
- Exact provenance, run-level validation splits, uncertainty diagnostics, and
  conservative fail-safe behavior.

### Out of scope for the first implementation

- Claiming that a learned latent coordinate equals radius, wall velocity, or
  another physical state.
- Recovering a unique continuous-time vector field from one sampled data set.
- Treating singular-value decay as proof of observability or finite Koopman
  closure.
- Requiring one latent scalar when validation supports a larger state.
- Replacing independent laboratory safety instrumentation with echo tracking.
- Publishing generated solves, screenshots, plots, or Viewer payloads in Git.

## Terminology: echo-only versus input-free

“Echo-only” describes the plant-state measurement available to estimation and
feedback. It does not mean that actuator commands may be discarded.

Let

\[
y_k=e(t_k)
\]

be the measured echo and let

\[
u_k
\]

be the command applied during the transition from sample \(k\) to sample
\(k+1\). If possible, a calibrated measurement of delivered acoustic pressure
should be retained separately from electrical command voltage. The autonomous
sampled dynamics may be studied from \(y_k\) alone. The input response cannot
be identified uniquely unless the applied input, or a sufficiently faithful
proxy for it, is recorded.

Accordingly, the final data-only contract is better summarized as:

\[
\text{echo-only state measurement} + \text{known actuation history}.
\]

## Why discrete time is the primary formulation

The present continuous-time observable pipeline evaluates Lie derivatives by
combining analytic observables with known equations of motion. That route is
powerful while the RPE state and vector field are available, but it is not the
right production dependency for a laboratory-only system.

A continuous-time model could still be estimated from data by smoothing and
differentiating the echo, differentiating delayed copies, or using weak and
integral formulations. None is logically impossible. They nevertheless add
filter phase, derivative uncertainty, boundary handling, and continuous-time
generator ambiguity to a problem that is natively sampled by the measurement
and control hardware.

The discrete relation

\[
\psi(z_{k+1})=\mu\psi(z_k)
\]

uses measured snapshot pairs directly. If a continuous eigenvalue \(\lambda\)
is prescribed for scientific comparison, the sampled multiplier is

\[
\mu=\exp(\lambda\Delta t).
\]

The implementation should retain \(\mu\) as the primary identified quantity.
Mapping it back to \(\lambda\) requires a logarithm branch and an explicit
anti-aliasing argument. Controller design should likewise be discrete rather
than applying a continuous Riccati design to a model known only at the sample
interval.

## Proposed mathematical pipeline

### 1. Experimental records

Each run must retain at least:

- echo samples \(y_k\);
- commanded input \(u_k\);
- calibrated delivered input when available;
- sample timestamps rather than only an assumed sample count;
- run and apparatus identity;
- sensor, amplifier, transducer, and filter configuration;
- equilibrium preparation and relevant environmental conditions;
- saturation, clipping, drop-out, and timing-quality flags.

Data must be split by complete run, not by randomly shuffling overlapping
delay windows. Random window splits would leak almost identical histories into
training and validation sets.

The initial experimental program should include distinct data roles:

1. unforced runs for autonomous evolution and all representation training;
2. safe persistently exciting runs for separate input-response identification;
3. held-out periodic targets;
4. held-out amplitudes, phases, and perturbations;
5. repeated runs acquired on different days for drift assessment.

The training policy is unforced-only. Delay-basis fitting, normalization,
rank selection, autonomous sampled dynamics, and Koopman-eigenfunction fitting
must use free-response runs with no applied driver or feedback. Forced target,
validation, and actuator-identification records are downstream assets; they
must not refit or otherwise leak into the learned echo-history representation.

### 2. Causal delay state

For a delay length \(q\), construct

\[
d^y_k=
\begin{bmatrix}
y_k & y_{k-1} & \cdots & y_{k-q+1}
\end{bmatrix}^{T}.
\]

The window ends at the current measurement. Future-looking or centered
windows are prohibited in the online-capable contract.

If actuator memory is material, construct an input-history vector

\[
d^u_k=
\begin{bmatrix}
u_{k-1} & u_{k-2} & \cdots & u_{k-p}
\end{bmatrix}^{T}.
\]

The current action \(u_k\) remains separate: it causes the transition to the
next sample. A measured actuation delay may shift this alignment by an integer
or fractional number of samples. That delay must be estimated and archived,
not silently absorbed by the regression.

No delay window may cross a run boundary, discontinuity, missing-data interval,
or sample-clock reset.

### 3. Training-only normalization

Using training windows only, calculate

\[
\bar d=\frac{1}{N}\sum_k d^y_k
\]

and a fixed scale matrix \(S\). The normalized current window is

\[
\widetilde d_k=S^{-1}(d^y_k-\bar d).
\]

The mean and scale are immutable parts of the learned observable. They must be
used identically for training, validation, target construction, replay, and
online control. Validation or target data must never influence them.

### 4. Hankel/POD coordinate extraction

Assemble normalized training windows as columns:

\[
\widetilde D=
\begin{bmatrix}
\widetilde d_{k_1} & \widetilde d_{k_2} & \cdots
\end{bmatrix}
=U\Sigma V^{*}.
\]

For a retained basis \(U_r\), define raw POD scores

\[
\xi_k=U_r^{*}\widetilde d_k.
\]

If the HAVOK-normalized convention is used instead, include the corresponding
\(\Sigma_r^{-1}\) factor. The choice is not scientifically decisive, but its
convention must be versioned because it changes all downstream coefficients.

The first candidate chart may be

\[
z_k=\begin{bmatrix}y_k&\xi_{1,k}\end{bmatrix}^{T}.
\]

However, the current echo is already represented inside every delay vector.
Blindly selecting the first POD score can therefore be redundant. Candidate
selection should compare at least:

- the leading ordinary POD coordinate;
- the leading delay coordinate after regressing out current echo;
- a coordinate selected for next-step prediction conditional on current echo;
- two or more retained coordinates when one coordinate fails closure tests.

POD orders directions by variance in the chosen norm. It does not order them
by observability, controllability, or suitability for feedback.

### 5. Output-only closure and observability diagnostics

The learned chart must be tested before fitting a controller. These tests can
be performed without radius or wall velocity.

For pairs satisfying

\[
z_i\approx z_j,\qquad u_i\approx u_j,
\]

compare their successors \(z_{i+1}\) and \(z_{j+1}\). Large successor
dispersion means the proposed chart has hidden memory and is not approximately
Markov at the requested resolution.

Required diagnostics include:

- nearest-neighbor forward-dispersion distributions;
- conditional variance of \(z_{k+1}\) given \((z_k,u_k)\);
- false-nearest-neighbor behavior as delay dimension and latent rank increase;
- one-step and multi-step prediction on held-out runs;
- sensitivity to window length, sample interval, and retained rank;
- local conditioning and noise amplification of the encoder;
- consistency across amplitudes, phases, transients, and acquisition days;
- comparison with direct delay-vector models that perform no severe
  compression.

Simulation may use hidden \((r,v)\) values as an audit channel during early
research, but those values must be withheld from training and online
evaluation. This produces a useful distinction between “the echo-only method
works” and “the simulation's hidden state explains why it works.”

### 5a. Progressive latent-dimension escalation policy

Echo-only experiments must expose the consequences of compression instead of
silently replacing a failed low-dimensional representation with a larger one.
The default experimental ladder is therefore

\[
r=2\;\longrightarrow\;r=3\;\longrightarrow\;r=4,
\]

starting with \((\xi_1,\xi_2)\). A study may proceed directly from two to four
coordinates only when the reason is documented, but it must still retain a
head-to-head two-versus-four comparison. A known or expected failure of the
two-coordinate chart is not permission to omit it: reproducing that failure
and tracing its downstream consequences are required scientific evidence.

The following rules govern every echo-only rank study:

1. **Begin at the lowest declared rank.** Fit and evaluate eigenfunctions,
   prediction, input response, and control first with \((\xi_1,\xi_2)\).
2. **Hold the parent representation fixed.** Rank comparisons must reuse the
   same training runs, causal delay window, training-only normalization,
   learned \(U\), singular-value convention, validation runs, and sample
   alignment. Promotion adds retained columns; it does not perform a new SVD
   or retune preprocessing unless a separate experiment explicitly studies
   that change.
3. **Record failures rather than overwriting them.** Each rank receives a
   separate immutable result with its coordinate-collision evidence,
   successor dispersion, eigenfunction residuals, prediction errors, control
   behavior, boundedness, and operating envelope.
4. **Promote for a stated reason.** Examples include non-injectivity,
   false-nearest-neighbor events, non-Markov successor dispersion,
   multivalued inverse behavior, unacceptable held-out eigenfunction residual,
   or failed closed-loop tracking.
5. **Measure the repair and its costs.** A higher rank must be judged both on
   whether it fixes the triggering failure and on new complications it may
   introduce, including conditioning, noise amplification, data and
   excitation requirements, polynomial-library growth, sparsity instability,
   overfitting, controller complexity, and new failure modes.
6. **Keep comparisons controlled.** Controller targets, horizons, costs,
   saturation limits, random seeds, training/validation partitions, and
   scoring rules remain identical across ranks unless the changed quantity is
   the explicit subject of the experiment.
7. **Use hidden state only for audit.** Simulation-only \((r,v)\) may diagnose
   folds and physical ambiguity, but it may not enter representation fitting,
   eigenfunction evaluation, or online control.
8. **Do not promote merely because energy capture is low.** Rank promotion is
   driven by held-out sufficiency, injectivity, prediction, and control
   evidence. Conversely, visually attractive phase curves or high captured
   energy do not excuse a failed collision or closure test.

The initial inviscid RPE audit motivating this policy used one fixed basis and
five unseen unforced trajectories. With \((\xi_1,\xi_2)\), median cycle-closure
error was approximately 0.45% for trajectories inside the training-state
envelope and 4.65% outside it. Proper representative-cycle intersections were
distributed as two inside–inside crossings, seven inside–outside crossings,
and ten outside–outside crossings. Adding \(\xi_3\) separated every detected
two-dimensional crossing, while retaining through \(\xi_4\) increased the
minimum and median separation margins. This evidence supports the escalation
policy, but it does not pre-authorize four coordinates as universally
sufficient.

### 5b. Delay-window and multi-trajectory training policy

The causal history length is a representation hyperparameter, not merely an
implementation detail. For each unforced training run (j), construct its own
Hankel block

\[
D^{(j)}=
\begin{bmatrix}
d^{(j)}_{q-1}&d^{(j)}_q&\cdots&d^{(j)}_{N_j-1}
\end{bmatrix},
\qquad
d^{(j)}_k=
\begin{bmatrix}
e^{(j)}_{k-q+1}&\cdots&e^{(j)}_k
\end{bmatrix}^{\mathsf T},
\]

and fit one shared basis from the horizontal concatenation

\[
D_{\mathrm{pool}}=
\begin{bmatrix}
w_1D^{(1)}&w_2D^{(2)}&\cdots&w_MD^{(M)}
\end{bmatrix}.
\]

No history window may cross a trajectory boundary. Runs must share the same
sample interval and echo convention, normalization statistics must be fitted
from training runs only, and complete held-out runs—not random windows from a
training run—must be used for validation. Equal-duration runs provide equal
block weight when (w_j=1); otherwise choose (w_j\propto N_j^{-1/2}) when
equal contribution per run is intended. Weighting, amplitude stratification,
and operating-envelope coverage must be recorded with the retained basis.

The first inviscid RPE delay study refit the top-two basis at windows from
0.125 to 4.0 reference periods and evaluated the same five unseen unforced
orbits. The former 0.885-period window had 19 proper cross-orbit intersections.
Two periods still had 19, whereas three and four periods had zero. The effect
was non-monotone: one- and 1.5-period windows were worse than the former
baseline. At four periods, median cycle-closure error was approximately 0.41%
inside and 1.50% outside the original training-state envelope. The top-two
energy fraction did not increase with the improvement, so energy capture is
not an acceptable proxy for injectivity.

A companion pooled fit used five equal-duration unforced training orbits and
the same held-out set. Pooling did not repair the short 0.885-period chart: it
also produced 19 intersections. It reduced the two-period result from 19 to
10 intersections and reduced four-period outside-envelope closure error from
about 1.50% to 1.28%; three and four periods again had zero intersections.
Thus broader training coverage improves generalization and separation margin,
but it does not substitute for enough causal history. The default sequence is:

1. sweep causal window length using held-out injectivity, closure, successor
   dispersion, conditioning, and noise sensitivity—not singular-value energy;
2. pool multiple unforced trajectories spanning the intended amplitude or
   action range using boundary-safe Hankel blocks;
3. freeze the chosen window, normalization, pooled basis, and data split;
4. only then apply the progressive (r=2\rightarrow3\rightarrow4) rank study.

The current noiseless inviscid evidence makes a roughly three-period window a
reasonable minimum candidate, not a universal constant. Laboratory noise,
transients, damping, forcing at evaluation time, and causal warm-up cost still
require separate validation.

### 6. Sparse sampled-Koopman eigenfunctions

Choose a dictionary evaluated only on online-accessible coordinates. For a
two-variable polynomial example,

\[
\Theta(y,\xi)=
\begin{bmatrix}
1&y&\xi&y^2&y\xi&\xi^2&\cdots
\end{bmatrix}.
\]

Represent one eigenfunction as

\[
\psi(z)=\Theta(z)c.
\]

For a prescribed sampled multiplier \(\mu\), form paired library matrices and
solve

\[
\left(\Theta_+-\mu\Theta_-\right)c\approx0
\]

subject to a nontrivial normalization, phase/sign convention, sparsity, and
held-out validation. The trivial vector \(c=0\) must be excluded explicitly.
Complex conjugate modes must be stored as a paired real representation or as a
versioned complex coefficient contract.

This stage may reuse SINDy-style polynomial libraries and thresholding, but it
must be named accurately in the public API: it is a sparse sampled-Koopman or
sparse EDMD eigenfunction fit, not a continuous-time derivative-based SINDy
fit.

When \(\mu\) is not prescribed, EDMD may first fit a finite sampled operator
and extract candidate multipliers. A prescribed value must be rejected or
reported as approximate when the data and dictionary do not support a
nontrivial low-residual eigenfunction.

#### First sampled-eigenfunction feasibility result

The first isolated experiment used the five-orbit, three-period causal Hankel
basis and a fifth-degree polynomial dictionary. Dense EDMD was used before
sparsification so that representability failure would not be confused with a
thresholding failure. State variables were withheld from fitting.

On the original single unforced orbit, the leading oscillatory pair recovered
an implied period of (1.440570\tau), compared with the echo-measured period
of (1.440531\tau). Its (0.01\tau) sampled-step residual was approximately
0.203%, its four-(\tau) semigroup error was approximately 0.358%, and its
phase regression was effectively linear. This establishes that a polynomial
of the causal latent coordinates can represent a useful *orbit-local*
eigenfunction.

The analogous five-orbit fit did not produce a global oscillatory
eigenfunction. It selected an implied period near (1.433\tau), while the
training periods ranged from approximately (1.020\tau) to (1.532\tau).
The fitted polynomial's observed phase speed followed each orbit's own period
rather than the single EDMD multiplier. The mean training sampled-step
residual was about 2.53%, but held-out four-(\tau) errors ranged from about
0.79 to 1.78. The small-step residual was therefore a near-identity illusion,
not evidence of a global semigroup eigenfunction. This is the expected
spectral obstruction for a nonlinear conservative oscillator whose frequency
depends on action: fixed nonzero point-spectrum eigenvalues exist on an
individual invariant orbit, but generally not over a continuum of
non-isochronous orbits.

A separate (mu=1) fit sought a nonconstant action-like invariant by
maximizing between-trajectory polynomial variation while minimizing
within-trajectory variation. This produced a useful approximate invariant,
but the two-coordinate version generalized unevenly to held-out orbits. A
controlled rank comparison retained the same data, three-period parent SVD,
normalization, dictionary degree, and validation set. Rank three improved the
held-out median invariant variation from about 0.061 to 0.055 and reduced the
held-out error of a quadratic (\omega(I)) diagnostic from about 16.3% to
11.9%. Rank four expanded the polynomial library to 126 terms and worsened
both results, exposing a conditioning/overfitting cost rather than providing a
monotone repair.

Consequently, the pooled oscillatory EDMD pair is rejected for a global KROM.
The viable next paths are deliberately distinct:

1. fit and validate an orbit-local oscillatory pair around the intended target
   action for local KLQR tracking;
2. retain a data-driven invariant and learn a state-dependent phase law
   (\omega(I)), recognizing that this is a nonlinear or parameter-varying
   lifted model rather than a fixed-eigenvalue KROM;
3. improve invariant regularization and trajectory coverage before deciding
   whether rank three should be promoted;
4. sparsify only after the dense representability and long-horizon tests pass.

### 7. Input-aware statistical control model

The autonomous relation is insufficient for control. Controlled data should
fit a residual such as

\[
\psi_{k+1}
=\mu\psi_k+widehat B_{\psi}(z_k,d^u_k)u_k+\varepsilon_k.
\]

For multiple eigenfunctions,

\[
\Psi_{k+1}
=A\Psi_k+B(z_k,d^u_k)u_k+\varepsilon_k.
\]

Candidate model classes should include, in increasing order of complexity:

1. constant input matrix;
2. affine state-dependent input matrix;
3. bilinear lift/input interaction;
4. sparse nonlinear input dictionary;
5. local or uncertainty-aware regression when a global map fails.

Input design must provide persistent excitation while respecting safe
amplitude, duty-cycle, and cavitation limits. A single perfectly periodic
driver is inadequate for separating autonomous evolution from actuation over a
broad region. Multi-sines, chirps, randomized phases, or bounded designed
perturbations should be considered.

Closed-loop identification requires additional care because \(u_k\) is
correlated with the measured state. Early control-map training should prefer
open-loop designed excitation or use an identification method that explicitly
accounts for feedback correlation.

### 8. Target construction

Targets must pass through exactly the same causal observable:

\[
y^{\mathrm{ref}}
\longrightarrow d^{y,\mathrm{ref}}_k
\longrightarrow \xi^{\mathrm{ref}}_k
\longrightarrow \psi^{\mathrm{ref}}_k.
\]

No separate target-only centering, basis, or polynomial coefficients are
allowed. A prerecorded target provides its prehistory, but the plant still
requires an online warm-up policy.

### 9. Discrete controller

The initial controller should operate on the same sample clock as the learned
model. A nominal linear lift may use

\[
\Psi_{k+1}=A\Psi_k+B_k u_k
\]

with a discrete Riccati equation. If \(B_k\) depends on the current chart, the
existing receding or gain-scheduled KLQR logic must be reformulated and
validated in sampled time rather than borrowing continuous-time gains.

Controller evaluation must separate:

- eigenfunction-target error;
- echo-target error;
- actuator effort and saturation;
- long-term boundedness;
- recovery after perturbations;
- sensitivity to timing and measurement noise;
- performance on unseen runs and operating conditions.

Recovering a known sine-wave driver is a valuable system-identification test,
but successful output tracking is the primary criterion. Exact driver recovery
is not expected when multiple inputs can produce observationally equivalent
echo histories.

### 10. Online evaluator

At each sample, the production evaluator will:

1. acquire and quality-check \(y_k\);
2. update the fixed-length echo ring buffer;
3. update the known-input history buffer;
4. center and scale the echo window;
5. project it through the saved Hankel/POD basis;
6. evaluate the saved polynomial dictionary and coefficients;
7. evaluate the statistical control map;
8. compute and saturate the next input;
9. archive timing, uncertainty, and safety diagnostics.

No SVD, sparse regression, derivative estimation, or physical-state
interpolation occurs online. The main costs are a small matrix-vector product,
a polynomial evaluation, and the controller update.

The startup policy must be explicit. Options include collecting \(q-1\)
uncontrolled samples, initializing from a recorded equilibrium prehistory, or
running an estimator with increasing confidence until the buffer is full.
Control must remain disabled or conservative while the history is invalid.

## Versioned artifact contracts

The first implementation should define independent, removable contracts rather
than embedding method-specific data in UI state.

### `dynsyslab.echo-delay-observable/v1`

- sample interval and timestamp policy;
- echo channel and units;
- delay count and lag spacing;
- input-history count and alignment;
- training mean and scale;
- projection basis and score convention;
- singular values and rank-selection evidence;
- training-run references and data-role provenance;
- causal-window and warm-up policy;
- supported operating envelope.

### `dynsyslab.sampled-koopman-eigenfunction/v1`

- sampled multiplier \(\mu\);
- optional continuous \(\lambda\) and branch provenance;
- dictionary definition and ordering;
- coefficient vector and normalization convention;
- complex-pair representation;
- sparsity settings;
- training and validation residuals;
- parent delay-observable identity;
- supported operating envelope.

### `dynsyslab.echo-control-map/v1`

- autonomous eigenfunction identity;
- input channel, units, delay, and calibration identity;
- regressor definition;
- constant, affine, bilinear, or nonlinear model class;
- uncertainty and extrapolation diagnostics;
- excitation-data provenance;
- held-out controlled-rollout evidence.

### `dynsyslab.echo-only-controller/v1`

- observable, eigenfunction, and control-map identities;
- discrete cost matrices and Riccati solution;
- sample time;
- reference construction;
- saturation, slew-rate, and fail-safe limits;
- controller-engagement and buffer-warm-up policy;
- validation and release status.

## Architecture placement

The implementation must follow the existing native-core boundaries.

- The independently compilable native core owns causal delay transforms,
  normalization, fixed projection, sparse sampled-Koopman fitting, and model
  evaluation.
- Control owns the discrete controller and statistical input-map interfaces.
- Workflow composes archive-backed jobs and publishes versioned results.
- The compiled CLI exposes validation, fitting, replay, and evaluation gates
  before any UI integration.
- Solver and Viewer consume public workflow contracts. They do not reimplement
  Hankel construction, regression, or controller mathematics.
- Archive records refer to immutable source runs instead of copying their
  arrays into every derived KROM or control attempt.
- Python scripts may support bounded research experiments, but production
  behavior must not depend on a parallel Python implementation.

Generated solves, plots, report screenshots, and Viewer payloads remain local
under the policy in [Local scientific and media artifacts](../LOCAL_ARTIFACT_POINTERS.md).
Git should contain job recipes, model contracts, checksums or archive IDs,
tests, and narrative conclusions rather than bulk numerical output.

## Validation ladder

### Gate A: analytic synthetic signals

- Pure sinusoid with a quarter-period delay must recover a quadrature chart.
- The causal evaluator must reproduce offline coordinates sample for sample.
- Mutating future samples must not alter the current latent coordinate.
- Saved and reloaded observable/eigenfunction packages must be bitwise or
  tolerance-equivalent.

### Gate B: known hidden-state benchmarks, learner blinded

- Use oscillators with known states and equations only as independent truth.
- Supply the learner only output and, where applicable, input histories.
- Compare one-latent and multi-latent closure, prediction, and eigenfunction
  residuals.
- Diagnose when a two-coordinate chart folds even though a larger delay state
  succeeds.
- Preserve the low-rank failure and compare the controlled rank ladder under
  identical data, targets, horizons, and scoring rules.

### Gate C: RPE simulation in model-blind mode

- Generate data with the existing trusted solver.
- Prevent the new learner and controller from reading \(r\), \(v\), the RPE
  equations, analytic echo expressions, or observable derivatives.
- Use hidden states only in a separate audit report.
- Reproduce mild and nonlinear sine targets with causal online replay.
- Compare against current echo/echo-dot and state-aware baselines.

### Gate D: controlled identification

- Collect deliberately excited echo/input pairs.
- Fit autonomous and input components without leakage from validation runs.
- Compare constant, affine, bilinear, and nonlinear statistical control maps.
- Reject maps that perform well only on their training driver.

### Gate E: closed-loop simulation

- Evaluate long horizons, transient recovery, unseen phases, amplitude shifts,
  timing jitter, measurement noise, saturation, and dropout.
- Require bounded behavior and documented failure detection.
- Test the same saved package through CLI replay and the online execution path.

### Gate F: laboratory shadow mode

- Run the estimator and proposed controller without applying its output.
- Compare proposed commands with safe operator-approved excitation.
- Verify timing, calibrated input alignment, drift, extrapolation alarms, and
  computational latency.

### Gate G: bounded laboratory engagement

- Begin within a conservative operating envelope.
- Retain independent safety limits and emergency disengagement.
- Expand amplitude or duration only after held-out and repeated-run evidence.

## Required tests

The code phase should add focused tests for:

- causal ring-buffer ordering and wraparound;
- absence of cross-run windows;
- training-only normalization;
- raw-POD versus HAVOK-normalized score conventions;
- offline/online projection parity;
- polynomial dictionary ordering and serialization;
- prescribed-multiplier eigenfunction residuals;
- complex conjugate pairs and normalization;
- input/output sample alignment and actuator delay;
- constant, bilinear, and nonlinear control regressions;
- discrete Riccati and controller saturation behavior;
- missing samples, clipping, jitter, and warm-up;
- archive lineage and derived-record deduplication;
- future-sample mutation as an explicit causality test.

## Principal scientific risks

### Insufficient observability

Distinct physical conditions may generate indistinguishable finite echo
histories. No algorithm can recover a state component absent from the measured
signal. Mitigation is to enlarge the history, retain more latent coordinates,
add a physically available sensor, or reduce the claimed operating envelope.

### Non-Markov compression

The full delay vector may predict well while \((y,\xi_1)\) does not. The
two-variable representation is therefore a candidate to test, not a required
architecture. Rank must follow validation evidence.

### Input/output confounding

If input records are absent, poorly aligned, or insufficiently exciting, the
autonomous and controlled effects are not separately identifiable. The
experiment must record commands and preferably delivered acoustic pressure.

### Apparatus dynamics

The learned system may include amplifier, transducer, propagation, sensor, and
filter dynamics. That can be acceptable for control of a fixed apparatus, but
it limits transfer. Apparatus identity and calibration must be part of model
provenance.

### Nonstationarity

Equilibrium radius, gas content, temperature, position, and boundary conditions
may drift. Models require repeated-run validation, drift detection, and either
retraining or scheduling rather than silent extrapolation.

### Hidden-state safety

Echo tracking alone does not prove that an unmeasured radius remains inside a
safe physical bound. Laboratory deployment must retain independent protective
limits or additional safety sensing until that relationship is established.

### Unsupported eigenvalues

A user-specified \(\lambda\) does not guarantee that the associated multiplier
has a nontrivial eigenfunction in the selected dictionary and data regime. The
fit must report failure or approximation error rather than returning an
arbitrary normalized near-null vector.

## Phased implementation plan for `echo-only`

### Phase 0: contracts and experiment harness

- Seal the schemas above.
- Add a model-blind experiment flag that prevents hidden-state access.
- Add CLI validation for causal windows and input alignment.
- Establish lightweight pointer-only evidence records.

### Phase 1: causal echo history

- Implement fixed-size echo and input ring buffers in the native core.
- Implement training-only centering/scaling and saved Hankel projections.
- Prove offline/streaming parity with synthetic tests.

### Phase 2: sampled eigenfunctions

- Add sparse sampled-Koopman fitting for fixed and discovered multipliers.
- Reuse canonical polynomial dictionaries without importing continuous
  derivative assumptions.
- Add residual, sparsity, conditioning, and held-out-run diagnostics.

### Phase 3: sufficiency study

- Run the explicit retained-coordinate ladder \(r=2\), then \(r=3\), then
  \(r=4\) as required by evidence; retain every result, including failures.
- Compare the ranked latent models with the uncompressed delay state.
- Select by held-out closure and prediction, not singular values alone.
- Document negative results, the failure that triggered each promotion, the
  improvement after promotion, newly introduced complications, and
  operating-envelope boundaries.

### Phase 4: statistical input response

- Fit constant, affine, bilinear, and sparse nonlinear input maps.
- Add persistent-excitation and closed-loop-bias diagnostics.
- Retain uncertainty and extrapolation indicators online.

### Phase 5: discrete control

- Implement discrete KLQR/LQR around the sampled lift.
- Add target-history construction, warm-up, saturation, and fail-safe logic.
- Re-run mild and nonlinear sine-driver recovery studies in model-blind mode.

### Phase 6: application integration

- Add archive jobs and CLI gates.
- Expose compact configuration and diagnostic views in Solver/Post-processor.
- Add Viewer panels for delay basis, closure, eigenfunction residual, input-map
  uncertainty, and controlled tracking.

### Phase 7: laboratory transition

- Define acquisition and calibration adapters.
- Run shadow-mode inference.
- Validate repeated-run robustness and timing.
- Proceed to bounded actuation only under independent safety controls.

Each phase should be independently reviewable and removable. A later phase may
not hide a failed earlier gate.

## Success criteria

The research program succeeds only if all of the following are demonstrated:

1. The current latent coordinate and eigenfunctions are exactly reproducible
   from current and past measurements; future data are unnecessary.
2. A compact delay representation gives low, stable held-out transition
   dispersion under matched inputs.
3. Sampled eigenfunction residuals remain low across complete held-out runs,
   not only overlapping windows from the training trajectory.
4. The statistical input model generalizes across driver phase, waveform, and
   amplitude within its declared envelope.
5. The discrete controller tracks echo targets for long horizons and recovers
   after bounded perturbations without hidden access to \(r\) or \(v\).
6. Timing, saturation, noise, dropout, and extrapolation failures are detected
   and produce safe behavior.
7. The same saved contracts drive offline replay, CLI evaluation, and online
   execution.
8. Every rank promotion retains a controlled low-rank baseline and documents
   both the repaired failure and any new complication caused by the larger
   representation.

Failure of the one-latent model is not failure of the overall program. It is a
scientific result directing the implementation toward additional latent
coordinates or another measured observable.

## Open decisions for the first experiment

- Which echo channel and preprocessing most closely match the eventual
  laboratory acquisition path?
- Should the first coordinate be ordinary Hankel POD, residualized POD, or a
  supervised simulation-only candidate selected for state reconstruction but
  evaluated from echo online?
- What lag spacing and total window duration best balance observability,
  latency, and drift?
- Should the first model prescribe the current successful KROM eigenvalues or
  discover sampled multipliers from echo history and match them afterward?
- Is a bilinear input map sufficient, or does the actuator/bubble combination
  require a nonlinear input dictionary?
- Which independent sensor or conservative limit will protect against hidden
  unsafe radius excursions in laboratory tests?

These questions should be answered by bounded experiments on this branch,
with no change to the validated mainline pipeline until the evidence supports
promotion.

## Source spine

- F. Takens, “Detecting Strange Attractors in Turbulence,” in *Dynamical
  Systems and Turbulence, Warwick 1980*, 1981.
  <https://doi.org/10.1007/BFb0091924>
- J. Stark, D. S. Broomhead, M. E. Davies, and J. Huke, “Takens Embedding
  Theorems for Forced and Stochastic Systems,” *Nonlinear Analysis*, 1997.
  <https://doi.org/10.1016/S0362-546X(96)00149-6>
- M. O. Williams, I. G. Kevrekidis, and C. W. Rowley, “A Data-Driven
  Approximation of the Koopman Operator: Extending Dynamic Mode
  Decomposition,” *Journal of Nonlinear Science*, 2015.
  <https://doi.org/10.1007/s00332-015-9258-5>
- J. L. Proctor, S. L. Brunton, and J. N. Kutz, “Dynamic Mode Decomposition
  with Control,” *SIAM Journal on Applied Dynamical Systems*, 2016.
  <https://doi.org/10.1137/15M1013857>
- H. Arbabi and I. Mezić, “Ergodic Theory, Dynamic Mode Decomposition, and
  Computation of Spectral Properties of the Koopman Operator,” *SIAM Journal
  on Applied Dynamical Systems*, 2017.
  <https://doi.org/10.1137/17M1125236>
- S. L. Brunton, B. W. Brunton, J. L. Proctor, E. Kaiser, and J. N. Kutz,
  “Chaos as an Intermittently Forced Linear System,” *Nature Communications*,
  2017. <https://doi.org/10.1038/s41467-017-00030-8>

## Decision summary

The future DynSysLab echo-only pathway will be model-free in its learned and
online dependencies, causal by construction, sampled-data first, and honest
about what output-only measurements can establish. It will learn latent
coordinates and eigenfunctions from echo histories, learn actuation response
from synchronized echo/input experiments, and reject rather than conceal
insufficient observability or closure. Physical models and hidden simulation
states may remain valuable independent truth sources during development, but
they will not be runtime inputs or training requirements for the final
laboratory-capable controller.
