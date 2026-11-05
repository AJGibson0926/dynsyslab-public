# Korda--Mezic Duffing control-matrix candidate study

Date: 2026-08-27  
Status: completed comparison experiment; **not promoted to the standard workflow**

## Purpose and safety boundary

This study implements removable candidate tooling for three related tasks:

1. projecting a declared observable onto a supplied full-complex Koopman basis,
   (Y \approx C\Psi+c_0);
2. identifying a state-independent lifted input matrix from randomly forced
   controlled trajectories using the Korda--Mezic multistep output regression;
3. identifying and evaluating the bilinear approximation
   (B_\Psi(x)\approx B_1\Psi(x)+B_2), including pointwise matrix-norm error
   fields over a bounded phase-space rectangle.

The candidate modules are isolated under `rom/candidates/` and are not imported
by `rom/__init__.py`. They do not attach a result to a KROM, do not alter a
KROM's default control representation, and do not change the accepted exact

\[
B_\Psi(x)=\nabla\Psi(h(x))\,J_h(x)\,B(x)
\]

path. Existing exact-
`B_psi` Duffing and RPE cases remain untouched. The separate Korda--Mezic and
MATLAB cleanup trees were treated as read-only sources; only a new wrapper in
DynSysLab was written.

The candidate implementation and subsequent revisions are now part of the
canonical Git history. The temporary out-of-tree rollback snapshots used
during the experiment have been retired; ordinary branches and commits are the
sole source-control and rollback mechanism.

## Source pathway reproduced

The parity source is the companion `Duffing.m` program for Korda and Mezić's
optimal Koopman-eigenfunction work. Its controlled equation simplifies to

\[
\dot x=v,\qquad
\dot v=-0.5v+x-4x^3+0.5u_{\mathrm{source}}.
\]

DynSysLab's canonical equivalent uses unit input gain,

\[
\dot v=-0.5v+x-4x^3+u_{\mathrm{canonical}},\qquad
u_{\mathrm{canonical}}=0.5u_{\mathrm{source}}.
\]

The MATLAB wrapper reproduces, without modifying the source:

- random seed `2141442`;
- 100 free trajectories;
- normalized two-component Gaussian initial conditions on the unit circle;
- 800 free samples at `dt=0.01`;
- the supplied 20 precomputed discrete eigenvalues;
- source boundary-value regression and `scatteredInterpolant` eigenfunctions;
- 100 controlled trajectories, 200 intervals each;
- independent source inputs uniformly distributed over `[-1,1]` and held over
  each interval;
- rejection when a state magnitude exceeds 3 or becomes NaN;
- the exact published multistep regression and open-loop prediction test.

All 100 controlled trajectories passed the source rejection test.

### Source RK detail retained only for parity

The published script labels its fixed step as RK4, but computes the fourth
stage from `x + dt*k1`, not `x + dt*k3`. The wrapper intentionally reproduces
that line. DynSysLab's standard controlled-data generator remains classical
RK4. With identical initial states and identical physical forcing, the
classical-RK4 trajectories differ from the MATLAB source-compatible data by:

| State | all-sample RMSE | maximum absolute difference | final-time RMSE |
|---|---:|---:|---:|
| position | `1.6546e-5` | `4.1055e-5` | `1.7467e-5` |
| velocity | `2.4022e-5` | `6.3144e-5` | `2.3911e-5` |

This difference is small but measurable. It is an expected integrator-form
difference, not an input-scaling bug, and is not a reason to copy the source bug
into the standard solver.

## Observable reconstruction C

The source uses the explicit matrix

\[
C_{\mathrm{official}}=
\begin{bmatrix}
\mathbf 1_{1\times10}&0\\
0&\mathbf 1_{1\times10}
\end{bmatrix}.
\]

The wrapper and Python candidate also compute the projection

\[
C_{\mathrm{projection}}=X\,\operatorname{pinv}(\Psi).
\]

Python reproduces MATLAB's projected matrix to `1.00e-14` relative Frobenius
error. The projected and official matrices themselves differ by `0.16428`
relative Frobenius norm, but their free-data reconstruction errors are close:

| Reconstruction | position RMSE | velocity RMSE |
|---|---:|---:|
| official block-sum C | `0.0027456` | `0.0084222` |
| projected C | `0.0026049` | `0.0084154` |

The candidate projection is therefore numerically correct, but this experiment
does not justify automatically replacing the official C. More generally, C is
basis-normalization dependent and must be validated after H-DMD-to-SINDy
eigenfunction rescaling. The existing standard KROM observable-reconstruction
gates remain authoritative.

## Constant lifted input matrix

For every trajectory and prediction horizon `j`, the source constructs

\[
x_j-C A_d^j z_0
=\sum_{k=0}^{j-1}
u_k\,C A_d^{j-k-1}B_d.
\]

Stacking the 100 trajectories and 200 horizons produces a complex design
matrix of shape `40000 x 20`. The Python candidate rebuilds this matrix rather
than trusting the MATLAB export and uses the MATLAB-equivalent SVD cutoff,
`max(Q.shape)*eps = 8.8818e-12`.

### Cross-language parity

| Quantity | relative Frobenius difference | maximum absolute entry difference |
|---|---:|---:|
| regression design Q | `4.1475e-16` | `2.6407e-14` |
| regression target | `1.9092e-14` | `1.7117e-15` |
| discrete source-input B | `2.1443e-9` | `3.0564e-7` |
| continuous source-input B | `2.1423e-9` | `3.0590e-5` |
| discrete canonical-input B versus `2*B_source` | `2.1443e-9` | `6.1127e-7` |
| continuous canonical-input B versus `2*B_source` | `2.1423e-9` | `6.1180e-5` |

The matrices agree. The remaining `~2e-9` relative difference is numerical:
the design has condition number `1.5814e7`, and MATLAB and NumPy call different
LAPACK pseudoinverse paths after receiving matrices that already agree nearly
to machine precision.

The discrete source matrix has Frobenius norm `259.7613`. It is converted to a
continuous matrix using the exact zero-order-hold identity

\[
B_d=\left(\int_0^{\Delta t}e^{A_c s}\,ds\right)B_c.
\]

The integral has condition number `1.0063`, full rank 20, and conversion
residual `3.05e-16`. The resulting source-input continuous norm is
`26014.6697`; canonical input coordinates double it to `52029.3394`.

Despite the regression target's relative residual `0.62015`, reconstructed
state behavior is useful:

| Diagnostic | position RMSE | velocity RMSE |
|---|---:|---:|
| one-step output | `0.0029741` | `0.0083138` |
| 200-step rollout output | `0.0070006` | `0.0145883` |
| published two-second open-loop test | `0.0107665` | `0.0148337` |

The large target-relative residual and good state reconstruction can coexist
because this is an output-weighted regression through C, not a pointwise fit of
the lifted vector field.

## Pointwise exact-B_Psi error

The free-training rectangle is

\[
x\in[-0.9997722,0.9999837],\qquad
v\in[-0.9998768,0.9996671].
\]

The wrapper evaluates the interpolated eigenfunctions on a `101 x 101` grid.
Centered spatial differences supply the velocity derivative, and the exact
source-coordinate input map is

\[
B_\Psi(x)=\frac{1}{2}\frac{\partial\Psi}{\partial v}(x).
\]

At every point the candidate computes

\[
e_{\mathrm{abs}}(x)=\|B_\Psi(x)-\widehat B_\Psi(x)\|_F,
\qquad
e_{\mathrm{rel}}(x)=
\frac{e_{\mathrm{abs}}(x)}
{\max(\|B_\Psi(x)\|_F,\epsilon_{\mathrm{rel}})}.
\]

The training convex hull is retained as a mask and boundary inside the full
requested rectangle. The exact-map norm on the rectangle has median `2.9426`,
mean `5.3825`, 95th percentile `19.3822`, and maximum `143.0125`.

### Constant output-regression matrix

Within the training convex hull, the constant matrix has absolute error near
`2.6015e4` everywhere and relative error:

- median `10198.83`;
- mean `14055.31`;
- 95th percentile `38764.06`;
- maximum `111514.36`.

This does **not** contradict the output-prediction result. Korda's constant
matrix minimizes multi-step error after projection through C. It can place very
large, mutually canceling components in weakly observed lifted directions.
It was never optimized to approximate local `B_psi` in Frobenius norm.

### Bilinear fits

Two candidate bilinear fits were retained because they answer different
questions.

1. **Continuous dynamics residual**, matching the cleanup experiment:

   \[
   \dot\Psi-A_c\Psi=(B_1\Psi+B_2)u.
   \]

   The design is `20000 x 21`, full rank, with condition number `18880.5`.
   Its relative lifted-residual error is `0.99552`.

2. **Direct control-map projection**:

   \[
   \min_{B_1,B_2}\sum_i
   \|B_\Psi(x_i)-(B_1\Psi(x_i)+B_2)\|_F^2.
   \]

   This uses the same 20,000 random-forced interval-start states and exact
   spatial-gradient `B_psi` samples. Its design condition number is `19079.7`
   and relative fit residual is `0.97889`.

Both are poor for this particular interpolated Korda basis. Inside the training
hull, direct-map bilinear error has median relative error `1.1844`, mean
`93.8900`, 95th percentile `472.393`, and maximum `5287.30`. The dynamics fit
has median `1.1822`, mean `101.724`, 95th percentile `508.705`, and maximum
`5680.42`. Heavy tails near high-gradient interpolation regions dominate the
means.

The most likely causes are structural:

- the source eigenfunctions are piecewise-linear scattered interpolants fitted
  to a discrete-time trajectory relation;
- they are not constrained to satisfy the continuous generator PDE
  `grad(Psi)*f = lambda*Psi` between samples;
- gradients of piecewise-linear interpolants change discontinuously at
  triangulation boundaries;
- an affine function of the 20 lift coordinates is insufficient to represent
  those spatial gradient changes over the whole domain.

This is also why direct control-map projection improves the tails only modestly
over dynamics-residual regression. The dominant failure is representation
capacity and eigenfunction regularity, not merely noisy time differentiation.

## Interpretation limits

The Frobenius error is meaningful only within one fixed eigenfunction basis and
normalization. Rescaling an eigenfunction rescales the corresponding row of
`B_psi`, the constant/bilinear matrix, and C. Therefore:

- these heat maps compare representations within this exact exported Korda
  basis;
- their absolute magnitudes must not be compared directly with another KROM
  unless eigenfunction normalization is aligned;
- output prediction and control performance remain necessary complementary
  metrics.

## Candidate implementation inventory

- `rom/candidates/observable_projection.py`: full-complex affine observable
  projection with deterministic or independent validation diagnostics.
- `rom/candidates/control_identification.py`: immutable pre-lifted datasets,
  multistep discrete constant-B fit, exact zero-order-hold conversion,
  continuous bilinear dynamics fit, and direct bilinear control-map projection.
- `rom/candidates/error_fields.py`: absolute and relative pointwise Frobenius
  fields for constant and bilinear representations, with training-hull metrics.
- `scripts/matlab_parity/export_korda_mezic_duffing_control_candidate.m`:
  non-mutating source-compatible MATLAB exporter.
- `scripts/run_korda_duffing_control_candidate_comparison.py`: independent
  Python reproduction, cross-language comparison, arrays, JSON, and heat map.
- `tester/tests/test_candidate_koopman_control_identification.py`: synthetic
  exact-recovery and error-field tests.

Generated experiment bundles live under ignored
`archives/experiments/korda_duffing_candidate/`.

## Recommendation

Do not promote either approximation from this interpolated Korda basis into a
default controller. Retain exact `B_psi` as the standard representation.

The next controlled experiment should use the standard Duffing H-DMD-to-SINDy
pipeline and its smooth polynomial SINDy eigenfunctions. For that case:

1. compute `B_psi` analytically from polynomial gradients;
2. fit constant output-regression, bilinear dynamics, and direct control-map
   projection candidates on the same controlled dataset;
3. hold out complete trajectories rather than random points;
4. cross-validate ridge strength and report design rank/condition number;
5. compare matrix error, lifted derivative error, reconstructed-state rollout,
   and actual closed-loop control separately;
6. promote nothing unless it passes explicit gates in its declared input,
   state, and time coordinates.

## Reversal and removal

The experiment can be removed without touching accepted frameworks by deleting
only:

- `rom/candidates/`;
- `scripts/matlab_parity/export_korda_mezic_duffing_control_candidate.m`;
- `scripts/run_korda_duffing_control_candidate_comparison.py`;
- `tester/tests/test_candidate_koopman_control_identification.py`;
- this report and generated ignored experiment artifacts;
- the corresponding README paragraph and `/archives/experiments/` ignore line.

No rollback of `rom/control_regression.py`, `rom/koopman.py`, controllers,
solvers, Viewer, accepted archives, or successful KLQR examples is required.
