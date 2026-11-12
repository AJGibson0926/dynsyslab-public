# Critical Koopman reference: exact unstable linear Duffing system

**Priority: critically important, explicitly designated by the user on
2026-11-05. Preserve this example for future Koopman analysis.**

Reference ID: `koopman-linear-duffing-exact`.

This is the exact, position-observed benchmark used to verify the new
KLQE-to-KLQR feedback connection. It provides independently checkable answers
for eigenfunctions, spectra, input maps, reconstruction, observability and
controllability. Preserve the baseline separately from future variations.

## Origin and purpose

The parameter choice was constructed for this project using the existing
DynSysLab Duffing implementation. It was not taken from an external paper or
discovered by fitting trajectory data.

DynSysLab uses the convention

$$
\dot x=v,\qquad \dot v=\alpha x-\beta x^3-\delta v+g u.
$$

The benchmark specifies

$$
\alpha=2,\quad \beta=0,\quad \delta=1,\quad g=1,
$$

so its dimensionless equations are

$$
\dot x=v,\qquad \dot v=2x-v+u,\qquad y=x.
$$

It is a linear special case of Duffing with an unstable saddle equilibrium at
(0,0). The measured quantity is position; velocity is unmeasured. There is no
echo observable or history transform in this baseline.

The system already exists in the library: select `duffing` with parameters
`alpha: 2, beta: 0, delta: 1, input_gain: 1`. No separate system implementation
or special controller is needed. The original jobs omit input_gain because its
resolved default is 1.

## Exact eigenfunctions: derivation and direct verification

With state s=(x,v)^T, the autonomous system is sdot=A s, where

$$
A=\begin{bmatrix}0&1\\2&-1\end{bmatrix},\qquad
B=\begin{bmatrix}0\\1\end{bmatrix}.
$$

The characteristic polynomial is

$$
\det(\lambda I-A)=\lambda^2+\lambda-2
 =(\lambda-1)(\lambda+2).
$$

For a linear system, a left eigenvector satisfying w^T A=lambda w^T defines
the Koopman eigenfunction psi(s)=w^T s, because

$$
\nabla\psi\cdot As=w^T As=\lambda w^T s=\lambda\psi(s).
$$

The two left eigenvectors [2,1] and [-1,1] therefore give

$$
\psi_1(x,v)=2x+v,\qquad \lambda_1=1,
$$

$$
\psi_2(x,v)=-x+v,\qquad \lambda_2=-2.
$$

Direct substitution with u=0 independently verifies both identities:

$$
\dot\psi_1=2v+(2x-v)=2x+v=\psi_1,
$$

$$
\dot\psi_2=-v+(2x-v)=2x-2v=-2\psi_2.
$$

Thus every autonomous trajectory satisfies

$$
\psi_1(t)=e^t\psi_1(0),\qquad
\psi_2(t)=e^{-2t}\psi_2(0).
$$

These are exact global polynomial identities. They do not depend on training
coverage, numerical differentiation, SVD truncation, approximate closure or
empirical error estimates. They describe Koopman eigenfunctions on the smooth
observable algebra; claims about a particular L2 operator space would require
specifying that space and its measure.

## Exact coordinate change, reconstruction and control map

The retained Koopman coordinates satisfy

$$
\Psi=T s,\qquad T=\begin{bmatrix}2&1\\-1&1\end{bmatrix},
\qquad \det T=3.
$$

This is a globally invertible change of coordinates:

$$
s=C_x\Psi,\qquad
C_x=T^{-1}=\frac13\begin{bmatrix}1&-1\\1&2\end{bmatrix}.
$$

In particular, x=(psi1-psi2)/3 and v=(psi1+2psi2)/3. The measurement map is

$$
y=C_y\Psi,\qquad C_y=\frac13\begin{bmatrix}1&-1\end{bmatrix}.
$$

Including the input in the chain rule gives

$$
\dot\psi_1=\psi_1+u,\qquad
\dot\psi_2=-2\psi_2+u,
$$

so the controlled lifted model closes exactly:

$$
\dot\Psi=\Lambda\Psi+B_\psi u,\qquad
\Lambda=\begin{bmatrix}1&0\\0&-2\end{bmatrix},\qquad
B_\psi=TB=\begin{bmatrix}1\\1\end{bmatrix}.
$$

The eigenfunction statement applies to the autonomous dynamics. Under input,
the extra Bpsi u term is part of the controlled evolution.

**The benchmark Bpsi was supplied analytically, not regressed.** Its stored
representation name `fit` is an implementation label, not evidence of a
fitting step. The polynomial coefficients, eigenvalues, Cx and Cy were also
supplied exactly through the standard KROM classes.

## Exact observability and controllability

For measurements y=x, the physical-state observability matrix is

$$
\mathcal O_x=\begin{bmatrix}C\\CA\end{bmatrix}
=\begin{bmatrix}1&0\\0&1\end{bmatrix}.
$$

The lifted observability and controllability matrices are

$$
\mathcal O_\Psi=
\begin{bmatrix}C_y\\C_y\Lambda\end{bmatrix}
=\frac13\begin{bmatrix}1&-1\\1&2\end{bmatrix},
\qquad \det\mathcal O_\Psi=\frac13,
$$

$$
\mathcal C_\Psi=
\begin{bmatrix}B_\psi&\Lambda B_\psi\end{bmatrix}
=\begin{bmatrix}1&1\\1&-2\end{bmatrix},
\qquad \det\mathcal C_\Psi=-3.
$$

Both have full rank. Position measurements over time can identify both modes;
the observer does not require a measured velocity or numerical position
derivative. Although instantaneous y alone does not identify s, the dynamical
model and repeated measurements make the pair observable.

## Continuous and discrete time

For sampling period h, the autonomous discrete eigenvalues are

$$
\mu_1=e^h,\qquad \mu_2=e^{-2h}.
$$

For input held constant over each interval, the exact sampled model is

$$
\Psi_{k+1}=\begin{bmatrix}e^h&0\\0&e^{-2h}\end{bmatrix}\Psi_k
+\begin{bmatrix}e^h-1\\(1-e^{-2h})/2\end{bmatrix}u_k.
$$

The sampled input map is not the continuous Bpsi, and it is not exactly h Bpsi.
This gives a useful reference for checking time-domain conversions. The saved
successful benchmark uses continuous-gain KLQE and KLQR with sampled,
zero-order-held signals; it is not a discrete Kalman measurement-update design.

## Saved successful run

| Setting | Value |
|---|---|
| Initial physical states | (+0.1,0) and (-0.1,0) |
| Initial estimated Koopman coordinates | (0,0), independent of plant truth |
| Measured quantity | y=x only |
| Target | Psi=0, equivalent to (x,v)=(0,0) |
| KLQR Q and R | Q=I2, R=0.1 |
| LQE process and measurement intensities | I2 and 0.001 |
| Injected measurement noise | None |
| Input bound | abs(u) <= 0.8 |
| Active-control interval | [0.1,12] |
| Measurement/controller interval | 0.005 |
| Archived plant-output interval | 0.02 |
| Final time | 12 |
| Plant solver | RK45; relative tolerance 1e-10, absolute tolerance 1e-12 |

The observer uses the actually applied input, including clipping and gating.
It passes estimated Psi directly to the usual KLQR law. Cx reconstructs the
physical-state estimate separately.

| Quantity at time 12 | Both symmetric initial conditions |
|---|---:|
| Controlled state norm | 1.1976463038e-8 |
| State-estimation error norm | 2.1852136878e-10 |
| Uncontrolled state norm | 15344.6688912 |
| Peak absolute applied input | 0.8 |

The uncontrolled comparator uses identical physics, initial conditions and
integration tolerances with control disabled. The norm curves for the two
initial conditions coincide by symmetry.

## Preserved evidence and reproduction

The definitive native constructor and integration checks are in
[koopman_feedback_tests.cpp](../../workflow/native/tests/koopman_feedback_tests.cpp).
The [example folder](../../examples/native/case-studies/koopman-feedback-benchmark-20261105/README.md)
contains the standard job files and Viewer plot helper.

The original KROM, four trajectories, two hierarchy manifests, original job/
result records, figures and metrics are preserved outside the disposable build
tree under:

`archives/experiments/koopman-feedback-benchmark-20261105/`

The [local asset manifest](../../examples/native/case-studies/koopman-feedback-benchmark-20261105/local-assets.json)
records original locations, preserved locations, byte sizes and SHA-256 hashes.
The copied scientific records are byte-identical; embedded provenance is not
rewritten. Use the manifest's preserved paths if an original build path has
subsequently been cleaned. Numerical archives and media remain local and
excluded from Git; source, documentation, job files and pointers are lightweight.

- [Full Viewer figure](../../archives/experiments/koopman-feedback-benchmark-20261105/figures/full.png)
- [Transient Viewer figure](../../archives/experiments/koopman-feedback-benchmark-20261105/figures/transient.png)
- [Measured result metrics](../../archives/experiments/koopman-feedback-benchmark-20261105/figures/metrics.json)
- [Controlled solve job](../../examples/native/case-studies/koopman-feedback-benchmark-20261105/controlled-job.json)
- [Uncontrolled solve job](../../examples/native/case-studies/koopman-feedback-benchmark-20261105/uncontrolled-job.json)
- [Observer/controller API and timing contract](../architecture/KOOPMAN_OUTPUT_FEEDBACK.md)

From the checkout root, the compiled CLI can plan/run either saved job using a
new archive output location. The controlled job refers to the preserved local
KROM. If the local artifacts are unavailable in another clone, build and run
`dynsyslab_koopman_feedback_tests.exe <new-output-directory>` to generate a
fresh exact KROM and benchmark through the same native code.

## How to use this reference later

Use the exact baseline to separate errors in coordinate evaluation, fitting,
input regression, mode selection and observer/controller assembly. Its exact
answers make it useful for controlled comparisons with learned models.

The polynomial observable algebra also supplies additional exact autonomous
eigenfunctions:

$$
\psi_1^a\psi_2^b,\qquad \lambda_{a,b}=a-2b,
\qquad a,b\in\{0,1,2,\ldots\}.
$$

This includes the constant function and nonconstant invariant psi1^2 psi2.
It also illustrates that an eigenvalue does not uniquely specify an
eigenfunction: psi1 and psi1^3 psi2 both have eigenvalue 1. These extra
eigenfunctions were not used in the successful two-mode run.

Future experiments may learn these quantities from data, change the measured
observable, add noise, or introduce nonzero cubic stiffness. Keep such variants
separate and clearly labeled. This baseline demonstrates successful output
feedback under exact finite-dimensional closure; nonlinear RPE or history-based
models require their own accuracy and control validation.
