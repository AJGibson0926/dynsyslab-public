# Requirements for regressing B_psi and extending it to echo-history coordinates

Date: 2026-10-05
Project: DynSysLab, echo-only experiments
Purpose: Preserve the complete discussion for later reading. Implementation observations below describe the code inspected for this discussion; this report does not introduce code changes.

## Question

> What are the requirements, what kind of information is needed, in order to regress a B_psi matrix that is accurate and useful for control, in both the case of eigenfunctions which are functions of the dynamical state (here r and r-dot) or other observables with explicit expressions (e.g. echo and echo-dot)? How does this extend -- or does it? -- to the case of data-driven observables like latent variables derived from the echo history?

## Answer

Yes, B_ψ can be learned entirely from measurements, including for eigenfunctions built from echo-history latent variables. But there are two separate requirements: the data must identify the effect of actuation, and the chosen coordinates must contain enough information to predict that effect. Regression can address the first; it cannot repair a failure of the second.

For history coordinates, there is an additional issue: the natural object to learn is usually a **sampled input-response map**, rather than automatically the same instantaneous, continuous-time B_ψ we use with (r, ṙ).

Let me make those distinctions precise.

### What B_psi represents

Suppose the bubble dynamics are control-affine:

$$
\dot x=f(x)+G(x)u,\qquad x=(r,\dot r),
$$

and collect the eigenfunctions into

$$
z=\Psi(x)=
\begin{bmatrix}
\psi_1(x)&\cdots&\psi_p(x)
\end{bmatrix}^{T}.
$$

If they are exact eigenfunctions of the **unforced** dynamics, then

$$
\dot z=\Lambda z+B_\psi(x)u,
\qquad
B_\psi(x)=D_x\Psi(x)G(x).
$$

Thus B_ψ measures how actuation changes the eigenfunction values. It is generally a **matrix-valued function**, not a constant matrix.

A constant regression,

$$
B_\psi(x)\approx B_0,
$$

is a further approximation. A bilinear model allows more variation:

$$
B_\psi(x)\approx B_0+\sum_{\ell}z_\ell B_\ell,
$$

possibly including conjugate coordinates or, equivalently, real and imaginary components. Whether that approximation is adequate depends on the region and inputs used for control. The distinction between constant, state-dependent, and bilinear input maps is developed explicitly in [Koopman form of nonlinear systems with inputs](https://doi.org/10.1016/j.automatica.2024.111525).

### Eigenfunctions of r and r-dot

There are two ways to obtain regression targets.

**With known dynamics and actuation**, we can compute

$$
B_\psi(x_i)=D_x\Psi(x_i)G(x_i)
$$

at many sample points and regress those values.

This requires:

- the eigenfunction expression and its derivatives;
- the actuation vector field G(x), including the units and sign of the input;
- state samples covering the intended operating region.

The full unforced vector field f is not needed to evaluate this formula once Ψ and G are known. Forced trajectories are not mathematically necessary for this route because the control sensitivity is supplied analytically.

**Without equations of motion**, we instead record controlled trajectories, evaluate the eigenfunctions along them, and form the residual

$$
q(t)=\dot z(t)-\Lambda z(t).
$$

Then regress

$$
q(t)\approx B_\psi(x(t))u(t).
$$

For a constant map, this is simply a regression of q against u. For a state-dependent map with features φ_ℓ,

$$
B_\psi(x)\approx\sum_\ell \phi_\ell(x)B_\ell,
$$

the regression columns are the products

$$
\phi_\ell(x_i)u_{i,j}.
$$

This route needs synchronized state measurements, applied inputs, timestamps, and enough temporal resolution to estimate ż. It does **not** require analytic equations.

Numerical differentiation can be avoided with an integral formulation:

$$
z(t_b)-z(t_a)
-\int_{t_a}^{t_b}\Lambda z(t)\,dt
\approx
\int_{t_a}^{t_b}B_\psi(x(t))u(t)\,dt.
$$

That is often preferable for noisy laboratory data.

A crucial complication is eigenfunction error. If the unforced relation is only approximate,

$$
\dot z=\Lambda z+\delta(x)+B_\psi(x)u,
$$

then the measured residual contains both δ(x) and the control response. If we omit δ, regression can incorrectly attribute autonomous model error to actuation.

That is why reliable free-dynamics validation, unforced baseline data, and input variation at comparable states matter.

### Eigenfunctions of explicit observables, such as echo and echo-dot

Suppose

$$
y=h(x),\qquad z=\Psi(y).
$$

Then the exact control map is

$$
\boxed{
B_\psi(x)=D_y\Psi(h(x))\,D_xh(x)\,G(x).
}
$$

These are precisely the additional chain-rule factors you previously identified.

If we know h, its Jacobian, and G, we can compute this map offline and regress it. If we only have measurements of y and the applied input, we can use the same residual or integral regression on z(t) = Ψ(y(t)); the measured evolution already includes the chain-rule effects.

However, to express the result as B_ψ(y), we need the observable coordinates to distinguish states with different control responses. If

$$
h(x_a)=h(x_b)
\quad\text{but}\quad
B_\psi(x_a)\ne B_\psi(x_b),
$$

there is no single exact B_ψ(y). Regression will average incompatible responses. More samples may estimate that average more precisely without making it useful for control.

There is also a specific caveat for echo-dot: **a derivative measured during forcing need not equal the unforced derivative observable we used to train the eigenfunction.**

For example, when e = g(x),

$$
\frac{de}{dt}
=Dg(x)f(x)+Dg(x)G(x)u.
$$

A laboratory derivative therefore includes actuation. Acoustic measurements can also have direct dependence on forcing or actuator/sensor dynamics. We must establish whether our deployed coordinates are h(x), h(x,u), or a filtered measurement history. The formula above assumes h(x); it cannot silently be applied to all three.

### What changes for echo-history latent variables

Let the causal delay vector be

$$
d_k=
\begin{bmatrix}
e_k&e_{k-1}&\cdots&e_{k-L+1}
\end{bmatrix}^{T},
$$

and let the fixed training transform produce

$$
\xi_k=P(d_k),\qquad z_k=\Psi(\xi_k).
$$

Here P includes the archived Hankel basis and its normalization.

We can evaluate zₖ from measured echo history exactly as we have implemented. We can also record applied inputs and regress their effect. No expression for echo, r, ṙ, or their dynamics is needed.

The main question becomes whether

$$
z_{k+1}=F(z_k,u_k)
$$

is sufficiently well defined.

Under unforced evolution, a history can sometimes be treated as a function of the current physical state. Under forced evolution, the same current physical state may have different histories because different inputs were applied previously. Consequently, a basis that captures free oscillations well does not automatically produce a closed coordinate system for controlled motion.

There are two possible failures:

1. **Observation failure:** the echo history itself does not distinguish states relevant to future response.
2. **Compression failure:** the history distinguishes them, but keeping only ξ₁, ξ₂ discards that distinction.

We must test controlled response as well as free evolution to distinguish these.

For the measured-history setting, I would initially regress

$$
z_{k+1}-Mz_k
\approx
B_d(\xi_k)u_k,
\qquad M=\operatorname{diag}(\mu_j),
$$

with constant and bilinear versions as competing approximations.

Here B_d is the response over one specified sampling interval. Even if the physical continuous-time dynamics are affine in u, their finite-time response need not be exactly affine in u. We may need smaller sampling intervals, nonlinear input features, or additional coordinates.

And B_d is not generally B_ψ Δt. For a constant continuous-time lifted input map and a held input,

$$
B_d=
\left(\int_0^{\Delta t}e^{\Lambda s}\,ds\right)B_\psi.
$$

That conversion requires the assumed continuous-time lifted model to be valid. For a history-based predictor, using the learned sampled model directly is often better justified.

Measurement-based prediction and control with delay coordinates have established precedents, including the input-output formulation in [Korda and Mezić](https://arxiv.org/abs/1611.03537) and the sparse-measurement control experiments in [Arbabi, Korda, and Mezić](https://arxiv.org/abs/1804.05291). They establish feasibility, not a guarantee that our particular two-coordinate compression is sufficient.

### What the identification data must contain

The requirements are substantially the same across these representations:

| Requirement | Why it matters |
|---|---|
| Synchronized observations and **actual applied input** | Establishes which actuation produced each response; commanded input may differ from delivered input. |
| Multiple independent runs and operating conditions | Separates state dependence from one trajectory's particular input pattern. |
| Inputs that excite the regression features | Prevents rank deficiency and ambiguous coefficients. |
| Coverage of target motion and likely deviations | A map fit near equilibrium may fail during nonlinear oscillations or recovery. |
| Fixed coordinate definitions and scaling | The learned map is tied to the eigenfunction coefficients, basis, normalization, and sample interval. |
| Unforced baseline accuracy | Prevents autonomous prediction error from being mistaken for control authority. |
| Independent forced validation | Tests whether the learned actuation response generalizes. |

For our bubble experiment, a single sine-driven run is unlikely to identify a broadly useful map: input, phase, amplitude, and state are strongly correlated. Multiple initial conditions, phases, amplitudes, and suitably varied inputs provide much better evidence.

Our policy of **unforced eigenfunction training remains appropriate**. Learning actuation from measurements nevertheless requires a separate collection of forced identification runs. Unforced data alone contain no evidence about the effect of an unknown actuator.

Input history may also be needed:

$$
\eta_k=
\begin{bmatrix}
\xi_k\\u_{k-1}\\ \vdots\\u_{k-q}
\end{bmatrix}.
$$

The controller knows its own past inputs, so this does not require measuring hidden bubble states. It supplies information about how the observed history was produced.

### How we judge whether the map is useful

A small training residual is insufficient. I would require:

- acceptable regression rank and conditioning;
- correct held-out response to changes in input, including sign and phase;
- accurate multistep predictions under new input sequences;
- little systematic residual dependence on discarded coordinates or past inputs;
- adequate controllability or stabilizability of the relevant lifted dynamics;
- successful tracking over the long horizons we care about, with realistic input limits.

A particularly informative experiment is to find histories that map to nearly the same (ξ₁, ξ₂), apply comparable probes, and compare their subsequent lifted responses. If those responses differ substantially, no sufficiently accurate map depending only on those two coordinates can exist over that region. That would give us a concrete reason to upgrade to four coordinates or include input memory.

Finally, I checked our current implementation. The existing native regression code supports exact-map fitting, controlled-residual fitting, and multistep output fitting. However, its regression helper still evaluates samples through `model.evaluate(state)`. The new history model requires `evaluateCoordinates(xi)`.

So my earlier statement that the approximate-map route was available needs qualification: **storage and evaluation of a learned map are supported, but fitting it from history-coordinate measurement records still needs an explicit native adapter and the sampled regression formulation.** That is an implementation gap, alongside the scientific closure tests—not a reason to abandon the approach.

For our next experiment, the strongest approach is to freeze the existing unforced-trained basis and eigenfunctions, collect separate input/echo identification records, fit a sampled actuation map using two latent coordinates, and test its failures before moving to four.
