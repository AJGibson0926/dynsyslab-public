# Colbrook spectral-method reference inventory and integration proposal

Checked 11 November 2026. This is the research and acquisition inventory in the miscellaneous computes side-car. All 14 accessible repositories have now been downloaded at the recorded commits. The native pipeline implementation, API, boundaries, and validation are described in [Spectral KROM](../../architecture/SPECTRAL_KROM.md).

The [manifest](repositories.json) records 14 publicly accessible repositories at exact commits, one published repository URL currently returning 404, entry points, licenses, paper links, external datasets, and bibliographic corrections. The scope includes Colbrook's relevant methods and coauthored implementations. It excludes unrelated mathematics repositories and unaffiliated reimplementations. A separate public repository for the 2025 kernelized Rigged DMD extension was not verified.

## Meaning of the computed objects

For a time-step map T, the Koopman operator acts on an observable by K g = g composed with T. It is linear in g even when T is nonlinear.

- An eigenfunction phi is a function of state with phi(T(x)) = lambda phi(x). Along a trajectory its value evolves as lambda^n phi(x_0). A continuous-time eigenfunction has exponential evolution with a generator eigenvalue. Discrete-to-continuous logarithms have frequency aliases.
- A Koopman mode is an output reconstruction coefficient or vector multiplying an eigenfunction. It is different from the eigenfunction itself. SVD/POD vectors describe variance or energy in the supplied data and need not be Koopman eigenfunctions.
- A pseudoeigenfunction has a small norm of (K-lambda)phi relative to phi. A residual landscape searches over lambda and dictionary coefficients. For nonnormal operators, small residual does not imply that lambda is comparably close to the spectrum.
- For an invertible measure-preserving system on L2 of an invariant measure, K is unitary. Its spectral projections E define the scalar spectral measure nu_g(B) = ||E(B)g||^2. This assigns an observable's squared norm to frequency sets B. The total mass is ||g||^2, or its variance when centered. Atoms represent projections onto eigenspaces. Continuous mass can have a density or be singular continuous.
- The correlations satisfy <g,K^k g> = integral exp(i k theta) d nu_g(theta), with inner product conjugate-linear in its first argument. This connects spectral measures to power spectra. A frequency histogram of finite-matrix eigenvalues is not automatically the underlying spectral measure.
- Generalized eigenfunctions may be distributions outside L2. Rigged DMD constructs finite-width wave packets approximating them. They allow a spectral description when ordinary L2 eigenfunctions do not span the observables.

These objects depend on the chosen function space and its measure or kernel. Changing from L2 to an RKHS changes the operator's spectral problem. See the [DMD review](https://arxiv.org/abs/2312.00137), [Rigged DMD](https://arxiv.org/abs/2405.00782), and [SpecRKHS](https://arxiv.org/abs/2506.15782).

## Method inventory

| Method and source | What is approximated and how | Proposed use |
|---|---|---|
| [Residual DMD / ResDMD](https://arxiv.org/abs/2111.14889), with the 2023 robust-and-verified applications | Dictionary eigenfunctions and residual-minimizing pseudoeigenfunctions. Uses an additional successor Gram matrix to estimate full-operator residuals rather than only residuals of a compressed square matrix. Searches can reveal spectral regions missed by EDMD eigenvalues. | First priority. Attach residuals and pseudospectral grids to ordinary, nonlinear-dictionary, and Hankel DMD. |
| [Kernel and dual ResDMD](https://arxiv.org/abs/2403.05891) | Kernel dictionaries allow nonlinear eigenfunction estimates. The fewer-snapshots-than-features regime needs special care because interpolated residuals can vanish spuriously. The dual least-squares formulation provides a different residual construction without the original training/quadrature split. | Support high-dimensional histories and nonlinear features. Record precisely which residual and operator interpretation were used. |
| [Stochastic ResDMD](https://arxiv.org/abs/2308.10697) | The stochastic Koopman operator advances conditional expectations. An observed squared residual contains both expectation-model error and conditional variance. Repeated successor samples can separate them; variance-pseudospectra measure statistical coherence. | Future stochastic ensemble analysis. Distinguish process noise from sensor noise and solver error. |
| [mpEDMD](https://doi.org/10.1137/22M1521407) | A Gram-weighted orthogonal Procrustes problem produces a finite measure-preserving approximation. Its eigendecomposition approximates spectral projections and measures under the convergence hypotheses. Finite unit-circle eigenvalues can discretize continuous spectrum. | Optional structured fitter using the same dictionaries and snapshot pairs. Require an explicit invariant-measure interpretation. |
| [Spectral measures via ResDMD](https://arxiv.org/abs/2111.14889) | Resolvent smoothing and filtered Fourier series of autocorrelation moments estimate observable-specific measures. Resolution, kernel order, sample error, and dictionary error remain separate. | Dedicated spectral-measure postprocessor and Viewer frequency plots. Useful even without a state-reconstruction model. |
| [Rigged DMD](https://arxiv.org/abs/2405.00782) | Combines a convergent unitary approximation with rational resolvent kernels, producing wave packets for generalized eigenfunctions and associated reconstruction modes. Covers discrete and continuous spectral components. | Downstream of mpEDMD and a shared resolvent/smoothing layer. Use a separate wave-packet result contract. |
| [Kernelized Rigged DMD](https://doi.org/10.1063/5.0271497), Hales, Colbrook and Jiang, 2025 | Introduces nonlinear kernel dictionaries into the generalized-eigenfunction and spectral-measure computation. | Combine frozen kernel features with the Rigged DMD stage. No separate implementation repository verified. |
| [Multiplicative DMD / MultDMD](https://arxiv.org/abs/2405.05334) | Preserves K(fg) = (Kf)(Kg) using a partition/indicator dictionary and a constrained finite transition map. Eigenvectors become piecewise-constant observable functions. | Alternative partition dictionary and structured fitter. Benchmark noise sensitivity, resolution, and predictions separately. |
| [Hermitian DMD convergence](https://arxiv.org/abs/2401.03192) | Preserves self-adjoint structure and obtains convergence of spectral measures under stated assumptions. Skew-Hermitian generators can be related by multiplication by i. | Specialized backend for justified self-adjoint problems. Generic dissipative systems do not satisfy that assumption. |
| [Weighted Birkhoff methods](https://arxiv.org/abs/2511.17772) | Tapered ergodic averages improve the moment estimates underlying wtDMD, wtEDMD, wtSINDy, spectral measures, and diffusion forecasting. Faster convergence depends on dynamics and regularity. | Shared explicit sample-weight policy. Apply independently within trajectories, with normalization and effective sample size recorded. |
| [SpecRKHS](https://arxiv.org/abs/2506.15782) | Kernel-space methods exploit the exact action of the adjoint on kernel sections to compute residuals, spectra, pseudospectra, and prediction bounds. Spectral measures require the relevant unitary/self-adjoint setting. | Dedicated kernel-space backend. Store Koopman versus adjoint/Perron-Frobenius orientation explicitly. |
| [TrustKoopman / Principal Angle Decomposition](https://arxiv.org/abs/2603.15091) | Measures angles between a feature subspace and its Koopman image. PAD chooses dynamically useful directions; further routines estimate multi-step mode-decomposition errors. | First priority for diagnosing the dense-Hankel rank studies. Compare against POD at equal ranks using one frozen dictionary and held-out data. |
| [EDMD + Personalized PageRank](https://arxiv.org/abs/2605.24666) | Uses the EDMD matrix as a graph to select approximately invariant sub-dictionaries around chosen seed observables. Results concern detection and leakage under stated assumptions. | Dictionary selection before fitting, seeded by echo/rate or requested outputs. Does not itself prove global nonlinear coordinate closure. |
| [Deep Embedded Multiplicative DMD](https://arxiv.org/abs/2606.05131) | Learns latent coordinates and partitions while alternating exact multiplicative-operator updates with differentiable clustering. | Later optional research adapter. Freeze/export the learned encoder and partition before native evaluation. Published code URL currently unavailable. |
| [Convergent algorithms and adversarial systems](https://www.nature.com/articles/s41467-026-74220-8) | Gives spectral search and spectral-type algorithms together with limits on learning across broad system classes. Some tasks require multiple ordered limits. | Benchmark suite and guarantee metadata. More data or modes alone is not a universal convergence certificate. |
| [SpecSolve](https://arxiv.org/abs/2006.01766) | General resolvent-based spectral-measure computations for self-adjoint operators. This is supporting spectral numerics, not a DMD fitter. | Reusable rational smoothing and shifted-solve facilities behind spectral postprocessors. |

The [DMD-Multiverse review code](https://github.com/MColbrook/DMD-Multiverse) also supplies exact, forward-backward, total-least-squares, randomized, and controlled DMD examples. These are useful baselines reviewed by Colbrook; the review does not claim he originated every algorithm.

The [ENSO application](https://arxiv.org/abs/2606.09369) demonstrates residual-minimum mode selection and harmonic organization with kernel EDMD. Treat it as a worked application of the residual framework.

## Shared numerical contract

Use M paired observations and N dictionary functions. Store evaluated dictionaries as P_X and P_Y with shape M by N. For our existing real Hankel archives these are the transposes of the stored feature-by-window matrices, possibly followed by a frozen dictionary transform. Let W contain explicit nonnegative sample weights. Then

G = P_X* W P_X,
A = P_X* W P_Y,
L = P_Y* W P_Y.

An EDMD eigenvector c defines phi(x) = sum_j c_j psi_j(x), with K_N = G^dagger A. Its empirical squared residual is

res(lambda,c)^2 =
c* [L - lambda A* - conjugate(lambda) A + |lambda|^2 G] c / (c* G c).

The asterisk denotes conjugate transpose. Null/ill-conditioned Gram directions must be handled before interpreting this quotient. Stable QR/SVD implementations are preferable to explicitly forming badly conditioned normal equations. Weighted sufficient statistics can be accumulated in batches.

Minimizing this quotient over c at a grid of lambda values estimates a residual landscape. In a nonnormal problem its interpretation is pseudospectral. Sampling or quadrature errors must be controlled before it becomes an operator-level certificate. Learned dictionaries must retain a frozen out-of-sample evaluator, and independent validation or the appropriate dual/RKHS construction must prevent interpolation from masquerading as zero operator error.

RPE flow derivatives can supply extra channels through the observable relationships already implemented. Generator methods would instead evaluate derivative dictionaries and need their own clock, domain, and theorem contract. Discrete-time proofs should not be silently applied to that different problem.

## Pipeline placement

These are proposed additions, not newly implemented stages.

1. **Measurement and Hankel preparation.** Keep the existing observable catalog, derivative links, channel scales, causal alignment, and one-analysis-tick X/Y pairing in [the Hankel contract](../../architecture/OBSERVABLE_HANKEL.md). Do not put fitting or spectral assumptions inside the matrix builder. Non-Hankel paired observations must also be supported.
2. **Frozen dictionary and weighting.** Add a reusable native Koopman input contract for paired feature evaluations, sample weights, trajectory membership, fit/validation roles, kernel/partition parameters, and frozen evaluation maps. POD, PAD, nonlinear feature dictionaries, and PageRank selection should compose here. Preserve original trajectory boundaries.
3. **Native Koopman numerics.** Small independent modules should own weighted EDMD statistics, residual/pseudospectral calculations, structured fits, invariance diagnostics, and resolvent/spectral-measure kernels. Use the existing compiled linear-algebra boundary. MATLAB/Python reference packages remain external comparison implementations.
4. **Neutral workflow.** Proposed methods include resdmd, koopman_invariance, mpedmd, multdmd, spectral_measure, and rigged_dmd. Resolve data through existing postprocessing/Archiver gateways. Produce immutable derived records with matrix/dictionary hashes, clock, metric, residual semantics, kernel, quadrature weights, regularization, rank, smoothing order/width, validation source, and assumptions.
5. **ROM.** Extend KROM only for compatible finite evolution/readout models. Spectral measures and generalized-eigenfunction packets deserve typed analysis results rather than being forced into a list of autonomous exponential modes.
6. **Viewer and workbench.** Compose those public methods. Viewer shows residual-colored eigenvalues, pseudospectral contours, spectral mass/density with smoothing resolution, principal-angle spectra, and error-versus-horizon plots. No algorithms belong in plotting handlers.
7. **Benchmarks.** Preserve the exact linear Duffing reference. Add rotation/periodic examples, mixing maps with continuous spectrum, a nonnormal linear example, a conservative pendulum, and then the same five-ring RPE train/held-out study. Keep observable set, physical window, sampled pairs, weights, ranks, selection budget, and error metric fixed across comparisons.

## RPE recommendation and guarantee limits

Start with ResDMD and principal-angle diagnostics on the exact shared Hankel pairs and the frozen POD rank prefixes. This tests whether high-energy retained directions remain within the feature span after evolution. Then compare PAD-selected directions and nonlinear dictionaries, always reporting measured-window error separately from free-rollout error.

This addresses finite **linear Koopman** closure. It does not decide whether a finite coordinate map admits some **nonlinear** autonomous law. An injective coordinate change can have nonlinear closed dynamics without its coordinate functions spanning a finite invariant linear space.

For a damped free RPE ringdown, an invariant measure supported at the final equilibrium cannot encode the transient oscillations. Forcing unitary structure on a transient sampling distribution would remove the desired decay interpretation. Measure-preserving methods may instead suit sustained invariant-measure dynamics, a conservative comparison system, or an appropriately extended driven system. Do not equate volume dissipation with the absence of any invariant measure.

Finite residuals and observed rank-sweep convergence do not guarantee global recovery. State the function space, norm, sampling assumptions, finite-data uncertainty, and limiting process. The adversarial results are about what can be guaranteed across specified classes, not a proof that this particular RPE problem is impossible.

## Download preparation

From the DynSysLab source root, print the plan.

~~~powershell
# (the fetch helper is in Git history; the pinned commits below identify each upstream)
~~~

The acquisition has been run with --download for all accessible pins. To reproduce it, use the same command with --download. Select one or more repositories with --repository owner/name. The default destination is the already ignored build/reference_sources/colbrook directory. --destination accepts an alternative research directory.

The helper fetches only the pinned revision, checks HEAD, skips unchanged matching checkouts, and refuses to reset or overwrite existing modified/mismatched checkouts. It runs no examples or installers and does not fetch external datasets or submodules. If a fetch is interrupted, the partial directory is retained for inspection.

Most implementations are MATLAB. PageRank has Python sources and notebooks; ENSO is supplied as notebooks. Some examples need Chebfun or separately hosted large datasets. Verified repository tree contents total about 225 MiB before Git overhead; external datasets are additional.

TrustKoopman, Adversarial-Dynamical-Systems, and ENSO-Fundamental-Modes have no license file in the inspected commit. Keep them as reference evidence until reuse terms are established. Other root licenses are MIT or BSD-2-Clause, with bundled dependencies requiring their own attribution review. DeepMDMD has no downloadable pin because its URL returned 404.

Two upstream citation errors are corrected in the manifest: the weighted-methods README points to an unrelated arXiv identifier, and the TrustKoopman BibTeX identifier differs from its correct paper link.
