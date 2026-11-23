# Toda lattice

The ODE registry key is `toda_lattice`, with `periodic` and `open` presets.
Parameters are `sites` (2–32, default 3), `alpha` (positive, default 1),
`damping` (nonnegative, default 0), and `input_gain` (default 1).
Mass is one and native coordinates are dimensionless,
\(x=(q_1,\ldots,q_N,p_1,\ldots,p_N)\).

\[
 H=\frac12\sum_jp_j^2+\alpha^2\sum_{\text{bonds }j}e^{q_j-q_{j+1}},
 \qquad \dot q_j=p_j.
\]
Each bond exerts opposite forces on its endpoints. Periodic boundaries include
the \(N\to1\) bond; open boundaries omit it and have no wall forces. Damping
adds \(-\gamma p_j\); the scalar external input acts on \(p_1\).
At zero damping/input, energy and total momentum are conserved.

Observables include every state and its linked zero-input derivative,
`total_energy`, `energy_dot`, `input_power_gain`, `total_momentum`,
`momentum_dot`, `center_of_mass`, `center_velocity`, `relative_energy`,
each `bond_j` and `relative_q_j`, and Flaschka coordinates `flaschka_a_j`
and `flaschka_b_j`. Indices in names begin at one.

## Analytic scope

The registered evaluator covers two sites. The periodic relative coordinate
\(R=q_1-q_2\) obeys \(\ddot R=-4\alpha^2\sinh R\), so it reuses the cosh
oscillator's physical-time elliptic evaluator. Center of mass advances at
\(P/2\), where \(P=p_1+p_2\). A full-state periodic torus is supplied only
when \(P=0\) and the relative motion is not an equilibrium.

For the open pair, let \(E=\dot R^2/4+\alpha^2e^R\). Its hyperbolic solution is
\[
 R(t)=\log(E/\alpha^2)-2\log\cosh(\sqrt E(t-t_*)).
\]
The implementation determines the phase from the initial state and evaluates
log-cosh stably. This is a scattering trajectory, not a periodic Fourier torus.

General periodic Toda has established Birkhoff/action-angle theory, but the
general theta-function evaluator is not implemented here. For more than two
sites, the analytic request returns an explicit inapplicability reason.
The ODE itself supports all registered site counts.

## Lax API

`core/Toda.hpp` exposes `todaLaxPair(q,p,alpha,periodic)`.
The same calculation is available through the `mathematics` CLI:

```json
{"schema":"dynsyslab.math-query/v1","operation":"inverse_scattering",
 "action":"toda_lax","positions":[0.3,-0.3,0],
 "momenta":[0.1,-0.2,0.1],"alpha":1,"periodic":true}
```

It returns row-major matrices \(L,P_L\) with \(\dot L=[P_L,L]\) and
\(a_j=\alpha e^{(q_j-q_{j+1})/2}/2,\ b_j=-p_j/2\).
The periodic matrix uses Bloch multiplier \(+1\); its eigenvalues alone do
not describe the complete periodic spectral curve. Two periodic sites have
two bonds, so their matrix entries are added, not overwritten.
Lax eigenvalues are conserved spatial spectral data, not Koopman eigenvalues.

Tests check the analytic solutions in both time directions, ODE/Jacobian
consistency, conservation, center-of-mass qualification and the Lax commutator.

Sources: [Teschl's Toda exposition](https://diana.mat.univie.ac.at/~gerald/ftp/book-jac/toda.html)
and [Henrici–Kappeler, global periodic Toda Birkhoff coordinates](https://arxiv.org/abs/0802.4252).
