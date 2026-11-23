# Inverse-scattering building blocks

`dynsyslab_math` exposes `InverseScattering.hpp`; the compiled `mathematics`
query and thin Python `dynsyslab_client.mathematics.query` expose the same
algorithms. These are spatial spectral and scattering tools. A Schrödinger/Lax
eigenfunction is not automatically a Koopman eigenfunction on the space of
potentials. All boundary and normalization conventions below are explicit.

## Conventions and scope

The spatial operator is \(L=-\partial_x^2+V\), with KdV convention
\(V_t-6VV_x+V_{xxx}=0\). A right-incident scattering state has
\(e^{-ikx}+R(k)e^{ikx}\) at \(+\infty\) and \(T(k)e^{-ikx}\) at
\(-\infty\). A bound state is normalized in spatial \(L^2\), with right tail
\(c_j e^{-\kappa_jx}\). The API stores \(\log(c_j^2)\), not \(c_j\).

`rightScattering` uses exact constant-cell transfer matrices, with zero
potential outside the supplied cells. `hillMonodromy` uses their periodic
extension instead. Its discriminant is the **full trace**; the real spectral
band condition is \(|\Delta|\le2\). Some papers use half the trace and threshold
one. Transfer products can become ill-conditioned in long forbidden regions;
overflow is rejected and determinant/flux defects are reported. They do not
bound discretization error or prove a spectral gap.

`dirichletNegativeSpectrum` uses centered finite differences at internal nodes
of a finite box and Sturm bisection for negative eigenvalues. These are finite
matrix eigenvalues. Whole-line bound states require box and grid convergence;
the routine does not compute their norming constants.

For reflectionless data, the inverse is the finite-rank Marchenko determinant:
\[
 V(x)=-2\partial_x^2\log\det(I+A(x)),\qquad
 A_{ij}(x)=\frac{c_i c_j e^{-(\kappa_i+\kappa_j)x}}{\kappa_i+\kappa_j}.
\]
The derivatives are evaluated analytically using Cholesky solves, not numerical
second differences. The single-soliton case has a stable explicit expression.
The multisoliton exponential chart rejects excessive growth rather than
silently producing nonfinite data. There is no general reflection-containing
Marchenko solver or periodic finite-gap inverse in this API yet.

Evolution is
\[
 R(k,t)=R(k,0)e^{8ik^3t},\qquad
 \log c_j^2(t)=\log c_j^2(0)+8\kappa_j^3t.
\]
Thus the rate for \(c_j\) itself is \(4\kappa_j^3\). These conventions are
important when comparing IST coordinates with Koopman rates.

## Query contract

`operation: "inverse_scattering"`, schema `dynsyslab.math-query/v1`:

| `action` | Required fields |
|---|---|
| `right_scattering` | `potential` cell values, `left`, `cell_width`, positive `wavenumber` |
| `hill_monodromy` | `potential` cell values, `cell_width`, spectral `energy` |
| `dirichlet_negative_spectrum` | Internal-node `potential`, box `length` |
| `evolve_reflection` | `reflection:[real,imag]`, `wavenumber`, `elapsed` |
| `evolve_norming` | `kappa`, `log_norming_squared`, optional `elapsed` |
| `reflectionless_potential` | `kappa`, `log_norming_squared`, `positions`, optional `elapsed` |
| `toda_lax` | `positions`, `momenta`, optional positive `alpha` and boolean `periodic`; finite Toda Lax matrices, see [Toda API](TODA_LATTICE.md) |

```python
from dynsyslab_client.mathematics import query
result = query({"operation":"inverse_scattering",
                "action":"reflectionless_potential",
                "kappa":[1], "log_norming_squared":[0.6931471805599453],
                "positions":[-1,0,1], "elapsed":0})
```

This gives the negative well \(-2\operatorname{sech}^2x\). Tests check free
transmission, reflectionless convergence, the finite-box bound eigenvalue,
norming evolution, Hill monodromy and a two-soliton KdV residual. These are
numerical checks rather than interval certificates.

## Algorithm references

- Gardner, Greene, Kruskal and Miura (1967): https://doi.org/10.1103/PhysRevLett.19.1095
- Drazin and Johnson, *Solitons: An Introduction*: https://doi.org/10.1017/CBO9781139172059
- Parker and Page, *Koopman Analysis of Isolated Fronts and Solitons*:
  https://doi.org/10.1137/19M1305033
- Parker and Valva, periodic finite-gap KdV (a distinct boundary problem):
  https://doi.org/10.1063/5.0137088
