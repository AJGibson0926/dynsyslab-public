# Exact Lie algebra and Koopman module analysis

## Central research question

For a control-affine system on a specified state manifold M,

`x_dot = f(x) + sum_k u_k g_k(x)`,

which finite-dimensional function spaces V simultaneously:

1. admit a basis of drift Koopman eigenfunctions;
2. are invariant under every input Lie derivative L_gk;
3. define a smooth embedding Psi of the physical state;
4. preserve useful controllability/stabilizability descriptions on Psi(M);
5. remain useful under quantified approximation, differentiation, and readout errors?

For constant input matrix B, full partial-derivative closure is sufficient but
stronger than closure along B's columns. The global question is on all of M,
not inferred from a bounded collection of samples. Local, regional, and global
statements must retain their domains.

The exact lifted law is generally bilinear:

`z_dot = A_0 z + sum_k u_k A_k z`.

An eigenfunction basis makes A_0 diagonal when such a basis exists. A finite
invariant space can instead have Jordan blocks. Constant input maps give the
special LTI case. Matrix dynamics outside the image Psi(M) are not extra
physical states.

## Ownership and public gates

`dynsyslab_lie` is an independent native analysis library. It depends only on
foundation JSON and private Boost rational arithmetic. It does not own ODE
integration, fit models, write archives, or depend on control, Viewer, Python,
Eigen, or the workflow library. Solver and Post-processor clients may compose it.
It is intentionally not part of the solver's time-stepping loop.

- C++: `dynsyslab/core/LieAlgebra.hpp`, namespace `dynsyslab::core::lie`.
- Workflow: `evaluateLieAlgebraQuery`.
- Stable C ABI: optional additive ABI-v4 symbol `dsl_workflow_lie_algebra_json`;
  output uses the existing result-handle lifetime and error contract.
- Compiled CLI: `dynsyslab-cli lie-algebra query.json`.
- Optional Python: `NativeWorkflow.lie_algebra(query)` or
  `python -m dynsyslab_client lie-algebra query.json`. Both delegate to C++.

No desktop dialog is added. These are read-only scientific queries, like the
existing analytic and observable-closure gates.

## Implemented computations

| Operation | Result |
|---|---|
| Polynomial differentiation/product | Canonical sparse exact polynomial |
| Lie derivative L_X h | Exact polynomial directional derivative |
| Vector-field bracket | `[X,Y]=DY X-DX Y` |
| Poisson bracket | Any supplied constant antisymmetric tensor |
| Flow Jacobi bracket | `a L_f b-b L_f a` |
| Generated Lie algebra | Independent polynomial fields and bracket-word provenance |
| Structure constants | Exact coordinates of every basis bracket |
| Derived/lower-central series | Dimensions, solvability, nilpotence |
| Center | Basis in Lie-basis coordinates |
| Generated ideal | Smallest ambient-adjoint-invariant subspace containing selected fields |
| Killing form | `trace(ad_i ad_j)` |
| Pointwise distribution rank | Exact rank at rational points; scope explicitly limited |
| Common invariant module | Bounded constructive completion under selected derivations |
| Operator matrices/residuals | Exact invariance identities or counterexamples |
| Linear state readout | Global polynomial identity `x=C Psi`, when it exists |
| Supplied-model defects | Exact `L_X Psi-A_X Psi` |
| Residual box bounds | Rigorous uniform and Lipschitz induced 1-norm bounds |
| Second-order infinite cascade | Analytic nonexistence certificate under checked hypotheses |
| Sinusoidal mechanical cascade | Infinite angular-harmonic certificate on an exact circle lift; damping-independent |

Sparse rational Gaussian elimination checks independence as functions, not
sampled values. The bracket search repeatedly applies generator adjoints to a
growing basis. If the complete span stops growing, Jacobi guarantees closure
under all brackets. The function-space search similarly applies all supplied
derivations. Its reduction coordinates are echelon coordinates, not a
least-squares or L2 projection.

Algorithms follow classical derivation, adjoint, and Lie-series definitions;
the finite representation criterion follows Katayama/Susuki [1]. The cascade
and deterministic coefficient bounds are project derivations documented in
the private case-study report, not copied implementations from those papers.

## JSON v1

The sinusoidal certificate recognizes
`F=(c*v,-s*v,A*s+B*c+C+D*v)` in coordinates `(s,c,v)`,
with a nonzero sine/cosine restoring coefficient, generated `d_v`, and
circle-preserving supplied generators. Its harmonic witnesses stay independent
on `s^2+c^2=1`, so it does not confuse ambient polynomial growth with a physical
obstruction. For `s=sin(x_1),c=cos(x_1),v=x_2`, the physical angular
partial is `R=(c,-s,0)`; do not substitute ambient `d_s,d_c`.
See the pendulum proof and executable cases.
The general exact query API below remains unchanged.

Requests require `schema: dynsyslab.lie-algebra-query/v1`, `dimension`, and
`action`. Results use `dynsyslab.lie-algebra-result/v1`.
An exact polynomial is a list of terms, for example x_2-x_1^2:

```json
[
  {"coefficient":"1", "powers":[0,1]},
  {"coefficient":"-1", "powers":[2,0]}
]
```

The empty list is zero. Coefficients and point coordinates are integer or p/q
strings; floating JSON numbers are rejected there. Exponents are nonnegative
JSON integers. A vector field is a list of d such polynomials. Rationalizing
measured floating coefficients is an explicit caller decision, not evidence
that the fitted model is exact. These proof objects are distinct from the
canonical floating-point `PolynomialEigenfunctionSet` runtime bank.

Actions and required fields:

- `algebra`: `generators`; optional rational `points`, `ideal_generators`.
  Ideal seeds must belong to the certified finite ambient algebra. Ideal results
  include a basis and solvability/nilpotence diagnostics; they do not establish
  the separate eigenfunction hypotheses of an embedding theorem.
- `module`: `observables`, `operators`; optional `include_partials`, `extend`
  (default true). Partials are appended after the requested operators.
- `defect`: `observables`, `field`, flat row-major `matrix`.
- `box_bound`: `residuals`, positive rational `radii`.
- `bracket`: vector fields `x`, `y`.
- `poisson`: polynomials `a`, `b`, row-major d-by-d `tensor`.
- `flow_bracket`: `field`, `a`, `b`.
- `lie_derivative`: `field`, `observable`.
- `derivative`: `observable`, zero-based `axis`.
- `product`: `a`, `b`.

Optional `coordinate_names` and `provenance` are carried through. Supply the
physical coordinates, units, exact parameter choices and source equations.
Examples in `examples/native/lie-algebra` are analytic query specifications,
not generated solution archives or default physics settings.

Budgets are `maximum_basis` (default 32, cap 64), `maximum_operations`
(4096, cap 100000), and `maximum_degree` (24, cap 64). Arithmetic has additional
term, coefficient-size and exponent guards. The search is deliberately bounded.

Statuses mean:

- `closed`: exact identities certify closure of the returned span.
- `not_invariant`: an unextended supplied observable span has a nonzero defect.
- `limit_reached`: incomplete; neither existence nor nonexistence follows.
- `proven_infinite`: the explicit cascade theorem applies, not merely that the
  search kept growing. `basis_dimension` then counts discovered independent
  fields, not the dimension of the infinite algebra.

For modules, `matrices[k][i*r+j]` implements `L_Xk Psi_i=sum_j A_k[i,j] Psi_j`.
Matrices on coefficient columns are transposes: `T_k=A_k^T`. Therefore
`A_[X,Y]=A_Y A_X-A_X A_Y`, while `[T_X,T_Y]=T_[X,Y]`.
This convention prevents a silent sign error between column features and
operator representations.

`global_linear_state_readout_certified=false` means no linear readout in that
span; it does not rule out a nonlinear inverse. The module never labels an
arbitrary invariant basis an eigenfunction basis. Numerical eigendecomposition
and runtime KROM construction remain owned by existing numerical/KROM modules.

## Proof boundaries and control

The infinite certificate recognizes a generated field `Y2=x1^2*d_x2` and a
generator `F=x2*d_x1+(A(x1)*x2+C(x1))*d_x2`. It checks these polynomial
hypotheses exactly. For `H=[F,Y2]`,

`[H,x1^n*d_x2]=-(n+2)*x1^(n+1)*d_x2`, for every integer n>=2.

These fields are independent on every open subset of R2. A finite smooth
immersive common-invariant lift would give a faithful finite matrix
representation, a contradiction. This excludes more than polynomial lifts.
It does not exclude scalar observables that lose state information.

Pointwise rank is a distribution-rank computation. It is not a global
controllability test, especially in the presence of drift/input constraints.
Finite algebra dimension is necessary for a faithful finite lift, not sufficient.

For known polynomial residual R and |x_j|<=rho_j, the box utility bounds
each Jacobian entry by the sum of absolute monomial coefficients times the
appropriate powers of rho. The maximum column sum bounds the induced 1-norm.
If R(0)=0, the mean-value integral proves `||R(x)||_1 <= L ||x||_1`.
Otherwise it returns a uniform bound but refuses the vanishing proportional
claim. This is a deterministic bound on supplied equations, not a finite-data
probability bound and not a controller stability certificate.

## Papers and connection to existing control work

1. [Katayama and Susuki, Lie-algebraic embedding conditions (2026 preprint)](https://arxiv.org/abs/2611.20438):
   Theorem 1 motivates the obstruction test. Theorem 2 additionally requires a
   suitable nilpotent ideal and a finite-algebraic-multiplicity embedding;
   finding a nilpotent ideal alone does not establish its hypotheses.
2. [Goswami and Paley, global bilinearization (2017)](https://cdcl.umd.edu/papers/cdc17.pdf):
   common invariant spaces and bilinear control; relevant to the operator matrices.
3. [Zhang and Li, Koopman bilinearization (2022 preprint)](https://arxiv.org/abs/2211.07112):
   operator/Lie-group perspective. No blanket global controllability conclusion
   from sampled ranks is implemented.
4. [Heeg and Worthmann, limitations of LTI modeling (2026 preprint)](https://arxiv.org/abs/2604.25290):
   retain coordinate-map inclusion and directional-excitability assumptions.
   Nonlinear readouts and unactuated nonlinear subsystems require separate treatment.
5. [Strasser et al., feedback with stability guarantees, TAC 2025](https://arxiv.org/abs/2312.01441):
   the invariant-dictionary assumption is a direct target of the module checker.
6. [Strasser et al., SafEDMD, Automatica 2026](https://arxiv.org/abs/2402.03145):
   proportional residuals connect to robust feedback; Corollary 3.2 explicitly
   accounts for projection error. Our box tool addresses known polynomial defects.
7. [Bold, Grune, Schaller and Worthmann, EDMD MPC](https://arxiv.org/abs/2308.00296):
   cost controllability and practical stability require more than Lie rank.
8. [Strasser et al., control overview](https://arxiv.org/abs/2509.02839):
   places approximation error between representation and closed-loop guarantees.
9. [Strasser et al., input-output feedback (2026 preprint)](https://arxiv.org/abs/2606.07758):
   quantitative injectivity of observation windows is stronger than pointwise
   Jacobian rank. This matters for the previous echo/history-coordinate studies.
10. [Applied Koopmanism](https://arxiv.org/abs/1206.3164) and
    [Jacobi structures and Spencer operators](https://archive.mpim-bonn.mpg.de/id/eprint/428/1/preprint_2014_46.pdf):
    distinguish multiplication of eigenfunctions from Lie/Poisson/Jacobi operations.

Detailed reading notes and downloaded papers remain local under
`research/literature/lie-algebra`, alongside the earlier Koopman-control notes.
The case-study discussion is indexed in the private study records.

## Validation and current limits

The native test target covers exact arithmetic, Heisenberg and sl2 algebras,
centers/series, polynomial embeddings, cascade certificates, finite search
limits, input rejection, matrix defects, and proportional residual bounds.
The CLI and optional Python/C ABI are checked against the same saved requests.

This first module supports real rational polynomial equations. It does not
automatically prove general nonpolynomial identities, classify all Lie
algebras, find arbitrary global inverses, synthesize LMI/SOS controllers, or
infer a global certificate from finite trajectory data. The existing solver,
linear-control, eigenfunction and spectral services retain their responsibilities.

## Inverse-radius mechanics

The pressure-controlled RPE family now has an exact infinite-algebra recognizer,
using the auxiliary variable w=1/r and a bracket recurrence valid on r*w=1.
See the explicit proof
and `examples/native/lie-algebra/rpe-pressure-algebra.json`.
For saved approximate polynomial KROMs, use [KROM audits](KROM_AUDIT.md)
to distinguish fitting-coordinate derivative completion from physical drift/control closure.
