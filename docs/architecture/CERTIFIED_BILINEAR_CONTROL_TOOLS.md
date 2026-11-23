# Certified bilinear control tools

Native C++ tools for fitting a bilinear K-ROM, bounding its residual, and designing a state-feedback
controller with a checked Lyapunov certificate. They use Eigen only and hold no scientific data.

| Piece | Header | Role |
|---|---|---|
| Dense SDP solver and verifier | `core/DenseSdp.hpp` | LMI feasibility and linear objectives over PSD blocks |
| Bilinear model | `core/BilinearKrom.hpp` | Fit (free or frozen A), simulate, linearize, residuals, JSON |
| K-ROM bridge | `core/KoopmanBilinearBridge.hpp` | Stored bilinear map to a discrete target-centred model |
| Residual sector | `core/ResidualSector.hpp` | `‖r‖ ≤ c_x‖z‖ + c_u‖u‖` from residual data (joint LP envelope) |
| Supported sector | `core/SupportedResidualSector.hpp` | Two-stage estimate with both channels identified; required for a data-based certificate |
| Certified design | `core/CertifiedDesign.hpp` | Robust LMI design and independent check |
| JSON operations | `math/ControlCertificateQuery.hpp` | `sdp_feasibility`, `bilinear_fit`, `residual_sector`, `certified_design` |
| Workflow step | `bilinear_certificate` post-processing | Derived archive with model, sector, status, `K`, `X` |
| Deployment | KLQR setting `bilinear_certificate` | Uses the checked gain; KLQE keeps the same map |

## Model and design

The model is discrete, with lifted coordinates that vanish at the target:

```text
z+ = A z + B0 u + sum_i u_i B_i z + r,    z(target) = 0,  u(target) = 0
‖r‖ <= c_x ‖z‖ + c_u ‖u‖
```

Each data pair must hold `u` constant over its sampling interval. The fit is rejected when the design
is rank deficient, unless the caller allows it.

The design searches for `X > 0`, `L` and `mu > 0`, with `K = L X^-1` and **`u = K z`** (the opposite
sign to the LQR convention `u = -K x`). For every vertex `theta` of the input box, with
`Phi = A X + sum_i theta_i B_i X + B0 L`:

```text
[ lambda X   Phi'       g_x X   g_u L' ]
[ Phi        X - mu I   0       0      ]  >= 0
[ g_x X      0          mu I    0      ]
[ g_u L      0          0       mu I   ]

[ ubar_i^2  L_i ]
[ L_i'      X   ]  >= 0,        X <= R^2 I
```

`g = sqrt(2) c` when both sector coefficients are nonzero (Young's inequality), and `g = c` otherwise.
The objective maximizes the ball inscribed in the certified set.

## What is and is not guaranteed

If `certified` is true, then **for the model with any residual in the stated sector**, on
`E = {z : z' X^-1 z <= 1}`:

- every input satisfies `|u_i| <= ubar_i`;
- `V(z) = z' X^-1 z` satisfies `V(z+) <= lambda V(z)`;
- `E` is forward invariant and lies inside `‖z‖ <= R`.

Limits:

- **The certificate is about the model and the sector.** It transfers to a physical plant only if
  the sector really holds for that plant on `‖z‖ <= R`. The tools do not check that the data covers
  this ball.
- **A fitted sector is a data envelope, not a bound on unseen states.** A probability is reported
  only when the caller declares the sampling (`sampling_assumption`) and fixes both LP weights before
  sampling. The result is then a scenario bound conditional on that declaration, which DynSysLab
  cannot check. Trajectory samples are generally not i.i.d.
- **Verification is floating-point.** The verifier rebuilds each LMI from the deployed `(K, X, mu)`
  and checks its minimum eigenvalue against a margin plus a rounding allowance. It is not interval
  arithmetic or a formal proof.
- **`infeasible` is a statement about this sufficient LMI** at the requested sector, box, radius,
  contraction and margin. It does not prove that no stabilizing controller exists.
- **Output feedback is not covered.** The design assumes `z` is known exactly.
- **Size.** The solver is dense and first-order, and the design enumerates `2^m` input vertices
  (at most 12 inputs). It is meant for small reduced models.

The design returns a gain only when the solver converged, `X` is positive definite, and the
independent verifier passed. Otherwise the status is `infeasible`, `inconclusive` or
`verification_failed`, with no controller.

## Frozen-A bilinear K-ROMs

`fitBilinearKrom` estimates `A`. For a K-ROM whose autonomous dynamics came from unforced training,
use the frozen path, which never changes them:

- `fitBilinearKromFrozenA` regresses only `B0` and `B_i` on `z_plus - A z` and copies `A` bit for bit.
- `bridgeKoopmanBilinear` converts a stored bilinear control map. With `w = split(psi)`,
  `z = w - w*` and `v = u - u*`:

```text
A    = exp(G dt)                                     exact
B0   = int_0^dt exp(G s) ds b0                       exact zero-order hold
B_i  = int_0^dt exp(G (dt-s)) N_i exp(G s) ds        exact first order in u (Van Loan block exponential)

z+   = A_eff z + B0_eff v + sum_i v_i B_i z + c
A_eff  = A + sum_i u*_i B_i       (equal to A, bit for bit, when u* = 0)
B0_eff = B0 + [B_1 w* ... B_m w*]
c      = A_eff w* + B0 u* - w*
```

The `O(u^2)` part of `exp((G + u N) dt)` is outside this model family. `conversion_defect` reports it
at the input-box vertices, and the residual sector must absorb it.

**A target that is not an equilibrium.** If `||c|| > equilibrium_tolerance * max(1, ||w*||)`, the
status is `target_not_equilibrium`, `constant_term` is stored, and no gain is produced. A nonzero
residual at `z = 0` cannot lie in the sector, so the constant is never dropped. Choosing a `target_input`
that cancels it is the caller's decision.

### Post-processing method `bilinear_certificate`

A pooled method over **forced** records. It reads the K-ROM and writes a new archive under
`derived/koopman-bilinear-certificate/`.

| Setting | Meaning |
|---|---|
| `krom` | K-ROM holding a stored bilinear map (usually the `input_map` output). |
| `control_representation` | Map name (default `bilinear`). |
| `target` | Exactly one of `{"state": [...]}`, `{"coordinates": [...]}`, `{"lifted": [...]}`. History K-ROMs need coordinates or lifted values. |
| `target_input` | Held `u*` (default 0). Not available with `operators: "fit"`. |
| `input_bounds`, `state_radius`, `contraction` | The design box on `v`, radius `R` on `z`, and `lambda`. |
| `operators` | `converted` (default: the stored map, so KLQE and the design share one input model) or `fit` (frozen-A regression; the design then differs from the map KLQE propagates, and the run metadata says so). |
| `inflation`, `zero_input_tolerance`, `zero_residual_tolerance`, `minimum_coverage` | Supported-sector options (below). |
| `one_sided_sector` | `{"channel": "state" or "input", "reason": "..."}`: declare `c_x = 0` or `c_u = 0`. The reason is recorded. |
| `equilibrium_tolerance`, `ridge`, `fit_start`, `state_window`, `history_transform`, `solver_margin`, `certificate_margin`, `tolerance`, `max_iterations` | As named. |

Only pairs with `||z|| <= R` and `|v| <= ubar` enter the fit and sector. Input holds follow the strict
rule in `INPUT_MAP_B_PSI.md`. The status is `certified`, `infeasible`, `inconclusive`,
`verification_failed`, `target_not_equilibrium`, `insufficient_data_in_region`, `operator_fit_failed`
or `sector_unsupported`. Every outcome writes an archive; `certificate_K` and `certificate_X` exist only
when certified.

### Supported sector

The joint LP returns a vertex, so it can put the whole residual on one channel. With records that never
have `u = 0` it can return `c_x = 0` even though the autonomous residual is not zero, and a certificate
built on that sector is false for the plant. `bilinear_certificate` therefore estimates in two stages:

```text
c_x = max over zero-input holds of ||r|| / ||z||
c_u = larger of
      max over forced pairs of max(||r|| - c_x ||z||, 0) / ||u||                 (envelope)
      max over pairs with ||u|| >= 0.1 max ||u|| of
          max(||r|| - largest hold residual in the same ||z|| bin, 0) / ||u||    (identification)
```

The identification term matters because the global `c_x` is set by the largest radius and can otherwise
hide a real input-channel residual at smaller radii.

Both are then multiplied by `inflation`. The status is `sector_unsupported`, with the reason and no gain,
unless all of these hold:

- at least 8 zero-input holds and 8 forced pairs;
- zero-input holds in every occupied radial bin of `||z||` (4 bins, at least 2 holds each);
- zero-input holds reach `minimum_coverage * state_radius` (default 0.5);
- forced pairs reach `minimum_coverage * input_bound`;
- a zero coefficient means that channel's residual is zero to tolerance.

A one-sided sector is otherwise accepted only through `one_sided_sector`, and a declaration that the
data contradicts is refused. In the `residual_sector` query, `"method": "zero_input_first"` selects this
estimate. The joint LP result is marked `supported: false`, and `certified_design` refuses a sector so
marked. A sector the caller writes by hand, without that field, remains the caller's assumption.

The two-stage estimate is a DynSysLab construction on the SafEDMD error structure, not a published
method. It is still a data envelope: passing these checks is evidence about the sampled region, not a
bound on unseen states.

### Using the gain

```json
{"kind": "klqr", "settings": {"krom": "<input_map archive>", "control_representation": "bilinear",
  "bilinear_certificate": "<certificate archive>", "Q": 1, "R": 1}}
```

The law is `u = u* + K (z - z*)` with no Riccati solve. Note the sign: the certificate uses `u = K z`,
while LQR uses `u = -K x`. The run is refused unless the archive is certified and its K-ROM identity,
representation, generator and map match bit for bit, the controller samples at the certified `dt`, any
`control_limit` contains the certified input box, and any controller target equals `z*`. It cannot be
combined with integral action, PID or balanced truncation. `Q` and `R` are still parsed by the
controller but do not affect this gain.

With `estimator`, KLQE predicts with the same frozen generator and stored map. **The certificate does
not cover this output-feedback loop**, because the design assumes `z` is known exactly. Run metadata
records the certificate's scope and this limit under `bilinear_certificate`.

## JSON operations

Requests go through `dynsyslab_cli mathematics query.json`, with schema `dynsyslab.math-query/v1`.

LMI with a known value (`min x` subject to `[x+0.1, 1; 1, x+0.1] >= 0`, so `x = 0.9`):

```json
{"schema": "dynsyslab.math-query/v1", "operation": "sdp_feasibility",
 "variables": 1, "objective": [1], "margin": 0,
 "blocks": [{"dimension": 2, "constant": [[0.1, 1], [1, 0.1]],
             "coefficients": [[[1, 0], [0, 1]]]}]}
```

Fit (`z`, `z_next` and `u` are arrays of sample rows):

```json
{"schema": "dynsyslab.math-query/v1", "operation": "bilinear_fit",
 "state_dimension": 2, "input_dimension": 1, "sample_time": 0.1,
 "z": [[0.1, 0.2], "..."], "z_next": [[0.12, 0.19], "..."], "u": [[0.5], "..."]}
```

Add `"frozen_A": [[...], [...]]` to regress only `B0` and `B_i` with that `A` held fixed.

Sector, from a fitted model and held-input pairs (optional `inflation`, `state_weight`,
`input_weight`, `sampling_assumption: {"declaration": "...", "confidence": 0.99}`):

```json
{"schema": "dynsyslab.math-query/v1", "operation": "residual_sector",
 "model": {"schema": "dynsyslab.bilinear-krom/v1", "...": "..."},
 "z": [["..."]], "z_next": [["..."]], "u": [["..."]], "inflation": 1.1}
```

Design (optional `solver_margin`, `certificate_margin`, `tolerance`, `max_iterations`):

```json
{"schema": "dynsyslab.math-query/v1", "operation": "certified_design",
 "model": {"schema": "dynsyslab.bilinear-krom/v1", "...": "..."},
 "sector": {"c_x": 0.004, "c_u": 0.002},
 "input_bounds": [1.0], "state_radius": 1.0, "contraction": 0.995}
```

The result holds `status`, `certified`, solver residuals and, only when certified, `K`, `X`, `mu`,
`inscribed_radius` and `certificate`.

## References

- B. O'Donoghue, E. Chu, N. Parikh, S. Boyd. Conic optimization via operator splitting and
  homogeneous self-dual embedding. *J. Optim. Theory Appl.* 169 (2016) 1042-1068. The solver is a
  compact dense reimplementation of this method, not the SCS library.
- Y. Ye, M. Todd, S. Mizuno. An O(sqrt(n)L)-iteration homogeneous and self-dual linear programming
  algorithm. *Math. Oper. Res.* 19 (1994) 53-67.
- R. Strässer, K. Worthmann, et al. SafEDMD. arXiv:2402.03145.
- C. F. Van Loan. Computing integrals involving the matrix exponential. *IEEE Trans. Automat.
  Control* 23 (1978) 395-404.
- D. Goswami, D. A. Paley. Bilinearization, reachability, and optimal control of control-affine
  nonlinear systems: a Koopman spectral approach. *IEEE Trans. Automat. Control* 67 (2022) 2715-2728.
- D. Bruder, X. Fu, R. Vasudevan. Advantages of bilinear Koopman realizations for the modeling and
  control of systems with unknown dynamics. *IEEE Robot. Autom. Lett.* 6 (2021) 4369-4376.
- J. Bernussou, P. Peres, J. Geromel. A linear programming oriented procedure for quadratic
  stabilization of uncertain systems. *Syst. Control Lett.* 13 (1989) 65-72.
- S. Boyd, L. El Ghaoui, E. Feron, V. Balakrishnan. *Linear Matrix Inequalities in System and
  Control Theory*. SIAM, 1994.
- M. Campi, S. Garatti. The exact feasibility of randomized solutions of uncertain convex programs.
  *SIAM J. Optim.* 19 (2008) 1211-1230.
