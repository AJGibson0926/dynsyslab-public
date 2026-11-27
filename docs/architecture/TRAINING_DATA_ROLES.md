# Training data roles: autonomous eigenfunctions, forced input maps

DynSysLab builds a controlled Koopman reduced-order model (KROM) in two
separate identification problems:

    psi' = Lambda psi + B_psi(psi) u,        y = C_y psi,   x = C_x psi

1. **The autonomous part** — the spectrum `Lambda`, the eigenfunctions `psi`,
   and the readouts `C` — describes the *unforced* flow. It is identified
   from **unforced (autonomous) training data only**.
2. **The input map** `B_psi` describes how the input moves `psi`. It is either
   computed analytically from `psi` and the input vector field (`exact`,
   `history_exact`), or **regressed from forced training data**.

The one deliberate exception is the DMDc / EDMDc family, whose purpose is to
identify the autonomous operator and the input matrix *jointly* from forced
data.

The pipeline enforces this. Every post-processing method has a data role
(`trainingDataRole`, `workflow/native/src/KromDatasetSupport.cpp`). Each
source is classified as `forced`, `unforced` or `unknown` before the method
runs (`datasetDynamics`):

- **forced:** a non-zero recorded input, a controller other than `none`, or a
  forcing type other than free.
- **unforced:** recorded zero input, `controller none`, or `forcing_type free`.
- **unknown:** anything else.

Coordinate changes carry this forcing provenance forward, so an H-DMD or
SINDy stage that reads a coordinate archive is checked against the record the
coordinates came from.

| Role | Methods | Accepted sources |
|---|---|---|
| `autonomous` | `dmd`, `hdmd`/`hankel_dmd`, `sindy`, `edmd`, `havok`, `koopman_spectral`, `resdmd`, `mpedmd`, `hermitian_dmd`, `tlsdmd`, `multdmd`, `stochastic_resdmd`, `specrkhs`, `boundary_eigenfunctions`, `eigenvalue_optimization`, `echo_history_koopman`, `history_basis`, `sindy_dynamics`, `koopman_validation` | unforced, or `unknown` (no forcing provenance: legacy imports, synthetic or measured records). **Forced is refused.** |
| `forced_allowed` | `dmdc`, `edmdc` | forced or unforced |
| `forced_required` | `control_regression`, `history_control_regression`, `input_map` with a regressed representation (`constant`, `bilinear`) | forced, with a recorded control history |
| `any` | decompositions (`svd`, `pca`, `pod`, `hankel_svd`, `hankel_pod`, `hankel_matrix`, …), `coordinate_transform`, `observable_reconstruction`, `krom_extract`, `linear_analysis`, `controllability`, `observability`, `klqe`, `input_map` with `exact`/`history_exact` | anything |

A refusal names the method, the source file, its classification and the
right alternative: `dmdc`/`edmdc` to identify from forced data, or
`input_map` to regress `B_psi` from it.

## Why each algorithm needs what it needs

### Autonomous eigenfunction methods

A Koopman eigenfunction satisfies `psi(x(t)) = exp(lambda t) psi(x(0))` along
trajectories of the **unforced** flow. Every method in the `autonomous` row
estimates that relation from consecutive samples (or from `dpsi/dt`). Along a
forced trajectory the measured evolution is instead
`psi' = lambda psi + (grad psi . g) u`. Fitting it as if it were free
therefore smears the input response into the spectrum and the eigenfunctions.

- **Hankel DMD / H-DMD** (Arbabi & Mezić 2017, *SIAM J. Appl. Dyn. Syst.*
  16(4); arXiv:1611.06664). Delay embedding of an ergodic, autonomous
  trajectory. The Hankel columns are time shifts of one free record, and the
  convergence to Koopman spectra relies on the record being a trajectory of
  that flow. A forced record does not satisfy the shift relation.
- **EDMD / DMD** (Williams, Kevrekidis & Rowley 2015, *J. Nonlinear Sci.*
  25; Schmid 2010). Least-squares approximation of the Koopman operator for
  the snapshot pairs `x_k -> x_{k+1}` of the autonomous map.
- **ResDMD, mpEDMD, Hermitian/TLS/multiplicative DMD, specRKHS**
  (e.g. Colbrook & Townsend 2024, *Comm. Pure Appl. Math.*,
  arXiv:2111.14889; Colbrook 2023, *SIAM J. Numer. Anal.*,
  arXiv:2209.02244). Residual and measure-preserving refinements of EDMD.
  They inherit its autonomous-snapshot assumption, and their residual
  certificates are only meaningful for the free operator.
- **HAVOK** (Brunton et al. 2017, *Nature Communications* 8:19). Treats the
  last delay coordinate as an *intrinsic* forcing of an autonomous chaotic
  system. External actuation would be indistinguishable from it.
- **SINDy eigenfunctions** (Kaiser, Kutz & Brunton 2021, *Mach. Learn.: Sci.
  Technol.* 2:035023; arXiv:1707.01146). They solve the PDE
  `grad psi . f(x) = lambda psi` for the **drift** `f`. Control enters
  afterwards through `B = grad psi . g`, which is exactly DynSysLab's
  `exact` / `history_exact` input map. Training on forced data (or with
  `include_input_terms`) puts `g u` into the PDE and is refused.
- **Causal-history (xi) bases** (`history_basis`, the history step of a
  coordinate chain, `echo_history_koopman`). Along a free trajectory a
  measurement history is a function of the current state, so `psi(xi)` can be
  an eigenfunction. Along a forced trajectory the same state can carry
  different histories, so a basis or `psi` fitted there need not close.

### DMDc and EDMDc (forced data by design)

DMDc (Proctor, Brunton & Kutz 2016, *SIAM J. Appl. Dyn. Syst.* 15(1);
arXiv:1409.6358) regresses `x_{k+1} = A x_k + B u_k` jointly. EDMDc (Korda &
Mezić 2018, *Automatica* 93; arXiv:1611.03537) does the same on a lifted
dictionary. With rich enough excitation the joint regression separates `A`
from `B`, which is why these two methods may, and usually must, use forced
data. Their eigenfunctions are those of the identified `A`. DynSysLab keeps
them as their own model family (`koopman-dmdc`, `koopman-edmdc`); they do not
replace the unforced eigenfunction bank.

Both accept `time_window` and `state_window` (`[[low, high], ...]` per snapshot component): only
pairs whose two snapshots lie inside both windows are fitted, which restricts the model to the
neighbourhood of an operating point. `edmdc` with `export_krom: true` also publishes the eigenfunction
bank as a standard K-ROM for `krom_extract` (including its `principal` mode rule: principal
eigenfunctions as in Mauroy & Mezić, IEEE TAC 61(11), 2016, with a ResDMD-style consistency filter
after Colbrook & Townsend, Comm. Pure Appl. Math. 77, 2024), `input_map`, KLQE and KLQR. Centering (`normalize: true`) without a constant term places the model's equilibrium at the
data mean, so use `normalize: false` (or a constant) when the target is not the data mean.

### Input-map regression (forced data only)

`input_map` (see `docs/architecture/INPUT_MAP_B_PSI.md`) keeps `psi` and
`Lambda` fixed and regresses only `B_psi`. That needs input excitation.
On unforced data `u = 0` and the regression is singular, so unforced sources
are refused.

## Practical rules

- Give each pipeline two inputs: `training` (unforced, for every autonomous
  stage) and `forced_training` (for `input_map` / control regression), as
  `examples/native/experiments/rpe-echo-jet-klqr-rediscovery/20kpa.json` does.
- To fit on part of a free record (for example the small-amplitude tail of a
  decay), use the method's `time_window`. Both `hdmd` and `sindy` accept
  `[start, stop]`. Do not splice forced data in.
- Validation of a free KROM (`koopman_validation`) is also autonomous. Forced
  hold-out data validates the controlled model, for example through the
  offline inverse-input test, not the free semigroup law.
