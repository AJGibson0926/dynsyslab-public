# The input map B_psi: one switchable pipeline step

A controlled KROM is `psi' = Lambda psi + B_psi(psi) u`. The eigenfunctions
`psi` and `Lambda` come from unforced training (see
`TRAINING_DATA_ROLES.md`). `B_psi` is chosen independently, so the same
experiment can be rerun with only the input map changed.

## Post-processing method `input_map`

`workflow/native/src/InputMap.cpp`. It is a pooled method: all sources are
read together and one KROM archive is written to
`derived/koopman-input-map/`.

| Setting | Meaning |
|---|---|
| `krom` | The KROM whose `psi` is used (required). |
| `representation` | `exact` (analytic `grad psi . g(x)`, state/observable KROMs), `history_exact` (analytic `grad_xi psi` · newest-sample weight · `L_g y`; causal-history KROMs; the controller reads the simulator state), `constant`, or `bilinear` (`B(psi) = b0 + sum_i b_i psi_i`). |
| `fit` | For regressed maps. `held_input` (every KROM kind; the only fit for causal-history KROMs and for bilinear maps on measured coordinates), `controlled_residual`, `multistep_output` (constant only), or `exact_map` (state/observable KROMs). |
| `conjugate_complete` | Bilinear only. Adds `conj(psi_i)` regressors, so a complex mode can have a general real-linear dependence on its own real and imaginary parts. Default `false`. |
| `strict_hold` | `held_input` only. Keep only pairs whose input was really held (below). Default `true`, except `false` for causal-history KROMs, which keeps their earlier behaviour. |
| `drift` | `held_input` only. `none`, or `affine`, which jointly fits a nuisance term `a(psi) = a0 + sum a_i psi_i` that absorbs zero-input closure error. The controller does not use `a`. |
| `ridge` | Tikhonov weight relative to the largest normal-equation diagonal. |
| `fit_start` | Skip forced samples before this time (start-up transients). Honoured by every regressed fit on state/observable KROMs and by `held_input`. |
| `state_window` | `[[low, high], ...]` per state component: regress only on forced samples inside this box (each run of consecutive inside samples is its own record). Use the window the eigenfunctions were fitted on, so `B_psi` describes the same operating region. |
| `map_name` | Name of the stored map (default: the representation). |
| `measurement_observable` | History measurement (default: the frozen transform's). |
| `map` | An object whose fields are merged into the settings. A pipeline passes its whole `B_psi` parameter here. |

The sources must be **forced** records with a recorded scalar control
history when the representation is regressed. Analytic representations ignore
them. The archive's `metadata_json.input_map` records the representation,
fit, forced sources and regression diagnostics: pairs, relative residual
(overall and by mode), and input RMS.

### `held_input` regression

`psi_k` comes from the state (state/observable KROMs), the measured coordinates,
or the measured observable pushed through the frozen history transform
(causal-history KROMs). The eigenvalues stay frozen. Per mode `j`:

    q_k = psi_{k+1,j} - mu_j psi_{k,j} = [a_j(psi_k)] + D_j(psi_k) u_k,   mu_j = exp(lambda_j dt)

This is solved by complex least squares. The stored continuous map is the
zero-order-hold inverse, exact to first order in `u` (Van Loan, *IEEE Trans.
Automat. Control* 23 (1978) 395-404). With

    phi(a, b) = int_0^dt exp(a (dt - s)) exp(b s) ds

the offset is `b0_j = D0_j / phi(lambda_j, 0)`, the operator on `psi_i` is
`N_ji = D_ji / phi(lambda_j, lambda_i)`, and the operator on `conj(psi_i)` uses
`phi(lambda_j, conj(lambda_i))`. Before 2026-11-23 every operator used
`phi(lambda_j, 0)`, which is wrong at order `dt` when `i != j`, so bilinear maps
refitted now can differ slightly from older archives. Constant maps are unchanged.

A pair is valid only if `u` was constant across it. Interval-aligned
`control_interval_values` and sampled control on a clock that is a whole multiple
of the record step qualify. For other sample-aligned records, strict mode keeps
only pairs with `u_{k+1} == u_k`, and a sampled clock that is not a whole multiple
of the step is refused. `diagnostics.pairs_dropped_at_input_change` reports the count.

### How the map reaches the estimator and controller

- **KLQE** propagates `z' = (G - L C + sum_a u_a N_a) z + L y + b0 u` by exact
  zero-order hold for each held input. `N_a` and `b0` are the exact real-split form
  of the stored map, conjugate operators included.
- **KLQR** without a certificate re-solves the Riccati equation on `(G, B(z))` at
  each evaluation. This is a state-dependent factorization (SDRE). It does not add
  `u N` to `G`, has no feedforward for a nonzero target, and carries no guarantee.
- **Certified KLQR** (`bilinear_certificate`) replaces that gain with a checked
  fixed gain. See `CERTIFIED_BILINEAR_CONTROL_TOOLS.md`.

## Selecting B_psi in a pipeline

The rediscovery pipeline
(`examples/native/pipelines/rpe-echo-jet-klqr-rediscovery/20kpa.json`) has:

- `B_psi`: the map, e.g. `{"representation": "exact"}` or
  `{"representation": "bilinear", "drift": "affine", "fit_start": 5}`.
- `forced_training`: a solve job or `{kind: archive}` producing the forced
  record(s). The default is a 60 τ random-pulse run (20 kPa, 2 τ holds) from
  rest.
- Stage `input_map` reads `${forced_training.archive}` with
  `{"krom": "${echo_krom.archive}", "map": "${params.B_psi}"}`.
- The controller uses `krom: ${input_map.archive}` and
  `control_representation: ${params.B_psi.representation}`. Nested parameter
  fields such as `${params.B_psi.representation}` are resolved by the
  pipeline.

Examples:

    dynsyslab-cli run .../20kpa.json --param 'B_psi={"representation":"constant","fit_start":5}'
    dynsyslab-cli run .../20kpa.json --param 'B_psi={"representation":"bilinear","drift":"affine"}' \
        --param 'forced_training={"kind":"archive","archives":["<a forced run you generated>"]}'

`KLQR` (like every feedback controller) accepts `control_limit`: a number for a
symmetric clip, `[low, high]` for an asymmetric one, and `null`/`false` for none.
In this pipeline it is the parameter `control_limit` (default `null`, uncapped),
e.g. `--param control_limit=0.5`, i.e. ±50 kPa.

The older `history_control_regression` method remains: it is a constant
`held_input` map without drift.

## Results: the 20 kPa, 100 kHz rediscovery in causal-history (xi) coordinates (2026-11-14)

The eigenfunctions came from the late, unforced fine baseline. H-DMD ran on
τ 50–80 and gave λ = −0.0601 ± 6.2829i. SINDy then fit ψ(ξ) on τ 40–80 from a
50-delay (0.05 τ) echo history with 2 latents and degree 1. All runs used
KLQR with R = 1000, a ±50 kPa limit and 300 τ, and all are scored over
τ ≥ 100 against the withheld sine.

| B_psi | Forced training | Recovered fundamental | Phase error | Correlation | Radius RMSE |
|---|---|---|---|---|---|
| `history_exact` (analytic) | none | 3.0 kPa | −132° | −0.03 | 50 (diverged) |
| constant, `held_input` | random pulses 5 kPa / 0.5 τ | 18.4 kPa | 1.2° | 0.9997 | 0.0030 |
| constant | 20 kPa / 0.5 τ | 18.4 kPa | 1.2° | 0.9997 | 0.0030 |
| constant | 20 kPa / 2 τ | 18.4 kPa | 1.3° | 0.9997 | 0.0031 |
| bilinear, ± affine drift | any of the above | 18.4 kPa | 1.2–1.4° | 0.9997 | 0.003 |
| any regressed form | 100 kPa / 0.5 τ | 1.5–6.6 kPa | wrong sign | ≈ 0 | 8–13 (diverged) |

- **The form of B_psi barely matters.** Constant, bilinear and the drift
  variants are within 0.1 kPa of each other. The forced training does matter:
  pulses beyond the near-linear regime (100 kPa) give the wrong sign.
- **R trades smoothness against drive amplitude.** At R = 300 the loop
  chatters at the actuator limit (correlation 0.35, fundamental still
  18.5 kPa). R = 1000 is best. R = 3000 is smooth but gives 17.2 kPa.
- **Only the short, 2-latent history works in closed loop.** 4 latents or
  longer histories (100–1275 delays) fail even where the offline inverse-driver
  gate looked good. The analytic `history_exact` map fails in every
  configuration.
- **The remaining ~8% amplitude shortfall is structural.** It is the steady
  tracking error that an LQR without feedforward needs in order to produce any
  drive.
