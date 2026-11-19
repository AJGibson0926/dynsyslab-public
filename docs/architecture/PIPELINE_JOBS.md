# Pipeline jobs: one submitted, deterministic Koopman study

A **pipeline** is a compute job (`"kind":"pipeline"`) that runs an ordered list of
ordinary compute jobs ("stages"). An agent or person writes it once (training
data, post-processing, KROM rule, controller, target) and submits it with
`dynsyslab-cli run`. No intermediate decisions or path copying are needed. It
is implemented in the native workflow layer (`workflow/native/src/Pipeline.cpp`)
and reuses `executeJob` for every stage, so stages keep their existing
contracts, validation, archives and index publication.

## Contract

```json
{"schema":"dynsyslab.compute-job/v1","kind":"pipeline","label":"study",
 "parameters":{"coordinates":{"kind":"identity"}},
 "stages":[
  {"id":"train","job":{"kind":"solve", ...}},
  {"id":"coords","job":{"kind":"postprocess","sources":["${train.archive}"],
     "methods":[{"name":"coordinate_transform","settings":"${params.coordinates}"}]}},
  {"id":"hdmd","job":{"kind":"postprocess","sources":["${coords.archive}"], ...}},
  {"id":"sindy","job":{"kind":"postprocess","sources":["${coords.archive}"],
     "methods":[{"name":"sindy","settings":{"spectrum_archive":"${hdmd.archive}", ...}}]}}]}
```

References are resolved when a stage starts:

| Reference | Meaning |
|---|---|
| `${id.archive}` | the stage's single published archive (error if it published several) |
| `${id.archives[N]}` | its N-th archive |
| `${id.output:METHOD}` | the archive a postprocess stage produced for METHOD |
| `${id.manifest}` | its hierarchy, campaign or pipeline manifest |
| `${params.NAME}` | a named parameter; as a whole string it is replaced by the JSON value |

`dynsyslab-cli plan` rejects references to later or unknown stages before any
work runs. It plans independent stages immediately and marks dependent ones
deferred.

## Running without editing files

Every input can be supplied on the command line (also honored by `queue submit`):

```powershell
dynsyslab-cli run base.json --param coordinates='{"kind":"observable","observable":"echo_jet2"}' `
  --set control.controller.settings.Q=1000 --set sindy.methods.0.settings.threshold=0.02
dynsyslab-cli run base.json --until krom          # stop after a stage
dynsyslab-cli run base.json --refresh sindy       # recompute a stage (and its dependents)
dynsyslab-cli run base.json --fresh               # recompute everything
```

`--param NAME=VALUE` sets a pipeline parameter. `--set STAGE.path=VALUE` sets
any field of a stage's job (numeric segments index arrays). VALUE is JSON,
otherwise text. The inputs are recorded in the manifest under `command_inputs`.

## Per-stage report

`python tools/pipeline_report.py run JOB [same flags]` runs the pipeline and
writes `viewer_outputs/pipeline-reports/<run>/report.html`. That is a
self-contained page with the job's dependency graph, then one section per stage:
the graph with that stage bolded, its method settings, measured results and its
native Viewer plot (training trajectory, H-DMD spectrum and spectral-measure stem
plots, SINDy/KROM eigenfunction dashboard, controlled trajectory against target).
Failed and not-run stages are marked. `pipeline_report.py report MANIFEST...`
rebuilds reports for earlier runs.

## Stages and what feeds them

A study has five kinds of stage, each chosen by a parameter of the job:

| Stage | Parameter | Examples |
|---|---|---|
| Training data | `training` | a solve job, or `{"kind":"archive","archives":[path]}` for an existing trajectory |
| Coordinate changes | `coordinates` | identity, `echo_jet2`, history, or a chain of them |
| Post-processing methods | `methods` | `["hdmd","sindy"]`, `["hdmd"]`, `["mpedmd"]` |
| Post-processing artefacts | `artefact` | K-ROM extraction rule |
| Control and estimation | `control_*`, `target` | KLQR, KMPC |

A stage `{"id":"methods","sequence":"${params.methods}","sources":[...]}` expands
into one stage per method (`methods_1`, `methods_2`, ...), each reading the same
sources and each reused or rerun on its own. A method named by a string takes its
settings from `parameters.method_settings.NAME`. A method that takes a spectrum
(SINDy, boundary eigenfunctions, eigenvalue optimization) is fed the previous
method's archive unless it names `spectrum_archive`. `${methods.archive}` is the
last method's output, so the artefact stage reads whichever method ran last: SINDy
eigenfunctions, or the interpolated eigenfunctions of H-DMD, mpEDMD or another
spectral method (all publish the same K-ROM contract: eigenfunction values,
physical training states and continuous rates). `--until methods` stops after the
sequence; `--refresh methods` recomputes all of it.

```powershell
dynsyslab-cli run base-r.json --param methods='["hdmd"]'           # K-ROM straight from H-DMD
dynsyslab-cli run base-r.json --param methods='["mpedmd"]'         # mpEDMD instead of H-DMD + SINDy
dynsyslab-cli run base-r.json --param training.initial_conditions.values='[[1.8,0]]'
dynsyslab-cli run base-r.json --param training='{"kind":"archive","archives":["archives/solutions/rpe/run.npz"]}'
dynsyslab-cli run base-r.json --param method_settings.sindy.threshold=0.02
```

`--param NAME` accepts a dotted path into the merged parameters.

## Coordinate system

Every study runs training data -> coordinate change -> post-processing. The
`coords` stage of the RPE example is one `coordinate_transform` whose settings
are the `coordinates` parameter; H-DMD and SINDy both read its archive:

| `coordinates` | Working coordinates |
|---|---|
| `{"kind":"identity"}` | the state (r, v) |
| `{"kind":"observable","observable":"echo_jet2"}` | (e, e'), exact Jacobian kept for SINDy/KLQR |
| `{"kind":"history","measurement":"far_field_acoustic_echo","delays":1275,...}` | xi, Hankel-SVD of the echo history, fitted in this stage |
| `{"kind":"chain","steps":[{"kind":"observable","observable":"echo_jet2"},{"kind":"history","measurement":"far_field_acoustic_echo",...}]}` | (r, v) -> (e, e') -> xi built from the echo column |


A chain runs any number of steps, each reading the previous step's output; the
stage publishes the last step's coordinates and records the intermediate
archives under `coordinate_chain`. Every step carries the physical states
forward. A catalog observable (`kind: observable`, `channels`) maps the
physical state, so it must be the first step; later steps select columns.
`hdmd_observable` selects which coordinate column H-DMD embeds (default 0).
See [Observable Hankel matrices](OBSERVABLE_HANKEL.md#module-boundary).

## Variants

`"extends":"base.json"` inherits a base pipeline (paths resolve against the
working directory, like all job paths). A derived pipeline may:

- override `parameters` (merged as an RFC 7386 patch),
- replace a stage by id or add one (`"before":"id"` inserts it ahead of a base stage),
- apply `"overrides":{"id":{...}}`, an RFC 7386 merge patch on a stage's job.

Arrays such as `methods` are replaced wholesale by merge patches. Expose a value
the study will vary as a parameter instead. A coordinate study is then a few
lines (see `examples/native/pipelines/rpe-viscous-klqr-multisine/`).

## Provenance

Each run publishes `archives/campaigns/pipelines/<label>__<utc>.pipeline.json`
(`dynsyslab.pipeline/v1`). It holds the extends chain, the fully expanded job,
and each stage's fully resolved job, archives and producing methods. A study is
reproducible from that one file.

## Stage features added with pipelines

- **SINDy row selection**: `sample_stride` and `time_window` choose the fitting
  rows explicitly. H-DMD's Hankel `sample_stride` does not propagate.
- **SINDy nullspace solver** (`nullspace_solver`, default `gram`): feature Gram
  blocks are formed once and each thresholded null vector comes from a
  feature-sized Hermitian eigenproblem. It matches the tall-SVD reference (`svd`)
  to about 1e-11 with identical sparsity on the standard RPE run.
- **SINDy `"reconstruction":"none"`** skips the all-mode least-squares C fit
  (samples x retained modes). That fit dominated SINDy runtime (6.7 min to 5 s
  on the standard RPE run) and is unused by state-Jacobian/coordinate-target control.
- **KROM rule selection**: `krom_extract` accepts
  `"mode_selection":{"rule":"oscillatory_pairs"|"leading","count":N}` over the
  parent's stored (energy) order, completing conjugates. An optional `"max_real_part": R`
  skips growing modes (continuous rate with real part above R).
- **H-DMD rank rule for noisy data**: `"rank_rule": "optimal_hard_threshold"` keeps singular values above
  the Gavish–Donoho threshold ω(β)·median(σ), so noise modes do not reach the spectrum (`rank` still caps it).
- **Noise-robust derivatives and measurement noise**: see `NUMERICAL_DIFFERENTIATION.md` (coordinate kind
  `measured_derivatives`, SINDy/`input_map` `derivative`, KLQE `measurement_noise` / `measurement_derivative`,
  and the `derivative_benchmark` method).
- **Target ramp toggle**: any state target (constant, waveform, trajectory,
  rpe_canonical) may carry `"ramp":{"duration":T,"from":[x0]}`. It blends from
  a rest state with a quintic smoothstep. For each `derivative_pairs` entry
  `[position, velocity]` (default `[[0,1]]` when the state has two or more
  components), the velocity stays the exact derivative of the ramped position.
  `from` should have zero velocity. `"ramp"` absent, `null`, `false` or
  `{"enabled":false}` uses the target as given, with no ramp.
- **Control cap toggle**: every feedback controller (KLQR, LQR, SDRE, PID,
  KLQE→KLQR, and prescribed inputs) accepts `"control_limit"`. A number `L`
  clips the command to `[-L, L]`, and `[low, high]` clips asymmetrically.
  Absent, `null` or `false` leaves it uncapped. MPC and KMPC treat
  `control_limit` as an optimizer constraint instead.
- **Pipeline toggles**: `rpe-viscous-klqr-multisine/base-r.json` exposes both toggles as
  parameters: `--param target_ramp=null` (no ramp) or `--param 'target_ramp={"duration":10,"from":[1,0]}'`,
  and `--param control_limit=null` (uncapped, the default) or `--param control_limit=1.5`.
  `rpe-echo-jet-klqr-rediscovery/20kpa.json` exposes `control_limit` and takes its target as
  `--param target_state=<state target>` (with its own optional `ramp`).
- **Offline target lifting (`control_target`)**: this post-processing method maps a state target
  into a KROM's coordinates offline, on a fixed clock. Measured coordinates such as (r, r') or
  (e, e') are evaluated along the target; causal-history coordinates give xi_ref through the
  frozen transform. It writes `time`, `values` and `state_target_values`. The controller then
  reads `{"space":"observable","kind":"trajectory","source":<archive>}`, so online control
  evaluates no physical model for its reference. `20kpa.json` runs it as stage `target_lift`.
  `target_step` (default 0.001, which is the history clock) sets the sampling;
  `continuous_rhs` evaluation interpolates it linearly. Data roles of each controller:
  `CONTROL_DATA_ROLES.md`.
- **History-coordinate control**: a `from_state` target is projected through the
  frozen history transform (`xi_ref`). `control_representation:"history_exact"`
  uses B = grad_xi psi * m_newest * L_g y(x), which reads the simulator state for
  L_g y. `history_control_regression` fits a measurement-only constant map from
  forced runs (held-input regression, ZOH-inverted to continuous B).
