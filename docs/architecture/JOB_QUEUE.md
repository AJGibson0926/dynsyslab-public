# Native job queue, campaigns and benchmarks

The native job store (`<archive_root>/job_manager/state.sqlite`, `native/src/JobStore.cpp`) is the
only job history. The queue (`workflow/native/src/Queue*.cpp`, library `dynsyslab_queue`) submits,
schedules and runs `dynsyslab.compute-job/v1` jobs; the CLI is `dynsyslab-cli queue ...`.

## Scheduling

Scheduling follows the batch-scheduler conventions of Slurm (Yoo, Jette & Grondona 2003,
"SLURM: Simple Linux Utility for Resource Management", JSSPP, LNCS 2862). A job is claimed only
when every rule below allows it; jobs are considered in submission order.

| Rule | Job key | Meaning |
|---|---|---|
| Retry | `queue.retry_limit` | A failed attempt is requeued while attempts <= retry_limit. |
| Dependencies | `queue.depends_on: [id, {"job": id, "condition": c}]` | `afterok` (default): the job completed. `afterany`: it ended (completed, failed, cancelled or quarantined). `afternotok`: it failed, was cancelled or quarantined. `aftercorr`: `job` names an array id; task *i* waits for task *i* of that array to complete. A dependency that can no longer be met cancels the job with the reason. |
| Arrays | campaign `systems_per_job`, `queue.array`, `queue.max_concurrent` | A campaign job becomes one queue task per chunk of points; at most `max_concurrent` tasks of one array run at once. |
| Spec claims | `queue.spec_claim` (default true) | Jobs with the same scientific content (the job minus `label`, `queue`, `execution` and the campaign id) never run at the same time. |
| Admission | `queue.resources.{threads, memory_mib}` (default: `execution` limits) | `queue serve --max-running N --pool-threads T --pool-memory-mib M` bounds what may run beside the running jobs. A single job is admitted when nothing else runs, so a limit never deadlocks the queue. |
| Side tray | `queue.hold` or `queue submit --hold` | The job is `reserved`: listed by `queue tray`, run after `queue release ID`, deleted by `queue remove ID`. |
| Reuse | `queue.reuse: "identical"` or `"extend"` | An identical completed spec is reused (the submit result names `reused_job`) instead of queued; `extend` continues a shorter completed run (see Reuse and extension). |
| Deadline | `queue.wall_time_seconds`, `queue.checkpoint_before_deadline_seconds`, `queue.max_requeues` | Past the wall time the job stops and is requeued. With a margin, a solve checkpoints and pauses that long before the deadline and is requeued to continue from the checkpoint. Requeues spend the retry budget, which `max_requeues` raises. |
| Staging | `queue.staging {inputs, cleanup_scratch, output_location, stage_out}` | Verified input copies (below). |

Statuses: `queued`, `running`, `completed`, `failed`, `cancelled`, `paused`, `quarantined`,
`reserved`.

## Execution limits

`execution.max_threads` and `execution.max_memory_mib` in a job, and `--max-threads N` /
`--max-memory-mib M` on `run`, `queue run-next` and `queue serve`, set `ExecutionLimits`; the
smaller value wins. The memory limit is checked before buffers are allocated (integrators,
lattices, PDE fields and the post-processing preflight below). The thread limit caps kernels that
spawn workers (`core::workerThreads`) and sets `OMP_NUM_THREADS` and the BLAS thread variables.

## Pause, watchdog and recovery

`queue pause ID` pauses a queued job at once. A running job is asked to pause: it stops at the next
cooperative checkpoint and becomes `paused`; `queue resume ID` requeues it. It restarts from its
beginning, except a solve with `history.pause_at_checkpoint`: the queue then writes the solve's
pause file (`job_manager/pause/<id>.pause`), the solve stops at its next checkpoint and archives
its samples, and the resumed job continues from there.
`queue watchdog [--grace-seconds S]` (also run by `queue serve`) checks running jobs with a pause
request: a worker process that has exited, or a worker whose heartbeat is older than the grace
period, is `quarantined` rather than requeued, so a live rogue worker is never duplicated.
`queue resume ID` returns a reviewed job to the queue. `queue recover AGE` requeues running jobs
without a pause request whose heartbeat is older than AGE seconds.

### Checkpoints and continuation (`QueueContinuation.hpp`)

A solve that pauses at a checkpoint (deadline or pause) leaves its outcome record's
`resume_jobs` (`NATIVE_SOLVER.md`) in `job_manager/continuations/<id>.json`
(`dynsyslab.queue-continuation/v1`). The next attempt of the same queue job runs those jobs with
the current `history` settings instead of the original job, so the job id stays the same. The
completed job's result lists every archive, the paused prefixes first, and a `continuation`
record with the number of rounds. Each continuation archive names its parent with a SHA-256.
The deadline policy ports `job_manager/scheduling.py` `DeadlinePolicy`.

### History policy (`HistoryPolicy.hpp`)

`history_policy` on a solve (or `campaign.checkpoint_policy`, set by `checkpoint_profile`) decides
before execution whether a run is checkpointed (port of `job_manager/history_policy.py`):
`enabled`, `estimated_seconds` (default `queue.estimated_seconds`, then `queue.wall_time_seconds`),
`if_estimated_seconds_over`, `when_estimate_unknown` (`disable` or `enable`), `mode`
(`checkpoint` or `telemetry`), `checkpoint_interval` (model time), `checkpoint_interval_seconds`
(wall clock, to `history.checkpoint_wall_seconds`) and `pause_at_checkpoint`. The decision is
recorded in `history.policy_resolution`.

### Input staging (`QueueStaging.hpp`)

`queue.staging.inputs` lists files copied into `job_manager/staging/<id>/` before the job runs.
Each copy's SHA-256 is checked against its source, and `staging-manifest.json`
(`dynsyslab.staging-manifest/v1`) lists source, staged path, digest and size. The manifest is
added to the job result under `staging`. A remote command gets the folder as `{staging_dir}`.
The folder is removed afterwards unless `cleanup_scratch` is false. Primary archives are only
read (port of `job_manager/staging.py`).

## Resident loop and remote execution

`queue serve [--workers N] [--poll-seconds S] [--max-jobs N] [--idle-exit S]` runs admissible jobs
until stopped, idle for S seconds, or N jobs have run. With
`--command-template "ssh host dynsyslab-cli run {job_file} --archive-root {archive_root}"` each
claimed job is written to `job_manager/remote/<id>.json` and the command runs instead of the
in-process workflow; `{job_id}` and `{result_file}` are also substituted. Exit status 0 completes
the job (with the JSON the command left in `{result_file}`, if any); any other status fails it
under its retry limit. Use the same template for scheduler front ends (`sbatch --wait ...`).

## Campaigns

A solve job with `campaign.enabled` runs a Cartesian product of axes (`Campaign.cpp`). Design keys
(`CampaignDesign.cpp`, `dynsyslab-cli campaign profiles` lists every table):

- Axis `{"path", "values"}` or `{"path", "start", "stop", "points"}` (linspace). A `driver.<field>`
  path addresses the periodic driver, the `waveform` controller's settings.
- `profile`: `low` 31 points / 10 per job, `medium` 101 / 25, `high` 301 / 50, `ultra_high` 1001 /
  100 (bifurcation resolution levels).
- `systems_per_job`: points per queue task. `point_range: [first, count]` runs one chunk.
- `execution_profile`: `fixed` (retry 0), `adaptive_local` (retry 1), `adaptive_portable` (retry 1,
  remote allowed).
- `checkpoint_profile`: `fast` or `recoverable`. `recoverable` checkpoints a run estimated to take
  over 300 s every 600 s of wall time (History policy above); a campaign array pauses and resumes
  at chunk boundaries.
- Cycle-based time: `time.transient_cycles`, `time.retained_cycles` and `time.samples_per_cycle`
  with the driver frequency (`time.frequency_hz` or the waveform's `frequency_hz`) give
  `start`, `stop = start + (transient + retained) / f`, `output_step = 1 / (f samples_per_cycle)`
  and `retained_start`. `time.cycle_profile` `lyapunov_screen` (2500 transient, 1024 retained) or
  `lyapunov_verification` (5000, 4096) fills the cycle counts for Lyapunov-exponent screens
  (Wolf, Swift, Swinney & Vastano 1985, "Determining Lyapunov exponents from a time series",
  Physica D 16:285-317). A swept frequency re-resolves the grid per point.

Every campaign run records its original plan and each chunk's manifest in the campaign registry
(`core/CampaignRegistry.hpp`, tables in the job store). `campaign list`, `campaign status ID`
(queue state and progress per chunk) and `campaign validate ID` (every planned point published
once, with its planned coordinates and existing solution archives; exit code 4 when invalid)
read it. Schemas: `dynsyslab.campaign-list/v1`, `dynsyslab.campaign-status/v1`,
`dynsyslab.campaign-validation/v1`.

`campaign repair ID` resubmits unpublished chunks whose job failed, was cancelled or is missing.
`campaign restructure ID --points-per-chunk N [--reason TEXT]` splits every unstarted (queued,
reserved or paused) unpublished chunk larger than N into smaller queue jobs over the same points.
It cancels the original job and records the split in `job_manager/restructure/`
(`dynsyslab.job-restructure/v1`). `campaign review ID` is a read-only coverage review: every planned
point is covered by exactly one chunk, every chunk is backed by a queue job or an existing
manifest, and the chunks needing repair are listed. It is written to `job_manager/reviews/`
(`dynsyslab.campaign-review/v1`, exit code 4 when not safe). These port the Python JobManager's
`restructure`, `controlled_parallel_transition`, `review_parallel_transition` and
`repair_transition_interruptions`. Native arrays run in parallel from submission, so the serial
to parallel transition is a re-chunking of unstarted work. A paused solve's checkpoint is
continued by its own job rather than replaced.

## Reuse and extension

`queue.reuse: "identical"` returns the newest completed job with the same spec key. A solve that
differs from a completed solve only in `time.stop` is reported in `extension_candidates`
(`relation` shorter or longer). With `queue.reuse: "extend"` the longest completed shorter run
whose archives cover every trajectory is continued: the queued job has `time.start` at its stop,
`initial_conditions {mode: continue_from, archives}` and an `extension_of` record (job, stop,
archives, the requested start and initial conditions). The submit result names
`extended_from`. Archives stay immutable: the parent and the continuation archives together hold
the full trajectory, linked by `continuation.parent_archive` and its SHA-256 (port of
`job_manager/reuse_services.py` `_try_single_solution_extension`). The restart is exact for
fixed-step methods and memoryless controllers (`NATIVE_SOLVER.md`).

## Telemetry

Progress files (`dynsyslab.compute-progress/v1`) carry `heartbeat_utc`, `elapsed_seconds`,
`estimated_remaining_seconds` (elapsed time scaled by the remaining fraction; -1 before any
progress) and `seconds_per_trajectory`. Queue heartbeats are the job's `updated_utc`.
Derived archives record `timing` (`dynsyslab.stage-timing/v1`: method and `stage_seconds`
from the start of the step to the archive's publication).

## Post-processing plan geometry and memory preflight

`dynsyslab-cli plan` on a post-processing job adds `geometry` (`dynsyslab.postprocess-geometry/v1`):
for each source archive (up to 256 MiB) its samples, state dimension and uniform clock; for Hankel
methods the matrix shape, the retained-rank bound and a dense working-set estimate (four
matrix-sized buffers for the SVD); for SINDy the polynomial feature count. Running the job refuses
it before any allocation (`ExecutionLimitExceeded`) when the estimate exceeds the smallest of the
execution memory limit, `execution.max_memory_mib` and any stage's `settings.max_memory_mb`.
Hankel analysis follows Arbabi & Mezic 2017 (SIAM J. Appl. Dyn. Syst. 16:2096-2126).

## Benchmarks and performance gates

`dynsyslab-cli benchmark [--iterations N] [--warmup N] [--case NAME]... [--output FILE]` times the
native kernel suite (`benchmark cases`: RK4, RK45, DOP853 and BDF Duffing solves, a 32x32 Ising
lattice, a degree-4 EDMD) and writes `dynsyslab.benchmark-report/v1` with every sample, the
min/median/max and the runtime fingerprint (version, commit, compiler, build type, machine).
`benchmark gate <reference> <candidate> [--required-speedup X] [--require CASE=X]...
[--any-runtime]` writes `dynsyslab.performance-gate/v1`: median speedup (reference / candidate)
per case, failing (exit code 4) when a requirement is missed or the runtimes differ. Medians resist
the outliers of a busy machine (Chen & Revels 2016, "Robust benchmarking in noisy environments",
arXiv:1608.04295).

## Templates

`dynsyslab-cli template solve|postprocess|estimate [--output FILE]` prints a starter job.
