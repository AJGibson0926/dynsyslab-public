# Native port ledger

Goal: every capability lives in the self-contained C++ core and its job pipeline, with only a thin
Python client (`dynsyslab_client`, the front door, the audit, the test runner). No external MATLAB
code: no toolbox adapters and no MATLAB parity scripts. A Python, MATLAB or side-car copy is
deleted only after its native replacement exists, is tested and is documented with citations.

The four ledgers below were compiled on 2026-11-25 by reading the Python and side-car code and
searching the native source for each capability. Every row cites a native file and symbol, or the
search that found nothing.

| Area | Ledger | Covered | Partial | Gap | Other |
|---|---|---|---|---|---|
| Post-processing, spectral, analysis data | [POSTPROCESSING.md](POSTPROCESSING.md) | 64 | 1 | 0 | 4 plumbing |
| ROM, control, estimation, campaigns, operations | [ROM_CONTROL_ESTIMATION.md](ROM_CONTROL_ESTIMATION.md) | 40 | 22 | 15 | 3 plumbing |
| Solver, systems, archives, jobs, performance | [SOLVER_SYSTEMS_JOBS.md](SOLVER_SYSTEMS_JOBS.md) | 54 | 25 | 33 | 3 plumbing |
| Python Viewer, tools, side-cars, scripts, case studies | [VIEWER_TOOLS_EXPERIMENTS.md](VIEWER_TOOLS_EXPERIMENTS.md) | 63 | 46 | 34 | 28 presentation, 9 plumbing |

Added after the audit:
- The LCS (LCStool) and SSM (SSMtool, SSMLearn) adapters held no algorithm. They become native
  implementations: finite-time Lyapunov exponent fields and coherent-structure extraction (Haller
  2015; Shadden, Lekien & Marsden 2005), and spectral-submanifold computation and data-driven
  reduced models (Haller & Ponsioen 2016; Jain & Haller 2022; Cenedese et al. 2022).
- The MATLAB parity scripts (`scripts/matlab_parity`, `tools/reference_checks`, the Colbrook
  `check_native_reference.m`) become native tests against frozen reference data (done in phase I).

## Phases

Each phase ends with its native methods tested and documented, then the Python copies it replaces
deleted, the build and both test suites green, and a commit.

| Phase | Scope |
|---|---|
| A. Archives and I/O | read `time_grid_ref` archives; LZMA writes and `archive recompress`; `DYNSYSLAB_ARCHIVE_DIR`; audit reports and campaign scope; storage locations, verified copy, relocation |
| B. Solver and systems | Marmottant shell laws; tabulated, command-tape, telegraph-multisine and symmetric-pulse inputs; RK23, DOP853, Radau, LSODA, CVODE-Adams; output sample times and PCHIP/nearest resampling; dimensional integration; integrator guards; preset metadata; vector and additive actuation; checkpoint, pause, resume and extend; delay rollout; PDE and DDE extensions |
| C. Jobs and queue | `retry_limit`; memory and thread limits from jobs and the CLI; dependency conditions, arrays, admission, pause watchdog; reuse and extension; campaign registry, axes and profiles; benchmarks and gates; telemetry; queue daemon and remote submit; CLI `template` |
| D. Identification | done: pooled H-DMD; multi-trajectory and domain-trained SINDy; SINDy options (rescale, sparsification, directional derivative, non-polynomial dictionaries) and vector-field SINDy ridge, weights and weak form; HAVOK derivative order; readout transfer; observable-spectrum comparison; eigenvalue summaries; spectral weight; SVD driver; archive coordinates; IC files; delay taps at physical-time lags and delay jets ([IDENTIFICATION.md](../IDENTIFICATION.md)). Stage timing and memory preflight came with phase C; POD amplitude interpolation goes with phase E's POD-Galerkin ROM; readout features (asinh k-NN) with phase F |
| E. K-ROM, control, estimation | held-out reconstruction gate; symmetric-pulse input maps; training-domain policy; classical estimators; mode-selection rules; control-error comparison; deeper control search; training-plan requirements; multi-input control; PCHIP, expression and lifted targets; SDRE factorisations; `krom_rollout`; KLQE design reuse; POD-Galerkin ROM; eigenfunction refit, gauge and pruning; validation windows; `control_score` metrics; KLQE subset assessment; control-map path diagnostics |
| F. Experiments and analyses | done: `sweep` and `select` experiment stages; cross-orbit injectivity and closure (`orbit_geometry`); multi-receiver echoes, neighbour/affine/linear readouts and ring errors (receiver side-car deleted); side-car kernels (affine A, c and exact propagation, prehistory, clock audit, run invariant, high-order differences, crossing refinement and visit matching, pooled and windowed EDMDc with held-out scores); model library families; case-study programs rewritten as experiment jobs (forced EDMDc, its validation, manual and retrained PDF, joint H-DMD/KLQE). Fractional lags and delay jets: phase D (`delay_taps`). Remaining (VIEWER_TOOLS_EXPERIMENTS.md): locked-basis extension, endpoint decoder, tangent readout, compact/local/paired/structured fields, waveform calibration (PARTIAL/GAP), echo_rate_fit (needs `sindy_dynamics` derivative options), and the presentation programs of koopman-feedback-benchmark and rpe-corrected-echo (phase H) |
| G. LCS and SSM | done: native `ftle` (FTLE fields, ridges, variational LCS) and `ssm` (parameterization method, backbone curves, data-driven SSMs) methods; `postprocessing/lcs` and `postprocessing/ssm` deleted ([LCS_AND_SSM.md](../LCS_AND_SSM.md)) |
| H. Viewer and reports | Viewer gaps and payload builders; `figure_dataset` operations; experiment report composer; SVG/HTML and multi-page report export |
| I. MATLAB removal | done: frozen-reference native tests for every parity script; all `.m` files deleted |
| J. Python retirement | delete each Python package once nothing it offers is missing; final thin client; docs |

## Phase I: MATLAB removal

The repository holds no MATLAB code. Each retired script was a wrapper that exported results from an
external MATLAB code for comparison. The native behaviour it checked is now held by a native test.

| Retired script | What it checked | Native test now |
|---|---|---|
| `tools/.../colbrook/check_native_reference.m` | native mpEDMD and Rigged DMD against the pinned upstream routines on a rotation | `dynsyslab_colbrook_reference_tests`: frozen upstream operator, order-1/4 densities and packets (`workflow/native/tests/data/colbrook_rotation_reference_v1.json`, 8 KB) and the exact two-atom answer |
| `tools/reference_checks/compare_rigged_reference.m` | native fit on whitened data versus upstream; whitening Gram audit | same test: whitening changes only coefficients; Gram defect below 1e-12 |
| `tools/reference_checks/rigged_history_reference.m`, `test_rigged_history_reference.m` | block-QR compression of weighted history rows before fitting | same test: compressed and raw weighted fits agree to 1e-10 |
| `scripts/matlab_parity/export_dynsyslab_rpe_hdmd_sindy.m` | RPE solve, Hankel construction, H-DMD spectrum, SINDy eigenfunctions | `dynsyslab_rpe_hdmd_reference_tests`: checkpoint parameters and natural-period scaling, energy balance on the canonical trajectory, analytic linear ringdown, exact Hankel entries, analytic H-DMD eigenvalues and forecast. SINDy eigenfunctions: `dynsyslab_sindy_tests` (analytic eigenfunctions) |
| `scripts/matlab_parity/run_matlab_rpe_klqr_control_experiment.m` | KLQR with a state-dependent input map, re-solved at every call, from the equilibrium IC | `dynsyslab_rpe_hdmd_reference_tests` (closed-form Riccati, controlled IC, stabilizing RPE gain); re-solve per changing map: `dynsyslab_coordinate_runtime_tests`; `dynsyslab_koopman_control_tests` |
| `scripts/matlab_parity/export_korda_mezic_duffing_control_candidate.m` | Korda–Mezić boundary eigenfunctions for given eigenvalues, multistep constant-B regression, zero-order-hold conversion | `dynsyslab_boundary_eigenfunction_tests`, `dynsyslab_koopman_control_tests` (multistep fit and exact ZOH factor). The source's RK4 slip and `0.5*u` input are not native conventions and are not reproduced |
| `scripts/matlab_parity/export_playground_duffing_lambdaopt_sindy_regression.m` | eigenvalue optimization, SINDy, regressed B matrices on 30 Duffing trajectories | `dynsyslab_eigenvalue_optimization_tests`, `dynsyslab_eigenvalue_optimization_postprocess_tests`, `dynsyslab_sindy_tests`, `dynsyslab_koopman_control_tests` |
| `scripts/matlab_parity/export_legacy_random_forced_rpe.m` | flattening a saved legacy MATLAB object into arrays (no computation) | none needed; legacy data enters through the native CSV importer, `dynsyslab_cli_legacy_import_contract` |

The pinned upstream repositories stay documented in
[`docs/provenance/colbrook/`](../../provenance/colbrook/README.md) (pinned manifest) and
`docs/provenance/`. The removed scripts remain in Git history.

Open: the SINDy, Korda–Mezić and KLQR chains are covered by analytic tests of their parts but have
no end-to-end frozen RPE test. Such a test would read Koopman model types and archives, which are
being merged into one `krom` object, so it should be written against that object once it lands.
