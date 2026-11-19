---
name: dynsyslab-spectral-analysis
description: Fit and diagnose DynSysLab spectral KROMs with native DMD/Koopman methods, spectral measures, frozen dictionaries, or existing-KROM observable references. Use for these computations and method comparisons in the DynSysLab pipeline.
---

# Spectral Koopman analysis in DynSysLab

This skill is project-scoped: the repository containing it is the DynSysLab
source root. (The Codex copy of this skill lives in `~/.codex/skills`; keep the
two in step when the guidance changes.)
Read `docs/architecture/SPECTRAL_KROM.md` in that root for the maintained API,
implemented-method table, assumptions, and examples. Do not infer implementation
from a method's paper title: some upstream workflows have only reusable primitives
ported, and those distinctions are recorded there.

## Run a computation

Use the compiled Post-Processor or `NativeWorkflow.koopman_spectral`; the latter
is a thin C-ABI adapter. The CLI is `dynsyslab-cli koopman-spectral query.json`.
Scientific fitting belongs to `dynsyslab_spectral`; result/evaluation/archive
contracts belong to `dynsyslab_koopman`. Workflow resolves inputs and publishes
results. Observables in foundation carry passive references and never load KROMs.
Do not recreate production numerical routines in a Python experiment script.

For comparisons, reuse one archived Hankel pair matrix or exactly the same
explicit X/Y pairs. State the observables, normalization, delays, pair shift,
rank, dictionary and training sources. Column stride is not the DMD pair shift.
Fit normalization and centers on training data only, and freeze them for
evaluation. Keep hidden priors, waveform calibration and data-label sources
explicit when comparing with older SINDy experiments.

Report measured-window readout error separately from an autonomous free rollout.
The spectral KROM rollout consumes only the initial coordinate vector and then
applies its learned operator. Initialization in history coordinates requires the
corresponding observed window; a single instantaneous echo/rate pair does not
determine an arbitrary delay vector.

Use held-out snapshot pairs for residual validation when available. Training
interpolation can have a near-zero residual without predictive value. Principal
angles and projection defect diagnose a finite linear trial space, not existence
of an arbitrary nonlinear closed state law. A small residual is not a spectral
distance certificate for a nonnormal operator.

## Choose a spectral interpretation

- Default unconstrained data regression does not assume measure preservation.
- Apply mpEDMD only when the invariant-measure/unitary assumption is justified;
  a damped RPE transient is not such evidence. The API requires explicit acknowledgement.
- Hermitian DMD requires a self-adjoint problem, a different assumption.
- Correlation measures require stationary/invariant sampling; retain total mass,
  grid convention and smoothing parameters. Frequencies are radians per sample.
- Rigged wave packets are finite-width approximations. Higher-order smoothing
  kernels may have negative lobes; do not claim certified distributional limits.
- SpecRKHS results describe the adjoint acting on kernel sections. Do not relabel
  its eigenvectors as a forward eigenfunction model or feed them into a rollout.
- Stochastic diagnostics need independent repeated successors at identical source
  points. Sample-mean error and innovation variation are distinct quantities.

For existing-KROM observables, use the stored model identity and named component.
Require the exact input-coordinate contract; dimension equality alone is insufficient.
A measure can be associated with an observable but is not a pointwise time signal.

## Maintain and verify

For changes, run the focused native spectral test and Python C-ABI/CLI/workflow
tests documented in `SPECTRAL_KROM.md`, then the affected suite and architecture
audit. Optional pinned MATLAB parity compares mpEDMD and Rigged DMD without making
MATLAB a runtime dependency. Source pins, acquisition helper and the reference
fixture live in `tools/experiments/miscellaneous_computes/references/colbrook`.
Keep upstream license notices and distinguish downloaded reference code from
ported, validated native functionality.
