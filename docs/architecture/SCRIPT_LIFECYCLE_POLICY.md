# Script and scientific-package lifecycle policy

DynSysLab keeps executable code only when it remains part of a supported
workflow. A completed experiment is evidence; it is not automatically a
permanent software feature.

## Permanent code

Reusable scientific capability belongs in an importable package with a small
typed API, tests, and—when useful for automation—a lightweight CLI. Permanent
scripts are limited to thin entry points for:

- standard training or control-identification data generation;
- general archive import/export and KROM error evaluation;
- bounded, repeatable research campaigns whose policies and evaluators live in
  reusable packages;
- report generation from a supported campaign schema; and
- repository operations such as the agent snapshot and Antenna deployment.

These entry points must not contain the scientific implementation. They parse
arguments, call a package API, and report paths or compact summaries.

## Evidence-only work

One-off MATLAB wrappers, language-parity checks, ad-hoc plotting commands,
single-target controller reruns, and temporary diagnostics are removed after
their conclusions are captured. Reproducibility is retained through:

- immutable numerical archives and experiment manifests;
- versioned JSON metadata and parameter snapshots;
- case-study reports containing the method, inputs, results, and conclusions;
- canonical API/CLI commands for the reusable operations; and
- tests for any behavior promoted into the supported suite.

Dated case studies may quote a historical script name. Such a reference is an
audit trail, not a current command. New documentation must point readers to the
supported UI/API/CLI instead of reviving a retired task script.

## Promotion and retirement gate

Before retaining experiment code, ask whether another study can configure and
reuse it without editing source. If yes, extract the capability into the
appropriate package and keep only a thin entry point. If no, preserve the
evidence and remove the executable. A retired script must also be removed from
console entry points, tests, current operating guides, and package exports.

This policy keeps the repository focused on the active scientific system while
preserving the provenance of important positive and negative experiments.
