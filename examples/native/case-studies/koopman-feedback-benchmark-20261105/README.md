# Measured-output feedback benchmark

**Critically important reference, designated by the user on 2026-11-05.**
Reference ID: `koopman-linear-duffing-exact`.
See the [complete derivation and retained reference](../../../../docs/case-studies/KOOPMAN_LINEAR_DUFFING_EXACT_REFERENCE.md).
The [local asset manifest](local-assets.json) records SHA-256-verified copies
of the original evidence under
`archives/experiments/koopman-feedback-benchmark-20261105/`, outside the
disposable build directory. Use those preserved paths when revisiting this run.

The standard [controlled](controlled-job.json) and
[uncontrolled](uncontrolled-job.json) jobs are saved alongside this README.
Run them from the checkout root with `dynsyslab-cli run <job> --archive-root <new-output-directory>`.
The controlled job points to the preserved exact KROM. The original build paths
below document how the first result was produced; they are not the sole retained
copies.

This is an exact linear special case of the existing Duffing system, not a
nonlinear RPE or data-driven fitting result. The standard C++ workflow runs
x'' = 2x - x' + u, with position-only measurements y=x. The continuous
eigenfunctions are psi1=2x+v (lambda1=1) and psi2=-x+v (lambda2=-2).
The constant input map [1,1]^T and output maps are supplied analytically.

The reproducible benchmark is implemented in
[the workflow test](../../../../workflow/native/tests/koopman_feedback_tests.cpp).
Passing an output directory to that test retains its fixture KROM, solve job and
result archives. The plot program only reads archives and assembles standard
Viewer figures; it does not implement another solver or controller.

Settings: initial conditions (+0.1,0) and (-0.1,0), zero observer initialization,
control Q=I2 and R=0.1, LQE process intensity I2 and measurement intensity 0.001,
sample period 0.005, output period 0.02, final time 12, input bound 0.8,
control enabled at time 0.1. No measurement noise is injected in this test.
The uncontrolled reference uses the same system, initial conditions and solver,
with the controller disabled.

Local artifacts (not distributed in Git):

- Controlled source manifest: build/feedback-verification/installed-cli-result.json
- Controlled job: build/feedback-verification/feedback-job.json
- Uncontrolled source manifest: build/feedback-verification/uncontrolled-result.json
- Uncontrolled job: build/feedback-verification/uncontrolled-job.json
- Viewer full figure: viewer_outputs/koopman-feedback-benchmark-20261105-full.png
- Viewer transient figure: viewer_outputs/koopman-feedback-benchmark-20261105-transient.png
- Metrics: viewer_outputs/koopman-feedback-benchmark-20261105-metrics.json

The source manifests point to the exact archives used. Both controlled cases
finish with state norm 1.19765e-8 and state-estimation error 2.18521e-10.
The uncontrolled state norm reaches 15344.6689 at time 12.
The two norm curves coincide because these initial conditions are symmetric.

Build the plot helper using its standalone CMake file, setting
DYNSYSLAB_SOURCE_DIR and DYNSYSLAB_SUITE_BUILD to the checkout and native build.
With the suite bin directory on PATH, run from the checkout root:

```powershell
koopman_feedback_plots.exe build/feedback-verification/installed-cli-result.json build/feedback-verification/uncontrolled-result.json viewer_outputs/koopman-feedback-benchmark-20261105
```

Existing images are not overwritten. Choose a new output prefix for another run.
