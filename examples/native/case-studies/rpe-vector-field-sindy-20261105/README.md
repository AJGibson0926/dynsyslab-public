# RPE dynamics SINDy in three two-dimensional coordinates

This case study uses the standard native pipeline; no standalone fitting,
coordinate-projection or integration script is required.

Read the rendered report at
[RPE dynamics discovery](../../../../docs/reports/RPE_TWO_COORDINATE_DYNAMICS_SINDY_2026-11-05.html)
and the [API/workflow guide](../../../../docs/architecture/COORDINATE_DYNAMICS_SINDY.md).

## Reproduction

Run each saved job from the checkout root:

    build/native-suite-fast/bin/dynsyslab-cli.exe plan JOB.json
    build/native-suite-fast/bin/dynsyslab-cli.exe run JOB.json --archive-root archives

The source archives are local scientific assets, not bundled sample data.
Paths in these jobs identify the author's original records. On a different
installation, supply equivalent archives and update paths explicitly.

1. 01-coordinates.json applies the existing identity/echo-jet/frozen-history
   mappings to training IC 0 and validation ICs 1 and 4. Exact reference
   derivatives are explicitly requested. No history basis is refitted.
2. 02-fit-state/echo/latent/inverse.json compare numerical and exact derivative
   targets with polynomial or compact inverse-radius libraries.
3. 03-broad-state/echo/latent.json record the original broad experiments, with
   threshold 0.01 and fitting stride 10. They are preserved as evidence, not
   silently replaced with improved results.
4. 04-broad-state-ablation.json holds fitting rows fixed (stride 4) and varies
   pruning and SVD rank tolerance.
5. 05-verify-embedding.json evaluates the compact model after its coefficients
   are embedded in the broader dictionary. This is coefficient reassignment,
   not a fit. Native evaluation checks its equivalence to the original.
6. 06-matched-threshold.json tests threshold 0.001 with the same SVD cutoff and
   samples, plus an unpruned broad fit using exact derivatives.

The broad dictionary contains the inverse-radius law. The original worse
result was mainly caused by sparse support selection: threshold 0.01 discarded
the viscous v/r^2 term. The matched 0.001 run reduces exact acceleration relative
RMS error to 1.45e-6; the compact fit is 3.49e-7. Unpruned exact-rate broad fitting
reaches 1.9e-11. The embedded compact model agrees to approximately 1e-16 on all
three records. Larger dictionaries may fit derivative noise and have unstable
coefficient decompositions; they need separate validation and rollout checks.

## Evidence

The local asset manifest lists archive/PNG locations and SHA-256 hashes.
Full metrics are stored in
archives/experiments/rpe-vector-field-sindy-20261105/metrics.json.
The scientific assets remain local and are not added to Git or installers.

Viewer views are sindy_dynamics and sindy_dynamics_transient, with separate pages
for the training record and two validation records. The selected plot shows the
first completed rollout and explicitly labels its initial time; earlier failed
rollouts remain in the metadata and report captions.

Both the identity derivative and echo/latent models are two-dimensional.
The saved history projection was reproduced bit for bit. All 115 native tests
passed, including original Koopman tests, the new dynamics module, coordinate
transforms, model save/load, validation and Viewer payloads. The final test log
is build/sindy-final-ctest.log. The architecture audit is build/sindy-audit.log.

An archive-index rebuild initially overlapped broad-run publication, causing
SQLite lock errors after the scientific files had been written. Those completed
files were recovered from their metadata and indexed; no result was fabricated
or recomputed to hide the publication issue. Index reconciliation reports zero
failed archives. A separate relative-root publication bug was fixed and covered
by a regression test.
