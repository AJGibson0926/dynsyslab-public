# Viscous RPE ring-down: corrected intrinsic monopole

[Ring-down, observable equations, SINDy and legacy campaign inventory](../../../../docs/reports/RPE_CORRECTED_MONOPOLE_RINGDOWN_2026-11-05.html) · [Joint HDMD and echo-only KLQE](../../../../docs/reports/RPE_CORRECTED_MONOPOLE_KLQE_2026-11-05.html).

The intrinsic observable is S(2 r v² + r² a₀), where S = R₀ P* and a₀ is the zero-input dimensionless acceleration. It equals physical distance times monopole pressure for the unforced motion studied here. Historical September/October assets named “echo” can use S(2 r v² + r a₀). An identical system parameter fingerprint does not make these sensors interchangeable. An explicit saved measurement or serialized expression keeps its original meaning.

All five input solutions contain physical state samples, independent of the echo formula. Keep the existing immutable solutions. All fresh scientific outputs and figures remain local; the case contains only specifications and pointers. The work date is November 5, 2026, local time; filenames also carry UTC timestamps. A standard installation does not bundle the named archives.

1. In the source checkout, run the compiled dynsyslab-cli with arguments **run JOB --archive-root archives**. Run 01-coordinates for state, analytic echo jet, and scalar echo (three outputs per IC). Run 02-history-basis on the corrected scalar of IC 0. Run 04-history-coordinates using that frozen basis on the five physical records; 05-projection independently checks projection of the training record.
2. Run the 03-fit, 03-broad and 03-matched-broad jobs against the corresponding new coordinate archives. Train only on IC 0 and validate on ICs 1 and 4. Compare finite-difference and stored derivative labels as specified. The 06-state jobs reproduce the original state-only ablations and the compact-dictionary embedding check. When replaying, replace dated output paths with those printed by the new upstream job; match roles and ICs, not filenames alone.
3. Run joint-large, joint-rank600 and joint-rank80 on the raw training solution with the corrected native library loaded. For diagnostics use the existing [joint experiment](../rpe-joint-hdmd-klqe-20261105/PROTOCOL.md) tools: inspect; modal (large); screen; full-real and assess (80); run with zero initial estimate, Q = 0.0001 I and the saved measurement noise setting. Each destination must be new. Keep an observer's candidate, Cx, Cy and measurements on the same sensor version. Use the separate validation record for the held-out observer.
4. For SINDy Viewer panels, use DynSysLabViewerCli render with --source FIT_ARCHIVE --view sindy_dynamics_transient. For the ring-down figure use the local CMake recipe with DYNSYSLAB_SOURCE_DIR and DYNSYSLAB_SUITE_BUILD, then ringdown_figure STATE ECHO_JET HISTORY_COORDINATES BASIS END_TIME OUTPUT.png. This assembles existing Viewer components; no fitting or state evolution occurs in the renderer.

[Local artifact manifest](local-assets.json). Additional machine-readable outputs (expanded coefficient tables, full fit metrics, closure witnesses, physical/projection checks, matched historical sensor probes and observer records) live in archives/experiments/rpe-corrected-echo-20261105. Do not overwrite older campaigns or publish archives with source.

## Two hydrophones

The [two-receiver job](09-two-receiver-monopole.json) exports the leading
monopole pressure for the training trajectory and held-out trajectory 4 at
5 and 6 mm. Channels share the emission clock. Convert to a common receiver
clock using the per-channel travel time recorded in `acoustic_receiver`.

## Echo derivatives from physical state

[Direct coefficient comparison](../../../../docs/reports/RPE_ECHO_RATES_FROM_PHYSICAL_STATE_2026-11-05.html). Build the case target **echo_rate_fit** with the existing CMake configuration. From the checkout root, place the suite bin directory on PATH and run **echo_rate_fit examples/native/case-studies/rpe-corrected-echo-20261105/07-state-to-echo-rates.json NEW_RESULT.json**. It calls the native SINDy fitter with (r,v) as predictors and the time derivatives of the corrected (e,e-prime) series as responses. The report identifies the scored output of each supplied-support fit, the per-output fits from the full library, and the local results.
