# RPE prior control cases: clean-pathway rerun

**Date:** 2026-08-27  
**Status:** completed experiment set  
**Controller contract:** unchanged target, equilibrium start, dimensionless KLQR, continuous-RHS state dependence, BDF, no feedforward, no cap, no startup shaping

## Results

| Case | Status | r RMSE (second half) | dr/dtau RMSE (second half) | Peak abs u* | Solve time (s) |
|---|---|---:|---:|---:|---:|
| Fixed target (2, 0) | complete | 0.000844141 | 2.89325e-06 | 44.6534 | 13.186 |
| Single sine, offset 2 | complete | 0.00995277 | 0.109097 | 2314.73 | 114.703 |
| Hard mixed multisine, offset 4 | complete | 0.00432634 | 0.0229518 | 63267.5 | 46.341 |
| Small equilibrium-centered multisine | complete | 0.000672826 | 0.00217474 | 0.707291 | 15.609 |
| Offset-1.5 feedback-only waveform | complete | 0.303291 | 0.557932 | 6.20816 | 37.167 |
| Periodic multisine, offset 3.5 | complete | 0.0101856 | 0.125355 | 8583.73 | 149.802 |
| Periodic multisine, offset 4 | complete | 0.00516589 | 0.0491299 | 41016.3 | 90.628 |
| Mixed multisine, offset 2.85 | complete | 0.0104565 | 0.120189 | 6428.67 | 55.830 |

## Conclusions

The corrected clean-isolation KROM/controller pathway resolves the three
explicit historical execution failures:

- Fixed `(2,0)` no longer leaves the bubble near equilibrium. Its second-half
  radius RMSE fell from approximately `1.000039` to `0.000844141`, and its
  final radius is `1.999156`.
- The offset-2 experimental single sine now produces an accepted ten-`tau`
  archive, with second-half radius RMSE `0.00995277`.
- The hard offset-4 mixed multisine now completes with second-half radius RMSE
  `0.00432634`; this reproduces the successful clean-isolation case.

The result is not universal success. The earlier offset-1.5 feedback-only
waveform completes but has second-half radius/velocity RMSE of `0.303291` and
`0.557932`, markedly worse than its prior specialized two-function KROM result.
It starts exactly at the target, so this is not an initialization spike problem:
pure feedback initially commands zero while the moving nonlinear reference
requires forcing, and this energy-selected KROM/error direction does not build
the needed tracking action accurately. KROM suitability remains target- and
control-objective-dependent.

The startup command remains strongly target-dependent. The equilibrium-centered
small multisine has only `0.707291 u*` peak control, while remote or high-
velocity targets range from `2314.73` to `63267.5 u*`. These are retained
results, not production acceptability claims. Several targets also leave the
free-training radius domain, and this comparison deliberately uses the same
explicit `target_domain_policy="ignore"` required by the clean-isolation
baseline.

## Alternative minimal-KROM follow-up

The remaining offset-1.5 failure was rerun with the next distinct low-frequency
complex pair from the same full SINDy archive. Only the eigenfunction pair
changed:

| Minimal KROM | Full-spectrum pair | Continuous eigenvalue | r RMSE, second half | dr/dtau RMSE, second half | Peak abs u* |
|---|---:|---:|---:|---:|---:|
| Current successful KROM | `1/2` | `-0.210616 +/- 4.818635i` | 0.303291 | 0.557932 | 6.20816 |
| Alternative KROM | `3/4` | `-0.125796 +/- 5.596275i` | 0.311481 | 0.580473 | 7.02241 |

The alternative completed but was slightly worse in every tracking measure;
therefore this target remains unresolved by a simple adjacent-pair swap. The
comparison changed neither target nor controller settings.

The same audit exposed a selection-description error: the current pair `1/2`
is retained energy-order position `999`, while alternative pair `3/4` is
position `1000`. The current archive is numerically correct, but earlier text
calling it the first energy-ranked pair was inaccurate.

The complete retained diagnosis, reproduction contract, target-path Jacobian
audit, inverse-pressure comparison, and rejected hypotheses are in
`RPE_OFFSET1P5_MINIMAL_KROM_FAILURE_CASE_STUDY_2026-08-27.md`. Exact
machine-readable values are in
`../notes/RPE_OFFSET1P5_EIGENFUNCTION_FAILURE_DIAGNOSTICS_2026-08-27.json`.

Viewer comparison:
`viewer_outputs/rpe-offset1p5-current-vs-alternative-minimal-krom.html`.

## Interpretation boundary

These runs test whether the corrected clean-isolation representation and controller pathway resolve earlier failures. They do not introduce the optional C2 target engagement, feedforward, clipping, or target-preserving gain modulation. Large initial transients are retained in the archives and handled only by Viewer's automatic detail inset.

Every case has an individual standard Viewer page. The combined dashboard is `viewer_outputs/rpe-prior-control-cases-clean-rerun.html`; machine-readable metrics are in `viewer_outputs/rpe-prior-control-cases-clean-rerun.json`.

## Case provenance

- **Fixed target (2, 0):** explicit prior failure: completed but remained near equilibrium
- **Single sine, offset 2:** explicit prior failure: no accepted solution archive
- **Hard mixed multisine, offset 4:** explicit prior failure: invalid radius/state or prohibitive runtime
- **Small equilibrium-centered multisine:** prior imperfect/chattering case retained for comparison
- **Offset-1.5 feedback-only waveform:** prior feedback-only case with nonzero tracking offset
- **Periodic multisine, offset 3.5:** canonical target coverage rerun
- **Periodic multisine, offset 4:** canonical target coverage rerun
- **Mixed multisine, offset 2.85:** canonical target coverage rerun
