# Suite architecture audit, 2026-11-22

This audit covers the tools added since the 2026-08-31 audit: the `dynsyslab_math` library, the
`dynsyslab_lie` library and the jet closure, the windowed EDMDc/K-ROM export, and the KLQG and
certified-control studies.

**Inputs:**
- `python -m dynsyslab audit --format markdown`;
- a review of `native/CMakeLists.txt`, `native/src`, `native/include` and the new docs.

Codex's uncommitted work in the tree was read, not changed.

## Healthy

- The audit reports no boundary violations and no package dependency cycles.
- No native file is over 1,000 lines. Three are over 500: `Pipeline.cpp` (608), `Postprocess.cpp`
  (629) and `NumericalDifferentiation.cpp` (669).
- `dynsyslab_math` and `dynsyslab_lie` are separate libraries. Each has its own test executable,
  versioned JSON queries and a CLI entry. `dynsyslab_math` depends only on foundation (plus private
  Boost multiprecision), which is the intended shape. Its documentation states every limit
  (`MATHEMATICS_TOOLKIT.md`).

## Problems worth fixing

1. **The audit counted git-ignored research. Fixed 2026-11-22.**
   - All 21 "hard size violations" were vendored mpmath files under
     `research/pendulum-koopman-eigenfunctions/_vendor/`.
   - `project_paths.source_files` and the audit now skip the side-research folder
     (`SIDE_RESEARCH_DIRECTORIES`, with a regression test).
   - Strict status went from DEBT REMAINS to **PASS**:
     - Python: 451 files and 57k lines, down from 560 and 118k;
     - preferred-size exceptions: 9, down from 45;
     - public-definition docstring coverage: 72.6%, up from 48.4%.
2. **73 docs are unindexed,** including every report since November and the KLQG and RPE
   certified-control records.
   - Either index them in `docs/README.md` or `docs/STUDY_RECORDS.md`, or let the audit accept a
     per-folder index.
3. **The closure family is spread across four places with no single map:**
   - `CoordinateClosure` (core);
   - `lie/` (its own library);
   - `jet/GlobalJetClosure` (a core public header, with sources under `native/src/jet/`);
   - `KoopmanDictionaryFit`/`SpectralDictionary`.

   Two changes would help:
   - Put the jet closure in `dynsyslab_lie`, or in its own `dynsyslab_jet`, so every certificate
     library sits at the same level.
   - Add one short architecture page saying which closure question each module answers (coordinate
     pairing, Lie-algebra cascade, jet/global certificate, dictionary invariance).
4. **`native/src` is still flat,** with 224 files.
   - The new `math/`, `lie/` and `jet/` subfolders are the right pattern.
   - Next groupings: `eigen_numerics/` (22 `EigenNumerics*` files), `analytic/` (`Analytic*`),
     `knowledge/` (`*Knowledge*`, `ResearchPapers*`), `history/` (`History*`, `Hankel*`) and
     `estimation/`.
   - This is a mechanical move. Do it only when no other agent has uncommitted native changes; at
     present Codex has 31 modified files.

## Sensible gaps

1. **`dynsyslab_math` is not connected to the spectral pipeline.** Only the CLI and its tests use
   it. Two small, generic connections would pay off immediately:
   - **Eigenvalue-lattice diagnostics for fitted K-ROMs.** Flag eigenvalues that are near-integer
     combinations of principal ones (λ ≈ mλ₁, using the existing detuning and candidate search).
     This is exactly the harmonic obstruction that breaks KLQR at unstable targets (KLQG report,
     finding 2). It would also make `krom_extract`'s `principal` rule less heuristic.
   - **Resonance reports on Hankel and Spectral-KROM frequencies,** with the same rule: candidates,
     not rationality claims.
2. **A generic period-function tool for one-degree-of-freedom Hamiltonians.**
   - It would give Ω(E), the action I(E), energy inversion Ω⁻¹ and separatrix asymptotics from a
     declared H.
   - The orbit-family research currently relies on a private Python elliptic oracle. A generic
     native quadrature version, with no private derivations, would let Viewer and the pipeline
     compute orbit families reproducibly.
   - The math library's `Distributions` (coarea, spectral density) already expects these roots
     from the caller.
3. **Reeb-graph or separatrix extraction.** `Topology.hpp` takes a caller-supplied graph. Building
   one from a sampled energy function, with lower-star persistence, is the missing step between H
   and the existing graph tools.
4. **There is no LMI/SDP solver in the native core.**
   - The KLQG margins, the SafEDMD-style certificate and the robust designs in
     the private research record run in Python research scripts, using cvxpy.
   - Policy puts algorithms in the C++ core.
   - The practical route is an optional dependency on a small C conic solver such as SCS (MIT
     licence), behind a `ControlDesign` interface.
   - Then: a residual-sector estimator next to `mode_consistency` in `edmdc`, and a certified-design
     query.
5. **Bilinear K-ROMs as a first-class model type.**
   - `edmdc` fits a bilinear term only for the KLQG harness and as a scheduled B_ψ.
   - The certified-control results show the guarantees need A, B₀ and Bᵢ learned from
     constant-input data with ψ(target) = 0.
   - A `bilinear` representation in `input_map`, with this structure, would serve KLQE, KLQR and
     certification together.

## Housekeeping

- **Python helpers.** With research excluded, C++ is the larger codebase: 82k lines against 57k
  of Python. The new Python additions (`dynsyslab_client/mathematics.py`) are thin clients, as
  policy requires.
- **Private research isolation.** `research/` is git-ignored and Codex marks its dossier private.
  Keep the vendored mpmath there, but exclude it from the audit (fix 1).
- **Uncommitted work.** There is a large uncommitted Codex change set across native, workflow and
  docs. Committing it in focused pieces before any reorganization (problem 4) avoids conflicts
  between the agents.
