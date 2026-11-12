# dynsyslab coordination-refactor evaluation — 2026-07-22

## Scope and evidence

This review covers the coordination path: JobManager state, archive access,
Viewer session handling, packaging exclusions, and regression verification.
Evidence was collected from the current source tree, a full automated test
run, a wheel build, and static boundary scans.

## Completed changes

- JobManager coordination state is now stored transactionally in
  `job_manager/state.sqlite` using SQLite WAL mode. Numerical arrays and Job
  manifests remain ordinary files. A one-time, non-destructive legacy
  `state.json` import is retained solely for migration.
- `ArchiveGateway` now supplies metadata-only reads, solution discovery,
  partial-checkpoint discovery, hierarchy discovery, manifest reads, and
  coordination-state reads. Viewer and JobManager access archive records
  through this facade.
- Archive inspection is isolated in `job_manager/archive_inspection.py`;
  campaign reporting and original-plan conformance checks are isolated in
  `job_manager/campaign_reporting.py`; Viewer session and passive-refresh
  coordination are isolated in `viewer/session.py`.
- Generated output is excluded through `.gitignore`, `MANIFEST.in`,
  setuptools package exclusion, and `project_paths.source_files()`. The
  existing `build/lib` directory is therefore not treated as source or
  included in built wheels.
- The Archiver boundary is used by Campaigner, Solver, Postprocessing, GUI,
  Viewer, JobManager, and Viewer’s private compute/session helpers for archive
  root resolution and archive operations.

## Verification

- Full suite: **118 tests passed** on 2026-07-22.
- Packaging: a wheel built successfully with the bundled Python runtime;
  inspection found **0** `build/` entries.
- Source discovery: 142 Python source files, with generated directories
  excluded.
- Static scan: no `ARCHIVE_ROOT` imports outside the `archives` and
  `archiver` ownership layers.

## Architecture assessment

The coordination design is materially safer: state changes are transactional,
archive policy is centralized, and Viewer remains read-only relative to the
primary archive. The metadata APIs avoid loading numerical arrays when
JobManager is only comparing compatibility or reporting progress. The
generated build tree is now inert to packaging and source-oriented tooling.

The two remaining large files are deliberate rendering/dispatch facades:
`viewer/core.py` still owns graphics construction, while
`job_manager/manager.py` still owns the public execution API. Future growth
should move additional bounded concerns (for example, transition/recovery
planning and result relocation) into services alongside the new reporter and
inspector, rather than adding more behavior to either facade. SQLite currently
stores one atomic coordination payload; if concurrent remote backends become
common, normalize Jobs, claims, and events into separate SQLite tables to
reduce whole-payload rewrite contention.
