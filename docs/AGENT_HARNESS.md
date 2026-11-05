# dynsyslab agent harness

This file is the durable operating context for an AI or human maintainer. It
is intentionally concise, version-controlled, and contains no credentials or
raw private conversation transcript.

## Mission

dynsyslab is a standalone numerical-dynamics suite. Work only on the user's
current goal; this document does not authorize a future scientific, UI, MCP,
or infrastructure phase by itself.

## Read first

1. `docs/START_HERE.md` — safe human/agent entry point.
2. `ARCHITECTURE.md` — canonical package relationships.
3. `docs/AI_CONTEXT.md` — architecture and invariants.
4. The feature guide or case study relevant to the active request.

Then run `dynsyslab-cli doctor`, `dynsyslab-cli capabilities`, and
`python -m quality --format markdown` when auditing the source tree.

## Working boundaries

- Use `archiver.ArchiveGateway` for archive I/O. Do not add direct application
  package I/O against `archives`.
- Viewer is read-only relative to the primary archive. Its private artifacts
  live only in `viewer_compute_archive`.
- Solver executes one solve request. Campaigner designs manifests and
  JobManager schedules them; neither owns solution arrays.
- Do not delete or overwrite archived numerical results. Extension must use the
  archive extension path and preserve provenance.
- The original RPE chaos-campaign archive is protected. It may be read for
  analysis, but must not be altered without the user's explicit agreed
  authorization in that conversation.
- Keep raw chat exports and credentials out of Git. The curated context files
  in this directory are the versioned memory record.
- MATLAB is read-only evidence unless the user explicitly authorizes a wrapper
  or adaptation. Never infer permission to copy the legacy suite wholesale.

## Active checkpoint

The historical RPE two-parameter campaign is presently paused. Its separate
checkpoint files show five completed screening points and no positive
finite-time Lyapunov screening result. Do not resume it solely because this
file exists; obtain current user authorization.

## Versioning and release practice

- `pyproject.toml` is the canonical package version.
- Use Semantic Versioning: patch for compatible fixes, minor for compatible
  features, major for incompatible API changes.
- Update `CHANGELOG.md`, run `python -m tester --all`, then commit and tag a
  release as `vX.Y.Z`.
- Conversation snapshots are created with
  `scripts/snapshot_chat_context.py`; see `docs/CHAT_SNAPSHOT_POLICY.md`.
