# Chat snapshot policy

## Purpose

dynsyslab preserves two forms of conversational continuity:

1. Versioned, curated context: `AI_CONTEXT.md` and `AGENT_HARNESS.md` describe
   the suite, its boundaries, and current development intent.
2. Private historical snapshots: source conversation exports are copied into a
   timestamped local archive under `history/chat/` and can later be mirrored to
   S3.

## Why exports are explicit

The chat host does not expose a reliable local raw-transcript file to the
project. The snapshot command therefore accepts a user-provided Markdown or
JSON export. It never fabricates missing turns. Raw exports are ignored by Git
because they can contain private information and grow rapidly.

## Create a local snapshot

The former `scripts/snapshot_chat_context.py` helper has been retired with the
Python application layer (it remains in Git history). Create a snapshot by
hand:

1. Make `history/chat/<UTC timestamp>-<label>/` (for example
   `history/chat/20261126T120000Z-dynsyslab-architecture/`).
2. Copy the chat export and the committed `docs/AI_CONTEXT.md` and
   `docs/AGENT_HARNESS.md` into it, for reproducibility.
3. Record a SHA-256 manifest, for example
   `Get-FileHash -Algorithm SHA256 * | Export-Csv manifest.csv`, and optionally
   compress the folder to a `.zip`.

Snapshot contents under `history/chat/` are git-ignored; never commit raw exports.

## Optional S3 mirror

After AWS credentials and a bucket are available, mirror the snapshot folder
explicitly, for example
`aws s3 cp history/chat/<folder> s3://your-bucket/dynsyslab/chat-history/<folder> --recursive`.
Nothing is uploaded unless you run that command.
