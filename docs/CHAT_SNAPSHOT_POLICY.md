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

```powershell
python scripts/snapshot_chat_context.py `
  --source C:\path\to\chat-export.json `
  --label dynsyslab-architecture
```

This creates `history/chat/<UTC timestamp>-<label>/` plus a compressed `.zip`
archive and SHA-256 manifest. The committed harness/context files are copied
into the snapshot for reproducibility.

## Optional S3 mirror

After AWS credentials and a bucket URI are available, either provide a target:

```powershell
python scripts/snapshot_chat_context.py --source C:\path\to\chat.md `
  --s3-uri s3://your-bucket/dynsyslab/chat-history
```

or copy `config/s3_snapshot.example.json` to the ignored
`config/s3_snapshot.json`. The command requires the AWS CLI and performs no
upload unless `--s3-uri` or configured `s3_uri` is supplied.

