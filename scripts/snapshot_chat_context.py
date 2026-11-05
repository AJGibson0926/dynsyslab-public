"""Create a private, reproducible local snapshot of a supplied chat export.

Raw transcript data is never added to Git. It may optionally be mirrored to
an S3 prefix when the AWS CLI and an explicit/configured URI are available.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

PROJECT_ROOT = Path(__file__).resolve().parents[1]
HISTORY_ROOT = PROJECT_ROOT / "history" / "chat"
CONFIG_PATH = PROJECT_ROOT / "config" / "s3_snapshot.json"
CONTEXT_FILES = tuple(
    PROJECT_ROOT / "docs" / name
    for name in ("AI_CONTEXT.md", "AGENT_HARNESS.md", "CHAT_SNAPSHOT_POLICY.md")
)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of one snapshot input."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configured_s3_uri() -> str | None:
    """Resolve the explicitly configured remote snapshot location."""
    if not CONFIG_PATH.exists():
        return None
    with CONFIG_PATH.open(encoding="utf-8") as stream:
        return json.load(stream).get("s3_uri")


def main() -> int:
    """Create a bounded, policy-filtered maintainer-context snapshot."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path, help="Chat export (.md, .json, or .txt).")
    parser.add_argument("--label", default="conversation", help="Short descriptive archive label.")
    parser.add_argument("--s3-uri", help="Optional S3 prefix, e.g. s3://bucket/dynsyslab/chat-history.")
    parser.add_argument("--dry-run", action="store_true", help="Describe the archive without writing or uploading.")
    args = parser.parse_args()
    source = args.source.expanduser().resolve()
    if not source.is_file():
        parser.error(f"Chat export not found: {source}")
    safe_label = "".join(ch if ch.isalnum() or ch in "-_" else "-" for ch in args.label).strip("-") or "conversation"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive_dir = HISTORY_ROOT / f"{stamp}-{safe_label}"
    zip_path = HISTORY_ROOT / f"{archive_dir.name}.zip"
    s3_uri = args.s3_uri or configured_s3_uri()
    if args.dry_run:
        print(json.dumps({"archive_dir": str(archive_dir), "zip": str(zip_path), "s3_uri": s3_uri}, indent=2))
        return 0

    archive_dir.mkdir(parents=True, exist_ok=False)
    copied_source = archive_dir / f"conversation{source.suffix.lower() or '.txt'}"
    shutil.copy2(source, copied_source)
    copied_context: list[Path] = []
    for context in CONTEXT_FILES:
        if context.is_file():
            destination = archive_dir / context.name
            shutil.copy2(context, destination)
            copied_context.append(destination)
    files = [copied_source, *copied_context]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_name": source.name,
        "label": safe_label,
        "files": [{"name": item.name, "sha256": sha256(item), "bytes": item.stat().st_size} for item in files],
    }
    (archive_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (archive_dir / "sha256sums.txt").write_text(
        "\n".join(f"{entry['sha256']}  {entry['name']}" for entry in manifest["files"]) + "\n", encoding="utf-8",
    )
    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED) as archive:
        for item in archive_dir.iterdir():
            archive.write(item, arcname=f"{archive_dir.name}/{item.name}")
    if s3_uri:
        aws = shutil.which("aws")
        if not aws:
            raise RuntimeError("S3 target specified, but the AWS CLI is not installed or not on PATH.")
        target = s3_uri.rstrip("/") + "/" + zip_path.name
        subprocess.run([aws, "s3", "cp", str(zip_path), target], check=True)
        manifest["s3_uri"] = target
        (archive_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"archive_dir": str(archive_dir), "zip": str(zip_path), "s3_uri": manifest.get("s3_uri")}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError, json.JSONDecodeError) as error:
        print(f"snapshot failed: {error}", file=sys.stderr)
        raise SystemExit(1)
