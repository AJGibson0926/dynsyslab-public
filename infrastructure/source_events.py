"""Passive source-change notifications without a dependency on Viewer.

Writers publish only tiny JSON events.  Viewer subscribes and polls them; the
primary numerical archive is never modified by this module.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any, Iterable

from project_paths import PROJECT_ROOT


ROOT = PROJECT_ROOT / "viewer_compute_archive" / "subscriptions"


def _stamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _normal(path: str | Path) -> str:
    return str(Path(path).expanduser().resolve())


def register_sources(session_id: str, sources: Iterable[dict[str, Any]]) -> Path:
    """Register the source records observed by one disposable UI session."""
    ROOT.mkdir(parents=True, exist_ok=True)
    target = ROOT / f"{session_id}.json"
    target.write_text(json.dumps({
        "session_id": session_id,
        "updated_utc": _stamp(),
        "sources": list(sources),
    }, sort_keys=True), encoding="utf-8")
    return target


def notify_source_update(source_path: str | Path, *, source_type: str,
                         reason: str) -> int:
    """Publish a refresh event to sessions that observe ``source_path``."""
    if not ROOT.is_dir():
        return 0
    updated = _normal(source_path)
    delivered = 0
    for subscription in ROOT.glob("*.json"):
        try:
            payload = json.loads(subscription.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        matches = [item for item in payload.get("sources", [])
                   if item.get("path") == updated or
                   (source_type == "solution" and
                    updated in item.get("solution_paths", []))]
        if not matches:
            continue
        event_path = ROOT / f"{payload['session_id']}.events.jsonl"
        with event_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({
                "event": "refresh",
                "time_utc": _stamp(),
                "source_path": updated,
                "source_type": source_type,
                "reason": reason,
                "source_names": [item.get("name") for item in matches],
            }) + "\n")
        delivered += 1
    return delivered


def poll_events(session_id: str, *, after: int = 0
                ) -> tuple[list[dict[str, Any]], int]:
    """Return events after an opaque line cursor and the next cursor."""
    path = ROOT / f"{session_id}.events.jsonl"
    if not path.is_file():
        return [], after
    lines = path.read_text(encoding="utf-8").splitlines()
    events: list[dict[str, Any]] = []
    for line in lines[after:]:
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    return events, len(lines)
