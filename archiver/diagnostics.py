"""Structured, durable non-fatal diagnostic events for dynsyslab."""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import json
from archives import ARCHIVE_ROOT

def warning_event(component: str, code: str, error: BaseException, **context: Any) -> dict[str, Any]:
    """Append a machine-readable warning without hiding its exception context."""
    row={"time_utc":datetime.now(timezone.utc).isoformat(),"severity":"warning","component":component,"code":code,
         "exception_type":type(error).__name__,"exception_message":str(error),"context":context}
    path=Path(ARCHIVE_ROOT)/"logs"/"warnings.jsonl"; path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("a",encoding="utf-8") as handle: handle.write(json.dumps(row,sort_keys=True,default=str)+"\n")
    return row
