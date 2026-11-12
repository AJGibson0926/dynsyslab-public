"""Read-only access to versioned cross-process coordination state."""
from __future__ import annotations

import json
from pathlib import Path
import sqlite3
from typing import Any


EMPTY_JOB_STATE={"format":"dynsyslab.job-manager-state/v2","jobs":{},"spec_claims":{}}


def read_job_coordination_state(database: str | Path) -> dict[str, Any]:
    """Read JobManager's public coordination DTO without importing its runtime."""
    path=Path(database)
    if not path.is_file():
        return {**EMPTY_JOB_STATE,"jobs":{},"spec_claims":{}}
    connection=sqlite3.connect(path)
    try:
        row=connection.execute("SELECT value FROM state WHERE key='payload'").fetchone()
    finally:
        connection.close()
    if row is None:
        return {**EMPTY_JOB_STATE,"jobs":{},"spec_claims":{}}
    payload=json.loads(row[0])
    if not isinstance(payload,dict):
        raise ValueError("Coordination payload must be a JSON object.")
    payload.setdefault("jobs",{})
    payload.setdefault("spec_claims",{})
    return payload


__all__=["EMPTY_JOB_STATE","read_job_coordination_state"]
