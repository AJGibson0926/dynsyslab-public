"""Read-only status aggregation for JobManager state."""
from __future__ import annotations
from typing import Any, Iterable

def progress(rows: dict[str,dict[str,Any]], job_ids: Iterable[str] | None=None) -> dict[str,Any]:
    selected=list(job_ids or rows.keys()); entries=[rows[item] for item in selected if item in rows]
    states=("pending","running","completed","failed","blocked","paused","interrupted")
    counts={state:sum(item.get("status")==state for item in entries) for state in states}
    complete=[item for item in entries if item.get("status")=="completed" and float(item.get("duration_seconds",0.))>0]
    units=[max(1,int(item.get("planned_system_count",1))) for item in complete]
    durations=[float(item["duration_seconds"]) for item in complete]
    mean=None if not durations else sum(durations)/sum(units)
    remaining_units=sum(max(1,int(item.get("planned_system_count",1))) for item in entries if item.get("status") in {"pending","running","blocked","interrupted"})
    return {"total":len(entries),"counts":counts,"fraction_complete":0. if not entries else counts["completed"]/len(entries),"mean_completed_system_seconds":mean,"remaining_planned_systems":remaining_units,"estimated_remaining_seconds":None if mean is None else mean*remaining_units}
