"""Opt-in checkpoint pause/resume helpers for long ODE trajectories."""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
from typing import Any
from archiver.gateway import SOLVER_ARCHIVE

def resume_ode_from_pause(paused_archive:str|Path, request:Any):
    """Continue a paused ODE from its archived terminal state and extend it.

    This deliberately restarts the numerical integrator at the checkpoint; it
    does not claim bitwise continuity of an adaptive solver's private state.
    The original archive is extended and annotated with the restart boundary.
    """
    path=Path(paused_archive).resolve(); series=SOLVER_ARCHIVE.read_time_series(path)
    with SOLVER_ARCHIVE.open_read(path) as data:
        import json
        metadata=json.loads(str(data["metadata_json"]))
    if metadata.get("status")!="paused": raise ValueError("resume_ode_from_pause requires an archive with status='paused'.")
    if request.kind!="ode": raise ValueError("Checkpoint continuation currently supports ODE requests only.")
    original_span=metadata.get("time_span") or request.time_span
    if original_span is None: raise ValueError("Paused archive has no original time span.")
    end=float(original_span[1]); start=float(series.time[-1])
    if start>=end: raise ValueError("Paused trajectory has no remaining integration interval.")
    history=replace(request.history,pause_at_checkpoint=False)
    continuation=replace(request,time_span=(start,end),initial_state=series.values[-1].copy(),history=history,label=f"{request.label}-resume")
    from . import solve
    result=solve(continuation)
    if result.status!="complete" or result.archive_path is None: return result
    SOLVER_ARCHIVE.extend_solution(path,result.archive_path)
    SOLVER_ARCHIVE.update_metadata(path,{"status":"complete","checkpoint_restart":{"continuation_archive":str(result.archive_path),"restart_time":start,"method":"integrator restart from durable archived state"}})
    result.archive_path=str(path); result.metadata["resumed_from_pause"]=str(path)
    return result
