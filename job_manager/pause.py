"""Durable pause/resume transitions, separate from scheduler orchestration."""
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
import json
import os
from time import sleep
from typing import Any, Callable, Iterable

def pause_jobs(ids:Iterable[str], status:Callable[[str],dict[str,Any]], snapshot:Callable[[str],str], set_status:Callable[...,None], event:Callable[...,None], stamp:Callable[[],str])->dict[str,list[str]]:
    """Request durable pause transitions at safe scheduler boundaries."""
    paused=[]; draining=[]
    for job_id in ids:
        row=status(job_id); immutable=snapshot(job_id)
        if row.get("status") in {"pending","blocked","failed"}:
            set_status(job_id,"paused",pause_requested_utc=stamp(),pause_snapshot_manifest=immutable,incomplete=True); event(job_id,"paused",snapshot_manifest=immutable); paused.append(job_id)
        elif row.get("status") in {"running","claiming"}:
            generation=int(row.get("pause_generation",0))+1
            set_status(job_id,row["status"],pause_after_current_solve=True,pause_generation=generation,pause_requested_utc=stamp(),pause_snapshot_manifest=immutable,incomplete=True); event(job_id,"pause_requested",mode="safe_boundary",pause_generation=generation,snapshot_manifest=immutable); draining.append(job_id)
    return {"paused":paused,"draining":draining}

def resume_jobs(ids:Iterable[str],status:Callable[[str],dict[str,Any]],set_status:Callable[...,None],event:Callable[...,None],stamp:Callable[[],str])->list[str]:
    """Return paused jobs to pending state without changing their plans."""
    resumed=[]
    for job_id in ids:
        if status(job_id).get("status")=="paused": set_status(job_id,"pending",resumed_utc=stamp()); event(job_id,"resumed"); resumed.append(job_id)
    return resumed

def _latest_telemetry(telemetry_root: Path, job_id: str) -> dict[str,Any] | None:
    newest=None
    for path in telemetry_root.glob("*.telemetry.jsonl") if telemetry_root.is_dir() else ():
        try:
            for line in path.read_text(encoding="utf-8").splitlines():
                event=json.loads(line)
                if event.get("job_id")!=job_id: continue
                if newest is None or str(event.get("time_utc",""))>str(newest.get("time_utc","")): newest=event
        except (OSError,json.JSONDecodeError):
            continue
    return newest

def _pid_alive(pid: Any) -> bool | None:
    if not isinstance(pid,int) or pid<=0: return None
    if os.name=="nt":
        # ``os.kill(pid, 0)`` is not a harmless existence probe on Windows.
        # Query the process handle instead, without signalling it.
        import ctypes
        kernel32=ctypes.windll.kernel32
        handle=kernel32.OpenProcess(0x1000,False,int(pid))  # PROCESS_QUERY_LIMITED_INFORMATION
        if not handle:
            return False
        kernel32.CloseHandle(handle)
        return True
    try: os.kill(pid,0)
    except ProcessLookupError: return False
    except PermissionError: return True
    return True

def _pid_cpu_seconds(pid: Any) -> float | None:
    """Read cumulative CPU time without signalling or suspending a process."""
    if not isinstance(pid,int) or pid<=0: return None
    if os.name=="nt":
        import ctypes
        from ctypes import wintypes
        kernel32=ctypes.windll.kernel32; handle=kernel32.OpenProcess(0x1000,False,int(pid))
        if not handle: return None
        try:
            creation=wintypes.FILETIME(); exit_time=wintypes.FILETIME(); kernel=wintypes.FILETIME(); user=wintypes.FILETIME()
            if not kernel32.GetProcessTimes(handle,ctypes.byref(creation),ctypes.byref(exit_time),ctypes.byref(kernel),ctypes.byref(user)): return None
            def ticks(value): return (int(value.dwHighDateTime)<<32)+int(value.dwLowDateTime)
            return (ticks(kernel)+ticks(user))*1e-7
        finally: kernel32.CloseHandle(handle)
    try:
        fields=Path(f"/proc/{int(pid)}/stat").read_text(encoding="utf-8").split()
        return (float(fields[13])+float(fields[14]))/float(os.sysconf("SC_CLK_TCK"))
    except (OSError,ValueError,IndexError): return None

def sample_worker_activity(ids:Iterable[str], *, status:Callable[[str],dict[str,Any]], sample_seconds:float=2.) -> dict[str,dict[str,Any]]:
    """Take a bounded, read-only CPU activity sample for local workers."""
    duration=max(0.,float(sample_seconds)); selected=list(ids)
    first={job_id:_pid_cpu_seconds(status(job_id).get("worker_pid")) for job_id in selected}
    if duration: sleep(duration)
    report={}
    for job_id in selected:
        row=status(job_id); second=_pid_cpu_seconds(row.get("worker_pid")); before=first[job_id]
        delta=None if before is None or second is None else max(0.,second-before)
        report[job_id]={"worker_pid":row.get("worker_pid"),"sample_seconds":duration,"cpu_seconds_before":before,"cpu_seconds_after":second,"cpu_seconds_delta":delta,"actively_computing":None if delta is None else delta>0.}
    return report

def watch_paused_jobs(ids:Iterable[str], *, status:Callable[[str],dict[str,Any]], set_status:Callable[...,None], event:Callable[...,None], stamp:Callable[[],str], telemetry_root:str|Path, grace_seconds:float=120.) -> dict[str,dict[str,Any]]:
    """Classify a requested pause without killing or duplicating active work.

    A live worker remains protected by its claim.  An unresponsive/exited
    worker is quarantined rather than released, so recovery requires explicit
    review and cannot duplicate a still-running scientific specification.
    """
    root=Path(telemetry_root); now=datetime.now(timezone.utc); report={}
    for job_id in ids:
        row=status(job_id)
        if not row.get("pause_after_current_solve"):
            continue
        latest=_latest_telemetry(root,job_id); alive=_pid_alive(row.get("worker_pid"))
        age=None
        if latest and latest.get("time_utc"):
            try: age=max(0.,(now-datetime.fromisoformat(str(latest["time_utc"]).replace("Z","+00:00"))).total_seconds())
            except ValueError: pass
        if row.get("status")=="paused": classification="paused_cleanly"
        elif alive is True and age is not None and age<=grace_seconds: classification="finishing_current_trajectory"
        elif alive is True: classification="rogue_suspected"
        elif alive is False: classification="exited_without_pause_acknowledgement"
        else: classification="awaiting_worker_heartbeat"
        details={"classification":classification,"worker_pid":row.get("worker_pid"),"worker_alive":alive,"worker_cpu_seconds":_pid_cpu_seconds(row.get("worker_pid")),"telemetry_age_seconds":age,"latest_telemetry":latest,"grace_seconds":float(grace_seconds),"pause_generation":row.get("pause_generation")}
        if classification in {"rogue_suspected","exited_without_pause_acknowledgement"}:
            set_status(job_id,"quarantined",pause_watchdog=details,incomplete=True,quarantined_utc=stamp())
            event(job_id,"pause_watchdog_quarantined",**details)
        else:
            set_status(job_id,row.get("status","running"),pause_watchdog=details)
            event(job_id,"pause_watchdog",**details)
        report[job_id]=details
    return report
