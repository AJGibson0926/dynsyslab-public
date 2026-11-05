"""Durable JobManager state, claims, events, and campaign progress."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
from threading import Event, Timer
from time import perf_counter
from typing import Any, Iterable
from uuid import uuid4
import json
import os
import re
import numpy as np

from hierarchy import Batch, Cluster, System, Trajectory
from .backends import ExecutionBackend, LocalProcessBackend, create_backend
from .claims import claim_specs, release_specs
from .history_policy import resolve as resolve_history_policy
from .manifests import ExecutionMode, JobCollection, JobManifest
from .parallel import evaluate_resources, write_resource_report
from .pause import (pause_jobs, resume_jobs, sample_worker_activity,
                    watch_paused_jobs)
from .status import progress as aggregate_progress
from ._runtime import MANAGER_ARCHIVE, _PauseAtBoundary, _safe, _stamp


class _JobStateMixin:
    """Durable JobManager state, claims, events, and campaign progress."""

    def _read_state(self) -> dict[str, Any]: return self.state_store.read()

    def _write_state(self, state: dict[str, Any]) -> None:
        self.state_store.write(state)

    @staticmethod
    def _spec_key(spec: dict[str,Any]) -> str:
        """Stable identity of a scientific solve, independent of job topology."""
        relevant={key:spec.get(key) for key in ("system_name","preset","parameters","driver","time_span","initial_state","solver","output")}
        return __import__("hashlib").sha256(json.dumps(relevant,sort_keys=True,separators=(",",":"),default=str).encode("utf-8")).hexdigest()

    def _claim_specs(self, job_id: str, specs: tuple[dict[str,Any],...]) -> tuple[bool,list[dict[str,Any]]]:
        """Atomically claim every new solve before a worker can integrate it.

        Claims are intentionally all-or-nothing for a Job: allowing a partial
        job to run while another worker owns one member would break the planned
        hierarchy. Claims owned by a non-running process are reclaimed here;
        crash recovery therefore cannot leave a permanent scheduling lock.
        """
        keys={self._spec_key(spec):str(spec.get("label",job_id)) for spec in specs}
        with self._lock,self._state_file_lock():
            state=self._read_state(); acquired,conflicts=claim_specs(state,job_id,keys,_stamp); self._write_state(state)
        return acquired,conflicts

    def _release_spec_claims(self, job_id: str) -> None:
        """Release completed/failed Job claims and wake dependent waiting jobs."""
        with self._lock,self._state_file_lock():
            state=self._read_state(); removed=release_specs(state,job_id,_stamp); self._write_state(state)
        if removed: self._event(job_id,"spec_claims_released",claim_count=len(removed))

    def _read_side_tray(self) -> dict[str,Any]:
        return json.loads(self.side_tray_path.read_text(encoding="utf-8"))

    def _write_side_tray(self, payload: dict[str,Any]) -> None:
        """Atomically publish reserved plans without placing work in the queue."""
        temporary=self.side_tray_path.with_suffix(self.side_tray_path.suffix+f".{uuid4().hex}.partial")
        with temporary.open("w",encoding="utf-8",newline="\n") as handle:
            handle.write(json.dumps(payload,indent=2,sort_keys=True)); handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary,self.side_tray_path)

    @contextmanager
    def _state_file_lock(self):
        """Cross-process lock for local worker state transitions."""
        self.state_lock_path.touch(exist_ok=True)
        with self.state_lock_path.open("a+b") as handle:
            if os.name == "nt":
                import msvcrt
                handle.seek(0); handle.write(b"0"); handle.flush()
                msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
                try: yield
                finally:
                    # Closing the descriptor releases the Windows byte-range
                    # lock reliably; explicit LK_UNLCK can fail on some
                    # networked/temporary Windows filesystems.
                    pass
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
                try: yield
                finally: fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def _event(self, job_id: str, event: str, **details: Any) -> None:
        row = {"time_utc": _stamp(), "job_id": job_id, "event": event, **details}
        with self._lock:
            with (self.log_dir / f"{job_id}.jsonl").open("a", encoding="utf-8") as handle: handle.write(json.dumps(row, sort_keys=True) + "\n")

    @staticmethod
    def _record_campaign_history(job: JobManifest, *, status: str, manifest_path: str | Path | None=None,
                                 solution_paths: list[str] | None=None) -> None:
        """Maintain the archive-owned campaign catalogue without solving.

        A failure to update a secondary discovery index must never invalidate a
        numerical result or a JobManager state transition, so this is invoked
        best-effort at durable job boundaries.
        """
        campaign_path=job.tags.get("campaign_manifest"); campaign_id=job.tags.get("campaign_id")
        if not campaign_path or not campaign_id or not Path(campaign_path).is_file(): return
        try:
            MANAGER_ARCHIVE.register_campaign(campaign_path,status="active")
            MANAGER_ARCHIVE.record_campaign_job(str(campaign_id),job.job_id,status=status,job_label=job.label,manifest_path=manifest_path)
            if solution_paths: MANAGER_ARCHIVE.record_campaign_solutions(str(campaign_id),job.job_id,solution_paths)
        except Exception as error:
            from archiver.diagnostics import warning_event
            warning_event("job_manager","campaign_history_update_failed",error,job_id=job.job_id,campaign_id=campaign_id,status=status)

    def _campaign_progress(self, job: JobManifest, *, force: bool=False) -> None:
        """Persist 10%-style aggregate updates for opted-in bifurcation campaigns."""
        campaign_path=job.tags.get("campaign_manifest")
        if not campaign_path or not Path(campaign_path).is_file(): return
        campaign=json.loads(Path(campaign_path).read_text(encoding="utf-8")); metadata=campaign.get("metadata",{}); reporting=metadata.get("progress_reporting",{})
        if not reporting.get("enabled") or "bifurcation" not in (str(campaign.get("goal",""))+str(metadata.get("campaign_type",""))).lower(): return
        plans=[Path(item) for item in metadata.get("job_plan_paths",[]) if Path(item).is_file()]
        if not plans: return
        state=self._read_state().get("jobs",{}); total=sum(len(JobManifest.read(path).systems) for path in plans); completed=0; generated=[]
        for path in plans:
            plan=JobManifest.read(path); row=state.get(plan.job_id,{})
            completed+=len(plan.systems) if row.get("status")=="completed" else int(row.get("completed_systems",0))
            generated.extend(row.get("generated_solution_archives",[])); generated.extend(row.get("solution_archives",[]))
        generated_paths=sorted({Path(path) for path in generated if Path(path).is_file() and "hierarchies" not in str(path)})
        generated_bytes=sum(path.stat().st_size for path in generated_paths)
        percent=int(100*completed/total) if total else 100; interval=int(reporting.get("interval_percent",10)); threshold=(percent//interval)*interval; log_path=Path(reporting["log_path"]); log_path.parent.mkdir(parents=True,exist_ok=True)
        running_estimates=[float(row.get("estimated_remaining_seconds")) for row in state.values() if row.get("estimated_remaining_seconds") is not None and row.get("status")=="running"]
        completed_rows=[row for row in state.values() if row.get("status")=="completed" and float(row.get("duration_seconds",0.))>0]
        mean_rate=None if not completed_rows else sum(float(row["duration_seconds"])/max(1,int(row.get("planned_system_count",1))) for row in completed_rows)/len(completed_rows)
        pending_units=sum(len(JobManifest.read(path).systems) for path in plans if state.get(JobManifest.read(path).job_id,{}).get("status") in {"pending","blocked"})
        campaign_remaining=(sum(running_estimates) if running_estimates else 0.) + (0. if mean_rate is None else mean_rate*pending_units)
        emitted=[]; checked=[]
        if log_path.exists():
            for line in log_path.read_text(encoding="utf-8").splitlines():
                try:
                    record=json.loads(line)
                    if record.get("event")=="progress": emitted.append(int(record.get("threshold_percent",0)))
                    if record.get("event")=="archive_progress_check": checked.append(int(record.get("attempted_threshold_percent",0)))
                except json.JSONDecodeError: pass
        # State is a cheap trigger only.  A persistent campaign progress value
        # must be counted from metadata-compatible archived solutions; this is
        # resilient to restarts, job repartitioning, and skipped reusable work.
        needs_check=(threshold and threshold not in emitted and threshold not in checked) or (force and percent==100 and 100 not in emitted)
        if needs_check:
            authoritative=self.campaign_status(campaign_path)
            completed=int(authoritative["archived_complete_trajectories"]); total=int(authoritative["planned_trajectories"])
            percent=int(100*completed/total) if total else 100; threshold=(percent//interval)*interval
            archive_paths=[Path(path) for path in MANAGER_ARCHIVE.campaign_paths(str(campaign.get("campaign_id"))) if Path(path).is_file()]
            generated_paths=sorted(set(archive_paths)); generated_bytes=sum(path.stat().st_size for path in generated_paths)
            if threshold and threshold not in emitted:
                event={"event":"progress","time_utc":_stamp(),"campaign_id":campaign.get("campaign_id"),"completed_trajectories":completed,"planned_trajectories":total,"percent":percent,"threshold_percent":threshold,"solution_file_count":len(generated_paths),"solution_bytes":generated_bytes,"estimated_remaining_seconds":authoritative.get("estimated_remaining_seconds")}
            else:
                event={"event":"archive_progress_check","time_utc":_stamp(),"campaign_id":campaign.get("campaign_id"),"attempted_threshold_percent":threshold,"completed_trajectories":completed,"planned_trajectories":total,"percent":percent,"solution_file_count":len(generated_paths),"solution_bytes":generated_bytes}
            with log_path.open("a",encoding="utf-8") as handle: handle.write(json.dumps(event)+"\n")
            if event["event"]=="progress":
                from infrastructure.source_events import notify_source_update
                notify_source_update(campaign_path,source_type="campaign",reason=f"campaign_progress_{threshold or percent}_percent")

    def _set_status(self, job_id: str, status: str, **details: Any) -> None:
        with self._lock, self._state_file_lock():
            state = self._read_state(); row = state["jobs"].setdefault(job_id, {}); row.update({"status": status, "updated_utc": _stamp(), **details}); self._write_state(state)

    def _record_interruption(self, job_id: str, *, reason: str="execution interrupted") -> dict[str,Any]:
        """Persist recoverable interruption state without discarding evidence."""
        current=self.status(job_id); partial=[item for item in current.get("generated_solution_archives",[]) if Path(item).is_file()]
        self._set_status(job_id,"interrupted",reason=reason,interrupted_utc=_stamp(),generated_solution_archives=partial,
                         generated_solution_bytes=sum(Path(item).stat().st_size for item in partial))
        self._event(job_id,"interrupted",reason=reason,solution_file_count=len(partial),solution_bytes=sum(Path(item).stat().st_size for item in partial))
        try:
            job=self._manifest(job_id); campaign_path=job.tags.get("campaign_manifest")
            if campaign_path and Path(campaign_path).is_file():
                campaign=json.loads(Path(campaign_path).read_text(encoding="utf-8")); reporting=(campaign.get("metadata") or {}).get("progress_reporting",{}); log_path=reporting.get("log_path")
                if log_path:
                    with Path(log_path).open("a",encoding="utf-8") as handle: handle.write(json.dumps({"event":"job_interrupted","time_utc":_stamp(),"campaign_id":campaign.get("campaign_id"),"job_id":job_id,"reason":reason,"solution_file_count":len(partial),"solution_bytes":sum(Path(item).stat().st_size for item in partial)})+"\n")
        except Exception as error:
            from archiver.diagnostics import warning_event
            warning_event("job_manager","interruption_campaign_log_failed",error,job_id=job_id,reason=reason)
        return self.status(job_id)

    def resume_interrupted(self, job_ids: Iterable[str] | None=None) -> list[str]:
        """Return interrupted Jobs to pending state; solver checkpoints remain intact."""
        ids=list(job_ids or self.status().keys()); resumed=[]
        for job_id in ids:
            if self.status(job_id).get("status") == "interrupted":
                self._set_status(job_id,"pending",resumed_utc=_stamp(),resume_of="interrupted")
                self._event(job_id,"resume_queued"); resumed.append(job_id)
        return resumed

    def recover_orphaned_jobs(self) -> list[str]:
        """Mark stale running Jobs recoverable after a manager/process crash."""
        recovered=[]
        for job_id,row in self.status().items():
            if row.get("status") in {"claiming","running"}:
                self._record_interruption(job_id,reason="manager restart detected a running Job without a live handle")
                self._release_spec_claims(job_id)
                recovered.append(job_id)
        return recovered
