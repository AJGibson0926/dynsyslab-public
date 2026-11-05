"""Manifest submission, reservation, hierarchy, and worker planning services."""
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


class _JobPlanningMixin:
    """Manifest submission, reservation, hierarchy, and worker planning services."""

    def submit(self, job: JobManifest) -> Path:
        from .scheduling import execution_provenance
        path = self.manifest_dir / f"{_safe(job.label)}__{job.job_id}.json"
        original = self.original_manifest_dir / f"{_safe(job.label)}__{job.job_id}.json"
        if not original.exists(): job.write(original)
        job.write(path); decisions=[]
        for spec in job.systems:
            candidate=MANAGER_ARCHIVE.efficiency_candidate(spec)
            if candidate: decisions.append({"label":spec.get("label",job.label),"action":candidate["kind"],"archive_path":candidate["archive_path"]})
        from .resource_limits import numerical_thread_policy
        provenance={**execution_provenance(),"numerical_thread_policy":numerical_thread_policy(int(job.resources.get("cores",1))),**job.provenance}
        self._set_status(job.job_id, "pending", label=job.label, manifest=str(path), original_manifest=str(original), dependencies=list(job.dependencies), dependency_conditions=job.dependency_conditions, resources=job.resources, staging=job.staging, deadline=job.deadline, provenance=provenance, tags=job.tags, plan_fingerprint=job.plan_fingerprint, planned_system_count=len(job.systems), efficiency_check=decisions)
        self._record_campaign_history(job,status="pending",manifest_path=path)
        self._event(job.job_id, "submitted", manifest=str(path),efficiency_check=decisions); return path

    def configure_scheduler_policy(self, **limits: Any) -> dict[str,Any]:
        from .scheduling import SchedulerPolicy
        policy=SchedulerPolicy(**limits); payload={key:value for key,value in policy.__dict__.items() if value is not None}
        self.policy_path.write_text(json.dumps(payload,indent=2,sort_keys=True),encoding="utf-8"); return payload

    def submit_array(self, array: Any) -> Path:
        """Submit indexed tasks and persist one immutable array manifest."""
        from .scheduling import JobArray
        if not isinstance(array,JobArray): raise TypeError("submit_array expects a JobArray.")
        for index,job in enumerate(array.jobs):
            tagged=replace(job,tags={**job.tags,"array_id":array.array_id,"array_index":index,"array_label":array.label,"array_max_concurrent":array.max_concurrent})
            self.submit(tagged)
        path=self.array_dir/f"{_safe(array.label)}__{array.array_id}.json"; path.write_text(json.dumps(array.payload(),indent=2,sort_keys=True),encoding="utf-8")
        return path

    def run_array(self, array_path: str|Path, *, backend: str|None=None) -> dict[str,Any]:
        payload=json.loads(Path(array_path).read_text(encoding="utf-8")); return self.run(payload["job_ids"],mode="parallel",max_workers=int(payload["max_concurrent"]),backend=backend)

    def snapshot_original_manifest(self, job_id: str) -> Path:
        """Create the one-time immutable reference copy for an existing Job."""
        row=self.status(job_id); source=Path(row["manifest"]); target=self.original_manifest_dir/source.name
        if not target.exists(): target.write_bytes(source.read_bytes())
        self._set_status(job_id,row.get("status","pending"),original_manifest=str(target))
        self._event(job_id,"original_manifest_snapshot",original_manifest=str(target))
        return target

    def submit_collection(self, collection: JobCollection) -> list[Path]:
        paths = [self.submit(item) for item in collection.jobs]
        payload = {"format":"dynsyslab.job-collection/v1", "collection_id":collection.collection_id, "label":collection.label, "tags":collection.tags, "jobs":[item.job_id for item in collection.jobs]}
        path = self.collection_dir / f"{_safe(collection.label)}__{collection.collection_id}.json"; path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        for item in collection.jobs:
            self._set_status(item.job_id, self.status(item.job_id).get("status", "pending"), collection_id=collection.collection_id, collection_manifest=str(path))
        return paths

    def reserve(self, manifest_path: str | Path, *, note: str="", tags: dict[str,Any] | None=None) -> dict[str,Any]:
        """Put an existing Job, collection, or Campaign manifest in the side tray.

        Reserving is deliberately not a scheduling operation: it creates no
        JobManager job state, submits no work, and does not call the solver. A
        byte-for-byte snapshot protects the intended future plan if the source
        manifest is later edited or moved.
        """
        source=Path(manifest_path).expanduser().resolve()
        if not source.is_file(): raise FileNotFoundError(source)
        payload=json.loads(source.read_text(encoding="utf-8"))
        if "campaign_id" in payload and "jobs" in payload: kind="campaign"; identifier=str(payload["campaign_id"])
        elif "collection_id" in payload and "jobs" in payload: kind="job_collection"; identifier=str(payload["collection_id"])
        elif "job_id" in payload and "systems" in payload: kind="job"; identifier=str(payload["job_id"])
        else: raise ValueError("Side-tray entries must be a dynsyslab campaign, job collection, or job manifest.")
        source_bytes=source.read_bytes(); fingerprint=__import__("hashlib").sha256(source_bytes).hexdigest()
        entry_id=uuid4().hex; snapshot=self.side_tray_snapshot_dir/f"{_safe(str(payload.get('label',kind)))}__{entry_id}.json"
        with snapshot.open("wb") as handle: handle.write(source_bytes); handle.flush(); os.fsync(handle.fileno())
        entry={"entry_id":entry_id,"kind":kind,"identifier":identifier,"label":str(payload.get("label",source.stem)),
               "source_manifest":str(source),"snapshot_manifest":str(snapshot),"source_sha256":fingerprint,
               "reserved_utc":_stamp(),"note":str(note),"tags":dict(tags or {}),"status":"reserved"}
        with self._state_file_lock():
            tray=self._read_side_tray()
            duplicate=next((item for item in tray.get("entries",[]) if item.get("kind")==kind and item.get("identifier")==identifier and item.get("status")=="reserved"),None)
            if duplicate is not None:
                snapshot.unlink(missing_ok=True)
                return duplicate
            tray.setdefault("entries",[]).append(entry); self._write_side_tray(tray)
        self._event(f"side-tray-{entry_id}","reserved",kind=kind,identifier=identifier,label=entry["label"])
        return entry

    def side_tray(self, *, include_removed: bool=False) -> list[dict[str,Any]]:
        """List reserved future work without querying or changing the queue."""
        entries=list(self._read_side_tray().get("entries",[]))
        return entries if include_removed else [entry for entry in entries if entry.get("status")=="reserved"]

    def remove_reservation(self, entry_id: str, *, delete_snapshot: bool=True) -> dict[str,Any]:
        """Remove a reservation only; it never cancels or deletes submitted work."""
        with self._state_file_lock():
            tray=self._read_side_tray(); entry=next((item for item in tray.get("entries",[]) if item.get("entry_id")==entry_id),None)
            if entry is None: raise KeyError(f"No side-tray entry {entry_id!r}.")
            if entry.get("status")!="reserved": return entry
            entry["status"]="removed"; entry["removed_utc"]=_stamp(); self._write_side_tray(tray)
        if delete_snapshot: Path(entry["snapshot_manifest"]).unlink(missing_ok=True)
        self._event(f"side-tray-{entry_id}","reservation_removed",kind=entry.get("kind"),identifier=entry.get("identifier"))
        return entry

    def _manifest(self, job_id: str) -> JobManifest:
        """Load a job manifest, including immutable plans removed from live state."""
        state = self._read_state(); row = state.get("jobs", {}).get(job_id, {})
        manifest = row.get("manifest") or row.get("original_manifest")
        if manifest and Path(manifest).is_file(): return JobManifest.read(manifest)
        snapshots = sorted(self.original_manifest_dir.glob(f"*__{job_id}.json"))
        if snapshots: return JobManifest.read(snapshots[-1])
        raise KeyError(f"Job manifest {job_id!r} is absent from live state and immutable originals.")

    def status(self, job_id: str | None = None) -> dict[str, Any]:
        state = self._read_state(); return state["jobs"] if job_id is None else state["jobs"].get(job_id, {})

    def progress(self, job_ids: Iterable[str] | None = None) -> dict[str, Any]:
        """Return inexpensive persisted progress and an empirical time estimate."""
        return aggregate_progress(self.status(),job_ids)

    def campaign_status(self, campaign_manifest: str | Path) -> dict[str, Any]:
        return self._campaign_reporter.status(campaign_manifest)

    def validate_campaign_against_original(self, campaign_manifest: str | Path) -> dict[str, Any]:
        return self._campaign_reporter.validate(campaign_manifest)

    def _postprocessors(self, plans: tuple[dict[str, Any], ...]):
        if not plans: return ()
        from postprocessing.koopman import dmd_archive, edmd_archive, hankel_dmd_archive
        lookup = {"dmd": dmd_archive, "edmd": edmd_archive, "hdmd": hankel_dmd_archive}
        output = []
        for plan in plans:
            name = str(plan.get("method", "")).lower()
            if name not in lookup: raise ValueError(f"Unsupported job postprocessor {name!r}; add an explicit adapter before scheduling it.")
            settings = dict(plan.get("settings", {})); function = lookup[name]
            output.append(lambda path, function=function, settings=settings: function(path, **settings))
        return tuple(output)

    def _build_target(self, job: JobManifest) -> System | Batch:
        from app import RunConfiguration
        storage_root=None
        location=str(job.schedule.get("result_location","primary"))
        if location!="primary":
            from archiver.storage import StorageRegistry
            record=StorageRegistry().locations().get(location)
            if not record or not record.get("writable",False): raise ValueError(f"Requested archive location {location!r} is not registered and writable.")
            storage_root=str(record["root"])
        systems = []
        for spec in job.systems:
            config = RunConfiguration(system_name=spec["system_name"], preset=spec.get("preset"), parameters=dict(spec.get("parameters", {})), driver=dict(spec.get("driver", {"kind":"zero"})), time_span=tuple(spec["time_span"]), label=spec.get("label", job.label))
            request = config.build_request()
            hierarchy = spec.get("hierarchy", {})
            clusters = []
            for cluster_spec in hierarchy.get("clusters", []):
                initial_conditions = cluster_spec["initial_conditions"]
                trajectories = []
                for index, value in enumerate(initial_conditions):
                    state = np.asarray(request.system.initial_state, dtype=float) if isinstance(value, str) and value == "__system_initial_state__" else np.asarray(value, dtype=float)
                    trajectories.append(Trajectory(state, label=f"{cluster_spec.get('label','cluster')}-t{index:03d}"))
                clusters.append(Cluster(trajectories, label=cluster_spec.get("label", "cluster")))
            if not clusters: raise ValueError("Every planned system needs at least one cluster of initial conditions.")
            output = spec.get("output", {})
            if "times" in output: times = np.asarray(output["times"], dtype=float); uniform = False; delta_t = None
            else:
                times = None; uniform = True; delta_t = float(output["delta_t"])
            systems.append(System(request.system, clusters, tuple(spec["time_span"]), times, uniform, delta_t, request.kind, label=spec.get("label", job.label), archive_root=storage_root))
        return systems[0] if len(systems) == 1 else Batch(systems, label=job.label)

    def _options(self, job: JobManifest) -> tuple[SolverOptions, HistoryOptions, str]:
        from solver import HistoryOptions, SolverOptions
        spec = job.systems[0].get("solver", {})
        options,history,resolved=resolve_history_policy(spec.get("history_policy"),resources=job.resources,options=dict(spec.get("options",{})),history=dict(spec.get("history",{})))
        deadline=dict(job.deadline)
        if deadline.get("requeue_on_timeout") and deadline.get("checkpoint_before_deadline_s"):
            history={**history,"mode":"checkpoint","track_during_solve":True,"checkpoint_strategy":"low_level_stepper","pause_at_checkpoint":True,
                     "checkpoint_interval":history.get("checkpoint_interval",max(1.,float(deadline["checkpoint_before_deadline_s"])/2))}
            options={**options,"adaptive_engine":"low_level_stepper"}
            resolved={**resolved,"deadline_checkpoint_enabled":True}
        self._event(job.job_id,"history_policy_resolved",**resolved)
        return SolverOptions(**options), HistoryOptions(**history), str(spec.get("interpolation_method", "pchip"))

    def _tag_archives(self, result: Any, job: JobManifest) -> list[str]:
        paths = [item.archive_path for item in result.trajectory_results if getattr(item, "archive_path", None)] + [result.archive_path]
        patch = {"job": {"job_id": job.job_id, "label": job.label, "tags": job.tags, "manifest_format": job.format, "plan_fingerprint":job.plan_fingerprint}}
        for path in paths: MANAGER_ARCHIVE.update_metadata(path, patch)
        return [str(path) for path in paths]

    def _hierarchy_matches(self, archive_path: str | Path, target: System | Batch) -> bool:
        """Manually verify that a hierarchy manifest matches the planned tree."""
        return self._archive_inspector.hierarchy_matches(archive_path,target)

    def _verify_or_restore_hierarchy(self, result: Any, target: System | Batch, job: JobManifest) -> tuple[str,bool]:
        original=str(result.archive_path)
        if self._hierarchy_matches(original,target):
            self._event(job.job_id,"hierarchy_verified",hierarchy_archive=original,restored=False)
            return original,False
        restored=str(MANAGER_ARCHIVE.save_hierarchy(target))
        self._event(job.job_id,"hierarchy_restored",original_hierarchy_archive=original,restored_hierarchy_archive=restored,
                    reason="saved hierarchy did not match the planned manifest")
        return restored,True

    def _dependencies_complete(self, job: JobManifest) -> bool:
        state = self._read_state()["jobs"]
        from .scheduling import dependency_satisfied
        for item in job.dependencies:
            condition=job.dependency_conditions.get(item,"afterok")
            parent=state.get(item,{})
            if condition!="aftercorr":
                if not dependency_satisfied(condition,parent.get("status")): return False
                continue
            index=job.tags.get("array_index")
            matches=[row for row in state.values() if row.get("tags",{}).get("array_id")==parent.get("tags",{}).get("array_id") and row.get("tags",{}).get("array_index")==index]
            if not matches or not dependency_satisfied("aftercorr",matches[0].get("status")): return False
        return True

    def _scheduled_ready(self, job: JobManifest) -> bool:
        start=job.schedule.get("start_at")
        if not start: return True
        scheduled=datetime.fromisoformat(str(start).replace("Z","+00:00"))
        if scheduled.tzinfo is None: scheduled=scheduled.replace(tzinfo=timezone.utc)
        if start and datetime.now(timezone.utc)<scheduled.astimezone(timezone.utc): return False
        campaigns={str(item) for item in job.schedule.get("after_campaign_ids",[])}
        if campaigns:
            rows=self.status().values()
            for campaign_id in campaigns:
                predecessors=[row for row in rows if row.get("tags",{}).get("campaign_id")==campaign_id]
                if not predecessors or any(row.get("status")!="completed" for row in predecessors if row.get("status") not in {"superseded","replaced","transitioned"}): return False
        return True

    def pause(self, job_ids: Iterable[str]) -> dict[str,Any]:
        """Pause dispatch safely; active solves finish their current atomic solve."""
        ids=list(job_ids)
        result=pause_jobs(ids,self.status,lambda job_id:str(self.snapshot_original_manifest(job_id)),self._set_status,self._event,_stamp)
        result["watchdog"]=self.watch_paused_workers(ids)
        return result

    def resume(self, job_ids: Iterable[str]) -> list[str]:
        return resume_jobs(job_ids,self.status,self._set_status,self._event,_stamp)

    def watch_paused_workers(self, job_ids: Iterable[str] | None=None, *, grace_seconds: float=120.) -> dict[str,dict[str,Any]]:
        """Refresh conservative post-pause worker classifications.

        This never terminates workers and never releases quarantined claims;
        it is safe to call repeatedly from a status UI, CLI, or scheduler.
        """
        ids=list(job_ids or self.status().keys())
        return watch_paused_jobs(ids,status=self.status,set_status=self._set_status,event=self._event,stamp=_stamp,
                                 telemetry_root=MANAGER_ARCHIVE.archive_root/"solver_logs",grace_seconds=grace_seconds)

    def probe_worker_activity(self, job_ids: Iterable[str] | None=None, *, sample_seconds: float=2.) -> dict[str,dict[str,Any]]:
        """Read-only local CPU probe for diagnosing apparently rogue workers."""
        return sample_worker_activity(list(job_ids or self.status().keys()),status=self.status,sample_seconds=sample_seconds)

    def register_worker(self, job_id: str, *, pid: int, backend: str="local_process") -> None:
        """Worker-side registration used only for durable watchdog evidence."""
        row=self.status(job_id)
        if not row: raise KeyError(f"Unknown job {job_id!r}.")
        self._set_status(job_id,row.get("status","running"),worker_pid=int(pid),worker_backend=backend,worker_started_utc=_stamp())
        self._event(job_id,"worker_registered",pid=int(pid),backend=backend)

    def record_worker_exit(self, job_id: str, *, pid: int, exit_code: int) -> None:
        row=self.status(job_id)
        if not row: return
        self._set_status(job_id,row.get("status","unknown"),worker_pid=int(pid),worker_exit_utc=_stamp(),worker_exit_code=int(exit_code))
        self._event(job_id,"worker_exit",pid=int(pid),exit_code=int(exit_code))

    def request_campaign_relocation(self, campaign_manifest: str | Path, *, destination_location: str) -> dict[str,Any]:
        """Arm a live-safe relocation: completed files copy now; new jobs route later."""
        campaign_path=Path(campaign_manifest).resolve(); campaign=json.loads(campaign_path.read_text(encoding="utf-8")); campaign_id=campaign["campaign_id"]
        from archiver.storage import StorageRegistry
        registry=StorageRegistry(); plan=registry.create_relocation_plan(campaign_id=campaign_id,destination=destination_location)
        changed=[]; active=[]; copied=[]; seen=set()
        for job_id,row in self.status().items():
            if row.get("tags",{}).get("campaign_id")!=campaign_id: continue
            for raw_path in list(row.get("generated_solution_archives",[]))+list(row.get("solution_archives",[])):
                item=Path(raw_path)
                # A campaign hierarchy manifest is part of the reproducible
                # result, not expendable bookkeeping.  It travels with its
                # trajectory archives and retains the intended tree.
                if item.is_file() and str(item) not in seen:
                    copied.append(registry.append_completed_copy(plan,item)); seen.add(str(item))
            if row.get("status")=="running": self._set_status(job_id,"running",relocation_plan=str(plan),relocation_destination=destination_location); active.append(job_id); continue
            if row.get("status") not in {"pending","blocked","failed","paused"}: continue
            job=self._manifest(job_id); updated=replace(job,schedule={**job.schedule,"result_location":destination_location})
            updated.write(row["manifest"]); self._set_status(job_id,row.get("status"),relocation_plan=str(plan),relocation_destination=destination_location); self._event(job_id,"relocation_routed",plan=str(plan),destination=destination_location); changed.append(job_id)
        payload=json.loads(plan.read_text(encoding="utf-8")); payload["pending_job_ids"]=changed; payload["active_job_ids"]=active; plan.write_text(json.dumps(payload,indent=2,sort_keys=True),encoding="utf-8")
        return {"plan_path":str(plan),"copied_completed_files":copied,"routed_pending_jobs":changed,"active_jobs_to_copy_on_completion":active}

    def _notify_relocation_completion(self, job_id: str, paths: list[str]) -> None:
        row=self.status(job_id); plan=row.get("relocation_plan")
        if not plan or not Path(plan).is_file(): return
        from archiver.storage import StorageRegistry
        registry=StorageRegistry(); copied=[]
        for path in paths:
            item=Path(path)
            if item.is_file(): copied.append(registry.append_completed_copy(plan,item))
        self._event(job_id,"relocation_completion_notified",plan=plan,copied=copied)
