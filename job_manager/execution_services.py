"""Local and backend execution, restructuring, and adaptive scheduling."""
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


class _JobExecutionMixin:
    """Local and backend execution, restructuring, and adaptive scheduling."""

    def execute(self, job_id: str) -> dict[str, Any]:
        from solver import solve
        job = self._manifest(job_id)
        previous = self._completed_plan_archive(job.plan_fingerprint)
        if previous is not None:
            self._set_status(job_id, "completed", reused=True, reused_from_archive=previous, hierarchy_archive=previous, solution_archives=[], duration_seconds=0.)
            self._event(job_id, "reused", existing_hierarchy_archive=previous, reason="identical plan already archived"); return self.status(job_id)
        if not self._dependencies_complete(job):
            self._set_status(job_id, "blocked", reason="uncompleted dependencies"); self._event(job_id, "blocked", dependencies=list(job.dependencies)); return self.status(job_id)
        # Claiming is a short durable state between eligibility and execution.
        # It prevents another worker from reclaiming the same scientific spec
        # during the otherwise unavoidable process-start race.
        self._set_status(job_id,"claiming",claim_started_utc=_stamp())
        claimed,conflicts=self._claim_specs(job_id,job.systems)
        if not claimed:
            self._set_status(job_id,"waiting_for_spec_claim",claim_conflicts=conflicts)
            self._event(job_id,"waiting_for_spec_claim",conflicts=conflicts)
            return self.status(job_id)
        previous_state=self.status(job_id)
        attempt_count=int(previous_state.get("attempt_count",0))+1
        retry_limit=max(0,int(job.execution.get("retry_limit",0)))
        self._set_status(job_id, "running", started_utc=_stamp(), attempt_count=attempt_count,
                         retry_limit=retry_limit, retry_eligible=False)
        self._record_campaign_history(job,status="running",manifest_path=self.status(job_id).get("manifest"))
        self._event(job_id, "started", attempt_count=attempt_count, retry_limit=retry_limit)
        started = perf_counter()
        try:
            from .staging import stage_inputs
            staging_manifest=stage_inputs(job,self.root/"staging")
            self._set_status(job_id,"running",staging=staging_manifest)
            self._event(job_id,"inputs_staged",input_count=len(staging_manifest["inputs"]),staging_directory=staging_manifest["staging_directory"])
            extension=self._try_single_solution_extension(job)
            if extension is not None:
                elapsed=perf_counter()-started; self._set_status(job_id,"completed",duration_seconds=elapsed,extended=True,solution_archives=[extension["extended_archive"]],continuation_archive=extension["continuation_archive"])
                self._record_campaign_history(job,status="completed",manifest_path=self.status(job_id).get("manifest"),solution_paths=[extension["extended_archive"]])
                self._event(job_id,"extended_existing_solution",duration_seconds=elapsed,**extension); self._release_spec_claims(job_id); return self.status(job_id)
            reused,remaining=self._partition_reusable_specs(job)
            reuse_manifest=None
            if reused:
                reuse_manifest=self._reuse_manifest(job,reused,remaining)
                self._event(job_id,"exact_archive_reuse",reused_count=len(reused),remaining_count=len(remaining),reuse_manifest=str(reuse_manifest),reused=reused)
            if not remaining:
                elapsed=perf_counter()-started; paths=[item["archive_path"] for item in reused]
                self._set_status(job_id,"completed",duration_seconds=elapsed,reused=True,reused_solution_archives=paths,
                                 solution_archives=paths,reuse_partition_manifest=None if reuse_manifest is None else str(reuse_manifest))
                self._record_campaign_history(job,status="completed",manifest_path=self.status(job_id).get("manifest"),solution_paths=paths)
                self._event(job_id,"completed_from_exact_reuse",duration_seconds=elapsed,solution_file_count=len(paths),reuse_manifest=None if reuse_manifest is None else str(reuse_manifest))
                self._release_spec_claims(job_id); return self.status(job_id)
            execution_job=job if len(remaining)==len(job.systems) else replace(job,systems=tuple(remaining))
            target = self._build_target(execution_job); options, history, interpolation_method = self._options(job)
            deadline_token=None; deadline_timer=None
            wall_time=job.resources.get("wall_time_s"); lead=job.deadline.get("checkpoint_before_deadline_s")
            if wall_time and lead and job.deadline.get("requeue_on_timeout"):
                delay=max(0.,float(wall_time)-float(lead)); deadline_token=Event(); deadline_timer=Timer(delay,deadline_token.set); deadline_timer.daemon=True; deadline_timer.start()
                self._event(job_id,"deadline_checkpoint_armed",wall_time_s=wall_time,checkpoint_lead_s=lead,trigger_after_s=delay)
            expected=sum(len(cluster.trajectories) for system in (target.systems if isinstance(target, Batch) else [target]) for cluster in system.clusters)
            def progress(event: dict[str, Any]) -> None:
                if event.get("event") != "hierarchy_progress": return
                completed=int(event["completed_trajectories"])
                current=self.status(job_id); generated=list(current.get("generated_solution_archives",[])); archive_path=event.get("archive_path")
                if archive_path and archive_path not in generated: generated.append(archive_path)
                elapsed_progress=perf_counter()-started; rate=None if not completed else elapsed_progress/completed
                self._set_status(job_id, "running", completed_systems=completed, planned_system_count=expected,
                                 fraction_complete=completed / expected if expected else 1., generated_solution_archives=generated,
                                 generated_solution_bytes=sum(Path(path).stat().st_size for path in generated if Path(path).is_file()),
                                 elapsed_seconds=elapsed_progress,seconds_per_trajectory=rate,estimated_remaining_seconds=None if rate is None else rate*(expected-completed))
                self._event(job_id, "progress", completed_systems=completed, planned_system_count=expected,
                            solution_file_count=len(generated),solution_bytes=sum(Path(path).stat().st_size for path in generated if Path(path).is_file()),estimated_remaining_seconds=None if rate is None else rate*(expected-completed))
                self._campaign_progress(job)
                if self.status(job_id).get("pause_after_current_solve"):
                    raise _PauseAtBoundary("pause requested after a completed trajectory")
            result = solve(target, options=options, history=history, interpolation_method=interpolation_method, postprocessors=self._postprocessors(job.postprocessing), on_event=progress, pause_token=deadline_token, job_id=job_id)
            if deadline_timer is not None: deadline_timer.cancel()
            if getattr(result,"status","complete")=="paused":
                elapsed=perf_counter()-started; partial=[item.archive_path for item in result.trajectory_results if getattr(item,"archive_path",None)]
                self._set_status(job_id,"paused",paused_utc=_stamp(),incomplete=True,pause_reason="deadline checkpoint",generated_solution_archives=partial,deadline_requeue_eligible=True)
                self._event(job_id,"deadline_checkpoint_paused",duration_seconds=elapsed,partial_solution_count=len(partial)); self._release_spec_claims(job_id); return self.status(job_id)
            hierarchy_archive, restored=self._verify_or_restore_hierarchy(result,target,execution_job)
            result.archive_path=hierarchy_archive
            paths = self._tag_archives(result, job); paths=list(dict.fromkeys([item["archive_path"] for item in reused]+paths)); elapsed = perf_counter() - started
            self._set_status(job_id, "completed", duration_seconds=elapsed, hierarchy_archive=result.archive_path, solution_archives=paths, hierarchy_restored=restored,
                             reused_solution_archives=[item["archive_path"] for item in reused],reuse_partition_manifest=None if reuse_manifest is None else str(reuse_manifest))
            self._record_campaign_history(job,status="completed",manifest_path=self.status(job_id).get("manifest"),solution_paths=paths)
            self._notify_relocation_completion(job_id,paths)
            self._campaign_progress(job,force=True)
            self._event(job_id, "completed", duration_seconds=elapsed, hierarchy_archive=result.archive_path, hierarchy_restored=restored); self._release_spec_claims(job_id); return self.status(job_id)
        except _PauseAtBoundary:
            current=self.status(job_id); partial=list(current.get("generated_solution_archives",[]))
            self._set_status(job_id,"paused",paused_utc=_stamp(),incomplete=True,pause_reason="safe trajectory boundary",generated_solution_archives=partial)
            self._event(job_id,"paused_at_safe_boundary",completed_systems=current.get("completed_systems",0),solution_file_count=len(partial))
            self._release_spec_claims(job_id); return self.status(job_id)
        except KeyboardInterrupt:
            interrupted=self._record_interruption(job_id,reason="KeyboardInterrupt during local solve"); self._record_campaign_history(job,status="interrupted",manifest_path=self.status(job_id).get("manifest")); self._release_spec_claims(job_id); return interrupted
        except Exception as error:
            elapsed=perf_counter()-started; partial=self.status(job_id).get("generated_solution_archives",[]); partial_paths=[Path(path) for path in partial if Path(path).is_file()]
            retry_eligible=attempt_count<=retry_limit
            self._set_status(job_id,"failed",duration_seconds=elapsed,reason=str(error),generated_solution_archives=partial,
                             generated_solution_bytes=sum(path.stat().st_size for path in partial_paths),attempt_count=attempt_count,
                             retry_limit=retry_limit,retry_eligible=retry_eligible)
            self._record_campaign_history(job,status="failed",manifest_path=self.status(job_id).get("manifest"),solution_paths=partial)
            self._event(job_id,"failed",duration_seconds=elapsed,reason=str(error),solution_file_count=len(partial_paths),
                        solution_bytes=sum(path.stat().st_size for path in partial_paths),attempt_count=attempt_count,
                        retry_limit=retry_limit,retry_eligible=retry_eligible)
            campaign_path=job.tags.get("campaign_manifest")
            if campaign_path and Path(campaign_path).is_file():
                campaign=json.loads(Path(campaign_path).read_text(encoding="utf-8")); reporting=(campaign.get("metadata") or {}).get("progress_reporting",{}); log_path=reporting.get("log_path")
                if log_path:
                    with Path(log_path).open("a",encoding="utf-8") as handle: handle.write(json.dumps({"event":"job_failed","time_utc":_stamp(),"campaign_id":campaign.get("campaign_id"),"job_id":job_id,"reason":str(error),"solution_file_count":len(partial_paths),"solution_bytes":sum(path.stat().st_size for path in partial_paths)})+"\n")
            self._release_spec_claims(job_id); return self.status(job_id)

    def restructure(self, job_id: str, *, systems_per_child: int, reason: str="manual adaptive restructuring") -> list[str]:
        """Split an unstarted Job into child Jobs without changing its science plan.

        The parent is retained as a logical execution node.  A durable
        restructuring manifest records the original hierarchy and all child
        archives, allowing later archive tools to reconstruct the intended
        parent collection without copying any solution data.
        """
        if systems_per_child < 1: raise ValueError("systems_per_child must be positive.")
        parent=self._manifest(job_id); state=self.status(job_id)
        if state.get("status") not in {"pending","blocked"}: raise ValueError("Only unstarted Jobs can be restructured.")
        if len(parent.systems) <= systems_per_child: return []
        children=[]
        for index in range(0,len(parent.systems),systems_per_child):
            systems=parent.systems[index:index+systems_per_child]
            child=replace(parent, label=f"{parent.label}-part-{index//systems_per_child:03d}", systems=tuple(systems),
                          job_id=uuid4().hex, tags={**parent.tags,"logical_parent_job":parent.job_id,"execution_restructured":True},
                          notes=f"{parent.notes}\nDerived from {parent.job_id}: {reason}".strip())
            self.submit(child); children.append(child)
        manifest=self.restructure_dir/f"{_safe(parent.label)}__{parent.job_id}.json"
        manifest.write_text(json.dumps({"format":"dynsyslab.job-restructure/v1","created_utc":_stamp(),"parent_job_id":parent.job_id,"reason":reason,"original_manifest":str(self.manifest_dir/f'{_safe(parent.label)}__{parent.job_id}.json'),"child_job_ids":[item.job_id for item in children],"original_hierarchy":"preserved in parent JobManifest; child hierarchy archives are linked after completion."},indent=2,sort_keys=True),encoding="utf-8")
        self._set_status(parent.job_id,"restructured",child_job_ids=[item.job_id for item in children],restructure_manifest=str(manifest),reason=reason)
        self._event(parent.job_id,"restructured",child_job_ids=[item.job_id for item in children],restructure_manifest=str(manifest),reason=reason)
        return [item.job_id for item in children]

    def _finalize_restructured(self, parent_id: str) -> bool:
        row=self.status(parent_id)
        if row.get("status")!="restructured": return False
        children=list(row.get("child_job_ids",[]))
        if not children or any(self.status(item).get("status")!="completed" for item in children): return False
        manifest=Path(row["restructure_manifest"]); payload=json.loads(manifest.read_text(encoding="utf-8"))
        payload["finalized_utc"]=_stamp(); payload["child_hierarchy_archives"]=[self.status(item).get("hierarchy_archive") for item in children]
        temporary=manifest.with_suffix(".partial.json"); temporary.write_text(json.dumps(payload,indent=2,sort_keys=True),encoding="utf-8"); temporary.replace(manifest)
        self._set_status(parent_id,"completed",execution_hierarchy_manifest=str(manifest),child_job_ids=children,restructured=True)
        self._event(parent_id,"restructure_completed",child_job_ids=children,execution_hierarchy_manifest=str(manifest))
        return True

    def _adapt_pending(self, ids: list[str], policy: str) -> list[str]:
        if policy == "fixed": return ids
        expanded=list(ids)
        for job_id in list(ids):
            row=self.status(job_id)
            if row.get("status") not in {"pending","blocked"}: continue
            job=self._manifest(job_id); limit=job.execution.get("max_systems_per_job")
            if limit and len(job.systems)>int(limit):
                children=self.restructure(job_id,systems_per_child=int(limit),reason=f"{policy} scheduler chunk limit")
                expanded.extend(children)
        return expanded

    def _adaptive_worker_cap(self, current: int) -> int:
        """Change a local worker cap only between fully completed waves."""
        rows=list(self.status().values()); cores=os.cpu_count() or 1
        if any(row.get("status")=="running" for row in rows): return current
        recent=[row for row in rows if row.get("status") in {"completed","failed","interrupted"}]
        failures=sum(row.get("status")!="completed" for row in recent[-max(1,current):])
        if failures: return max(1,current-1)
        return min(max(1,cores-1),current+1)

    def run(self, job_ids: Iterable[str] | None = None, *, mode: ExecutionMode = "series", max_workers: int | None = None,
            backend: str | ExecutionBackend | None = None, adaptive_policy: str="fixed") -> dict[str, Any]:
        """Schedule a dependency DAG through a JobManager-owned backend.

        ``parallel`` now uses isolated local worker processes.  Remote backend
        submission is explicit and never occurs merely because a backend class
        exists.  Campaigner remains a planner and the Solver remains a single
        solve engine throughout.
        """
        ids=list(job_ids or self.status().keys())
        if backend is None:
            requested={str(self._manifest(item).execution.get("backend","serial")) for item in ids if item in self.status()}
            backend_name=requested.pop() if len(requested)==1 else ("local_process" if mode=="parallel" else "serial")
        else: backend_name=backend
        if adaptive_policy=="fixed":
            requested_policies={str(self._manifest(item).execution.get("adaptive_policy","fixed")) for item in ids if item in self.status()}
            if len(requested_policies)==1: adaptive_policy=requested_policies.pop()
        ids=self._adapt_pending(ids,adaptive_policy)
        if isinstance(backend_name,str):
            selected=create_backend(backend_name,max_workers=max_workers) if backend_name in {"local_process","local_worker"} else create_backend(backend_name)
        else: selected=backend_name
        context={"manager_root":self.root,"project_root":Path(__file__).resolve().parents[1],
                 "resources_by_job":{item:self._manifest(item).resources for item in ids}}
        while True:
            for parent in list(ids): self._finalize_restructured(parent)
            rows={item:self.status(item) for item in ids}
            ready=[]
            for item,row in rows.items():
                status=row.get("status")
                retryable=status=="failed" and bool(row.get("retry_eligible",False))
                if status not in {"pending","blocked"} and not retryable:
                    continue
                if self._dependencies_complete(self._manifest(item)) and self._scheduled_ready(self._manifest(item)):
                    ready.append(item)
            # Apply local admission limits before dispatch.  The policy is
            # advisory for remote backends, which translate the same resource
            # request to their native scheduler/QoS declarations.
            from .scheduling import ResourceRequest, SchedulerPolicy, admission_reason
            limits=json.loads(self.policy_path.read_text(encoding="utf-8")) if self.policy_path.is_file() else {}
            policy=SchedulerPolicy(**limits); admitted=[]
            ready.sort(key=lambda item:(-int(self._manifest(item).schedule.get("priority",0)),str(item)))
            for item in ready:
                manifest=self._manifest(item); request=ResourceRequest(**{key:value for key,value in manifest.resources.items() if key in ResourceRequest.__dataclass_fields__})
                qos=str(manifest.execution.get("qos",""))
                if policy.allowed_qos and qos not in policy.allowed_qos:
                    self._set_status(item,"blocked",reason="scheduler admission: qos"); self._event(item,"admission_deferred",reason="qos"); continue
                reason=admission_reason(request,running=list(rows.values())+[{**rows[chosen],"status":"running"} for chosen in admitted],policy=policy)
                if reason is None: admitted.append(item)
                else: self._set_status(item,"blocked",reason=f"scheduler admission: {reason}"); self._event(item,"admission_deferred",reason=reason)
            ready=admitted
            if not ready: break
            if isinstance(selected,LocalProcessBackend) and (mode=="parallel" or selected.capabilities.supports_parallel):
                wave=ready
                if adaptive_policy in {"adaptive_local","adaptive_portable"}:
                    prior_workers=selected.max_workers or max(1,min(4,(os.cpu_count() or 2)-1))
                    selected.max_workers=self._adaptive_worker_cap(prior_workers)
                    wave=ready[:selected.max_workers]
                    if selected.max_workers != prior_workers:
                        change={
                            "timestamp":_stamp(), "backend":selected.capabilities.name,
                            "prior_workers":prior_workers, "new_workers":selected.max_workers,
                            "safety":"the preceding worker wave was quiescent; no running Jobs were rescheduled",
                        }
                        self._append_reparallelization_report(self._campaign_path_for(wave),change)
                        self._event(wave[0],"parallel_rebalanced",**change)
                    self._event(wave[0],"parallel_wave_started",backend=selected.capabilities.name,worker_cap=selected.max_workers,wave_job_ids=wave)
                try: selected.execute_many(wave,context=context)
                except KeyboardInterrupt:
                    for item in wave: self._record_interruption(item,reason="JobManager interrupted while local worker processes were active")
                    break
            elif selected.capabilities.supports_remote:
                for item in ready:
                    handle=selected.submit(item,context=context); self._set_status(item,"submitted_remote",backend=selected.capabilities.name,remote_handle=handle.payload()); self._event(item,"remote_submitted",backend=selected.capabilities.name,remote_handle=handle.payload())
                break
            else:
                try:
                    for item in ready: selected.execute(item,execute_local=self.execute,context=context)
                except KeyboardInterrupt:
                    for item in ready:
                        if self.status(item).get("status") == "running": self._record_interruption(item,reason="JobManager interrupted during scheduling")
                    break
        snapshot={item:self.status(item) for item in ids}; finished=[row for row in snapshot.values() if row.get("status")=="completed"]
        mean=float(np.mean([row.get("duration_seconds",0.) for row in finished])) if finished else None; remaining=sum(row.get("status")!="completed" for row in snapshot.values())
        return {"jobs":snapshot,"completed":len(finished),"total":len(ids),"backend":selected.capabilities.name,"adaptive_policy":adaptive_policy,"estimated_remaining_seconds":None if mean is None else mean*remaining}
