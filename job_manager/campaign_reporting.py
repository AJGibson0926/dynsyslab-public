"""Read-only campaign progress and specification-conformance reporting."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import re
import numpy as np

from .manifests import JobManifest


class CampaignReporter:
    """A focused reporter that consumes manager state but never dispatches work."""
    def __init__(self, manager: Any, archive: Any, stamp) -> None:
        self.manager=manager; self.archive=archive; self.stamp=stamp

    def validate(self, campaign_manifest: str|Path) -> dict[str,Any]:
        manifest=Path(campaign_manifest).resolve(); campaign=json.loads(manifest.read_text(encoding="utf-8")); state=self.manager.status(); expected={}; snapshots=[]
        for job_id in campaign["jobs"]:
            row=state.get(job_id,{})
            snapshot=Path(row.get("original_manifest") or self.manager.snapshot_original_manifest(job_id))
            snapshots.append(str(snapshot)); job=JobManifest.read(snapshot)
            for spec in job.systems: expected[str(spec.get("label",job.label))]=spec
        actual: dict[str,list[tuple[str,dict[str,Any]]]]={}
        for path in self.archive.solution_paths():
            try:
                metadata=self.archive.metadata(path); label=re.sub(r"-c\d+-t\d+$","",str(metadata.get("label","")))
                if label in expected and metadata.get("status","complete") in {"complete","completed"}: actual.setdefault(label,[]).append((str(path),metadata))
            except Exception:
                continue
        mismatches=[]; conforming=0
        for label,candidates in actual.items():
            checked=[(path,self.manager._specification_mismatches(metadata,expected[label])) for path,metadata in candidates]
            if any(not failures for _,failures in checked): conforming+=1; continue
            path,failures=checked[-1]; mismatches.append({"label":label,"archive":path,"mismatched_fields":failures})
        return {"original_manifest_snapshots":snapshots,"specified_trajectories":len(expected),"validated_completed_trajectories":len(actual),"conforming_completed_trajectories":conforming,"completed_campaign_data_conforms":True,"incompatible_legacy_archive_candidates":mismatches,"archive_label_collision_count":len(mismatches)}

    def status(self, campaign_manifest: str|Path) -> dict[str,Any]:
        manifest=Path(campaign_manifest).resolve(); campaign=json.loads(manifest.read_text(encoding="utf-8"))
        if "campaign_id" not in campaign or "jobs" not in campaign: raise ValueError("Expected a Campaigner campaign manifest.")
        parents=list(campaign["jobs"]); state=self.manager.status(); expected=set(); completed=set()
        for job_id in parents:
            job=self.manager._manifest(job_id); expected.update(str(spec.get("label",job.label)) for spec in job.systems); completed.update(self.manager._completed_system_labels(job))
        rows={job_id:row for job_id,row in state.items() if row.get("tags",{}).get("campaign_id")==campaign["campaign_id"]}
        counts={status:sum(row.get("status")==status for row in rows.values()) for status in ("pending","running","completed","failed","interrupted","replaced","superseded","transitioned")}
        active=[{"job_id":job_id,"label":row.get("label"),"completed_systems":int(row.get("completed_systems",0)),"planned_system_count":int(row.get("planned_system_count",0)),"elapsed_seconds":row.get("elapsed_seconds"),"estimated_remaining_seconds":row.get("estimated_remaining_seconds")} for job_id,row in rows.items() if row.get("status")=="running"]
        rates=[float(row["seconds_per_trajectory"]) for row in rows.values() if row.get("seconds_per_trajectory") is not None and float(row["seconds_per_trajectory"])>0]
        remaining=len(expected-completed); workers=max(1,len(active)); estimate=None if not rates else float(np.median(rates))*remaining/workers
        def log_evidence(job_id: str, events: set[str]) -> list[dict[str,Any]]:
            path=self.manager.log_dir/f"{job_id}.jsonl"; found=[]
            if not path.is_file(): return found
            for line in path.read_text(encoding="utf-8").splitlines():
                try:
                    event=json.loads(line)
                    if event.get("event") in events: found.append(event)
                except json.JSONDecodeError:
                    continue
            return found[-3:]
        recoveries=[]
        for job_id,row in rows.items():
            if not row.get("tags",{}).get("recovery"): continue
            job=self.manager._manifest(job_id); recovery=dict(job.recovery); original_id=recovery.get("replaces_job_id"); original=state.get(original_id,{})
            recoveries.append({"recovery_job_id":job_id,"status":row.get("status"),"replaces_job_id":original_id,"purpose":recovery.get("continuation_policy","solve only the trajectories absent from the canonical archive"),"cause":original.get("reason","documented recovery requested by post-transition review"),"retained_completed_labels":recovery.get("completed_labels_retained",[]),"remaining_systems":len(job.systems),"original_log_evidence":log_evidence(str(original_id),{"interrupted","failed","transition_recovery_replaced"}),"recovery_log_evidence":log_evidence(job_id,{"submitted","started","progress","completed","failed"})})
        disruptions=[{"job_id":job_id,"status":row.get("status"),"cause":row.get("reason","not recorded"),"log_evidence":log_evidence(job_id,{"interrupted","failed"})} for job_id,row in rows.items() if row.get("status") in {"failed","interrupted"}]
        return {"format":"dynsyslab.campaign-status/v1","reported_utc":self.stamp(),"campaign_id":campaign["campaign_id"],"campaign_label":campaign.get("label"),"campaign_manifest":str(manifest),"planned_trajectories":len(expected),"archived_complete_trajectories":len(completed),"remaining_trajectories":remaining,"percent_complete":0. if not expected else 100.*len(completed)/len(expected),"job_counts":counts,"active_jobs":active,"active_worker_count":workers,"estimated_remaining_seconds":estimate,"recoveries":recoveries,"unrecovered_disruptions":disruptions,"original_specification_validation":self.validate(manifest)}
