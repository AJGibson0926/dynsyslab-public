"""Parallel transition planning, review, and recovery services."""
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


class _JobTransitionMixin:
    """Parallel transition planning, review, and recovery services."""

    def evaluate_parallelization(self, job_ids: Iterable[str] | None=None, *, benchmark: bool=True,
                                 candidate_backends: tuple[str,...]=("local_process","local_worker","ssh","slurm","pbs","cloud_batch")) -> dict[str,Any]:
        """Assess current resources without touching an active numerical solve."""
        ids=list(job_ids or self.status().keys()); report=evaluate_resources(self,ids,candidate_backends=candidate_backends,benchmark=benchmark)
        path=write_resource_report(self,report); payload=report.payload(); payload["report_path"]=str(path)
        return payload

    def evaluate_proposal(self, collection: JobCollection, *, campaign_manifest: str | Path | None=None,
                          seconds_per_trajectory: float | None=None) -> dict[str, Any]:
        """Return a strictly read-only estimate for an unsubmitted job proposal.

        This deliberately does not call :meth:`submit`, create a state row,
        schedule a backend, or invoke the solver.  It may read the archive
        index to identify exact reusable work, but a proposal is never altered
        on that basis.  Campaigner can supply a declared planning-rate
        assumption through its manifest; callers can override it explicitly.
        """
        campaign={}
        if campaign_manifest is not None:
            path=Path(campaign_manifest).resolve()
            campaign=json.loads(path.read_text(encoding="utf-8"))
        metadata=dict(campaign.get("metadata",{}))
        specs=[spec for job in collection.jobs for spec in job.systems]
        if not specs: raise ValueError("A proposed collection must contain at least one planned system.")
        systems=[]; final_samples=[]; raw_bytes=0; reusable=[]; extendible=[]
        for spec in specs:
            start,stop=map(float,spec["time_span"]); output=dict(spec.get("output",{}))
            if "times" in output:
                samples=len(output["times"])
            else:
                delta=float(output["delta_t"]); samples=int(np.ceil((stop-start)/delta))+1
            # Loading a library object is configuration validation only; it
            # performs neither integration nor archive mutation.
            from app import RunConfiguration
            request=RunConfiguration(system_name=spec["system_name"],preset=spec.get("preset"),
                                     parameters=dict(spec.get("parameters",{})),driver=dict(spec.get("driver",{"kind":"zero"})),
                                     time_span=(start,stop),label=str(spec.get("label","proposal"))).build_request()
            dimension=int(np.asarray(request.system.initial_state,dtype=float).size)
            final_samples.append(samples); raw_bytes+=samples*(dimension+1)*8
            candidate=MANAGER_ARCHIVE.efficiency_candidate(spec)
            if candidate and candidate.get("kind")=="reuse": reusable.append({"label":spec.get("label"),"archive_path":candidate["archive_path"]})
            elif candidate and candidate.get("kind")=="extend": extendible.append({"label":spec.get("label"),"archive_path":candidate["archive_path"]})
            systems.append({"system_name":spec["system_name"],"state_dimension":dimension,"time_span":[start,stop],
                            "driver":dict(spec.get("driver",{})),"solver":dict(spec.get("solver",{}))})
        assumed_rate=seconds_per_trajectory
        if assumed_rate is None: assumed_rate=metadata.get("planning_seconds_per_trajectory")
        assumed_rate=None if assumed_rate is None else float(assumed_rate)
        new_count=len(specs)-len(reusable)
        # NPZ compression is data-dependent.  A deliberately conservative
        # no-calibration estimate preserves room for metadata and raw adaptive
        # history without pretending it is a measured result.
        compressed_estimate=int(raw_bytes*float(metadata.get("planning_compression_factor",0.70)))
        return {"format":"dynsyslab.campaign-planning-report/v1","generated_utc":_stamp(),"read_only":True,
                "scheduled":False,"executed":False,"campaign_manifest":None if campaign_manifest is None else str(Path(campaign_manifest).resolve()),
                "campaign_id":campaign.get("campaign_id"),"label":campaign.get("label",collection.label),"goal":campaign.get("goal"),
                "campaign_type":metadata.get("campaign_type"),"system":metadata.get("system_name"),"preset":metadata.get("preset"),
                "hierarchy":metadata.get("hierarchy"),"job_count":len(collection.jobs),"planned_trajectory_count":len(specs),"new_trajectory_count":new_count,
                "exact_reuse_candidates":reusable,"extension_candidates":extendible,
                "solver_options":metadata.get("solver_options",{}),"time_span":metadata.get("time_span"),"systems":systems,
                "sampling":{"samples_per_cycle":metadata.get("samples_per_cycle"),"delta_t":metadata.get("delta_t"),
                            "transient_cycles":metadata.get("transient_cycles"),"retained_cycles":metadata.get("retained_cycles"),
                            "estimated_output_samples_per_trajectory":max(final_samples),"minimum_samples_per_trajectory":min(final_samples),
                            "total_final_samples":sum(final_samples)},
                "parameter_space":metadata.get("parameter_space",{}),
                "storage":{"final_time_state_bytes_per_trajectory":raw_bytes//len(specs),"final_time_state_bytes_total":raw_bytes,
                           "compressed_projection":{"mean_total_bytes":compressed_estimate,"p90_total_bytes":compressed_estimate,
                                                    "method":"Conservative planning factor; no proposed solve was run."},
                           "note":"Planning projection only; adaptive raw histories and NPZ compressibility may change the final size."},
                "runtime":{"seconds_per_trajectory_assumption":assumed_rate,
                           "assumption_note":"Declared campaign planning assumption; no proposed solve or timing benchmark was run.",
                           "estimated_total_seconds":None if assumed_rate is None else assumed_rate*new_count,
                           "estimated_remaining_seconds_at_start":None if assumed_rate is None else assumed_rate*new_count,
                           "estimated_parallel_seconds_at_2_workers":None if assumed_rate is None else assumed_rate*new_count/2,
                           "basis":"Declared campaign planning assumption; no proposed solve or timing benchmark was run."}}

    def _completed_system_labels(self, job: JobManifest) -> set[str]:
        """Discover compatible completed systems through the archive inspector."""
        return self._archive_inspector.completed_system_labels(job)

    @staticmethod
    def _specification_mismatches(metadata: dict[str, Any], spec: dict[str, Any]) -> list[str]:
        """Fields that must agree before an archive may satisfy a planned solve."""
        failures=[]; driver=spec.get("driver",{}); source=((metadata.get("forcing",{}).get("sources") or [{}])[0].get("parameters",{}))
        for key in ("amplitude","frequency_hz"):
            if key in driver and not np.isclose(float(driver[key]),float(source.get(key,float("nan"))),rtol=1e-12,atol=1e-12): failures.append(f"driver.{key}")
        if tuple(float(x) for x in metadata.get("time_span",[])) != tuple(float(x) for x in spec.get("time_span",[])): failures.append("time_span")
        expected_dim=bool(spec.get("solver",{}).get("options",{}).get("use_dimensionless_form",False))
        if bool(metadata.get("options",{}).get("use_dimensionless_form",False)) != expected_dim: failures.append("solver.use_dimensionless_form")
        return failures

    def controlled_parallel_transition(self, job_ids: Iterable[str] | None=None, *, backend: str="local_process",
                                       strategy: str="fixed", max_workers: int | None=None, systems_per_child: int | None=None,
                                       apply: bool=False, campaign_path: str | Path | None=None,
                                       allow_descendant_repartition: bool=False) -> dict[str,Any]:
        """Prepare or apply a duplication-safe serial-to-parallel transition.

        ``apply=False`` returns a complete, read-only plan.  Applying a
        transition is refused while any selected serial Job is marked running,
        or when a partial archive is present.  Completed systems are detected
        directly from primary archives and excluded from every child Job.
        """
        if strategy not in {"fixed","adaptive_local","adaptive_portable"}: raise ValueError("strategy must be fixed, adaptive_local, or adaptive_portable.")
        if backend not in {"local_process","local_worker","ssh","slurm","pbs","cloud_batch"}: raise ValueError("Choose an explicit parallel-capable backend.")
        ids=list(job_ids or self.status().keys()); rows={item:self.status(item) for item in ids}
        running=[item for item,row in rows.items() if row.get("status")=="running"]
        partial=MANAGER_ARCHIVE.partial_solution_paths()
        if running or partial:
            return {"safe":False,"applied":False,"reason":"serialized execution has not demonstrably halted","running_job_ids":running,"partial_archives":[str(item) for item in partial],"new_job_ids":[]}
        already_transitioned=[item for item,row in rows.items() if row.get("parallel_transition_manifest") or row.get("status")=="transitioned"]
        descendants=[item for item,row in rows.items() if row.get("tags",{}).get("parallel_transition") or row.get("tags",{}).get("logical_parent_job")]
        if already_transitioned:
            return {"safe":False,"applied":False,"reason":"one or more selected Jobs already have an immutable parallel-transition lineage","already_transitioned_job_ids":already_transitioned,"new_job_ids":[]}
        if descendants and not allow_descendant_repartition:
            return {"safe":False,"applied":False,"reason":"refusing to re-partition transition descendants without explicit allow_descendant_repartition=True","descendant_job_ids":descendants,"new_job_ids":[]}
        workers=int(max_workers or max(1,min(4,(os.cpu_count() or 2)-1))); plans=[]; all_new=[]; globally_planned=set()
        for job_id in ids:
            row=rows[job_id]
            if row.get("status")=="completed": continue
            job=self._manifest(job_id); completed=self._completed_system_labels(job)
            labels=[str(spec.get("label",job.label)) for spec in job.systems]
            if len(labels)!=len(set(labels)): return {"safe":False,"applied":False,"reason":f"planned Job {job_id} contains non-unique system labels","new_job_ids":[]}
            remaining=[spec for spec in job.systems if str(spec.get("label",job.label)) not in completed]
            overlap=globally_planned.intersection(str(spec.get("label",job.label)) for spec in remaining)
            if overlap:
                return {"safe":False,"applied":False,"reason":"selected parent Jobs overlap scientifically identical trajectory labels","duplicate_labels":sorted(overlap),"new_job_ids":[]}
            globally_planned.update(str(spec.get("label",job.label)) for spec in remaining)
            chunk=int(systems_per_child or max(1,(len(remaining)+workers*2-1)//(workers*2)))
            plans.append({"job_id":job_id,"completed_system_labels":sorted(completed),"remaining_systems":len(remaining),"systems_per_child":chunk})
            if not apply: continue
            children=[]
            for index in range(0,len(remaining),chunk):
                child=replace(job,label=f"{job.label}-parallel-{index//chunk:03d}",systems=tuple(remaining[index:index+chunk]),job_id=uuid4().hex,
                              execution={**job.execution,"backend":backend,"adaptive_policy":strategy,"max_systems_per_job":chunk},
                              tags={**job.tags,"logical_parent_job":job.job_id,"parallel_transition":True,"parallel_backend":backend,"parallel_strategy":strategy},
                              notes=f"{job.notes}\nParallel transition from {job.job_id}; completed system labels excluded.".strip())
                self.submit(child); children.append(child.job_id); all_new.append(child.job_id)
            manifest=self.restructure_dir/f"parallel-transition__{job.job_id}.json"
            transition_payload={"format":"dynsyslab.parallel-transition/v1","created_utc":_stamp(),"parent_job_id":job.job_id,"backend":backend,"strategy":strategy,"completed_system_labels":sorted(completed),"child_job_ids":children,"safety":"primary archive labels were checked before creating child manifests; descendant re-partitioning is refused by default."}
            temporary=manifest.with_suffix(manifest.suffix+f".{uuid4().hex}.partial")
            with temporary.open("w",encoding="utf-8",newline="\n") as handle:
                handle.write(json.dumps(transition_payload,indent=2,sort_keys=True)); handle.flush(); os.fsync(handle.fileno())
            os.replace(temporary,manifest)
            self._set_status(job.job_id,"transitioned",parallel_transition_manifest=str(manifest),parallel_child_job_ids=children,completed_system_labels=sorted(completed))
            self._event(job.job_id,"parallel_transitioned",backend=backend,strategy=strategy,child_job_ids=children,completed_system_count=len(completed))
        result={"safe":True,"applied":apply,"backend":backend,"strategy":strategy,"recommended_workers":workers,"plans":plans,"new_job_ids":all_new}
        if apply:
            resolved_campaign=campaign_path or self._campaign_path_for(ids)
            self._append_transition_report(resolved_campaign,result)
            review=self.review_parallel_transition(ids,campaign_path=resolved_campaign)
            result["post_transition_review"]=review
            if not review["safe_to_execute"]:
                for child_id in all_new:
                    self._set_status(child_id,"held_for_review",reason="post-transition review found a safety issue",review_path=review["report_path"])
                    self._event(child_id,"held_for_review",review_path=review["report_path"])
        return result

    def review_parallel_transition(self, parent_job_ids: Iterable[str], *, campaign_path: str | Path | None=None) -> dict[str, Any]:
        """Perform the required post-transition, archive-first safety review.

        This review is deliberately read-only with respect to primary solution
        data.  It validates the transition ledger, hierarchy manifests, child
        plan labels, and canonical archive coverage before any child job runs.
        """
        parents=list(parent_job_ids); state=self.status(); children=[]; issues=[]; expected=set(); planned=set(); completed=set()
        for parent_id in parents:
            row=state.get(parent_id,{})
            transition_path=row.get("parallel_transition_manifest")
            if not transition_path or not Path(transition_path).is_file():
                issues.append({"parent_job_id":parent_id,"issue":"missing_transition_manifest"}); continue
            ledger=json.loads(Path(transition_path).read_text(encoding="utf-8")); child_ids=list(ledger.get("child_job_ids",[])); children.extend(child_ids)
            parent=self._manifest(parent_id)
            parent_labels={str(spec.get("label",parent.label)) for spec in parent.systems}; expected.update(parent_labels)
            ledger_completed=set(ledger.get("completed_system_labels",[])); archive_completed=self._completed_system_labels(parent)
            if not ledger_completed.issubset(archive_completed):
                issues.append({"parent_job_id":parent_id,"issue":"ledger_completed_labels_not_in_archive","labels":sorted(ledger_completed-archive_completed)})
            completed.update(archive_completed)
            for child_id in child_ids:
                child_row=state.get(child_id,{})
                if not child_row or not child_row.get("manifest"):
                    issues.append({"parent_job_id":parent_id,"child_job_id":child_id,"issue":"missing_child_manifest"}); continue
                child=self._manifest(child_id)
                labels=[str(spec.get("label",child.label)) for spec in child.systems]
                overlap=planned.intersection(labels)
                if overlap: issues.append({"parent_job_id":parent_id,"child_job_id":child_id,"issue":"duplicate_child_labels","labels":sorted(overlap)})
                planned.update(labels)
                if any(label in archive_completed for label in labels):
                    issues.append({"parent_job_id":parent_id,"child_job_id":child_id,"issue":"child_repeats_archived_solution","labels":sorted(set(labels)&archive_completed)})
                if child_row.get("status")=="completed" and child_row.get("hierarchy_archive"):
                    try:
                        if not self._hierarchy_matches(child_row["hierarchy_archive"],self._build_target(child)):
                            issues.append({"parent_job_id":parent_id,"child_job_id":child_id,"issue":"completed_child_hierarchy_mismatch"})
                    except Exception as error:
                        issues.append({"parent_job_id":parent_id,"child_job_id":child_id,"issue":"completed_child_hierarchy_unreadable","reason":str(error)})
        missing=expected-(completed|planned)
        excess=planned-expected
        if missing: issues.append({"issue":"uncovered_parent_labels","labels":sorted(missing)})
        if excess: issues.append({"issue":"child_labels_not_in_parent","labels":sorted(excess)})
        report={"format":"dynsyslab.parallel-transition-review/v1","reviewed_utc":_stamp(),"parent_job_ids":parents,"child_job_ids":children,"expected_label_count":len(expected),"archived_complete_label_count":len(completed),"planned_child_label_count":len(planned),"issues":issues,"safe_to_execute":not issues}
        path=self.review_dir/f"parallel-transition-review__{uuid4().hex}.json"; path.write_text(json.dumps(report,indent=2,sort_keys=True),encoding="utf-8"); report["report_path"]=str(path)
        for parent_id in parents: self._event(parent_id,"parallel_transition_review",safe_to_execute=report["safe_to_execute"],issue_count=len(issues),review_path=str(path))
        self._append_transition_review_report(campaign_path or self._campaign_path_for(parents),report)
        return report

    def _append_transition_review_report(self, campaign_path: str | Path | None, review: dict[str, Any]) -> None:
        if not campaign_path or not Path(campaign_path).is_file(): return
        campaign=json.loads(Path(campaign_path).read_text(encoding="utf-8")); base=(MANAGER_ARCHIVE.archive_root/"campaigns"/"reports")/f"{_safe(campaign['label'])}__{campaign['campaign_id']}"; markdown=base.with_suffix(".md"); payload=base.with_suffix(".json")
        line=(f"\n### Post-transition review - {review['reviewed_utc']}\n"
              f"- Parents: {len(review['parent_job_ids'])}; child Jobs: {len(review['child_job_ids'])}.\n"
              f"- Archive-complete labels: {review['archived_complete_label_count']}; planned remaining labels: {review['planned_child_label_count']}.\n"
              f"- Result: {'safe to execute' if review['safe_to_execute'] else str(len(review['issues'])) + ' issue(s); execution held'}.\n")
        if markdown.is_file(): markdown.write_text(markdown.read_text(encoding="utf-8")+line,encoding="utf-8")
        if payload.is_file():
            data=json.loads(payload.read_text(encoding="utf-8")); data.setdefault("post_transition_reviews",[]).append(review); payload.write_text(json.dumps(data,indent=2,sort_keys=True),encoding="utf-8")

    def repair_transition_interruptions(self, parent_job_ids: Iterable[str], *, campaign_path: str | Path | None=None,
                                        reason: str="post-transition safety review") -> dict[str, Any]:
        """Replace interrupted child plans with archive-aware recovery manifests.

        Finished trajectory archives are retained and excluded.  A replacement
        manifest contains only absent trajectories and points back to the
        interrupted child and original logical parent, preserving the intended
        hierarchy without copying or overwriting any solution content.  A true
        solver checkpoint (``*.partial.npz``) is never overwritten; it is held
        for the solver's continuation pathway instead of being restarted.
        """
        parents=list(parent_job_ids); state=self.status(); replacements=[]; held=[]
        partials={str(path) for path in MANAGER_ARCHIVE.partial_solution_paths()}
        for parent_id in parents:
            row=state.get(parent_id,{}); transition=row.get("parallel_transition_manifest")
            if not transition or not Path(transition).is_file(): continue
            ledger=json.loads(Path(transition).read_text(encoding="utf-8"))
            for child_id in ledger.get("child_job_ids",[]):
                child_row=self.status(child_id); status=child_row.get("status")
                if status not in {"running","interrupted","failed"}: continue
                child=self._manifest(child_id)
                checkpoint_paths=[item for item in child_row.get("generated_solution_archives",[]) if str(item) in partials]
                if checkpoint_paths:
                    held.append({"child_job_id":child_id,"checkpoint_paths":checkpoint_paths,"reason":"partial solver checkpoint requires continuation-aware solver execution"})
                    self._set_status(child_id,"interrupted",reason="held for continuation from partial solver checkpoint",checkpoint_paths=checkpoint_paths)
                    self._event(child_id,"transition_recovery_held",checkpoint_paths=checkpoint_paths); continue
                completed=self._completed_system_labels(child)
                remaining=tuple(spec for spec in child.systems if str(spec.get("label",child.label)) not in completed)
                if not remaining:
                    self._set_status(child_id,"completed",recovered_without_resolve=True,completed_systems=len(child.systems))
                    self._event(child_id,"transition_recovery_not_needed",completed_system_count=len(child.systems)); continue
                replacement=replace(child,job_id=uuid4().hex,label=f"{child.label}-recovery",systems=remaining,
                    recovery={"kind":"archive-aware-recovery","reason":reason,"replaces_job_id":child_id,"logical_parent_job":parent_id,
                              "completed_labels_retained":sorted(completed),"continuation_policy":"extend_existing_checkpoint_else_solve_only_missing"},
                    tags={**child.tags,"recovery":True,"replaces_job_id":child_id,"logical_parent_job":parent_id},
                    notes=f"{child.notes}\nRecovery manifest: completed child trajectories are retained; only absent labels are solved.".strip())
                self.submit(replacement)
                self._set_status(child_id,"replaced",reason=reason,recovery_job_id=replacement.job_id,completed_system_labels=sorted(completed))
                self._event(child_id,"transition_recovery_replaced",recovery_job_id=replacement.job_id,retained_completed_labels=sorted(completed),remaining_systems=len(remaining))
                replacements.append({"interrupted_job_id":child_id,"recovery_job_id":replacement.job_id,"retained_completed_labels":sorted(completed),"remaining_systems":len(remaining)})
        report={"format":"dynsyslab.transition-recovery/v1","created_utc":_stamp(),"parent_job_ids":parents,"reason":reason,"replacements":replacements,"held_checkpoints":held}
        path=self.review_dir/f"transition-recovery__{uuid4().hex}.json"; path.write_text(json.dumps(report,indent=2,sort_keys=True),encoding="utf-8"); report["report_path"]=str(path)
        for parent_id in parents: self._event(parent_id,"transition_recovery_plan",recovery_path=str(path),replacement_count=len(replacements),held_checkpoint_count=len(held))
        return report

    def repair_campaign_coverage(self, campaign_manifest: str | Path) -> dict[str, Any]:
        """Create only the Jobs needed to cover specs not met by conforming archives."""
        manifest=Path(campaign_manifest).resolve(); campaign=json.loads(manifest.read_text(encoding="utf-8")); state=self.status(); campaign_id=campaign["campaign_id"]
        assigned=set()
        for job_id,row in state.items():
            if row.get("tags",{}).get("campaign_id") != campaign_id or row.get("status") in {"superseded","transitioned"}: continue
            try: assigned.update(str(spec.get("label")) for spec in self._manifest(job_id).systems)
            except Exception: continue
        repairs=[]
        for parent_id in campaign["jobs"]:
            parent=JobManifest.read(state[parent_id].get("original_manifest") or self.snapshot_original_manifest(parent_id))
            completed=self._completed_system_labels(parent)
            missing=tuple(spec for spec in parent.systems if str(spec.get("label")) not in completed|assigned)
            if not missing: continue
            repair=replace(parent,job_id=uuid4().hex,label=f"{parent.label}-specification-repair",systems=missing,
                           recovery={"kind":"specification-conformity-repair","logical_parent_job":parent_id,"purpose":"replace label-matched archives that fail original solver specification validation"},
                           tags={**parent.tags,"recovery":True,"coverage_repair":True,"logical_parent_job":parent_id},
                           notes=f"Specification repair of {parent_id}: only original labels not covered by conforming archive data are solved.")
            self.submit(repair); repairs.append({"parent_job_id":parent_id,"repair_job_id":repair.job_id,"missing_systems":len(missing)})
            self._event(parent_id,"specification_coverage_repair",repair_job_id=repair.job_id,missing_systems=len(missing))
        report={"format":"dynsyslab.specification-coverage-repair/v1","created_utc":_stamp(),"campaign_id":campaign_id,"repairs":repairs}
        path=self.review_dir/f"specification-coverage-repair__{uuid4().hex}.json"; path.write_text(json.dumps(report,indent=2,sort_keys=True),encoding="utf-8"); report["report_path"]=str(path)
        return report

    def _campaign_path_for(self, job_ids: list[str]) -> str | None:
        for job_id in job_ids:
            path=self.status(job_id).get("tags",{}).get("campaign_manifest")
            if path and Path(path).is_file(): return str(path)
        return None

    def _append_transition_report(self, campaign_path: str | Path | None, transition: dict[str,Any]) -> None:
        if not campaign_path or not Path(campaign_path).is_file(): return
        campaign=json.loads(Path(campaign_path).read_text(encoding="utf-8")); reports=MANAGER_ARCHIVE.archive_root/"campaigns"/"reports"; base=reports/f"{_safe(campaign['label'])}__{campaign['campaign_id']}"; markdown=base.with_suffix(".md"); payload=base.with_suffix(".json")
        stamp=_stamp(); section=f"\n## Parallelization transition — {stamp}\n\n- Backend: {transition['backend']}\n- Strategy: {transition['strategy']}\n- Requested worker processes: {transition['recommended_workers']}\n- New child Jobs: {len(transition['new_job_ids'])}\n- Safety: serialized execution was halted; completed archive labels were excluded before child jobs were created.\n"
        if markdown.is_file(): markdown.write_text(markdown.read_text(encoding="utf-8")+section,encoding="utf-8")
        if payload.is_file():
            data=json.loads(payload.read_text(encoding="utf-8")); data.setdefault("parallelization_transitions",[]).append({"timestamp":stamp,**transition}); payload.write_text(json.dumps(data,indent=2,sort_keys=True),encoding="utf-8")

    def _append_reparallelization_report(self, campaign_path: str | Path | None, event: dict[str, Any]) -> None:
        """Append a deliberately compact record for a later adaptive worker change."""
        if not campaign_path or not Path(campaign_path).is_file(): return
        campaign=json.loads(Path(campaign_path).read_text(encoding="utf-8")); reports=MANAGER_ARCHIVE.archive_root/"campaigns"/"reports"; base=reports/f"{_safe(campaign['label'])}__{campaign['campaign_id']}"; markdown=base.with_suffix(".md"); payload=base.with_suffix(".json")
        section=(f"\n### Re-parallelization - {event['timestamp']}\n"
                 f"- {event['backend']}: {event['prior_workers']} -> {event['new_workers']} workers.\n"
                 f"- Safety: {event['safety']}.\n")
        if markdown.is_file(): markdown.write_text(markdown.read_text(encoding="utf-8")+section,encoding="utf-8")
        if payload.is_file():
            data=json.loads(payload.read_text(encoding="utf-8")); data.setdefault("adaptive_reparallelizations",[]).append(event)
            payload.write_text(json.dumps(data,indent=2,sort_keys=True),encoding="utf-8")
