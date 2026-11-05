"""Exact archive reuse and compatible single-solution extension services."""
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


class _JobReuseMixin:
    """Exact archive reuse and compatible single-solution extension services."""

    def _completed_plan_archive(self, fingerprint: str) -> str | None:
        """Find a completed hierarchy archive produced by the identical plan."""
        for row in self._read_state()["jobs"].values():
            if row.get("status") != "completed" or not row.get("hierarchy_archive"): continue
            candidate=row.get("plan_fingerprint")
            if candidate is None and row.get("manifest"):
                try: candidate=JobManifest.read(row["manifest"]).plan_fingerprint
                except Exception: continue
            if candidate == fingerprint: return str(row["hierarchy_archive"])
        return self._archive_inspector.completed_plan_archive(fingerprint,self._read_state())

    def _try_single_solution_extension(self, job: JobManifest) -> dict[str, Any] | None:
        """Extend one fully compatible archived trajectory instead of restarting it.

        Multi-system Jobs are intentionally left to normal hierarchy execution;
        Campaigner/JobManager can split them into singleton extension Jobs when
        that optimization is worthwhile.
        """
        from app import RunConfiguration
        from solver import HistoryOptions, InterpolationOptions, SolverOptions, solve
        if len(job.systems)!=1: return None
        spec=job.systems[0]; candidate=MANAGER_ARCHIVE.efficiency_candidate(spec)
        if not candidate or candidate["kind"]!="extend": return None
        prior=MANAGER_ARCHIVE.read_time_series(candidate["archive_path"]); prior_end=float(prior.time[-1]); final=float(spec["time_span"][-1])
        if prior_end >= final: return None
        options=SolverOptions(**dict(spec.get("solver",{}).get("options",{}))); history=HistoryOptions(**dict(spec.get("solver",{}).get("history",{})))
        output=spec.get("output",{}); step=float(output.get("delta_t",prior.delta_t or 0.))
        if step<=0: return None
        count=int(np.floor((final-prior_end)/step)); times=prior_end+step*np.arange(count+1,dtype=float)
        if not np.isclose(times[-1],final,rtol=1e-12,atol=1e-15): times=np.append(times,final)
        config=RunConfiguration(system_name=spec["system_name"],preset=spec.get("preset"),parameters=dict(spec.get("parameters",{})),driver=dict(spec.get("driver",{"kind":"zero"})),time_span=(prior_end,final),initial_state=tuple(np.asarray(prior.values[-1],dtype=float)),label=str(spec.get("label",job.label)))
        request=config.build_request()
        request=replace(request,options=options,history=history,interpolation=InterpolationOptions(times=times,method=str(spec.get("solver",{}).get("interpolation_method","pchip"))),label=str(spec.get("label",job.label)))
        result=solve(request)
        if result.status!="complete" or not result.archive_path: raise RuntimeError(f"continuation solve failed: {result.message}")
        extended=MANAGER_ARCHIVE.extend_solution(candidate["archive_path"],result.archive_path)
        MANAGER_ARCHIVE.update_metadata(extended,{"job":{"job_id":job.job_id,"label":job.label,"tags":job.tags,"extension_of":candidate["archive_path"]},"extension_execution":{"continuation_archive":result.archive_path,"requested_final_time":final}})
        return {"extended_archive":str(extended),"continuation_archive":str(result.archive_path),"prior_end":prior_end,"final_time":final}

    def _partition_reusable_specs(self, job: JobManifest) -> tuple[list[dict[str,Any]],list[dict[str,Any]]]:
        """Split a planned hierarchy into exact archive reuse and new solves.

        This check occurs immediately before execution, not only at submission:
        a parallel sibling or a recovered parent may have completed a matching
        trajectory after this Job was originally queued.  Only exact,
        scientifically compatible matches returned by Archiver are reused.
        """
        reused=[]; remaining=[]
        for spec in job.systems:
            candidate=MANAGER_ARCHIVE.efficiency_candidate(spec)
            if candidate is not None and candidate.get("kind")=="reuse":
                reused.append({"label":str(spec.get("label",job.label)),"archive_path":str(candidate["archive_path"])})
            else:
                remaining.append(spec)
        return reused,remaining

    def _reuse_manifest(self, job: JobManifest, reused: list[dict[str,Any]], remaining: list[dict[str,Any]]) -> Path:
        """Persist intended hierarchy membership without copying solution data."""
        path=self.restructure_dir/f"reuse-partition__{_safe(job.label)}__{job.job_id}.json"
        payload={"format":"dynsyslab.job-reuse-partition/v1","created_utc":_stamp(),"job_id":job.job_id,
                 "original_manifest":str(self.status(job.job_id).get("original_manifest",self.status(job.job_id).get("manifest",""))),
                 "intended_system_labels":[str(spec.get("label",job.label)) for spec in job.systems],
                 "reused":reused,"remaining_system_labels":[str(spec.get("label",job.label)) for spec in remaining],
                 "note":"Canonical trajectories are referenced by path; no numerical arrays are duplicated."}
        temporary=path.with_suffix(path.suffix+f".{uuid4().hex}.partial")
        with temporary.open("w",encoding="utf-8",newline="\n") as handle:
            handle.write(json.dumps(payload,indent=2,sort_keys=True)); handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary,path)
        return path
