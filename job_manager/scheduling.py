"""Backend-neutral scheduler contracts: resources, arrays, dependencies, and provenance."""
from __future__ import annotations
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Literal
from uuid import uuid4
import os, platform, sys

DependencyKind=Literal["afterok","afterany","afternotok","aftercorr"]

@dataclass(frozen=True)
class ResourceRequest:
    """Portable resource requirements for one scheduled job."""
    cores: int=1; memory_mb: int|None=None; wall_time_s: float|None=None
    gpus: int=0; gpu_type: str|None=None; nodes: int|None=None; scratch_mb: int|None=None
    features: tuple[str,...]=()
    def __post_init__(self):
        if self.cores<1 or self.gpus<0 or (self.nodes is not None and self.nodes<1): raise ValueError("Resource counts must be positive (GPUs may be zero).")
        if any(value is not None and value<=0 for value in (self.memory_mb,self.wall_time_s,self.scratch_mb)): raise ValueError("Resource limits must be positive when supplied.")
    def payload(self)->dict[str,Any]:
        return {"cores":self.cores,"memory_mb":self.memory_mb,"wall_time_s":self.wall_time_s,"gpus":self.gpus,"gpu_type":self.gpu_type,"nodes":self.nodes,"scratch_mb":self.scratch_mb,"features":list(self.features)}

@dataclass(frozen=True)
class Dependency:
    """One scheduler dependency and its completion condition."""
    job_id: str; condition: DependencyKind="afterok"
    def __post_init__(self):
        if self.condition not in {"afterok","afterany","afternotok","aftercorr"}: raise ValueError(f"Unknown dependency condition {self.condition!r}.")

@dataclass(frozen=True)
class DeadlinePolicy:
    """Checkpoint and requeue behavior near a scheduler deadline."""
    checkpoint_before_deadline_s: float|None=None; requeue_on_timeout: bool=False; max_requeues: int=0
    def __post_init__(self):
        if self.checkpoint_before_deadline_s is not None and self.checkpoint_before_deadline_s<=0: raise ValueError("checkpoint_before_deadline_s must be positive.")
        if self.max_requeues<0: raise ValueError("max_requeues cannot be negative.")

@dataclass(frozen=True)
class StagingPlan:
    """Input/output movement policy for local or remote execution."""
    inputs: tuple[str,...]=(); output_location: str="primary"; stage_out: bool=True; cleanup_scratch: bool=True

@dataclass(frozen=True)
class SchedulerPolicy:
    """Local admission limits; remote backends may translate these to QoS."""
    max_running_jobs: int|None=None; max_cores: int|None=None; max_memory_mb: int|None=None; max_gpus: int|None=None; allowed_qos: tuple[str,...]=(); reservation_windows: dict[str,tuple[str,str]]=field(default_factory=dict)
    def __post_init__(self):
        if any(value is not None and value<1 for value in (self.max_running_jobs,self.max_cores,self.max_memory_mb,self.max_gpus)): raise ValueError("Scheduler limits must be positive when supplied.")

@dataclass(frozen=True)
class JobArray:
    """Durable logical array with stable task indices and a concurrency cap."""
    label: str; jobs: tuple[Any,...]; max_concurrent: int=1; array_id: str=field(default_factory=lambda:uuid4().hex)
    def __post_init__(self):
        if not self.jobs or self.max_concurrent<1: raise ValueError("A JobArray needs tasks and a positive concurrency limit.")
        ids=[item.job_id for item in self.jobs]
        if len(ids)!=len(set(ids)): raise ValueError("JobArray task job IDs must be unique.")
    def payload(self)->dict[str,Any]:
        return {"format":"dynsyslab.job-array/v1","array_id":self.array_id,"label":self.label,"max_concurrent":self.max_concurrent,"job_ids":[item.job_id for item in self.jobs]}

def execution_provenance() -> dict[str,Any]:
    """Compact reproducibility facts captured at dispatch, not numerical data."""
    return {"captured_utc":datetime.now(timezone.utc).isoformat(),"python":sys.version,"executable":sys.executable,"platform":platform.platform(),"cwd":os.getcwd()}

def dependency_satisfied(condition: str, status: str|None) -> bool:
    """Evaluate a scheduler dependency condition against one job status."""
    if condition=="afterok": return status=="completed"
    if condition=="afterany": return status in {"completed","failed","interrupted","quarantined","paused"}
    if condition=="afternotok": return status in {"failed","interrupted","quarantined"}
    if condition=="aftercorr": return status=="completed"
    raise ValueError(f"Unknown dependency condition {condition!r}.")

def admission_reason(request: ResourceRequest, *, running: list[dict[str,Any]], policy: SchedulerPolicy) -> str|None:
    """Return a readable denial reason without allocating or modifying work."""
    active=[row for row in running if row.get("status") in {"running","claiming","submitted_remote"}]
    if policy.max_running_jobs is not None and len(active)>=policy.max_running_jobs: return "max_running_jobs"
    def used(name): return sum(int(row.get("resources",{}).get(name) or 0) for row in active)
    if policy.max_cores is not None and used("cores")+request.cores>policy.max_cores: return "max_cores"
    if policy.max_memory_mb is not None and used("memory_mb")+(request.memory_mb or 0)>policy.max_memory_mb: return "max_memory_mb"
    if policy.max_gpus is not None and used("gpus")+request.gpus>policy.max_gpus: return "max_gpus"
    return None
