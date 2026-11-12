"""Lightweight, serializable solve-plan manifests owned by :mod:`job_manager`."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal
from uuid import uuid4
import hashlib
import json
import os
from infrastructure.persistence import require_format, require_keys

ExecutionMode = Literal["series", "parallel"]
AdaptivePolicy = Literal["fixed", "adaptive_local", "adaptive_portable"]

@dataclass(frozen=True)
class JobManifest:
    """A solve plan only; numerical outputs always remain in ``archives``.

    ``systems`` preserves the solver's Batch > System > Cluster > Trajectory
    structure using serializable system construction and initial-condition data.
    It never stores trajectories produced by a solve, ROMs, or analyses.
    """
    label: str
    systems: tuple[dict[str, Any], ...]
    job_id: str = field(default_factory=lambda: uuid4().hex)
    dependencies: tuple[str, ...] = ()
    dependency_conditions: dict[str, str] = field(default_factory=dict)
    input_bindings: dict[str, str] = field(default_factory=dict)
    postprocessing: tuple[dict[str, Any], ...] = ()
    recovery: dict[str, Any] = field(default_factory=dict)
    schedule: dict[str, Any] = field(default_factory=lambda: {"start_at": None, "after_campaign_ids": [], "result_location": "primary"})
    resources: dict[str, Any] = field(default_factory=lambda: {"cores": 1, "memory_mb": None, "wall_time_s": None})
    staging: dict[str, Any] = field(default_factory=dict)
    deadline: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    execution: dict[str, Any] = field(default_factory=lambda: {"backend": "serial", "adaptive_policy": "fixed", "retry_limit": 0})
    tags: dict[str, Any] = field(default_factory=dict)
    notes: str = ""
    format: str = "dynsyslab.job/v1"

    def __post_init__(self):
        if not self.label or not self.systems:
            raise ValueError("A job manifest needs a label and at least one planned system.")
        if any(not isinstance(item, dict) for item in self.systems):
            raise TypeError("Job systems must be serializable dictionaries.")
        from .scheduling import Dependency, ResourceRequest, DeadlinePolicy, StagingPlan
        ResourceRequest(**{key:value for key,value in self.resources.items() if key in ResourceRequest.__dataclass_fields__})
        DeadlinePolicy(**{key:value for key,value in self.deadline.items() if key in DeadlinePolicy.__dataclass_fields__})
        StagingPlan(**{key:tuple(value) if key=="inputs" else value for key,value in self.staging.items() if key in StagingPlan.__dataclass_fields__})
        for item in self.dependencies: Dependency(item,self.dependency_conditions.get(item,"afterok"))

    def payload(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def plan_fingerprint(self) -> str:
        """Stable identity for repeated-work detection, independent of labels/tags."""
        def clean(value: Any):
            if isinstance(value, dict): return {key: clean(item) for key, item in value.items() if key not in {"label", "job_id", "tags"}}
            if isinstance(value, (list, tuple)): return [clean(item) for item in value]
            return value
        plan={"systems":self.systems,"dependencies":self.dependencies,"dependency_conditions":self.dependency_conditions,"input_bindings":self.input_bindings,"postprocessing":self.postprocessing,"recovery":self.recovery,"schedule":self.schedule,"resources":self.resources,"staging":self.staging,"deadline":self.deadline,"execution":self.execution}
        return hashlib.sha256(json.dumps(clean(plan),sort_keys=True,separators=(",",":"),default=str).encode("utf-8")).hexdigest()

    def write(self, path: str | Path) -> Path:
        target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
        # A manifest is the durable contract for a job.  Never expose a
        # half-written contract to a concurrent scheduler or status tool.
        temporary=target.with_suffix(target.suffix+f".{uuid4().hex}.partial")
        with temporary.open("w",encoding="utf-8",newline="\n") as handle:
            handle.write(json.dumps(self.payload(), indent=2, sort_keys=True))
            handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary,target)
        return target

    @classmethod
    def read(cls, path: str | Path) -> "JobManifest":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        require_format(data,"dynsyslab.job/"); require_keys(data,"job_id","label","systems")
        for name in ("systems", "dependencies", "postprocessing"):
            if name in data and isinstance(data[name], list): data[name] = tuple(data[name])
        return cls(**data)

@dataclass(frozen=True)
class JobCollection:
    """A named, dependency-aware collection of lightweight job manifests."""
    label: str
    jobs: tuple[JobManifest, ...]
    collection_id: str = field(default_factory=lambda: uuid4().hex)
    tags: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        ids = [item.job_id for item in self.jobs]
        if len(ids) != len(set(ids)): raise ValueError("A job collection cannot contain duplicate job IDs.")
        unknown = {dependency for item in self.jobs for dependency in item.dependencies if dependency not in ids}
        if unknown: raise ValueError(f"Collection dependencies must point to a collection job: {sorted(unknown)}.")
