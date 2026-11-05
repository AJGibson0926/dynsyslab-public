"""Resource assessment and immutable-plan parallel transition helpers."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Any
import json, os, shutil
import numpy as np

def _stamp() -> str: return datetime.now(timezone.utc).isoformat()

@dataclass(frozen=True)
class ParallelResourceReport:
    """Low-impact resource assessment for independent campaign work."""
    created_utc: str
    cpu_logical_cores: int
    archive_free_bytes: int
    benchmark_seconds: float | None
    benchmark_kind: str
    active_serial_jobs: int
    recommended_local_workers: int
    candidate_backends: tuple[str,...]
    estimated_independent_speedup: float
    notes: tuple[str,...]
    report_path: str | None=None
    def payload(self) -> dict[str,Any]: return asdict(self)

def evaluate_resources(manager: Any, job_ids: list[str], *, candidate_backends: tuple[str,...], benchmark: bool=True) -> ParallelResourceReport:
    """Low-impact local assessment; never starts a numerical solve."""
    cores=os.cpu_count() or 1; active=sum(manager.status(item).get("status")=="running" for item in job_ids)
    # Avoid a competing heavy integration while a campaign is running.  This
    # vector kernel measures basic array throughput in well under a second.
    timing=None
    if benchmark:
        values=np.linspace(0.,1.,250_000,dtype=float); start=perf_counter()
        for _ in range(4): values=np.sin(values)+np.cos(values)
        timing=perf_counter()-start
    workers=max(1,min(4,cores-1))
    if active: workers=max(1,min(workers, max(1,cores//2)))
    free=shutil.disk_usage(manager.root).free
    efficiency=.75 if workers>1 else 1.; speedup=workers*efficiency
    notes=("Synthetic vector benchmark only; it does not interrupt or duplicate active solves.",
           "Independent parameter trajectories are expected to scale near linearly until CPU, memory, or numerical-library threading becomes limiting.",
           "Remote adapters require explicit staging, scheduler commands, and credentials before they can be selected.")
    return ParallelResourceReport(_stamp(),cores,free,timing,"low-impact NumPy vector throughput" if benchmark else "not run",active,workers,candidate_backends,speedup,notes)

def write_resource_report(manager: Any, report: ParallelResourceReport, *, label: str="parallel-resource-evaluation") -> Path:
    """Persist one resource assessment as machine and human-readable reports."""
    directory=manager.root/"reports"; directory.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"); path=directory/f"{label}__{stamp}.json"
    path.write_text(json.dumps(report.payload(),indent=2,sort_keys=True),encoding="utf-8")
    markdown=path.with_suffix(".md")
    markdown.write_text(f"# Parallel resource evaluation\n\n- Logical CPU cores: {report.cpu_logical_cores}\n- Active serial Jobs: {report.active_serial_jobs}\n- Recommended local worker processes after controlled transition: {report.recommended_local_workers}\n- Estimated independent-work speedup: {report.estimated_independent_speedup:.2f}x\n- Free filesystem space: {report.archive_free_bytes/2**30:.2f} GiB\n- Benchmark: {report.benchmark_kind}; {report.benchmark_seconds!s} s\n- Candidate backends: {', '.join(report.candidate_backends)}\n\n"+"\n".join(f"- {note}" for note in report.notes)+"\n",encoding="utf-8")
    return markdown
