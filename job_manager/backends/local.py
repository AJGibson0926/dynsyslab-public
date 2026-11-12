"""Local execution adapters: deterministic serial and isolated processes."""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
import os, subprocess, sys
from typing import Any, Callable, Iterable
from .base import BackendCapabilities, BackendResult, ExecutionBackend, RemoteHandle
from ..resource_limits import numerical_thread_limit, numerical_thread_policy, worker_thread_environment


def _allocated_cores(job_id: str, context: dict[str, Any]) -> int:
    resources=context.get("resources_by_job",{}).get(job_id,{})
    return int(resources.get("cores",1))

@dataclass
class LocalSerialBackend(ExecutionBackend):
    capabilities=BackendCapabilities("serial",supports_parallel=False)
    def execute(self, job_id: str, *, execute_local: Callable[[str],dict[str,Any]], context: dict[str,Any]) -> BackendResult:
        cores=_allocated_cores(job_id,context)
        with numerical_thread_limit(cores) as policy:
            row=execute_local(job_id)
        return BackendResult("serial",job_id,str(row.get("status","unknown")),details={"in_process":True,"numerical_thread_policy":policy})

@dataclass
class LocalProcessBackend(ExecutionBackend):
    """One isolated Python worker per job; suitable for CPU-bound solves."""
    max_workers: int | None=None
    python_executable: str=sys.executable
    capabilities=BackendCapabilities("local_process",supports_parallel=True)
    def _command(self, job_id: str, context: dict[str,Any]) -> list[str]:
        return [self.python_executable,"-m","job_manager.worker","--manager-root",str(context["manager_root"]),"--job-id",job_id]
    def execute(self, job_id: str, *, execute_local: Callable[[str],dict[str,Any]], context: dict[str,Any]) -> BackendResult:
        # A lone process uses the same worker protocol as a many-job pool.
        cores=_allocated_cores(job_id,context); environment=worker_thread_environment(cores)
        completed=subprocess.run(self._command(job_id,context),cwd=str(context["project_root"]),capture_output=True,text=True,env=environment)
        return BackendResult("local_process",job_id,"submitted" if completed.returncode==0 else "failed",details={"returncode":completed.returncode,"stdout":completed.stdout[-4000:],"stderr":completed.stderr[-4000:],"numerical_thread_policy":numerical_thread_policy(cores)})
    def execute_many(self, job_ids: Iterable[str], *, context: dict[str,Any]) -> list[BackendResult]:
        workers=self.max_workers or max(1,(os.cpu_count() or 2)-1); results=[]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures={pool.submit(self.execute,job,execute_local=lambda _: {},context=context):job for job in job_ids}
            for future in as_completed(futures): results.append(future.result())
        return results

@dataclass
class LocalWorkerBackend(LocalProcessBackend):
    """Persistent-queue compatible local backend; currently process-isolated."""
    capabilities=BackendCapabilities("local_worker",supports_parallel=True,supports_checkpoint_resume=True)
