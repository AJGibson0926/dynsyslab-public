"""Remote scheduler adapters. They build/submit jobs only when configured."""
from __future__ import annotations
from dataclasses import dataclass, field
import subprocess
from typing import Any, Callable
from uuid import uuid4
from .base import BackendCapabilities, BackendConfigurationError, BackendResult, ExecutionBackend, RemoteHandle

@dataclass
class _CommandBackend(ExecutionBackend):
    submit_command: tuple[str,...] | None=None
    capabilities=BackendCapabilities("command",supports_parallel=True,supports_remote=True)
    def execute(self, job_id: str, *, execute_local: Callable[[str],dict[str,Any]], context: dict[str,Any]) -> BackendResult:
        return BackendResult(self.capabilities.name,job_id,"submitted",handle=self.submit(job_id,context=context))
    def submit(self, job_id: str, *, context: dict[str,Any]) -> RemoteHandle:
        if not self.submit_command: raise BackendConfigurationError(f"Configure {self.capabilities.name}.submit_command before remote submission.")
        command=[part.format(job_id=job_id,manager_root=context["manager_root"],project_root=context["project_root"]) for part in self.submit_command]
        completed=subprocess.run(command,cwd=str(context["project_root"]),capture_output=True,text=True)
        if completed.returncode: raise BackendConfigurationError(f"{self.capabilities.name} submission failed: {completed.stderr.strip()}")
        external=completed.stdout.strip().splitlines()[-1] if completed.stdout.strip() else uuid4().hex
        return RemoteHandle(self.capabilities.name,external,metadata={"command":command,"job_id":job_id})
    def status(self, handle: RemoteHandle) -> str: return "submitted"

@dataclass
class SSHBackend(_CommandBackend):
    host: str | None=None
    remote_project_root: str | None=None
    capabilities=BackendCapabilities("ssh",supports_parallel=True,supports_remote=True)

@dataclass
class SlurmBackend(_CommandBackend):
    capabilities=BackendCapabilities("slurm",supports_parallel=True,supports_remote=True)

@dataclass
class PBSBackend(_CommandBackend):
    capabilities=BackendCapabilities("pbs",supports_parallel=True,supports_remote=True)

@dataclass
class CloudBatchBackend(_CommandBackend):
    """Generic AWS/Azure/GCP batch adapter configured with a provider command."""
    provider: str | None=None
    capabilities=BackendCapabilities("cloud_batch",supports_parallel=True,supports_remote=True)
    def __post_init__(self):
        if self.provider and not self.submit_command:
            raise BackendConfigurationError("CloudBatchBackend requires an explicit provider-specific submit_command; credentials are never stored by dynsyslab.")
