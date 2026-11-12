"""Portable execution-backend protocol; no solver implementation belongs here."""
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable
from uuid import uuid4

def _stamp() -> str: return datetime.now(timezone.utc).isoformat()

class BackendConfigurationError(RuntimeError):
    """Raised when a remote backend has not been explicitly configured."""

@dataclass(frozen=True)
class BackendCapabilities:
    name: str
    supports_parallel: bool
    supports_remote: bool=False
    supports_dependency_graph: bool=True
    supports_checkpoint_resume: bool=True
    supports_resource_requests: bool=True

@dataclass(frozen=True)
class RemoteHandle:
    backend: str
    external_id: str
    submitted_utc: str=field(default_factory=_stamp)
    metadata: dict[str,Any]=field(default_factory=dict)
    def payload(self) -> dict[str,Any]: return asdict(self)

@dataclass(frozen=True)
class BackendResult:
    backend: str
    job_id: str
    status: str
    handle: RemoteHandle | None=None
    details: dict[str,Any]=field(default_factory=dict)

class ExecutionBackend(ABC):
    """JobManager-owned adapter. Solver receives no backend knowledge."""
    capabilities: BackendCapabilities
    @abstractmethod
    def execute(self, job_id: str, *, execute_local: Callable[[str],dict[str,Any]], context: dict[str,Any]) -> BackendResult: ...
    def submit(self, job_id: str, *, context: dict[str,Any]) -> RemoteHandle:
        raise BackendConfigurationError(f"{self.capabilities.name} does not expose remote submission.")
    def status(self, handle: RemoteHandle) -> str: return "unknown"
    def cancel(self, handle: RemoteHandle) -> None: raise BackendConfigurationError(f"{self.capabilities.name} does not expose cancellation.")
