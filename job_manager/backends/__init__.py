"""Execution backends used exclusively by :mod:`job_manager` scheduling."""
from .base import BackendCapabilities, BackendConfigurationError, BackendResult, ExecutionBackend, RemoteHandle
from .local import LocalProcessBackend, LocalSerialBackend, LocalWorkerBackend
from .remote import CloudBatchBackend, PBSBackend, SSHBackend, SlurmBackend
from .scripts import pbs_script, slurm_script

BACKEND_TYPES = {
    "serial": LocalSerialBackend,
    "local_process": LocalProcessBackend,
    "local_worker": LocalWorkerBackend,
    "ssh": SSHBackend,
    "slurm": SlurmBackend,
    "pbs": PBSBackend,
    "cloud_batch": CloudBatchBackend,
}

def create_backend(kind: str, **configuration):
    try: return BACKEND_TYPES[str(kind).lower()](**configuration)
    except KeyError as error: raise ValueError(f"Unknown JobManager backend {kind!r}; choose from {sorted(BACKEND_TYPES)}.") from error

__all__=["BackendCapabilities","BackendConfigurationError","BackendResult","ExecutionBackend","RemoteHandle","LocalSerialBackend","LocalProcessBackend","LocalWorkerBackend","SSHBackend","SlurmBackend","PBSBackend","CloudBatchBackend","slurm_script","pbs_script","create_backend","BACKEND_TYPES"]
