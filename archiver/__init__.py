"""Archive-management façade with lazy capability and management exports."""
from __future__ import annotations

from importlib import import_module


_LAZY_EXPORTS={
    "ArchiveAuditReport":(".core","ArchiveAuditReport"),
    "ArchiveScope":(".core","ArchiveScope"),
    "Archiver":(".core","Archiver"),
    "DeduplicationResult":(".core","DeduplicationResult"),
    "StorageLocation":(".storage","StorageLocation"),
    "StorageRegistry":(".storage","StorageRegistry"),
    "ArchiveAccessError":(".gateway","ArchiveAccessError"),
    "ArchiveGateway":(".gateway","ArchiveGateway"),
    "MANAGER_ARCHIVE":(".gateway","MANAGER_ARCHIVE"),
    "POSTPROCESSING_ARCHIVE":(".gateway","POSTPROCESSING_ARCHIVE"),
    "SOLVER_ARCHIVE":(".gateway","SOLVER_ARCHIVE"),
    "VIEWER_ARCHIVE":(".gateway","VIEWER_ARCHIVE"),
}


def __getattr__(name: str):
    target=_LAZY_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value=getattr(import_module(target[0],__name__),target[1])
    globals()[name]=value
    return value


def __dir__():
    return sorted(set(globals())|set(_LAZY_EXPORTS))


__all__=["ArchiveScope","ArchiveAuditReport","DeduplicationResult","Archiver","StorageLocation","StorageRegistry","ArchiveAccessError","ArchiveGateway","VIEWER_ARCHIVE","SOLVER_ARCHIVE","MANAGER_ARCHIVE","POSTPROCESSING_ARCHIVE"]
