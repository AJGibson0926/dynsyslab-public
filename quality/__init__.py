"""Repeatable quality and architecture checks for the DynSysLab source tree."""
from __future__ import annotations


def __getattr__(name: str):
    if name != "audit_project":
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from .architecture import audit_project
    globals()[name] = audit_project
    return audit_project


__all__ = ["audit_project"]
