"""Small lazy runtime primitives shared by JobManager services."""
from __future__ import annotations

from datetime import datetime, timezone


def _stamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe(value: str) -> str:
    return "".join(item if item.isalnum() or item in "-_" else "-"
                   for item in value)


class _ArchiveGatewayProxy:
    """Defer solver/SciPy-only archive imports for plan-only management work."""

    def __getattr__(self, name: str):
        from archiver.gateway import MANAGER_ARCHIVE
        return getattr(MANAGER_ARCHIVE, name)


MANAGER_ARCHIVE = _ArchiveGatewayProxy()


class _PauseAtBoundary(RuntimeError):
    """Internal control signal: current trajectory is safely archived."""
