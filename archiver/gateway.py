"""Capability-oriented archive access boundary for dynsyslab packages.

The gateway deliberately delegates to the storage kernel without buffering or
copying numerical arrays.  It is therefore a policy boundary, not an I/O
proxy: one gateway call maps to the existing direct archive operation.
"""
from __future__ import annotations

from dataclasses import dataclass
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Literal
import numpy as np
from infrastructure.coordination import read_job_coordination_state

from archives import (ARCHIVE_ROOT, ArchiveSession, campaign_solution_paths,
                      extend_solution_from_archive, find_efficiency_candidate,
                      load_control_target_time_series, load_control_time_series, load_time_series, record_campaign_job,
                      record_campaign_solutions, register_campaign_run,
                      save_hierarchy, update_archive_metadata)

ArchiveCapability = Literal["read", "write", "manage"]

class ArchiveAccessError(PermissionError):
    """Raised when a package crosses its declared archive boundary."""

@dataclass(frozen=True)
class ArchiveGateway:
    """The only public archive-I/O facade for application packages.

    ``read`` opens only existing solution data. ``write`` can create/extend
    solver outputs and add provenance. ``manage`` additionally permits archive
    management operations.  The class holds no solution arrays and adds no
    serialization, so its steady-state overhead is a capability check.
    """
    capability: ArchiveCapability
    root: Path = ARCHIVE_ROOT

    def _require(self, *allowed: ArchiveCapability) -> None:
        if self.capability not in allowed:
            raise ArchiveAccessError(f"Archive capability {self.capability!r} cannot perform this operation.")

    @property
    def archive_root(self) -> Path:
        return Path(self.root).resolve()

    def read_time_series(self, path: str | Path):
        self._require("read", "write", "manage")
        return load_time_series(path)

    def read_control_target_time_series(self, path: str | Path):
        self._require("read", "write", "manage")
        return load_control_target_time_series(path)

    def read_control_time_series(self, path: str | Path):
        self._require("read", "write", "manage")
        return load_control_time_series(path)

    @contextmanager
    def open_read(self, path: str | Path):
        """Open a primary archive read-only without copying its arrays."""
        self._require("read", "write", "manage")
        target=Path(path).expanduser().resolve()
        if not target.is_file(): raise FileNotFoundError(target)
        with np.load(target,allow_pickle=False) as archive:
            yield archive

    def iter_paths(self, section: str, pattern: str="*.npz") -> list[Path]:
        """Discover primary archive files through the policy boundary."""
        self._require("read", "write", "manage")
        return list((self.archive_root/section).rglob(pattern)) if (self.archive_root/section).exists() else []

    def solution_paths(self, *, include_partial: bool=False) -> list[Path]:
        """Return solution records without exposing archive layout to callers."""
        paths=self.iter_paths("solutions", "*.npz")
        return paths if include_partial else [path for path in paths if not path.name.endswith(".partial.npz")]

    def partial_solution_paths(self) -> list[Path]:
        """Return incomplete solution checkpoints through the read boundary."""
        return self.iter_paths("solutions", "*.partial.npz")

    def hierarchy_paths(self) -> list[Path]:
        """Return persisted hierarchy manifests through the read boundary."""
        return self.iter_paths("hierarchies", "*.npz")

    def metadata(self, path: str | Path) -> dict[str, Any]:
        """Read compact archive metadata without materializing numerical arrays."""
        self._require("read", "write", "manage")
        with self.open_read(path) as archive:
            if "metadata_json" not in archive:
                return {}
            import json
            return json.loads(str(archive["metadata_json"]))

    def read_json(self, path: str | Path) -> dict[str, Any]:
        """Read a small manifest through the archive policy boundary."""
        self._require("read", "write", "manage")
        import json
        target=Path(path).expanduser().resolve()
        if not target.is_file():
            raise FileNotFoundError(target)
        return json.loads(target.read_text(encoding="utf-8"))

    def refresh_index(self) -> dict[str, int]:
        """Reconcile the selected archive through the management boundary."""
        self._require("manage")
        from archives import build_index
        return build_index(reconcile=True, archive_root=self.archive_root)

    def new_session(self, request: Any) -> ArchiveSession:
        self._require("write", "manage")
        return ArchiveSession(request)

    def save_hierarchy(self, target: Any) -> Path:
        self._require("write", "manage")
        return save_hierarchy(target)

    def update_metadata(self, path: str | Path, patch: dict[str, Any]) -> Path:
        self._require("write", "manage")
        return update_archive_metadata(path, patch)

    def extend_solution(self, archive_path: str | Path, continuation_archive: str | Path) -> Path:
        self._require("write", "manage")
        return extend_solution_from_archive(archive_path, continuation_archive)

    def efficiency_candidate(self, spec: dict[str, Any]) -> dict[str, Any] | None:
        self._require("read", "write", "manage")
        return find_efficiency_candidate(spec)

    def campaign_paths(self, campaign_id: str) -> list[Path]:
        self._require("read", "write", "manage")
        return campaign_solution_paths(campaign_id)

    def coordination_state(self) -> dict[str, Any]:
        """Read JobManager coordination metadata without exposing its store."""
        self._require("read", "write", "manage")
        return read_job_coordination_state(self.archive_root/"job_manager"/"state.sqlite")

    def register_campaign(self, manifest_path: str | Path, *, status: str="submitted") -> dict[str, Any]:
        self._require("write", "manage")
        return register_campaign_run(manifest_path, status=status)

    def record_campaign_job(self, campaign_id: str, job_id: str, **kwargs: Any) -> None:
        self._require("write", "manage")
        record_campaign_job(campaign_id, job_id, **kwargs)

    def record_campaign_solutions(self, campaign_id: str, job_id: str, paths: list[str | Path]) -> None:
        self._require("write", "manage")
        record_campaign_solutions(campaign_id, job_id, paths)


VIEWER_ARCHIVE = ArchiveGateway("read")
SOLVER_ARCHIVE = ArchiveGateway("write")
MANAGER_ARCHIVE = ArchiveGateway("manage")
POSTPROCESSING_ARCHIVE = ArchiveGateway("write")
