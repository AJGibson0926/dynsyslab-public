"""Durable, archive-oriented execution of lightweight dynsyslab job manifests.

The public :class:`JobManager` remains a small facade. Focused internal
services own state/events, transition recovery, planning, archive reuse, and
execution without changing the manifest or state-store contracts.
"""
from __future__ import annotations

from pathlib import Path
from threading import Lock
import json

from ._runtime import (MANAGER_ARCHIVE, _ArchiveGatewayProxy,
                       _PauseAtBoundary, _safe, _stamp)
from .execution_services import _JobExecutionMixin
from .planning_services import _JobPlanningMixin
from .reuse_services import _JobReuseMixin
from .state_services import _JobStateMixin
from .transition_services import _JobTransitionMixin


class JobManager(_JobStateMixin, _JobTransitionMixin, _JobPlanningMixin,
                 _JobReuseMixin, _JobExecutionMixin):
    """Persist and schedule Jobs through focused, behavior-preserving services."""

    def __init__(self, root: str | Path | None = None):
        self.root = Path(root or MANAGER_ARCHIVE.archive_root / "job_manager").resolve()
        self.manifest_dir = self.root / "manifests"; self.original_manifest_dir=self.root / "original_manifests"; self.collection_dir = self.root / "collections"
        self.side_tray_dir=self.root/"side_tray"; self.side_tray_snapshot_dir=self.side_tray_dir/"snapshots"; self.side_tray_path=self.side_tray_dir/"reserved.json"; self.array_dir=self.root/"arrays"; self.policy_path=self.root/"scheduler_policy.json"
        self.log_dir = self.root / "logs"; self.restructure_dir=self.root / "restructures"; self.review_dir=self.root / "reviews"; self.state_path = self.root / "state.sqlite"; self.state_lock_path=self.root / ".state.lock"
        for directory in (self.manifest_dir, self.original_manifest_dir, self.collection_dir, self.side_tray_dir, self.side_tray_snapshot_dir, self.log_dir, self.restructure_dir, self.review_dir, self.array_dir): directory.mkdir(parents=True, exist_ok=True)
        self._lock = Lock()
        from .state_store import JobStateStore
        from .archive_inspection import ArchiveInspector
        from .campaign_reporting import CampaignReporter
        self.state_store=JobStateStore(self.state_path)
        self._archive_inspector=ArchiveInspector(MANAGER_ARCHIVE,self._specification_mismatches)
        self._campaign_reporter=CampaignReporter(self,MANAGER_ARCHIVE,_stamp)
        legacy=self.root/"state.json"
        if legacy.is_file() and not self.state_store.read().get("jobs"):
            self.state_store.write(json.loads(legacy.read_text(encoding="utf-8")))
        if not self.side_tray_path.exists(): self._write_side_tray({"format":"dynsyslab.job-manager-side-tray/v1","entries":[]})
        if not self.policy_path.exists(): self.policy_path.write_text(json.dumps({}),encoding="utf-8")
