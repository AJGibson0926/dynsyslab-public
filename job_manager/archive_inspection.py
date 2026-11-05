"""Read-only archive inspection used by JobManager planning and recovery.

The inspector is intentionally separated from dispatch.  Every primary
archive read goes through ``MANAGER_ARCHIVE``; it keeps the scheduler focused
on manifests and state transitions rather than archive file layout.
"""
from __future__ import annotations

from pathlib import Path
import re
from typing import Any, Callable

from .manifests import JobManifest


class ArchiveInspector:
    """Read-only archive reuse and hierarchy-completeness queries."""
    def __init__(self, gateway: Any, mismatches: Callable[[dict[str,Any],dict[str,Any]],list[str]]):
        self.gateway=gateway
        self.mismatches=mismatches

    def completed_system_labels(self, job: JobManifest) -> set[str]:
        wanted={str(spec.get("label",job.label)):spec for spec in job.systems}; found=set()
        for label,spec in wanted.items():
            candidate=self.gateway.efficiency_candidate(spec)
            if candidate and candidate["kind"]=="reuse": found.add(label)
        if len(found)==len(wanted): return found
        for path in self.gateway.solution_paths():
            try:
                metadata=self.gateway.metadata(path); label=str(metadata.get("label", ""))
                canonical_label=re.sub(r"-c\d+-t\d+$", "", label)
                if (canonical_label in wanted and metadata.get("status", "complete") in {"complete","completed"}
                        and not self.mismatches(metadata,wanted[canonical_label])):
                    found.add(canonical_label)
            except Exception:
                continue
        return found

    def completed_plan_archive(self, fingerprint: str, state: dict[str,Any]) -> str | None:
        for row in state.get("jobs",{}).values():
            if row.get("status")=="completed" and row.get("hierarchy_archive") and row.get("plan_fingerprint")==fingerprint:
                return str(row["hierarchy_archive"])
        for path in self.gateway.hierarchy_paths():
            try:
                if self.gateway.metadata(path).get("job",{}).get("plan_fingerprint")==fingerprint:
                    return str(path)
            except Exception:
                continue
        return None

    def hierarchy_matches(self, archive_path: str | Path, target: Any) -> bool:
        systems=target.systems if hasattr(target,"systems") else [target]
        try:
            saved=self.gateway.metadata(archive_path).get("systems",[])
            if len(saved)!=len(systems): return False
            for expected,actual in zip(systems,saved):
                if actual.get("label")!=expected.label or len(actual.get("clusters",[]))!=len(expected.clusters): return False
                for expected_cluster,actual_cluster in zip(expected.clusters,actual["clusters"]):
                    children=actual_cluster.get("trajectories",[])
                    if actual_cluster.get("label")!=expected_cluster.label or len(children)!=len(expected_cluster.trajectories): return False
                    if any(not child.get("solution_archive") or not Path(child["solution_archive"]).is_file() for child in children): return False
            return True
        except Exception:
            return False
