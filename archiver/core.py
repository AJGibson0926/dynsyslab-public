"""Read-mostly, provenance-aware management of dynsyslab archives."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timedelta, timezone
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable
import json
import os

import numpy as np

from archives import ARCHIVE_ROOT, rewrite_archive_metadata, campaign_record, campaign_solution_paths, list_campaign_runs, register_campaign_run, record_campaign_job, record_campaign_solutions


def _stamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe(value: str) -> str:
    return "".join(char if char.isalnum() or char in "-_" else "-" for char in value)


def _hash(array: np.ndarray | None) -> str | None:
    if array is None:
        return None
    item = np.asarray(array)
    digest = sha256()
    digest.update(str(item.shape).encode())
    digest.update(item.dtype.str.encode())
    digest.update(np.ascontiguousarray(item).tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class ArchiveScope:
    """Restrict an archive-management operation without touching old data."""
    days: float | None = 2.0
    sections: tuple[str, ...] = ("solutions",)
    systems: tuple[str, ...] = ()
    campaign_id: str | None = None
    include_derived: bool = False


@dataclass
class ArchiveAuditReport:
    created_utc: str
    root: str
    scope: dict[str, Any]
    scanned_files: int
    total_bytes: int
    exact_duplicate_groups: list[dict[str, Any]] = field(default_factory=list)
    same_plan_groups: list[dict[str, Any]] = field(default_factory=list)
    unreadable_files: list[dict[str, str]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    report_json: str | None = None
    report_markdown: str | None = None

    @property
    def redundant_file_count(self) -> int:
        return sum(max(0, len(group["files"]) - 1) for group in self.exact_duplicate_groups)

    @property
    def redundant_bytes(self) -> int:
        return sum(int(group.get("redundant_bytes", 0)) for group in self.exact_duplicate_groups)

    def payload(self) -> dict[str, Any]:
        result = asdict(self)
        result["redundant_file_count"] = self.redundant_file_count
        result["redundant_bytes"] = self.redundant_bytes
        return result


@dataclass
class DeduplicationResult:
    created_utc: str
    dry_run: bool
    canonical_files: list[str]
    removed_files: list[str]
    skipped_groups: list[dict[str, Any]]
    reclaimed_bytes: int
    audit_report: str | None = None
    rewritten_reference_files: list[str] = field(default_factory=list)


class Archiver:
    """Dedicated archive-management facade.

    The class deliberately has no solver dependency.  Its optional adapters
    only read campaign/job metadata and publish a source-change event after a
    successful mutation, allowing Campaigner, JobManager, Solver,
    Postprocessing, and Viewer to coexist without a circular dependency.
    """
    def __init__(self, root: str | Path | None = None):
        self.root = Path(root or ARCHIVE_ROOT).resolve()
        self.management_root = self.root / "management"
        self.report_root = self.management_root / "reports"
        self.report_root.mkdir(parents=True, exist_ok=True)

    def integration_status(self) -> dict[str, dict[str, Any]]:
        """Describe lightweight integration points without starting services."""
        return {
            "solver": {"role": "archives completed solutions; Archiver never invokes it", "available": True},
            "postprocessing": {"role": "may be audited as derived data when requested", "available": True},
            "job_manager": {"role": "its manifests/logs provide provenance for duplicate causes", "available": (self.root / "job_manager").exists()},
            "campaigner": {"role": "campaign manifests can scope archive searches", "available": (self.root / "campaigns").exists()},
            "viewer": {"role": "receives passive archive-change events after mutations", "available": True},
        }

    def campaign_history(self) -> list[dict[str, Any]]:
        """List JobManager-run campaigns from the indexed archive catalogue."""
        return list_campaign_runs()

    def campaign_lookup(self, campaign_id: str) -> dict[str, Any] | None:
        """Fast campaign-to-solution lookup for management and Viewer adapters."""
        record=campaign_record(campaign_id)
        if record is None: return None
        record["existing_solution_paths"]=[str(path) for path in campaign_solution_paths(campaign_id)]
        return record

    def rebuild_campaign_catalog(self) -> dict[str, int]:
        """Backfill archive history from already-submitted JobManager records.

        This touches only the SQLite discovery catalogue.  It neither changes
        campaign manifests nor numerical archives, so it is safe while a
        campaign is live.
        """
        from infrastructure.coordination import read_job_coordination_state
        database=self.root/"job_manager"/"state.sqlite"
        if not database.is_file(): return {"campaigns":0,"jobs":0,"solution_paths":0}
        state=read_job_coordination_state(database).get("jobs",{}); campaigns=jobs=paths=0; seen=set()
        for job_id,row in state.items():
            tags=row.get("tags",{}); campaign_id=tags.get("campaign_id"); manifest=tags.get("campaign_manifest")
            if not campaign_id or not manifest or not Path(manifest).is_file(): continue
            if campaign_id not in seen: register_campaign_run(manifest,status="historical"); seen.add(campaign_id); campaigns+=1
            record_campaign_job(str(campaign_id),str(job_id),status=str(row.get("status","unknown")),job_label=row.get("label"),manifest_path=row.get("manifest")); jobs+=1
            items=[item for key in ("solution_archives","generated_solution_archives") for item in row.get(key,[]) or []]
            usable=[item for item in items if Path(item).is_file()]
            if usable: record_campaign_solutions(str(campaign_id),str(job_id),usable); paths+=len(usable)
        # Older JobManager states did not always retain every completed child
        # path.  Solution metadata is the independent provenance source, so
        # backfill those paths once into the same fast catalogue.
        for archive_path in (self.root/"solutions").rglob("*.npz"):
            if archive_path.name.endswith(".partial.npz"): continue
            try:
                with np.load(archive_path,allow_pickle=False) as data: metadata=json.loads(str(data["metadata_json"]))
                job=metadata.get("job",{}); tags=job.get("tags",{}); campaign_id=tags.get("campaign_id")
                if not campaign_id: continue
                record_campaign_solutions(str(campaign_id),str(job.get("job_id","archived-provenance")),[archive_path]); paths+=1
            except (OSError,ValueError,KeyError,json.JSONDecodeError):
                continue
        return {"campaigns":campaigns,"jobs":jobs,"solution_paths":paths}

    def reconcile_index(self) -> dict[str,int]:
        """Reconcile the shared archive summary without changing scientific data."""
        from archives import build_index
        return build_index(reconcile=True,archive_root=self.root)

    def audit(self, scope: ArchiveScope | None = None, *, write_report: bool = True) -> ArchiveAuditReport:
        """Audit exact numerical duplicates in a bounded archive scope.

        By default this scans completed solution files modified within two
        days.  An exact duplicate has identical final time and state arrays,
        not merely matching forcing parameters.
        """
        scope = scope or ArchiveScope()
        records: list[dict[str, Any]] = []
        unreadable: list[dict[str, str]] = []
        for path in self._paths(scope):
            try:
                record = self._record(path)
            except Exception as error:  # An audit must continue past bad/partial files.
                unreadable.append({"path": str(path), "error": f"{type(error).__name__}: {error}"})
            else:
                records.append(record)
        exact: dict[str, list[dict[str, Any]]] = {}
        plans: dict[str, list[dict[str, Any]]] = {}
        for record in records:
            exact.setdefault(record["exact_key"], []).append(record)
            plans.setdefault(record["plan_key"], []).append(record)
        exact_groups = [self._group("exact_solution", value) for value in exact.values() if len(value) > 1]
        plan_groups = [self._group("same_plan", value) for value in plans.values() if len(value) > 1]
        report = ArchiveAuditReport(
            created_utc=_stamp(), root=str(self.root), scope=asdict(scope),
            scanned_files=len(records), total_bytes=sum(item["size_bytes"] for item in records),
            exact_duplicate_groups=exact_groups, same_plan_groups=plan_groups,
            unreadable_files=unreadable,
            notes=[
                "Exact duplicates require byte-identical final time and state arrays.",
                "Derived results are excluded unless scope.include_derived is enabled.",
                "Deduplication is dry-run by default and never alters the legacy RPE campaign outside this archive root.",
            ],
        )
        if write_report:
            json_path, markdown_path = self._write_report(report)
            report.report_json, report.report_markdown = str(json_path), str(markdown_path)
            # Persist links after paths are known.
            json_path.write_text(json.dumps(report.payload(), indent=2, sort_keys=True), encoding="utf-8")
        return report

    def deduplicate(self, scope: ArchiveScope | None = None, *, dry_run: bool = True,
                    relink_references: bool = False) -> DeduplicationResult:
        """Remove exact duplicate *solution* files without any Codex service.

        A canonical file is selected deterministically (newest mtime, then
        pathname).  Files referenced by a hierarchy, derived result, campaign,
        or JobManager record are skipped rather than leaving dangling
        provenance.  Set ``relink_references=True`` to atomically repoint
        provenance to the canonical file before deletion.  Pass
        ``dry_run=False`` to execute the independently runnable cleanup.
        """
        report = self.audit(scope, write_report=True)
        references = self._external_references()
        canonical: list[str] = []; removed: list[str] = []; skipped: list[dict[str, Any]] = []; reclaimed = 0
        replacements: dict[str, str] = {}
        for group in report.exact_duplicate_groups:
            candidates = sorted((Path(item) for item in group["files"]), key=lambda item: (item.stat().st_mtime, str(item)), reverse=True)
            keep = candidates[0]; canonical.append(str(keep)); deletable = []
            for candidate in candidates[1:]:
                if str(candidate) in references and not relink_references:
                    skipped.append({"canonical": str(keep), "candidate": str(candidate), "reason": "referenced by archive provenance"})
                else:
                    deletable.append(candidate)
                    if str(candidate) in references:
                        replacements[str(candidate)] = str(keep)
        rewritten: list[str] = []
        if replacements and not dry_run:
            rewritten = self._relink_references(replacements)
        for group in report.exact_duplicate_groups:
            candidates = sorted((Path(item) for item in group["files"]), key=lambda item: (item.stat().st_mtime, str(item)), reverse=True)
            keep = candidates[0]
            for candidate in candidates[1:]:
                if str(candidate) in references and not relink_references:
                    continue
                if not dry_run and candidate.is_file():
                    size = candidate.stat().st_size
                    candidate.unlink()
                    removed.append(str(candidate)); reclaimed += size
        result = DeduplicationResult(_stamp(), dry_run, canonical, removed, skipped, reclaimed, report.report_json, rewritten)
        if not dry_run and removed:
            self.reconcile_index()
            self._notify_archive_change(removed)
        return result

    def _relink_references(self, replacements: dict[str, str]) -> list[str]:
        """Replace source paths in NPZ metadata and JSON manifests atomically."""
        changed: list[str] = []
        for path in self.root.rglob("*.npz"):
            if path.name.endswith(".partial.npz"):
                continue
            try:
                with np.load(path, allow_pickle=False) as archive:
                    if "metadata_json" not in archive.files:
                        continue
                    metadata = json.loads(str(archive["metadata_json"]))
                updated = self._replace_paths(metadata, replacements)
                if updated != metadata:
                    rewrite_archive_metadata(path, lambda _: updated)
                    changed.append(str(path))
            except (OSError, ValueError, KeyError, json.JSONDecodeError):
                continue
        for path in self.root.rglob("*.json"):
            try:
                payload = json.loads(path.read_text(encoding="utf-8")); updated = self._replace_paths(payload, replacements)
                if updated != payload:
                    temporary = path.with_suffix(path.suffix + ".partial")
                    temporary.write_text(json.dumps(updated, indent=2, sort_keys=True), encoding="utf-8")
                    os.replace(temporary, path); changed.append(str(path))
            except (OSError, json.JSONDecodeError):
                continue
        return changed

    @staticmethod
    def _replace_paths(value: Any, replacements: dict[str, str]) -> Any:
        if isinstance(value, dict):
            return {key: Archiver._replace_paths(item, replacements) for key, item in value.items()}
        if isinstance(value, list):
            return [Archiver._replace_paths(item, replacements) for item in value]
        if isinstance(value, str):
            for old, new in replacements.items():
                value = value.replace(old, new)
        return value

    def _paths(self, scope: ArchiveScope) -> Iterable[Path]:
        cutoff = None if scope.days is None else datetime.now(timezone.utc) - timedelta(days=float(scope.days))
        sections = scope.sections or ("solutions",)
        campaign_paths = self._campaign_solution_paths(scope.campaign_id) if scope.campaign_id else None
        for section in sections:
            base = self.root / section
            if not base.exists():
                continue
            for path in base.rglob("*.npz"):
                # An active solve writes atomically through a transient
                # ``.partial.npz`` path.  Never audit it as completed data.
                if path.name.endswith(".partial.npz"):
                    continue
                try:
                    modified = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
                except FileNotFoundError:
                    # A completed archive can be renamed between discovery
                    # and inspection while a campaign is actively writing.
                    continue
                if cutoff is not None and modified < cutoff:
                    continue
                if scope.systems and not any(system in path.parts for system in scope.systems):
                    continue
                if campaign_paths is not None and path.resolve() not in campaign_paths:
                    continue
                if not scope.include_derived and "derived" in path.parts:
                    continue
                yield path

    def _campaign_solution_paths(self, campaign_id: str) -> set[Path]:
        """Resolve campaign-owned solution paths through JobManager state.

        Campaign manifests are intentionally lightweight plans; generated
        solution locations are recorded by JobManager, which is the canonical
        campaign-to-data bridge.
        """
        indexed=set(campaign_solution_paths(campaign_id))
        if indexed:
            return indexed
        database = self.root / "job_manager" / "state.sqlite"
        if not database.is_file():
            return set()
        try:
            from infrastructure.coordination import read_job_coordination_state
            jobs = read_job_coordination_state(database).get("jobs", {})
        except (OSError, ValueError):
            return set()
        paths: set[Path] = set()
        for row in jobs.values():
            if row.get("tags", {}).get("campaign_id") != campaign_id:
                continue
            for key in ("solution_archives", "generated_solution_archives"):
                for item in row.get(key, []) or []:
                    candidate = Path(item)
                    if candidate.is_file():
                        paths.add(candidate.resolve())
        # A campaign can be running under an earlier JobManager version that
        # has not yet recorded per-trajectory paths.  Its standalone job plans
        # still contain deterministic system labels, so use those as a
        # read-only fallback rather than declaring the campaign empty.
        manifest = None
        for candidate in (self.root / "campaigns").glob("*.json"):
            try:
                payload = json.loads(candidate.read_text(encoding="utf-8"))
                if payload.get("campaign_id") == campaign_id:
                    manifest = payload; break
            except (OSError, json.JSONDecodeError):
                continue
        if manifest is None:
            return paths
        for plan_name in manifest.get("metadata", {}).get("job_plan_paths", []):
            try:
                plan = json.loads(Path(plan_name).read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            for system in plan.get("systems", []):
                label = _safe(str(system.get("label", "")))
                if not label:
                    continue
                for solution in (self.root / "solutions").rglob(f"{label}*.npz"):
                    if solution.is_file() and not solution.name.endswith(".partial.npz"):
                        paths.add(solution.resolve())
        return paths

    def _record(self, path: Path) -> dict[str, Any]:
        with np.load(path, allow_pickle=False) as archive:
            metadata = json.loads(str(archive["metadata_json"])) if "metadata_json" in archive else {}
            values = np.asarray(archive["values"]) if "values" in archive else None
            time = np.asarray(archive["time"]) if "time" in archive else None
        if values is None:
            raise ValueError("not a solution-like archive with a values array")
        driver_sources = metadata.get("forcing", {}).get("sources") or []
        driver = driver_sources[0].get("parameters", {}) if driver_sources else {}
        system_parameters = metadata.get("system_parameters", {})
        plan = {
            "system": metadata.get("system"), "preset": system_parameters.get("preset", {}).get("key"),
            "initial_state": metadata.get("initial_state_override"), "integration_form": metadata.get("integration_form"),
            "time_span": metadata.get("time_span"), "solver_options": metadata.get("options"),
            "driver_amplitude": driver.get("amplitude"), "driver_frequency_hz": driver.get("frequency_hz"),
            "driver_phase_rad": driver.get("phase_rad"), "state_shape": list(values.shape),
            "time_shape": None if time is None else list(time.shape),
        }
        plan_key = json.dumps(plan, sort_keys=True, default=str)
        exact_key = json.dumps({**plan, "time_hash": _hash(time), "values_hash": _hash(values)}, sort_keys=True)
        return {"path": str(path), "size_bytes": path.stat().st_size, "plan_key": plan_key, "exact_key": exact_key, "plan": plan}

    def _group(self, kind: str, members: list[dict[str, Any]]) -> dict[str, Any]:
        ordered = sorted(members, key=lambda item: item["path"])
        # All but one file in an exact group are physically redundant.
        redundant = sum(item["size_bytes"] for item in ordered[:-1]) if kind == "exact_solution" else 0
        return {"kind": kind, "count": len(ordered), "files": [item["path"] for item in ordered],
                "plan": ordered[0]["plan"], "redundant_bytes": redundant,
                "likely_cause": self._likely_cause(ordered)}

    @staticmethod
    def _likely_cause(members: list[dict[str, Any]]) -> str:
        names = [Path(item["path"]).name for item in members]
        if any("job-manager-test" in name for name in names):
            return "Repeated JobManager test/workflow execution."
        if len({Path(item["path"]).parent for item in members}) == 1:
            return "Repeated direct solve or demonstration run with the same numerical plan."
        return "Matching numerical data; inspect provenance logs for the initiating run."

    def _external_references(self) -> set[str]:
        """Conservatively find file-path references outside a candidate itself."""
        refs: set[str] = set()
        for path in self.root.rglob("*.npz"):
            try:
                with np.load(path, allow_pickle=False) as archive:
                    if "metadata_json" in archive:
                        refs.update(self._paths_in_value(json.loads(str(archive["metadata_json"]))))
            except Exception as error:
                from .diagnostics import warning_event
                warning_event("archiver","reference_scan_npz_failed",error,path=str(path))
                continue
        for path in (self.root / "job_manager").rglob("*.json") if (self.root / "job_manager").exists() else ():
            try:
                refs.update(self._paths_in_value(json.loads(path.read_text(encoding="utf-8"))))
            except Exception as error:
                from .diagnostics import warning_event
                warning_event("archiver","reference_scan_manifest_failed",error,path=str(path))
                continue
        return refs

    @staticmethod
    def _paths_in_value(value: Any) -> set[str]:
        if isinstance(value, dict):
            return set().union(*(Archiver._paths_in_value(item) for item in value.values())) if value else set()
        if isinstance(value, list):
            return set().union(*(Archiver._paths_in_value(item) for item in value)) if value else set()
        if isinstance(value, str) and value.lower().endswith(".npz"):
            return {value}
        return set()

    def _write_report(self, report: ArchiveAuditReport) -> tuple[Path, Path]:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        base = self.report_root / f"archive-audit__{stamp}"
        json_path, markdown_path = base.with_suffix(".json"), base.with_suffix(".md")
        text = ["# dynsyslab archive audit", "", f"- Created: {report.created_utc}", f"- Scope: `{report.scope}`", f"- Scanned files: {report.scanned_files:,}", f"- Exact duplicate groups: {len(report.exact_duplicate_groups):,}", f"- Redundant files: {report.redundant_file_count:,}", f"- Potential reclaimable storage: {report.redundant_bytes / 2**20:.2f} MiB", ""]
        if report.exact_duplicate_groups:
            text.extend(["## Exact duplicate solutions", ""])
            for index, group in enumerate(report.exact_duplicate_groups, 1):
                text.extend([f"### Group {index}: {group['count']} files", "", f"Cause: {group['likely_cause']}", ""])
                text.extend([f"- `{item}`" for item in group["files"]]); text.append("")
        if report.same_plan_groups:
            text.extend(["## Repeated plans", "", f"{len(report.same_plan_groups)} groups have the same nominal plan; exact groups above are a subset.", ""])
        if report.unreadable_files:
            text.extend(["## Skipped files", "", f"{len(report.unreadable_files)} files were not readable at audit time (often an in-progress write).", ""])
        markdown_path.write_text("\n".join(text), encoding="utf-8")
        json_path.write_text(json.dumps(report.payload(), indent=2, sort_keys=True), encoding="utf-8")
        return json_path, markdown_path

    @staticmethod
    def _notify_archive_change(paths: list[str]) -> None:
        try:
            from infrastructure.source_events import notify_source_update
            for path in paths:
                notify_source_update(path, source_type="solution", reason="archiver_deduplicated")
        except Exception as error:
            # Viewer notifications are optional and must never make archive
            # management unavailable outside an interactive session.
            from .diagnostics import warning_event
            warning_event("archiver","viewer_notification_failed",error,changed_paths=paths)
