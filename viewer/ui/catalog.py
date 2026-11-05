"""Fast metadata catalog for the native Viewer workbench.

Only archive headers and compact metadata are read. Numerical payloads remain
lazy and are requested through :mod:`viewer.ui.payloads` after user selection.
"""
from __future__ import annotations

from hashlib import sha1
import json
from pathlib import Path
import re
from typing import Any, Iterable
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import numpy as np

from .contracts import CATALOG_SCHEMA, CatalogItem, Relation, VIEW_LABELS
from .grouping import archive_grouping, campaign_grouping, campaign_memberships, is_virtual_group
from .workflow import controlled_kind, prepare_workflow, workflow_metadata

_SHAPE_KEYS = {"values", "x_values", "y_values", "A", "C", "state_estimate",
               "eigenvalues", "eigenvalues_continuous", "eigenfunction_values",
               "discrete_eigenvalues", "continuous_eigenvalues", "mode_counts",
               "modes", "basis", "singular_values"}
_GENERATOR_REVISION = 9

_SYSTEM_LABELS = {
    "duffing": "Duffing",
    "lorenz": "Lorenz",
    "marmottant": "Marmottant",
    "rayleigh-plesset": "Rayleigh-Plesset",
    "rayleighplesset": "Rayleigh-Plesset",
    "reddy-szeri": "Reddy-Szeri",
    "rossler": "Rossler",
    "rpe": "Rayleigh-Plesset",
}

def _identifier(root: Path, path: Path) -> str:
    try:
        key = path.resolve().relative_to(root.resolve()).as_posix().lower()
    except ValueError:
        key = str(path.resolve()).lower()
    return sha1(key.encode("utf-8")).hexdigest()[:16]
def _metadata(data: Any) -> dict[str, Any]:
    if "metadata_json" not in data.files:
        return {}
    try:
        value = json.loads(str(data["metadata_json"]))
        return value if isinstance(value, dict) else {}
    except (TypeError, ValueError, json.JSONDecodeError):
        return {}
def _scalar(value: Any, fallback: str = "") -> str:
    if value is None:
        return fallback
    if isinstance(value, (str, int, float, bool)):
        return str(value)
    return fallback
def _label(path: Path, metadata: dict[str, Any]) -> str:
    analysis = metadata.get("analysis_metadata") or {}
    reservation = analysis.get("reservation_label") if isinstance(analysis, dict) else None
    return _scalar(metadata.get("label") or metadata.get("name") or reservation, path.stem)
def _system_label(value: str) -> str:
    key = value.strip().casefold().replace("_", "-").replace(" ", "-")
    return _SYSTEM_LABELS.get(key, value.strip())
def _system(metadata: dict[str, Any]) -> str:
    direct = _scalar(metadata.get("system") or
                     (metadata.get("metadata") or {}).get("system"))
    if direct:
        return _system_label(direct)
    source = metadata.get("source_archive")
    if isinstance(source, str) and source:
        parent = Path(source).parent.name
        if parent:
            return _system_label(parent)
    analysis = metadata.get("analysis_metadata") or {}
    krom = analysis.get("krom") if isinstance(analysis, dict) else {}
    descriptor = krom.get("system") if isinstance(krom, dict) else None
    match = re.search(r"system_name=['\"]([^'\"]+)", str(descriptor or ""))
    return _system_label(match.group(1)) if match else ""


def _krom_role(analysis: str) -> str:
    if analysis == "koopman-krom-extract":
        return "extracted KROM"
    if analysis == "koopman-krom-control-representations":
        return "control-map variant"
    return "base KROM"


def _model_label(path: Path, metadata: dict[str, Any], kind: str,
                 krom: dict[str, Any]) -> str:
    explicit = _label(path, metadata)
    if explicit != path.stem:
        return explicit
    source = metadata.get("source_archive")
    source_name = Path(source).stem.split("__", 1)[0] if isinstance(source, str) else "Koopman"
    method = _scalar(krom.get("method") or metadata.get("analysis"), "Koopman")
    if "sindy" in method.lower():
        method = "SINDy"
    elif "hankel" in method.lower() or "hdmd" in method.lower():
        method = "H-DMD"
    count = krom.get("eigenfunction_count") or krom.get("full_spectrum_eigenfunction_count")
    object_name = "KROM" if kind == "krom" else "spectrum"
    count_label = f" · {count} modes" if count else ""
    return f"{source_name} · {method} {object_name}{count_label}"


def _array_shapes(path: Path, keys: set[str]) -> dict[str, tuple[int, ...]]:
    """Read NPY headers inside an NPZ without inflating numerical arrays."""
    result: dict[str, tuple[int, ...]] = {}
    try:
        with zipfile.ZipFile(path) as archive:
            for key in keys:
                member = f"{key}.npy"
                if member not in archive.namelist():
                    continue
                with archive.open(member) as stream:
                    version = np.lib.format.read_magic(stream)
                    reader = (np.lib.format.read_array_header_1_0 if version == (1, 0)
                              else np.lib.format.read_array_header_2_0)
                    shape, _, _ = reader(stream)
                    result[key] = tuple(int(item) for item in shape)
    except (OSError, ValueError, zipfile.BadZipFile):
        pass
    return result


def _array_shape(shapes: dict[str, tuple[int, ...]], key: str) -> str:
    return " x ".join(str(item) for item in shapes.get(key, ()))


def _control_present(keys: set[str], metadata: dict[str, Any]) -> bool:
    history = metadata.get("control_history") or {}
    return bool(metadata.get("controller") or history.get("present") or
                {"control_values", "controls", "controlled_controls"}.intersection(keys))


def _infer_npz(root: Path, path: Path) -> CatalogItem | None:
    try:
        with np.load(path, allow_pickle=False) as data:
            keys = set(data.files)
            metadata = _metadata(data)
            shapes = _array_shapes(path, keys.intersection(_SHAPE_KEYS))
            fmt = _scalar(metadata.get("format"))
            system = _system(metadata)
            controlled = _control_present(keys, metadata)
            analysis_metadata = metadata.get("analysis_metadata") or {}
            if not isinstance(analysis_metadata, dict):
                analysis_metadata = {}
            krom_metadata = analysis_metadata.get("krom") or {}
            if not isinstance(krom_metadata, dict):
                krom_metadata = {}
            analysis_name = _scalar(metadata.get("analysis"))
            role = "primary"
            lineage_parent_path = ""
            # A compact archive may reference one shared output time grid.
            has_time_series = ("values" in keys and
                               ("time" in keys or bool(metadata.get("time_grid_ref"))))
            if has_time_series:
                kind = controlled_kind(metadata) if controlled else "trajectory"
                category = ("Control attempts" if kind == "control_attempt" else
                            "Control data" if kind == "control_data" else "Training data")
                dimension = shapes.get("values", (0, 0))[-1]
                views = ["overview", "time"]
                if dimension >= 2:
                    views.append("phase")
                if controlled:
                    views.append("control")
                if dimension >= 3:
                    views.append("three_d")
                views.extend(("fourier", "observable_time"))
                if dimension == 2:
                    views.append("observable_heatmap")
                state_shape = _array_shape(shapes, "values")
                properties = {"state samples": state_shape}
            elif {"observable_names", "mode_counts", "discrete_eigenvalues",
                  "continuous_eigenvalues"}.issubset(keys):
                kind, category = "observable_spectrum_comparison", "Koopman spectra"
                views, properties = ["observable_spectrum_comparison"], {
                    "observables": _array_shape(shapes, "mode_counts"),
                    "spectra": _array_shape(shapes, "discrete_eigenvalues"),
                }
            elif "control-approximation-error" in fmt or {
                    "x_values", "y_values", "field_0_relative_error"}.issubset(keys):
                kind, category = "control_error", "Control-map errors"
                views, properties = ["control_error"], {
                    "phase grid": f"{_array_shape(shapes, 'x_values')} x {_array_shape(shapes, 'y_values')}"
                }
            elif {"A", "C", "gramian", "singular_values", "directions"}.issubset(keys):
                kind, category = "observability", "Estimation"
                views, properties = ["observability"], {
                    "state dimension": _array_shape(shapes, "A"),
                    "measurement map": _array_shape(shapes, "C"),
                }
            elif "state_estimate" in keys and "measurement" in keys:
                kind, category = "klqe", "Estimation"
                views, properties = ["klqe"], {
                    "estimated state": _array_shape(shapes, "state_estimate")
                }
            elif (analysis_name != "koopman-krom-control-representations" and
                  any(name.startswith("bilinear") or name.startswith("constant_B")
                      for name in keys)):
                kind, category = "control_map", "Control maps"
                views, properties = ["data_overview"], {"arrays": str(len(keys))}
            elif ({"eigenvalues", "eigenvalues_continuous"}.intersection(keys) or
                  (krom_metadata and ("coefficients" in keys or
                   ("eigenfunction_values" in keys and
                    {"training_states", "eigenfunction_state_samples"}.intersection(keys))))):
                retained = krom_metadata.get("retained_indices") or []
                analytic = "coefficients" in keys
                sampled = ("eigenfunction_values" in keys and
                           bool({"training_states", "eigenfunction_state_samples"}.intersection(keys)))
                executable = bool(krom_metadata) and (analytic or sampled)
                kind = "krom" if executable else "koopman_spectrum"
                category = "KROMs" if executable else "Koopman spectra"
                views = ["koopman_spectrum"]
                reconstructable = "eigenfunction_values" in keys or (
                    "coefficients" in keys and bool(retained))
                if executable and retained and reconstructable:
                    views.append("krom")
                eigenfunction_shape = _array_shape(shapes, "eigenfunction_values")
                if not eigenfunction_shape and "coefficients" in keys and retained:
                    eigenfunction_shape = "deferred from coefficients"
                if not eigenfunction_shape:
                    eigenfunction_shape = "not materialized"
                properties = {
                    "eigenvalues": _array_shape(shapes, "eigenvalues") or
                                   _array_shape(shapes, "eigenvalues_continuous") or
                                   str(krom_metadata.get("eigenfunction_count") or "linked source"),
                    "eigenfunctions": eigenfunction_shape,
                }
                if executable:
                    role = _krom_role(analysis_name)
                    parent = analysis_metadata.get("source_krom_archive") or metadata.get("parent_archive")
                    if isinstance(parent, str):
                        lineage_parent_path = parent
                else:
                    role = "spectrum"
            elif "singular_values" in keys or "modes" in keys or "basis" in keys:
                kind, category = "pod", "POD and reduced models"
                views, properties = ["data_overview"], {
                    "singular values": _array_shape(shapes, "singular_values"),
                    "modes": _array_shape(shapes, "modes") or _array_shape(shapes, "basis"),
                }
            elif controlled and any("state" in key for key in keys):
                kind, category = "control_data", "Control data"
                views, properties = ["data_overview"], {
                    "arrays": str(len(keys)), "controlled": "yes"
                }
            elif path.parts and "experiments" in path.parts:
                kind, category = "experiment_data", "Experiment data"
                views, properties = ["data_overview"], {"arrays": str(len(keys))}
            else:
                return None

            properties.update({
                "format": fmt or "npz",
                "size": f"{path.stat().st_size / (1024 * 1024):.2f} MiB",
            })
            integration = _scalar(metadata.get("integration_form"))
            tags = [value for value in (system, integration,
                    "controlled" if controlled else "") if value]
            if role != "primary":
                tags.append(role)
                properties["catalog role"] = role
            relation_paths = []
            for field, relation in (
                ("source_archive", "derived from"),
                ("parent_archive", "extracted from"),
                ("training_archive", "trained on"),
                ("control_data_archive", "regressed from"),
            ):
                value = metadata.get(field)
                if isinstance(value, str) and value:
                    relation_paths.append((value, relation))
            value = analysis_metadata.get("eigenvalue_source_archive")
            if isinstance(value, str) and value:
                relation_paths.append((value, "uses spectrum from"))
            source_krom = analysis_metadata.get("source_krom_archive")
            if isinstance(source_krom, str) and source_krom:
                relation = "extracted from" if role == "extracted KROM" else "variant of"
                relation_paths.append((source_krom, relation))
            workflow_relations, workflow_tags, workflow_properties = workflow_metadata(metadata)
            relation_paths.extend(workflow_relations)
            tags.extend(value for value in workflow_tags if value not in tags)
            properties.update(workflow_properties)
            label = (_model_label(path, metadata, kind, krom_metadata)
                     if kind in {"krom", "koopman_spectrum"} else _label(path, metadata))
            item = CatalogItem(
                _identifier(root, path), kind, category, label, path.resolve(),
                system, fmt or kind.replace("_", " "), path.stat().st_mtime,
                path.stat().st_size, views, tags, properties, relation_paths,
                role=role, lineage_parent_path=lineage_parent_path,
            )
            item.grouping = archive_grouping(metadata)
            if item.grouping.get("physical_parameter_id"):
                item.properties["physical parameter identity"] = item.grouping.get(
                    "physical_parameter_id") or "not recorded"
            if kind in {"trajectory", "control_attempt"}:
                item.properties["dynamics formulation"] = item.grouping.get("formulation", "not recorded")
            return item
    except (OSError, ValueError, EOFError):
        return None


def _report_item(root: Path, path: Path) -> CatalogItem:
    kind = "report"
    label = path.stem.replace("-", " ").replace("_", " ")
    return CatalogItem(
        _identifier(root, path), kind, "Reports", label, path.resolve(), "",
        path.suffix.lstrip(".").upper() + " Viewer artifact", path.stat().st_mtime,
        path.stat().st_size, ["open"], [path.suffix.lstrip(".")],
        {"format": path.suffix.lstrip(".").upper(),
         "size": f"{path.stat().st_size / (1024 * 1024):.2f} MiB"},
    )


def _campaign_item(root: Path, path: Path) -> CatalogItem | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return None
    if not isinstance(payload, dict):
        return None
    grouping = campaign_grouping(payload)
    metadata = payload.get("metadata") if isinstance(payload.get("metadata"), dict) else {}
    kind = _scalar(payload.get("campaign_type") or metadata.get("campaign_type") or
                   payload.get("kind") or payload.get("analysis_type"), "campaign")
    views = ["time", "phase", "fourier"]
    if "bifurcation" in kind.lower() or "parameter" in kind.lower():
        views.append("bifurcation")
    label = _scalar(payload.get("label") or payload.get("name"), path.stem)
    system = _system_label(_scalar(payload.get("system") or payload.get("system_name") or
                                   metadata.get("system") or metadata.get("system_name")))
    item = CatalogItem(
        _identifier(root, path), "campaign", "Campaigns", label, path.resolve(), system,
        kind.replace("_", " "), path.stat().st_mtime, path.stat().st_size, views,
        [value for value in (system, kind) if value],
        {"format": _scalar(payload.get("format"), "JSON"),
         "campaign id": grouping.get("campaign_id") or "not recorded",
         "campaign type": grouping.get("campaign_type") or kind,
         "jobs": str(len(payload.get("jobs") or [])),
         "size": f"{path.stat().st_size / 1024:.1f} KiB"},
    )
    item.grouping = grouping
    return item


def _candidate_paths(root: Path) -> Iterable[Path]:
    for section in ("solutions", "derived", "experiments"):
        base = root / "archives" / section
        if base.is_dir():
            yield from base.rglob("*.npz")
    reports = root / "viewer_outputs"
    if reports.is_dir():
        for suffix in ("*.html", "*.pdf"):
            yield from reports.glob(suffix)
    campaigns = root / "archives" / "campaigns"
    if campaigns.is_dir():
        yield from campaigns.glob("*.json")


def _resolve_relations(root: Path, items: list[CatalogItem]) -> None:
    by_path = {str(item.source.resolve()).lower(): item for item in items}
    for item in items:
        for raw, kind in item.relation_paths:
            path = Path(raw).expanduser()
            if not path.is_absolute():
                path = root / path
            target = by_path.get(str(path.resolve()).lower())
            if target is not None:
                item.relations.append(Relation(target.identifier, kind, target.label))
                target.relations.append(Relation(item.identifier, "related result", item.label))


def _resolve_lineages(root: Path, items: list[CatalogItem]) -> None:
    """Build KROM parent/child lineages without changing immutable archives."""
    by_path = {str(item.source.resolve()).lower(): item for item in items}
    by_id = {item.identifier: item for item in items}
    for item in items:
        if item.kind != "krom" or not item.lineage_parent_path:
            continue
        path = Path(item.lineage_parent_path).expanduser()
        if not path.is_absolute():
            path = root / path
        parent = by_path.get(str(path.resolve()).lower())
        if parent is not None and parent.kind == "krom" and parent is not item:
            item.parent_identifier = parent.identifier
            item.primary = False

    def lineage(identifier: str) -> str:
        seen: set[str] = set()
        current = by_id[identifier]
        while current.parent_identifier and current.parent_identifier not in seen:
            seen.add(current.identifier)
            parent = by_id.get(current.parent_identifier)
            if parent is None:
                break
            current = parent
        return current.identifier

    children: dict[str, list[str]] = {}
    for item in items:
        item.lineage_identifier = lineage(item.identifier)
        if item.parent_identifier:
            children.setdefault(item.parent_identifier, []).append(item.identifier)

    def descendants(identifier: str, seen: set[str] | None = None) -> int:
        visited = set() if seen is None else seen
        if identifier in visited:
            return 0
        visited.add(identifier)
        direct = children.get(identifier, [])
        return len(direct) + sum(descendants(child, visited) for child in direct)

    for item in items:
        item.variant_count = descendants(item.identifier)
        if item.variant_count:
            item.properties["related variants"] = str(item.variant_count)


def _disambiguate_labels(items: list[CatalogItem]) -> None:
    groups: dict[tuple[str, str], list[CatalogItem]] = {}
    for item in items:
        groups.setdefault((item.category, item.label.casefold()), []).append(item)
    for duplicates in groups.values():
        if len(duplicates) < 2:
            continue
        qualifiers = [item.source.parent.name for item in duplicates]
        if len(set(value.casefold() for value in qualifiers)) < len(qualifiers):
            qualifiers = [datetime.fromtimestamp(item.modified, timezone.utc).strftime(
                "%Y-%m-%d %H:%M:%S UTC") for item in duplicates]
        if len(set(value.casefold() for value in qualifiers)) < len(qualifiers):
            qualifiers = [item.identifier[:8] for item in duplicates]
        for item, qualifier in zip(duplicates, qualifiers):
            item.label = f"{item.label} — {qualifier}"


def build_catalog(project_root: str | Path) -> dict[str, Any]:
    root = Path(project_root).expanduser().resolve()
    paths = list(_candidate_paths(root))
    reports = [_report_item(root, path) for path in paths
               if path.suffix.lower() in {".html", ".pdf"}]
    campaigns = [_campaign_item(root, path) for path in paths if path.suffix.lower() == ".json"]
    archives = [path for path in paths if path.suffix.lower() == ".npz"]
    with ThreadPoolExecutor(max_workers=min(12, max(2, len(archives)))) as workers:
        inferred = list(workers.map(lambda path: _infer_npz(root, path), archives))
    items = reports + [item for item in campaigns if item is not None] + [item for item in inferred if item is not None]
    _disambiguate_labels(items)
    items.extend(prepare_workflow(root, items, campaign_memberships(root)))
    _resolve_lineages(root, items)
    _resolve_relations(root, items)
    items.sort(key=lambda item: (-item.modified, item.category, item.label.lower()))
    archive_counts: dict[str, int] = {}
    counts: dict[str, int] = {}
    member_counts: dict[str, int] = {}
    for item in items:
        if not is_virtual_group(item):
            archive_counts[item.category] = archive_counts.get(item.category, 0) + 1
        if item.primary:
            counts[item.category] = counts.get(item.category, 0) + 1
        if item.kind == "trajectory":
            member_counts[item.category] = member_counts.get(item.category, 0) + 1
    return {
        "schema": CATALOG_SCHEMA,
        "generator_revision": _GENERATOR_REVISION,
        "project": {"name": "DynSysLab", "root": str(root)},
        "scan": {"candidate_count": len(paths),
                 "latest_mtime": max((path.stat().st_mtime for path in paths), default=0.0)},
        "view_labels": VIEW_LABELS,
        "counts": counts,
        "archive_counts": archive_counts,
        "member_counts": member_counts,
        "items": [item.as_dict(root) for item in items],
    }


def write_catalog(project_root: str | Path, output: str | Path) -> Path:
    target = Path(output).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    root = Path(project_root).expanduser().resolve()
    if target.is_file():
        try:
            previous = json.loads(target.read_text(encoding="utf-8"))
            paths = list(_candidate_paths(root))
            scan = previous.get("scan") or {}
            latest = max((path.stat().st_mtime for path in paths), default=0.0)
            if (previous.get("schema") == CATALOG_SCHEMA and
                    int(previous.get("generator_revision", -1)) == _GENERATOR_REVISION and
                    (previous.get("project") or {}).get("root") == str(root) and
                    int(scan.get("candidate_count", -1)) == len(paths) and
                    float(scan.get("latest_mtime", -1.0)) == latest):
                return target
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            pass
    target.write_text(json.dumps(build_catalog(root), indent=2), encoding="utf-8")
    return target
