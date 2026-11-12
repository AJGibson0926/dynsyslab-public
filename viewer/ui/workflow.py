"""System-centric workflow organization for the native Viewer catalog.

Physical realization is the primary scientific identity.  Data roles (free,
forced, campaign, identification, or control test) remain tags and provenance;
they never create a second copy of the underlying dynamics.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha1
import json
from pathlib import Path
from typing import Any, Iterable

from .contracts import CatalogItem
from .grouping import archive_grouping, organize_training_data


_PHYSICAL_KEYS = ("physical_parameter_id", "physical_parameters", "scientific_name",
                  "dimensionless_groups", "nondimensionalization")
_STAGES = {
    "trajectory": "Data",
    "control_data": "Data",
    "experiment_data": "Data",
    "campaign": "Data",
    "krom": "KROMs and spectra",
    "koopman_spectrum": "KROMs and spectra",
    "observable_spectrum_comparison": "KROMs and spectra",
    "pod": "KROMs and spectra",
    "control_map": "Control identification",
    "control_error": "Control identification",
    "control_attempt": "Control attempts",
    "observability": "Estimation",
    "klqe": "Estimation",
}
_COLLECTION_KINDS = {"forcing_campaign", "training_campaign", "trajectory_ensemble",
                     "formulation_group"}


def controlled_kind(metadata: dict[str, Any]) -> str:
    """Separate targeted controller evaluations from identification forcing."""
    controller = metadata.get("controller") if isinstance(metadata.get("controller"), dict) else {}
    name = str(controller.get("controller") or controller.get("type") or
               controller.get("_dynsyslab_class") or "").casefold()
    target = controller.get("target") or metadata.get("control_target_present")
    target_history = metadata.get("control_target_history")
    if isinstance(target_history, dict):
        target = target or target_history.get("present")
    return "control_attempt" if target or "lqr" in name else "control_data"


def workflow_metadata(metadata: dict[str, Any]) -> tuple[list[tuple[str, str]],
                                                         list[str], dict[str, str]]:
    """Extract workflow links and display modifiers embedded in one archive."""
    relations: list[tuple[str, str]] = []
    tags: list[str] = []
    properties: dict[str, str] = {}
    controller = metadata.get("controller") if isinstance(metadata.get("controller"), dict) else {}
    rom = controller.get("rom") if isinstance(controller.get("rom"), dict) else {}
    source = rom.get("source_archive")
    if isinstance(source, str) and source:
        relations.append((source, "controlled with"))
    name = str(controller.get("controller") or controller.get("type") or "")
    if name:
        properties["controller"] = name
    if "koopman" in name.casefold() and "lqr" in name.casefold():
        tags.append("KLQR")
    target = controller.get("target")
    if isinstance(target, str) and target:
        properties["control target"] = target
    representation = controller.get("control_representation")
    if isinstance(representation, str) and representation:
        properties["control representation"] = representation
        tags.append(f"Bψ: {representation}")
    encoding = metadata.get("encoding") or metadata.get("schema")
    if isinstance(encoding, str) and encoding:
        properties["data encoding"] = encoding
    purpose = metadata.get("identification_purpose")
    if isinstance(purpose, str) and purpose:
        properties["identification purpose"] = purpose
        tags.extend(("forced", "Bψ regression data"))
    elif isinstance(encoding, str) and "identification" in encoding.casefold():
        tags.extend(("forced", "Bψ regression data"))
    for container in (metadata, metadata.get("analysis_metadata") or {}):
        if not isinstance(container, dict):
            continue
        for field, relation in (("training_archive", "trained on"),
                                ("control_data_archive", "regressed from"),
                                ("controlled_data_archive", "regressed from"),
                                ("source_krom_archive", "uses KROM")):
            value = container.get(field)
            if isinstance(value, str) and value:
                relations.append((value, relation))
    return relations, tags, properties


def _resolved(root: Path, raw: str | Path) -> Path:
    path = Path(raw).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def _copy_identity(item: CatalogItem, source: CatalogItem, reason: str) -> bool:
    identity = source.grouping.get("physical_parameter_id")
    if not identity or item.grouping.get("physical_parameter_id"):
        return False
    for key in _PHYSICAL_KEYS:
        item.grouping[key] = source.grouping.get(key)
    if not item.system:
        item.system = source.system
    item.properties["physical identity provenance"] = reason
    return True


def _infer_system(item: CatalogItem) -> None:
    if item.system:
        return
    parameters = item.grouping.get("physical_parameters") or {}
    keys = {str(key).casefold() for key in parameters}
    source = str(item.source).casefold()
    if {"equilibrium_radius_m", "ambient_pressure_pa"}.issubset(keys) or "rpe" in source:
        item.system = "Rayleigh-Plesset"
    elif {"alpha", "beta", "delta"}.issubset(keys) or "duffing" in source:
        item.system = "Duffing"


def _metadata_edges(root: Path, items: list[CatalogItem]) -> list[tuple[CatalogItem, CatalogItem, str]]:
    by_path = {str(item.source.resolve()).casefold(): item for item in items}
    edges: list[tuple[CatalogItem, CatalogItem, str]] = []
    for item in items:
        paths = list(item.relation_paths)
        if item.lineage_parent_path:
            paths.append((item.lineage_parent_path, "KROM lineage"))
        for raw, relation in paths:
            target = by_path.get(str(_resolved(root, raw)).casefold())
            if target is not None and target is not item:
                edges.append((item, target, relation))
    return edges


def _archive_references(value: Any, prefix: str = "") -> Iterable[tuple[str, str]]:
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{prefix}.{key}" if prefix else str(key)
            yield from _archive_references(child, path)
    elif isinstance(value, list):
        for child in value:
            yield from _archive_references(child, prefix + "[]")
    elif isinstance(value, str) and value.casefold().endswith(".npz"):
        yield prefix, value


def _embedded_identities(value: Any) -> Iterable[tuple[str, dict[str, Any]]]:
    if isinstance(value, dict):
        parameters = value.get("physical_parameters") or value.get("system_parameters")
        system = value.get("system_name") or value.get("system")
        if isinstance(parameters, dict) and parameters and isinstance(system, str):
            grouping = archive_grouping({"system_parameters": parameters})
            if grouping.get("physical_parameter_id"):
                yield system, grouping
        for child in value.values():
            yield from _embedded_identities(child)
    elif isinstance(value, list):
        for child in value:
            yield from _embedded_identities(child)


def _manifest_edges(root: Path, items: list[CatalogItem]) -> list[tuple[CatalogItem, CatalogItem, str]]:
    """Recover explicitly recorded legacy experiment context.

    Some old numerical NPZ files have no metadata header.  Their adjacent JSON
    manifests do name the KROM, solution, and error archives.  These links are
    deliberately marked inferred instead of being presented as embedded facts.
    """
    experiments = root / "archives" / "experiments"
    if not experiments.is_dir():
        return []
    by_path = {str(item.source.resolve()).casefold(): item for item in items}
    by_directory: dict[Path, list[CatalogItem]] = defaultdict(list)
    for item in items:
        if experiments in item.source.parents:
            by_directory[item.source.parent.resolve()].append(item)
    edges: list[tuple[CatalogItem, CatalogItem, str]] = []
    seen: set[tuple[str, str]] = set()
    manifest_identities: dict[Path, list[tuple[str, dict[str, Any]]]] = defaultdict(list)
    for manifest in experiments.rglob("*.json"):
        try:
            payload = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        manifest_identities[manifest.parent.resolve()].extend(_embedded_identities(payload))
        references = []
        for key, raw in _archive_references(payload):
            target = by_path.get(str(_resolved(root, raw)).casefold())
            if target is not None:
                references.append((key.casefold(), target))
        data_refs = [(key, item) for key, item in references
                     if item.kind == "control_data" or "dataset" in key]
        map_refs = [(key, item) for key, item in references
                    if item.kind == "control_map" or "candidate" in key]
        error_refs = [(key, item) for key, item in references
                      if item.kind == "control_error" or "error" in key]
        for _, result in map_refs:
            for _, data in data_refs:
                if result is not data:
                    result.relation_paths.append((str(data.source), "regressed from"))
                    edges.append((result, data, "regressed from"))
        for _, result in error_refs:
            for _, control_map in map_refs:
                if result is not control_map:
                    result.relation_paths.append((str(control_map.source), "diagnoses"))
                    edges.append((result, control_map, "diagnoses"))
        anchors = [(key, item) for key, item in references
                   if item.kind in {"krom", "trajectory"} or "krom" in key or "training" in key]
        if not anchors:
            continue
        anchor = next((item for key, item in anchors if item.kind == "krom"), anchors[0][1])
        candidates = [item for _, item in references]
        candidates.extend(by_directory.get(manifest.parent.resolve(), []))
        for item in candidates:
            edge = (item.identifier, anchor.identifier)
            if item is anchor or edge in seen:
                continue
            seen.add(edge)
            edges.append((item, anchor, "inferred from experiment manifest"))
            if not any(str(_resolved(root, raw)).casefold() == str(anchor.source.resolve()).casefold()
                       for raw, _ in item.relation_paths):
                item.relation_paths.append((str(anchor.source), "experiment uses"))
    for directory, candidates in manifest_identities.items():
        unique = {(system.casefold(), grouping["physical_parameter_id"])
                  for system, grouping in candidates}
        if len(unique) != 1:
            continue
        system, grouping = candidates[0]
        for item in items:
            if directory == item.source.parent.resolve() or directory in item.source.parents:
                if not item.grouping.get("physical_parameter_id"):
                    for key in _PHYSICAL_KEYS:
                        item.grouping[key] = grouping.get(key)
                    item.system = item.system or ("Duffing" if system.casefold() == "duffing" else system)
                    item.properties["physical identity provenance"] = "experiment manifest parameters"
    folders: dict[Path, list[CatalogItem]] = defaultdict(list)
    for item in items:
        try:
            relative = item.source.resolve().relative_to(experiments.resolve())
        except ValueError:
            continue
        if relative.parts:
            folders[experiments / relative.parts[0]].append(item)
    for members in folders.values():
        anchors = [item for item in members if item.grouping.get("physical_parameter_id")]
        identities = {(item.system, item.grouping.get("physical_parameter_id")) for item in anchors}
        if len(identities) != 1:
            continue
        anchor = anchors[0]
        for item in members:
            if item is not anchor and not item.grouping.get("physical_parameter_id"):
                edges.append((item, anchor, "inferred from single-realization experiment folder"))
    return edges


def _propagate_identity(root: Path, items: list[CatalogItem]) -> None:
    for item in items:
        _infer_system(item)
    edges = _metadata_edges(root, items) + _manifest_edges(root, items)
    for _ in range(max(2, len(items))):
        changed = False
        for item, target, relation in edges:
            changed |= _copy_identity(item, target, relation)
            if not target.grouping.get("physical_parameter_id") and item.grouping.get("physical_parameter_id"):
                # A manifest can name a result before its source.  Propagate only
                # within the same system to avoid conflating parameter variants.
                if not target.system or not item.system or target.system == item.system:
                    changed |= _copy_identity(target, item, f"related result: {relation}")
        if not changed:
            break


def _virtual_id(kind: str, key: str) -> str:
    return sha1(f"viewer:{kind}:{key}".encode("utf-8")).hexdigest()[:16]


def _flatten(value: Any, prefix: str = "") -> Iterable[tuple[str, Any]]:
    if isinstance(value, dict):
        for key in sorted(value):
            yield from _flatten(value[key], f"{prefix}.{key}" if prefix else str(key))
    elif not isinstance(value, (list, tuple)):
        yield prefix, value


def _fallback_parameter_group(root: Path, item: CatalogItem) -> CatalogItem:
    identity = str(item.grouping.get("physical_parameter_id") or item.identifier)
    scientific = str(item.grouping.get("scientific_name") or "")
    label = scientific or f"Parameter set {identity[:8]}"
    properties = {"parameter identity": identity}
    for heading, key in (("parameter", "physical_parameters"),
                         ("dimensionless group", "dimensionless_groups"),
                         ("scaling", "nondimensionalization")):
        for name, value in _flatten(item.grouping.get(key) or {}):
            properties[f"{heading} · {name}"] = str(value)
    return CatalogItem(
        _virtual_id("physical_parameter_set", f"{item.system}:{identity}"),
        "physical_parameter_set", "Systems", label, root, item.system,
        "physical realization", item.modified, 0, [],
        [item.system, "physical realization"], properties,
        role="physical realization", grouping={"virtual": True, **item.grouping})


def _stage_for(item: CatalogItem) -> str:
    if item.kind in _COLLECTION_KINDS:
        return "Data"
    return _STAGES.get(item.kind, "")


def _role_tags(item: CatalogItem) -> list[str]:
    tags: list[str] = []
    if item.kind == "trajectory":
        tags.append("forced" if item.grouping.get("forcing_parameters") else "free")
        if item.grouping.get("campaign_id"):
            tags.append("campaign member")
    elif item.kind == "control_attempt":
        tags.extend(("forced", "control attempt"))
    elif item.kind == "control_data":
        tags.extend(("forced", "control-identification data"))
    elif item.kind == "control_map":
        tags.append("Bψ regression")
    elif item.kind == "control_error":
        tags.append("Bψ diagnostic")
    return tags


def _stage_item(parameter: CatalogItem, stage: str, members: list[CatalogItem]) -> CatalogItem:
    kinds = Counter(item.kind for item in members)
    properties = {"workflow object count": str(len(members))}
    properties.update({f"objects · {kind.replace('_', ' ')}": str(count)
                       for kind, count in sorted(kinds.items())})
    return CatalogItem(
        _virtual_id("workflow_stage", f"{parameter.identifier}:{stage}"),
        "workflow_stage", "Systems", stage, parameter.source, parameter.system,
        f"workflow stage · {len(members)} object{'s' if len(members) != 1 else ''}",
        max((item.modified for item in members), default=parameter.modified),
        sum(item.byte_size for item in members), [],
        [parameter.system, "workflow stage"], properties,
        role="workflow stage", parent_identifier=parameter.identifier,
        primary=False, grouping={"virtual": True, "workflow_stage": stage})


def prepare_workflow(root: Path, items: list[CatalogItem],
                     memberships: dict[str, str]) -> list[CatalogItem]:
    """Create the system → physical realization → workflow-stage hierarchy."""
    _propagate_identity(root, items)
    virtual = organize_training_data(root, items, memberships)
    combined = items + virtual
    parameters = {(item.system, str(item.grouping.get("physical_parameter_id") or
                                    item.properties.get("parameter identity") or "")): item
                  for item in virtual if item.kind == "physical_parameter_set"}
    for item in items:
        identity = str(item.grouping.get("physical_parameter_id") or "")
        key = (item.system, identity)
        if identity and key not in parameters:
            parameter = _fallback_parameter_group(root, item)
            parameters[key] = parameter
            virtual.append(parameter)
            combined.append(parameter)

    by_id = {item.identifier: item for item in combined}
    stage_members: dict[tuple[str, str], list[CatalogItem]] = defaultdict(list)
    for item in combined:
        if item.kind in {"physical_parameter_set", "workflow_stage"}:
            continue
        identity = str(item.grouping.get("physical_parameter_id") or "")
        parameter = parameters.get((item.system, identity))
        if parameter is None:
            continue
        stage = _stage_for(item)
        if not stage:
            continue
        item.properties["workflow stage"] = stage
        item.properties.setdefault("physical parameter identity", identity)
        for tag in _role_tags(item):
            if tag not in item.tags:
                item.tags.append(tag)
        stage_members[(parameter.identifier, stage)].append(item)

    stages: dict[tuple[str, str], CatalogItem] = {}
    for (parameter_id, stage), members in stage_members.items():
        parameter = by_id[parameter_id]
        group = _stage_item(parameter, stage, members)
        stages[(parameter_id, stage)] = group
        virtual.append(group)
        combined.append(group)

    collection_ids = {item.identifier for item in virtual if item.kind in _COLLECTION_KINDS}
    for (parameter_id, stage), members in stage_members.items():
        stage_id = stages[(parameter_id, stage)].identifier
        for item in members:
            if item.parent_identifier in collection_ids:
                continue
            item.parent_identifier = stage_id
            item.primary = False

    for item in virtual:
        if item.kind not in _COLLECTION_KINDS:
            continue
        parent = by_id.get(item.parent_identifier)
        if parent is None or parent.kind != "physical_parameter_set":
            continue
        stage = stages.get((parent.identifier, "Data"))
        if stage is not None:
            item.parent_identifier = stage.identifier
            item.primary = False

    for parameter in parameters.values():
        related = [item for (parameter_id, _), members in stage_members.items()
                   if parameter_id == parameter.identifier for item in members]
        counts = Counter(item.kind for item in related)
        parameter.category = "Systems"
        parameter.role = "physical realization"
        parameter.tags = [tag for tag in (parameter.system,
                          str(parameter.grouping.get("scientific_name") or ""),
                          "physical realization") if tag]
        parameter.properties["workflow object count"] = str(len(related))
        scientific = str(parameter.grouping.get("scientific_name") or "")
        if scientific:
            parameter.properties["scientific shorthand"] = scientific
        for item in related:
            for heading, key in (("dimensionless group", "dimensionless_groups"),
                                 ("scaling", "nondimensionalization")):
                for name, value in _flatten(item.grouping.get(key) or {}):
                    parameter.properties.setdefault(f"{heading} · {name}", str(value))
        for label, kinds in (("training datasets", {"trajectory"}),
                             ("KROMs and spectra", {"krom", "koopman_spectrum"}),
                             ("control maps", {"control_map"}),
                             ("control attempts", {"control_attempt"})):
            parameter.properties[label] = str(sum(counts[kind] for kind in kinds))
        parameter.subtitle = (f"physical realization · {len(related)} linked workflow "
                              f"object{'s' if len(related) != 1 else ''}")
        parameter.primary = True
    return virtual
