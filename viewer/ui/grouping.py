"""Scientific identity and hierarchy assembly for Viewer data.

Forcing, initial conditions, solver settings, and formulation are provenance,
not part of the physical equations fingerprint.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha1
import json
from pathlib import Path
import re
import sqlite3
from typing import Any, Iterable

import numpy as np

from .contracts import CatalogItem, Relation
_NON_PHYSICAL_KEYS = {
    "caution", "control_affine", "description", "driver", "initial_state",
    "initial_state_override", "key", "name", "prefer_dimensionless_form",
    "pressure_input_convention", "control_input_convention",
    "use_dimensionless_form", "_dynsyslab_class", "_dynsyslab_module",
}
_GROUP_KINDS = {"physical_parameter_set", "forcing_campaign", "training_campaign",
                "trajectory_ensemble", "formulation_group", "workflow_stage"}
_MEMBER_SUFFIX = re.compile(r"-c\d+-t\d+$", re.IGNORECASE)

def _plain(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if str(key) in _NON_PHYSICAL_KEYS:
                continue
            cleaned = _plain(item)
            if cleaned not in (None, {}, []):
                result[str(key)] = cleaned
        return result
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)
def _physical_parameters(metadata: dict[str, Any]) -> dict[str, Any]:
    for key in ("system_parameters", "physical_parameters", "model_parameters",
                "parameters"):
        value = metadata.get(key)
        if isinstance(value, dict) and value:
            return _plain(value)
    nested = metadata.get("metadata")
    if isinstance(nested, dict):
        return _physical_parameters(nested)
    return {}
def _raw_parameters(metadata: dict[str, Any]) -> dict[str, Any]:
    for key in ("system_parameters", "physical_parameters", "model_parameters",
                "parameters"):
        value = metadata.get(key)
        if isinstance(value, dict) and value:
            return value
    nested = metadata.get("metadata")
    return _raw_parameters(nested) if isinstance(nested, dict) else {}
def _scientific_name(metadata: dict[str, Any]) -> str:
    parameters = _raw_parameters(metadata)
    preset = parameters.get("preset") if isinstance(parameters.get("preset"), dict) else {}
    description = preset.get("description") or parameters.get("description")
    if isinstance(description, str) and description.strip():
        return description.strip()
    gas = preset.get("gas")
    return str(gas).strip().capitalize() if isinstance(gas, str) and gas.strip() else ""
def _dimensionless_metadata(metadata: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    groups: dict[str, Any] = {}
    for key in ("dimensionless_groups", "dimensionless_parameters"):
        value = metadata.get(key)
        if isinstance(value, dict):
            groups.update(_plain(value))
    scaling = metadata.get("nondimensionalization")
    if not isinstance(scaling, dict):
        scaling = {}
    nested_groups = scaling.get("dimensionless_groups")
    if isinstance(nested_groups, dict):
        groups.update(_plain(nested_groups))
    scales: dict[str, Any] = {}
    for key, value in _plain(scaling).items():
        if key == "dimensionless_groups" or key in groups:
            continue
        if any(token in key.casefold() for token in ("scale", "frequency", "convention")):
            scales[key] = value
        else:
            groups[key] = value
    return groups, scales
def _fingerprint(parameters: dict[str, Any]) -> str:
    if not parameters:
        return ""
    canonical = parameters
    if {"equilibrium_radius_m", "ambient_pressure_pa"}.issubset(parameters):
        preset = parameters.get("preset") if isinstance(parameters.get("preset"), dict) else {}
        fields = ("equilibrium_radius_m", "ambient_pressure_pa", "kappa", "mu_pa_s",
                  "rho_kg_m3", "sigma_n_m", "vapor_pressure_pa")
        canonical = {key: parameters.get(key, preset.get(key)) for key in fields
                     if parameters.get(key, preset.get(key)) is not None}
        radius = canonical.get("equilibrium_radius_m")
        pressure = canonical.get("ambient_pressure_pa")
        vapor = canonical.get("vapor_pressure_pa")
        sigma = canonical.get("sigma_n_m")
        gas_pressure = parameters.get("initial_gas_pressure_pa")
        if gas_pressure is None and all(isinstance(value, (int, float))
                                        for value in (radius, pressure, vapor, sigma)):
            gas_pressure = pressure - vapor + 2.0 * sigma / radius
        if gas_pressure is not None:
            canonical["initial_gas_pressure_pa"] = gas_pressure
    payload = json.dumps(canonical, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=False)
    return sha1(payload.encode("utf-8")).hexdigest()[:16]
def _campaign_id(metadata: dict[str, Any]) -> str:
    candidates = [metadata.get("campaign_id")]
    for key in ("campaign", "campaign_attribution"):
        value = metadata.get(key)
        if isinstance(value, dict):
            candidates.append(value.get("campaign_id") or value.get("id"))
    job = metadata.get("job")
    if isinstance(job, dict):
        tags = job.get("tags")
        if isinstance(tags, dict):
            candidates.append(tags.get("campaign_id"))
    for value in candidates:
        if isinstance(value, str) and value:
            return value
    return ""


def _forcing_parameters(metadata: dict[str, Any]) -> dict[str, Any]:
    forcing = metadata.get("forcing")
    sources = forcing.get("sources") if isinstance(forcing, dict) else None
    result: dict[str, Any] = {}
    if not isinstance(sources, list):
        return result
    for source_index, source in enumerate(sources):
        if not isinstance(source, dict):
            continue
        parameters = source.get("parameters")
        if not isinstance(parameters, dict):
            continue
        prefix = str(source.get("source") or source.get("name") or
                     ("forcing" if len(sources) == 1 else f"forcing_{source_index}"))
        for key, value in parameters.items():
            if str(key).startswith("_dynsyslab_") or isinstance(value, (dict, list, tuple)):
                continue
            if isinstance(value, (str, int, float, bool)):
                result[f"{prefix}.{key}"] = value
    return result


def _formulation(metadata: dict[str, Any]) -> str:
    direct = str(metadata.get("integration_form") or "").strip().casefold()
    if direct in {"dimensionless", "nondimensional", "non-dimensional"}:
        return "dimensionless"
    if direct == "dimensional":
        return "dimensional"
    options = metadata.get("options")
    if isinstance(options, dict) and options.get("use_dimensionless_form") is not None:
        return "dimensionless" if options.get("use_dimensionless_form") else "dimensional"
    return "not recorded"


def archive_grouping(metadata: dict[str, Any]) -> dict[str, Any]:
    """Extract immutable grouping facts from one archive's compact metadata."""
    physical = _physical_parameters(metadata)
    label = str(metadata.get("label") or metadata.get("name") or "")
    collection = _MEMBER_SUFFIX.sub("", label) if _MEMBER_SUFFIX.search(label) else ""
    initial = metadata.get("initial_state_override")
    if initial is None and isinstance(metadata.get("system_parameters"), dict):
        initial = metadata["system_parameters"].get("initial_state_override")
    job = metadata.get("job") if isinstance(metadata.get("job"), dict) else {}
    groups, scaling = _dimensionless_metadata(metadata)
    return {
        "physical_parameter_id": _fingerprint(physical),
        "physical_parameters": physical,
        "scientific_name": _scientific_name(metadata),
        "dimensionless_groups": groups,
        "nondimensionalization": scaling,
        "campaign_id": _campaign_id(metadata),
        "campaign_label": str(job.get("label") or ""),
        "collection_id": collection,
        "formulation": _formulation(metadata),
        "forcing_parameters": _forcing_parameters(metadata),
        "initial_state": _plain(initial),
    }


def campaign_grouping(payload: dict[str, Any]) -> dict[str, Any]:
    metadata = payload.get("metadata") if isinstance(payload.get("metadata"), dict) else {}
    return {
        "campaign_id": str(payload.get("campaign_id") or metadata.get("campaign_id") or ""),
        "campaign_type": str(payload.get("campaign_type") or metadata.get("campaign_type") or
                             payload.get("kind") or "campaign"),
        "campaign_metadata": _plain(metadata),
        "campaign_goal": str(payload.get("goal") or ""),
    }


def campaign_memberships(root: Path) -> dict[str, str]:
    """Read the optional campaign-to-solution index without mutating it."""
    database = root / "archives" / "management" / "archive_index.sqlite"
    if not database.is_file():
        return {}
    try:
        connection = sqlite3.connect(f"file:{database.as_posix()}?mode=ro", uri=True)
        try:
            rows = connection.execute(
                "SELECT campaign_id, solution_path FROM campaign_solutions").fetchall()
        finally:
            connection.close()
    except sqlite3.Error:
        return {}
    result = {}
    for campaign, path in rows:
        candidate = Path(str(path))
        if candidate.is_file():
            result[str(candidate.resolve()).casefold()] = str(campaign)
    return result


def _virtual_id(kind: str, key: str) -> str:
    return sha1(f"viewer:{kind}:{key}".encode("utf-8")).hexdigest()[:16]


def _number(value: Any) -> str:
    return f"{float(value):.6g}" if isinstance(value, (int, float)) else str(value)


def _parameter_summary(system: str, parameters: dict[str, Any], identity: str,
                       scientific_name: str = "") -> str:
    if not parameters:
        return f"Parameters not recorded · {identity[:8]}"
    preset = parameters.get("preset") if isinstance(parameters.get("preset"), dict) else {}
    if system == "Rayleigh-Plesset":
        parts = [scientific_name] if scientific_name else []
        radius = parameters.get("equilibrium_radius_m")
        pressure = parameters.get("ambient_pressure_pa")
        kappa = preset.get("kappa")
        if isinstance(radius, (int, float)):
            parts.append(f"R₀ {float(radius) * 1e6:.6g} µm")
        if isinstance(pressure, (int, float)):
            parts.append(f"p∞ {float(pressure) / 1e3:.6g} kPa")
        if isinstance(kappa, (int, float)):
            parts.append(f"κ {float(kappa):.8g}")
        if parts:
            return " · ".join(parts)
    preferred = ("alpha", "beta", "delta", "sigma", "rho", "a", "b", "c",
                 "input_gain")
    parts = [f"{key}={_number(parameters[key])}" for key in preferred
             if key in parameters and isinstance(parameters[key], (int, float, str))]
    return " · ".join(parts[:5]) if parts else f"Parameter set {identity[:8]}"


def _flatten(value: Any, prefix: str = "") -> Iterable[tuple[str, Any]]:
    if isinstance(value, dict):
        for key in sorted(value):
            name = f"{prefix}.{key}" if prefix else str(key)
            yield from _flatten(value[key], name)
    elif not isinstance(value, (list, tuple)):
        yield prefix, value


def _parameter_properties(parameters: dict[str, Any], identity: str,
                          groups: dict[str, Any] | None = None,
                          scaling: dict[str, Any] | None = None) -> dict[str, str]:
    result = {"parameter identity": identity or "not recorded"}
    for key, value in _flatten(parameters):
        if not key or key.endswith("._dynsyslab_class"):
            continue
        result[f"parameter · {key}"] = _number(value)
    for key, value in _flatten(groups or {}):
        result[f"dimensionless group · {key}"] = _number(value)
    for key, value in _flatten(scaling or {}):
        result[f"scaling · {key}"] = _number(value)
    return result


def _scaled(value: float, key: str, system: str) -> tuple[float, str]:
    lowered = key.casefold()
    if lowered.endswith("frequency_hz") and abs(value) >= 1e3:
        return value / 1e3, "kHz"
    if lowered.endswith("frequency_hz"):
        return value, "Hz"
    if lowered.endswith("phase_rad"):
        return value, "rad"
    if lowered.endswith("amplitude") and system == "Rayleigh-Plesset":
        return value / 1e3, "kPa"
    return value, ""


def _forcing_properties(members: list[CatalogItem], system: str) -> dict[str, str]:
    fields: dict[str, list[Any]] = defaultdict(list)
    for item in members:
        for key, value in (item.grouping.get("forcing_parameters") or {}).items():
            fields[key].append(value)
    result = {}
    for key, values in sorted(fields.items()):
        unique = sorted(set(values), key=str)
        label = re.sub(r" (hz|rad|pa|kpa)$", "",
                       key.split(".", 1)[-1].replace("_", " "), flags=re.I)
        if all(isinstance(value, (int, float)) for value in unique):
            low, unit = _scaled(float(min(unique)), key, system)
            high, _ = _scaled(float(max(unique)), key, system)
            text = _number(low) if np.isclose(low, high) else f"{_number(low)}–{_number(high)}"
            if unit:
                text += f" {unit}"
            if len(unique) > 1:
                text += f" · {len(unique)} values"
        else:
            text = ", ".join(str(value) for value in unique[:6])
        result[f"forcing · {label}"] = text
    return result


def _planned_campaign_properties(metadata: dict[str, Any]) -> dict[str, str]:
    result = {}
    units = {"pa": "Pa", "kpa": "kPa", "hz": "Hz", "khz": "kHz", "rad": "rad"}
    for key, value in metadata.items():
        match = re.fullmatch(r"(.+)_range_(pa|kpa|hz|khz|rad)", str(key), re.I)
        if match and isinstance(value, list) and len(value) >= 2:
            label = match.group(1).replace("_", " ") + " range"
            result[f"planned · {label}"] = (
                f"{_number(value[0])}–{_number(value[-1])} {units[match.group(2).casefold()]}")
            continue
        match = re.fullmatch(r"(frequency|pressure|amplitude)_(pa|kpa|hz|khz)", str(key), re.I)
        if match and isinstance(value, (int, float)):
            scaled, unit = _scaled(float(value), str(key), "")
            result[f"planned · {match.group(1)}"] = f"{_number(scaled)} {unit or units[match.group(2).casefold()]}"
    return result


def _initial_condition_properties(members: list[CatalogItem]) -> dict[str, str]:
    states = []
    for item in members:
        value = item.grouping.get("initial_state")
        if isinstance(value, list) and value and all(isinstance(entry, (int, float)) for entry in value):
            states.append(np.asarray(value, dtype=float))
    if not states or len({state.size for state in states}) != 1:
        return {}
    matrix = np.vstack(states)
    return {f"initial state x{index}": f"{_number(matrix[:, index].min())}–{_number(matrix[:, index].max())}"
            for index in range(matrix.shape[1])}


def _group_item(kind: str, key: str, label: str, source: Path, system: str,
                parent: str, members: list[CatalogItem], properties: dict[str, str],
                role: str) -> CatalogItem:
    count = len(members)
    properties = dict(properties)
    properties["dataset count"] = str(count)
    return CatalogItem(
        identifier=_virtual_id(kind, key), kind=kind, category="Training data",
        label=label, source=source, system=system,
        subtitle=f"{role} · {count} dataset{'s' if count != 1 else ''}",
        modified=max((item.modified for item in members), default=0.0),
        byte_size=sum(item.byte_size for item in members), views=[],
        tags=[system, role], properties=properties, role=role,
        parent_identifier=parent, primary=not bool(parent),
        grouping={"virtual": True, "member_count": count},
    )


def organize_training_data(root: Path, items: list[CatalogItem],
                           memberships: dict[str, str]) -> list[CatalogItem]:
    """Nest trajectory archives without changing any archive payload."""
    training = [item for item in items if item.kind == "trajectory" and
                item.category == "Training data"]
    if not training:
        return []
    campaign_records = {item.grouping.get("campaign_id"): item for item in items
                        if item.kind == "campaign" and item.grouping.get("campaign_id")}
    parameter_sets: dict[tuple[str, str], list[CatalogItem]] = defaultdict(list)
    for item in training:
        campaign = item.grouping.get("campaign_id") or memberships.get(
            str(item.source.resolve()).casefold(), "")
        item.grouping["campaign_id"] = campaign
        identity = item.grouping.get("physical_parameter_id") or ""
        if not identity:
            identity = (f"unknown-campaign-{campaign}" if campaign else
                        f"unknown-collection-{item.grouping.get('collection_id')}" if item.grouping.get("collection_id")
                        else f"unknown-archive-{item.identifier}")
            item.grouping["physical_parameter_id"] = identity
        parameter_sets[(item.system, identity)].append(item)

    virtual: list[CatalogItem] = []
    for (system, identity), members in parameter_sets.items():
        parameters = members[0].grouping.get("physical_parameters") or {}
        scientific_name = str(members[0].grouping.get("scientific_name") or "")
        groups = members[0].grouping.get("dimensionless_groups") or {}
        scaling = members[0].grouping.get("nondimensionalization") or {}
        parameter_id = _virtual_id("physical_parameter_set", f"{system}:{identity}")
        parameter_source = members[0].source.parent
        parameter = _group_item(
            "physical_parameter_set", f"{system}:{identity}",
            _parameter_summary(system, parameters, identity, scientific_name), parameter_source,
            system, "", members, _parameter_properties(parameters, identity, groups, scaling),
            "physical parameter set")
        parameter.grouping.update({"physical_parameter_id": identity,
                                   "physical_parameters": parameters,
                                   "scientific_name": scientific_name,
                                   "dimensionless_groups": groups,
                                   "nondimensionalization": scaling})
        parameter.identifier = parameter_id
        virtual.append(parameter)

        owners: list[tuple[CatalogItem | None, list[CatalogItem], str]] = []
        by_campaign: dict[str, list[CatalogItem]] = defaultdict(list)
        remainder: list[CatalogItem] = []
        for member in members:
            campaign = member.grouping.get("campaign_id") or ""
            (by_campaign[campaign] if campaign else remainder).append(member)
        for campaign, campaign_members in by_campaign.items():
            record = campaign_records.get(campaign)
            campaign_metadata = (record.grouping.get("campaign_metadata") if record else {}) or {}
            labels = [str(item.grouping.get("campaign_label") or "") for item in campaign_members]
            label = record.label if record else next((value for value in labels if value),
                                                      f"Campaign {campaign[:8]}")
            forced = any(item.grouping.get("forcing_parameters") for item in campaign_members)
            group_kind = "forcing_campaign" if forced else "training_campaign"
            group_role = "forcing campaign" if forced else "training campaign"
            properties = {
                "campaign id": campaign,
                "campaign type": (record.grouping.get("campaign_type") if record else "") or "campaign",
                **_planned_campaign_properties(campaign_metadata),
                **_forcing_properties(campaign_members, system),
            }
            owner = _group_item(
                group_kind, f"{system}:{identity}:{campaign}", label,
                (record.source if record else campaign_members[0].source.parent), system,
                parameter_id, campaign_members, properties, group_role)
            if record:
                owner.relations.append(Relation(record.identifier, "campaign record", record.label))
                record.relations.append(Relation(owner.identifier, "contains training data", owner.label))
            virtual.append(owner)
            owners.append((owner, campaign_members, "campaign member"))

        by_collection: dict[str, list[CatalogItem]] = defaultdict(list)
        standalone: list[CatalogItem] = []
        for member in remainder:
            collection = member.grouping.get("collection_id") or ""
            (by_collection[collection] if collection else standalone).append(member)
        for collection, collection_members in by_collection.items():
            if len(collection_members) < 2:
                standalone.extend(collection_members)
                continue
            owner = _group_item(
                "trajectory_ensemble", f"{system}:{identity}:{collection}", collection,
                collection_members[0].source.parent, system, parameter_id,
                collection_members, _initial_condition_properties(collection_members),
                "trajectory ensemble")
            virtual.append(owner)
            owners.append((owner, collection_members, "ensemble member"))
        if standalone:
            owners.append((None, standalone, "training dataset"))

        for owner, owner_members, member_role in owners:
            parent = owner.identifier if owner else parameter_id
            formulations: dict[str, list[CatalogItem]] = defaultdict(list)
            for member in owner_members:
                formulations[str(member.grouping.get("formulation") or "not recorded")].append(member)
            use_formulation_groups = len(formulations) > 1
            for formulation, formulation_members in formulations.items():
                member_parent = parent
                if use_formulation_groups:
                    form = _group_item(
                        "formulation_group", f"{parent}:{formulation}",
                        f"{formulation.capitalize()} formulation",
                        formulation_members[0].source.parent, system, parent,
                        formulation_members, {"dynamics formulation": formulation},
                        "formulation variant")
                    virtual.append(form)
                    member_parent = form.identifier
                for member in formulation_members:
                    member.parent_identifier = member_parent
                    member.primary = False
                    member.role = member_role
                    member.properties["dynamics formulation"] = formulation
                    if formulation not in member.tags:
                        member.tags.append(formulation)
    physical_groups = [item for item in virtual if item.kind == "physical_parameter_set"]
    duplicate_labels: dict[tuple[str, str], list[CatalogItem]] = defaultdict(list)
    for item in physical_groups:
        duplicate_labels[(item.system, item.label.casefold())].append(item)
    for duplicates in duplicate_labels.values():
        if len(duplicates) > 1:
            for item in duplicates:
                identity = item.properties.get("parameter identity", item.identifier)
                item.label += f" · {identity[:8]}"
    return virtual


def is_virtual_group(item: CatalogItem) -> bool:
    return item.kind in _GROUP_KINDS
