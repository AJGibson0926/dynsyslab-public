"""Discover systems, presets, observables, targets, methods, and KROMs."""
from __future__ import annotations

from dataclasses import fields, is_dataclass
import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from control import RPE_CONTROL_TARGETS
from dynsyslib import ODE_LIBRARY, SYSTEM_PRESETS, load_system
from provenance import serializable

from .contracts import REGISTRY_SCHEMA


def _value(value: Any) -> Any:
    cleaned = serializable(value)
    return cleaned if isinstance(cleaned, (dict, list, str, int, float, bool)) or cleaned is None else repr(value)


def _derivative_count(value: dict[str, Any]) -> int:
    derivatives = value.get("derivative_observables") or []
    return len(derivatives)


def _parameters(system: Any) -> list[dict[str, Any]]:
    excluded = {"driver", "name", "preset", "initial_state_override"}
    if not is_dataclass(system):
        return []
    result = []
    for field in fields(system):
        if field.name in excluded or not field.init:
            continue
        value = getattr(system, field.name)
        result.append({"name": field.name, "value": _value(value),
                       "type": type(value).__name__})
    physical = getattr(system, "preset", None)
    if is_dataclass(physical):
        known = {item["name"] for item in result}
        for field in fields(physical):
            if field.name in known or not field.init:
                continue
            value = getattr(physical, field.name)
            result.append({"name": field.name, "value": _value(value),
                           "type": type(value).__name__, "scope": "physical_preset"})
    return result


def _system_entry(name: str) -> dict[str, Any]:
    presets = SYSTEM_PRESETS.get(name, {})
    default_preset = next(iter(presets), None)
    system = load_system(name, **({"preset": default_preset} if default_preset else {}))
    coordinates = "dimensionless" if (hasattr(system, "dimensionless_rhs") and
        bool(getattr(system, "prefer_dimensionless_form", True))) else "dimensional"
    catalog = system.observable_catalog(coordinates).summary()
    return {
        "name": name,
        "label": str(getattr(system, "name", name)).replace("_", " ").title(),
        "state_dimension": int(np.asarray(system.initial_state).size),
        "state_names": list(system.state_component_names(coordinates)),
        "state_labels": list(system.state_component_labels(coordinates)),
        "initial_state": np.asarray(system.initial_state, dtype=float).tolist(),
        "dimensionless_initial_state": (np.asarray(system.dimensionless_initial_state, dtype=float).tolist()
            if hasattr(system, "dimensionless_initial_state") else None),
        "preferred_coordinates": coordinates,
        "parameters": _parameters(system),
        "presets": [{"key": key, "description": preset.description,
                     "parameters": _value(getattr(preset, "parameters", preset)),
                     "initial_state": _value(getattr(preset, "initial_state", None)),
                     "processing": _value(getattr(preset, "processing_metadata", {}))}
                    for key, preset in presets.items()],
        "observables": [{"key": key,
                         "name": value.get("name", key),
                         "output_dimension": value.get("output_dimension"),
                         "state_dimension": value.get("state_dimension"),
                         "value_type": value.get("value_type"),
                         "analytic_jacobian": value.get("analytic_jacobian", False),
                         "descriptor": {item: (value.get("descriptor") or {}).get(item)
                                        for item in ("label", "units", "kind")
                                        if (value.get("descriptor") or {}).get(item) is not None},
                         "derivative_observable_count": _derivative_count(value)}
                        for key, value in catalog["observables"].items()],
    }


def _kroms(project_root: Path, catalog: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    if catalog is None:
        from viewer.ui.catalog import write_catalog
        base = Path(os.environ.get("LOCALAPPDATA", str(project_root / "build")))
        cache = base / "DynSysLab" / "Viewer" / "catalog.json"
        write_catalog(project_root, cache)
        catalog = json.loads(cache.read_text(encoding="utf-8"))
    return [{"id": item.get("id"), "label": item.get("label"),
             "source": item.get("source"), "system": item.get("system"),
             "subtitle": item.get("subtitle", "")}
            for item in catalog.get("items", []) if item.get("kind") == "krom"]


def build_registry(project_root: str | Path, *, catalog: dict[str, Any] | None = None) -> dict[str, Any]:
    """Build the versioned scientific choices exposed by the native workbench."""
    root = Path(project_root).expanduser().resolve()
    from postprocessing.processor import PostprocessingPipeline
    return {
        "schema": REGISTRY_SCHEMA,
        "project_root": str(root),
        "systems": [_system_entry(name) for name in ODE_LIBRARY],
        "controllers": [
            {"kind": "none", "label": "Free dynamics"},
            {"kind": "waveform", "label": "Sine-wave mixture"},
            {"kind": "random_pulse", "label": "Random pulsing"},
            {"kind": "klqr", "label": "Koopman LQR"},
            {"kind": "kmpc", "label": "Koopman MPC"},
        ],
        "targets": [{"name": key, "system": "rpe", "description": value.description,
                     "offset": value.offset, "amplitudes": list(value.amplitudes),
                     "angular_frequencies": list(value.angular_frequencies)}
                    for key, value in RPE_CONTROL_TARGETS.items()],
        "postprocessing_methods": sorted(PostprocessingPipeline.supported_methods),
        "kroms": _kroms(root, catalog),
    }
