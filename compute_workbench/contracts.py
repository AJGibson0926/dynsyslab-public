"""Versioned, JSON-safe compute job contracts."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

JOB_SCHEMA = "dynsyslab.compute-job/v1"
RESULT_SCHEMA = "dynsyslab.compute-result/v1"
REGISTRY_SCHEMA = "dynsyslab.compute-registry/v1"


def solve_template() -> dict[str, Any]:
    """Return the versioned, non-executing default Solver job document."""
    return {
        "schema": JOB_SCHEMA,
        "kind": "solve",
        "label": "interactive-solve",
        "system": {"name": "duffing", "preset": "double_well_reference",
                   "parameters": {}},
        "time": {"start": 0.0, "stop": 10.0, "output_step": 0.01,
                 "coordinate_system": "automatic"},
        "initial_conditions": {"mode": "default", "count": 1, "clusters": 1,
                               "seed": 0, "values": [], "ranges": []},
        "controller": {"kind": "none", "settings": {},
                       "target": {"kind": "constant", "state": [0.0, 0.0]}},
        "campaign": {"enabled": False, "mode": "cartesian", "axes": []},
        "solver": {"method": "RK45", "rtol": 1e-8, "atol": 1e-10,
                   "max_step": None, "step_size": None,
                   "control_sample_time": None,
                   "control_evaluation": "automatic",
                   "execution_backend": "auto"},
        "archive": {"enabled": True},
    }


def postprocess_template() -> dict[str, Any]:
    """Return the versioned, non-executing default post-processing job."""
    return {
        "schema": JOB_SCHEMA,
        "kind": "postprocess",
        "label": "interactive-postprocess",
        "sources": [],
        "methods": [
            {"name": "hdmd", "enabled": True,
             "settings": {"observable": 0, "matrix_size": 200,
                          "sampling_window": 1000}},
            {"name": "sindy", "enabled": True,
             "settings": {"archive_storage": "coefficient-first"}},
        ],
        "archive": {"enabled": True},
    }


def load_job(source: str | Path | dict[str, Any]) -> dict[str, Any]:
    """Load and defensively copy one compute job mapping or JSON document."""
    if isinstance(source, dict):
        return json.loads(json.dumps(source))
    path = Path(source).expanduser().resolve()
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError("A compute job document must contain one JSON object.")
    return payload


def _positive_integer(value: Any, field: str) -> int:
    if isinstance(value, bool) or int(value) != value or int(value) < 1:
        raise ValueError(f"{field} must be a positive integer.")
    return int(value)


def validate_job(source: str | Path | dict[str, Any]) -> dict[str, Any]:
    """Validate structural and safety constraints without executing the job."""
    job = load_job(source)
    if job.get("schema") != JOB_SCHEMA:
        raise ValueError(f"Compute jobs must declare schema {JOB_SCHEMA!r}.")
    kind = job.get("kind")
    if kind not in {"solve", "postprocess"}:
        raise ValueError("Compute job kind must be 'solve' or 'postprocess'.")
    if kind == "solve":
        system = job.get("system") or {}
        if not system.get("name"):
            raise ValueError("A solve job requires system.name.")
        time = job.get("time") or {}
        start, stop = float(time.get("start", 0.0)), float(time.get("stop", 0.0))
        step = float(time.get("output_step", 0.0))
        if not start < stop or step <= 0:
            raise ValueError("Solve time requires start < stop and positive output_step.")
        ratio = (stop - start) / step
        if abs(ratio - round(ratio)) > 1e-9 * max(1.0, abs(ratio)):
            raise ValueError("The solve interval must be an integer multiple of output_step.")
        initial = job.get("initial_conditions") or {}
        count = _positive_integer(initial.get("count", 1), "initial_conditions.count")
        clusters = _positive_integer(initial.get("clusters", 1), "initial_conditions.clusters")
        if clusters > count:
            raise ValueError("initial_conditions.clusters cannot exceed trajectory count.")
        if initial.get("mode", "default") not in {"default", "explicit", "uniform", "normal"}:
            raise ValueError("Initial-condition mode must be default, explicit, uniform, or normal.")
        controller = (job.get("controller") or {}).get("kind", "none")
        if controller not in {"none", "waveform", "random_pulse", "klqr", "kmpc"}:
            raise ValueError("Unknown controller kind.")
        campaign = job.get("campaign") or {}
        if campaign.get("enabled", False):
            axes = campaign.get("axes") or []
            if not axes:
                raise ValueError("An enabled campaign requires at least one sweep axis.")
            if campaign.get("mode", "cartesian") != "cartesian":
                raise ValueError("The first-pass campaign contract supports cartesian sweeps.")
            variants = 1
            for axis in axes:
                path, values = str(axis.get("path", "")), axis.get("values") or []
                if not path.startswith("controller.") or not values:
                    raise ValueError("Campaign axes require a controller.* path and nonempty values.")
                variants *= len(values)
            if variants > 10000:
                raise ValueError("A UI campaign is limited to 10,000 System variants per job.")
    else:
        sources = job.get("sources") or []
        if not sources:
            raise ValueError("A postprocess job requires at least one archived source.")
        methods = [item for item in job.get("methods", []) if item.get("enabled", True)]
        if not methods:
            raise ValueError("A postprocess job requires at least one enabled method.")
        # Method names are checked by the native planner.
    return job
