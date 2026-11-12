"""Plan and execute compute jobs through the canonical DynSysLab services."""
from __future__ import annotations

from dataclasses import fields, is_dataclass, replace
from datetime import datetime, timezone
import itertools
import json
import os
from pathlib import Path
from typing import Any
import uuid

import numpy as np

from dynsyslib import load_system
from hierarchy import Batch, System
from postprocessing import PostprocessingPipeline, load_solution
from postprocessing.planning import summarize_processing_runs
from solver.hierarchy import solve_hierarchy
from solver.models import HistoryOptions, SolverOptions
from provenance import fingerprint, serializable

from .contracts import RESULT_SCHEMA, load_job, validate_job
from .controllers import build_controller
from .initial_conditions import generate_initial_conditions, make_clusters


def _system(job: dict[str, Any]):
    spec = job["system"]
    name, preset = spec["name"], spec.get("preset")
    parameters = dict(spec.get("parameters") or {})
    prototype = load_system(name, **({"preset": preset} if preset else {}))
    system_fields = {item.name for item in fields(prototype) if item.init} if is_dataclass(prototype) else set()
    direct = {key: value for key, value in parameters.items() if key in system_fields}
    system = load_system(name, **({"preset": preset} if preset else {}), **direct)
    remaining = {key: value for key, value in parameters.items() if key not in direct}
    physical = getattr(system, "preset", None)
    if remaining and is_dataclass(physical):
        physical_fields = {item.name for item in fields(physical) if item.init}
        physical_values = {key: value for key, value in remaining.items() if key in physical_fields}
        if physical_values:
            system = replace(system, preset=replace(physical, **physical_values))
            remaining = {key: value for key, value in remaining.items() if key not in physical_values}
    if remaining:
        raise ValueError(f"Unknown {name} parameters: {sorted(remaining)}")
    return system


def _coordinate_policy(job: dict[str, Any], system: Any) -> tuple[bool | None, str]:
    selected = (job.get("time") or {}).get("coordinate_system", "automatic")
    if selected == "automatic":
        dimensionless = None
        active = ("dimensionless" if hasattr(system, "dimensionless_rhs") and
                  bool(getattr(system, "prefer_dimensionless_form", True)) else "dimensional")
    elif selected == "dimensionless":
        if not hasattr(system, "dimensionless_rhs"):
            raise ValueError(f"{system.name} does not implement dimensionless dynamics.")
        dimensionless, active = True, "dimensionless"
    elif selected == "dimensional":
        dimensionless, active = False, "dimensional"
    else:
        raise ValueError("time.coordinate_system must be automatic, dimensionless, or dimensional.")
    return dimensionless, active


def _physical_identity(system: Any) -> str:
    excluded = {"driver", "name", "initial_state_override", "control_affine",
                "prefer_dimensionless_form", "description", "caution", "key"}
    def clean(value):
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items() if key not in excluded}
        if isinstance(value, list): return [clean(item) for item in value]
        return value
    return fingerprint(clean(serializable(system)))


def _controller_system(system: Any, controller: Any) -> Any:
    krom = getattr(controller, "rom", None)
    if krom is None: return system
    if _physical_identity(system) != _physical_identity(krom.system):
        raise ValueError("The selected KROM and Solver system use different physical parameters; choose the KROM's matching physical realization.")
    # Koopman controllers deliberately bind to the exact executable system
    # retained by the standalone KROM, not an equivalent parallel instance.
    return krom.system


def _physical_mesh(job: dict[str, Any], system: Any, active: str) -> tuple[float, float, float]:
    time = job["time"]
    start, stop, step = (float(time["start"]), float(time["stop"]),
                         float(time["output_step"]))
    if active == "dimensionless":
        convert = getattr(system, "to_dimensional_time", None)
        if not callable(convert):
            raise TypeError("Dimensionless jobs require system.to_dimensional_time().")
        start, stop = float(np.asarray(convert(start))), float(np.asarray(convert(stop)))
        physical_zero = float(np.asarray(convert(0.0)))
        step = float(np.asarray(convert(float(time["output_step"])))) - physical_zero
    return start, stop, step


def _physical_initial_conditions(job: dict[str, Any], system: Any,
                                 active: str, controlled: bool) -> np.ndarray:
    spec = job["initial_conditions"]
    values = generate_initial_conditions(system, spec, controlled=controlled)
    if active == "dimensionless" and spec.get("mode", "default") != "default":
        convert = getattr(system, "to_dimensional_state", None)
        if not callable(convert):
            raise TypeError("Dimensionless initial conditions require system.to_dimensional_state().")
        values = np.stack([np.asarray(convert(row), dtype=float) for row in values])
    return values


def _solver_options(job: dict[str, Any], dimensionless: bool | None) -> SolverOptions:
    settings = dict(job.get("solver") or {})
    allowed = {item.name for item in fields(SolverOptions)}
    settings = {key: value for key, value in settings.items() if key in allowed and value is not None}
    settings["use_dimensionless_form"] = dimensionless
    return SolverOptions(**settings)


def _source_path(root: Path, source: str | Path) -> Path:
    path = Path(source).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def _set_path(payload: dict[str, Any], path: str, value: Any) -> None:
    tokens = path.split(".")
    current: Any = payload
    for token in tokens[:-1]:
        if isinstance(current, list):
            current = current[int(token)]
        else:
            current = current.setdefault(token, {})
    last = tokens[-1]
    if isinstance(current, list): current[int(last)] = value
    else: current[last] = value


def _variants(job: dict[str, Any]) -> list[tuple[dict[str, Any], dict[str, Any]]]:
    campaign = job.get("campaign") or {}
    if not campaign.get("enabled", False):
        return [(json.loads(json.dumps(job)), {})]
    axes = campaign["axes"]
    variants = []
    for values in itertools.product(*(axis["values"] for axis in axes)):
        variant = json.loads(json.dumps(job)); coordinates = {}
        for axis, value in zip(axes, values):
            _set_path(variant, axis["path"], value); coordinates[axis["path"]] = value
        variant["campaign"]["enabled"] = False
        variants.append((variant, coordinates))
    return variants


def _campaign_manifest(root: Path, job: dict[str, Any], campaign_id: str,
                       variants: list[tuple[dict[str, Any], dict[str, Any]]]) -> Path:
    directory = root / "archives" / "campaigns"; directory.mkdir(parents=True, exist_ok=True)
    label = job.get("label", "compute-campaign")
    safe = "".join(value if value.isalnum() or value in "-_" else "-" for value in label)
    path = directory / f"{safe}__{campaign_id}.json"
    payload = {"format": "dynsyslab.campaign/v1", "campaign_id": campaign_id,
               "campaign_type": "control_sweep", "label": label,
               "goal": (job.get("campaign") or {}).get("goal", "Interactive control sweep"),
               "created_utc": datetime.now(timezone.utc).isoformat(),
               "jobs": [f"compute-{campaign_id}"],
               "system": job["system"]["name"],
               "metadata": {"campaign_id": campaign_id, "campaign_type": "control_sweep",
                            "system_name": job["system"]["name"],
                            "preset": job["system"].get("preset"),
                            "hierarchy": "Batch > swept System > Cluster > Trajectory",
                            "axes": job["campaign"]["axes"],
                            "variant_count": len(variants),
                            "trajectories_per_system": job["initial_conditions"].get("count", 1)}}
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def plan_reference_job(source: dict[str, Any], project_root: str | Path) -> dict[str, Any]:
    """Run the retained Python planning path for explicit parity studies."""
    job = validate_job(source)
    root = Path(project_root).expanduser().resolve()
    if job["kind"] == "solve":
        system = _system(job)
        controller = build_controller(job["controller"], system)
        system = _controller_system(system, controller)
        dimensionless, active = _coordinate_policy(job, system)
        start, stop, step = _physical_mesh(job, system, active)
        variants = _variants(job)
        count = int(job["initial_conditions"].get("count", 1))
        values = _physical_initial_conditions(job, system, active,
                                             controller is not None)
        return {"schema": RESULT_SCHEMA, "status": "planned", "kind": "solve",
                "system": system.name, "coordinate_system": active,
                "trajectory_count": count * len(variants),
                "system_count": len(variants),
                "campaign": len(variants) > 1,
                "cluster_count": int(job["initial_conditions"].get("clusters", 1)) * len(variants),
                "clusters_per_system": int(job["initial_conditions"].get("clusters", 1)),
                "state_dimension": int(values.shape[1]),
                "physical_time_span": [start, stop], "physical_output_step": step,
                "output_samples_per_trajectory": int(round((stop-start)/step))+1,
                "controller": job["controller"].get("kind", "none"),
                "archive_policy": "canonical solver archive (mandatory)"}
    details = []
    for value in job["sources"]:
        path = _source_path(root, value)
        trajectory = load_solution(path)
        details.append({"source": str(path), "samples": trajectory.time_series.sample_count,
                        "state_shape": list(trajectory.values.shape[1:]),
                        "system": trajectory.metadata.get("system")})
    methods = [{"name": item["name"].lower(), "settings": item.get("settings", {})}
               for item in job["methods"] if item.get("enabled", True)]
    return {"schema": RESULT_SCHEMA, "status": "planned", "kind": "postprocess",
            "sources": details, "methods": methods,
            "archive": bool((job.get("archive") or {}).get("enabled", True))}


def _refresh_catalog(root: Path, output: str | Path | None) -> str | None:
    if output is None:
        return None
    from viewer.ui.catalog import write_catalog
    return str(write_catalog(root, output))


def execute_reference_job(source: dict[str, Any], project_root: str | Path, *,
                          catalog_output: str | Path | None = None,
                          progress_output: str | Path | None = None) -> dict[str, Any]:
    """Run the retained Python executor for explicit parity studies only."""
    job = validate_job(source)
    root = Path(project_root).expanduser().resolve()
    reporter = None
    if progress_output is not None:
        from .progress import ProgressReporter
        planned = (len(_variants(job)) * int(job["initial_conditions"].get("count", 1))
                   if job["kind"] == "solve" else 0)
        reporter = ProgressReporter(progress_output, kind=job["kind"],
                                    planned_trajectories=planned)
    try:
        payload = _execute_job(job, root, catalog_output=catalog_output,
                               reporter=reporter)
    except Exception as error:
        if reporter is not None:
            reporter.fail(str(error))
        raise
    if reporter is not None:
        reporter.complete()
    return payload


def plan_job(source: dict[str, Any], project_root: str | Path) -> dict[str, Any]:
    """Plan through the canonical compiled workflow service."""
    from .native_workflow import plan_job as native_plan
    return native_plan(load_job(source), project_root=project_root)


def execute_job(source: dict[str, Any], project_root: str | Path, *,
                catalog_output: str | Path | None = None,
                progress_output: str | Path | None = None) -> dict[str, Any]:
    """Execute through the canonical compiled workflow service."""
    from .native_workflow import run_job
    root = Path(project_root).expanduser().resolve()
    archive = Path(os.environ.get("DYNSYSLAB_ARCHIVE_DIR", str(root / "archives")))
    payload = run_job(load_job(source), archive, project_root=root,
                      progress_path=progress_output)
    catalog = _refresh_catalog(root, catalog_output)
    if catalog is not None:
        payload["catalog"] = catalog
    return payload


def _execute_job(job: dict[str, Any], root: Path, *,
                 catalog_output: str | Path | None,
                 reporter: Any = None) -> dict[str, Any]:
    """Execute one validated job; ``execute_job`` owns lifecycle reporting."""
    if job["kind"] == "solve":
        variants = _variants(job)
        campaign_id = uuid.uuid4().hex if len(variants) > 1 else None
        manifest = _campaign_manifest(root, job, campaign_id, variants) if campaign_id else None
        systems, controllers = [], []
        dimensionless = None; active = "dimensional"; start = stop = step = 0.0
        for index, (variant, coordinates) in enumerate(variants):
            system = _system(variant)
            selected_controller = build_controller(variant["controller"], system)
            system = _controller_system(system, selected_controller)
            dimensionless, active = _coordinate_policy(variant, system)
            start, stop, step = _physical_mesh(variant, system, active)
            values = _physical_initial_conditions(variant, system, active, selected_controller is not None)
            clusters = make_clusters(values, int(variant["initial_conditions"].get("clusters", 1)))
            metadata = {"compute_job_schema": variant["schema"], "sweep_coordinates": coordinates,
                        "controller_configuration": variant["controller"]}
            if campaign_id:
                metadata.update({"campaign_id": campaign_id,
                                 "campaign_attribution": {"campaign_id": campaign_id,
                                     "manifest": str(manifest), "system_index": index},
                                 "job": {"label": job.get("label", "control-sweep"),
                                         "tags": {"campaign_id": campaign_id}}})
            systems.append(System(system, clusters=clusters, time_span=(start, stop),
                                  uniform_time_steps=True, delta_t=step,
                                  label=f"{job.get('label', 'interactive-solve')}-s{index:03d}",
                                  metadata=metadata))
            controllers.append(selected_controller)
        target = (systems[0] if len(systems) == 1 else
                  Batch(systems, label=job.get("label", "control-sweep"),
                        metadata={"campaign_id": campaign_id,
                                  "campaign_manifest": str(manifest),
                                  "sweep_axes": job["campaign"]["axes"]}))
        if manifest:
            from archiver.gateway import MANAGER_ARCHIVE
            MANAGER_ARCHIVE.register_campaign(manifest, status="running")
            MANAGER_ARCHIVE.record_campaign_job(campaign_id, f"compute-{campaign_id}",
                status="running", job_label=job.get("label"), manifest_path=manifest)
        history = HistoryOptions(telemetry_interval_seconds=1.0) if reporter else HistoryOptions()
        result = solve_hierarchy(target, options=_solver_options(job, dimensionless),
                                 history=history, controllers=controllers,
                                 on_event=None if reporter is None else reporter.event)
        archives = [item.archive_path for item in result.trajectory_results]
        if manifest:
            MANAGER_ARCHIVE.record_campaign_solutions(campaign_id, f"compute-{campaign_id}", archives)
            MANAGER_ARCHIVE.record_campaign_job(campaign_id, f"compute-{campaign_id}",
                status="complete", job_label=job.get("label"), manifest_path=manifest)
        payload = {"schema": RESULT_SCHEMA, "status": result.status, "kind": "solve",
                   "hierarchy_archive": result.archive_path,
                   "solution_archives": archives,
                   "trajectory_count": len(archives),
                   "system_count": len(systems),
                   "cluster_count": sum(len(item.clusters) for item in systems),
                   "coordinate_system": active,
                   "campaign_id": campaign_id,
                   "campaign_manifest": None if manifest is None else str(manifest)}
    else:
        archive = bool((job.get("archive") or {}).get("enabled", True))
        pipeline = PostprocessingPipeline(
            [_source_path(root, item) for item in job["sources"]], archive=archive)
        for method in job["methods"]:
            if method.get("enabled", True):
                pipeline.add_method(method["name"], **dict(method.get("settings") or {}))
        runs = pipeline.run()
        payload = {"schema": RESULT_SCHEMA, "status": "complete", "kind": "postprocess",
                   "processing": summarize_processing_runs(runs),
                   "derived_archives": [str(item.archive) for run in runs for item in run.results
                                        if item.archive is not None]}
    payload["catalog"] = _refresh_catalog(root, catalog_output)
    return payload
