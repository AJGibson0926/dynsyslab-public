"""Read-only archive records and dimension-aware display transforms."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import numpy as np

from archiver.gateway import VIEWER_ARCHIVE
from time_series import TimeSeries
from . import sources as _sources
from .contracts import ArchiveSource, _ViewerArtifact

def _manifest(path: Path) -> dict[str,Any]:
    with VIEWER_ARCHIVE.open_read(path) as data:
        if "metadata_json" not in data.files: raise ValueError(f"{path} is not a dynsyslab archive.")
        return json.loads(str(data["metadata_json"]))

def _trajectory_paths(manifest: dict[str,Any], system_index: int, cluster_index: int | None=None, trajectory_index: int | None=None) -> list[Path]:
    systems=manifest.get("systems")
    if systems is None: raise ValueError("The requested level requires a tiered hierarchy manifest.")
    system=systems[system_index]; clusters=system["clusters"] if cluster_index is None else [system["clusters"][cluster_index]]
    paths=[]
    for cluster in clusters:
        trajectories=cluster["trajectories"] if trajectory_index is None else [cluster["trajectories"][trajectory_index]]
        paths.extend(Path(item["solution_archive"]) for item in trajectories)
    if not paths or any(not item.is_file() for item in paths): raise FileNotFoundError("A referenced child trajectory archive is unavailable.")
    return paths

def _load_archive_source(source: ArchiveSource) -> list[tuple[str,TimeSeries,dict[str,Any]]]:
    return _sources.load_archive_source(source)

def _solution_record(path: Path) -> tuple[TimeSeries,dict[str,Any]]:
    with VIEWER_ARCHIVE.open_read(path) as data: metadata=json.loads(str(data["metadata_json"]))
    metadata["source_archive"]=str(path)
    return VIEWER_ARCHIVE.read_time_series(path),metadata

def _viewer_compute_record(source: _ViewerArtifact) -> tuple[dict[str,np.ndarray],dict[str,Any]]:
    """Read a display artifact; Viewer never creates or changes this file."""
    path=Path(source.path).expanduser().resolve()
    with np.load(path,allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        if metadata.get("format")!="dynsyslab.viewer-compute/v1": raise ValueError(f"{path} is not a viewer-compute artifact.")
        arrays={name:np.asarray(data[name]) for name in data.files if name!="metadata_json"}
    metadata["source_viewer_compute_archive"]=str(path)
    return arrays,metadata

def _axis_context(metadata: dict[str,Any]) -> dict[str,Any]:
    """Return display transforms and physical/state-aware labels without a solver import."""
    return _sources.axis_context(metadata)
    system=metadata.get("system",""); form=metadata.get("integration_form","dimensional")
    parameters=metadata.get("system_parameters",{})
    if system=="rayleigh-plesset":
        if form=="dimensionless":
            preset=parameters.get("preset",{}); radius=float(parameters["equilibrium_radius_m"]); pressure=float(parameters["ambient_pressure_pa"]); density=float(preset["rho_kg_m3"])
            scale_time=radius*np.sqrt(density/pressure)
            return {"time_scale":1/scale_time,"state_scale":np.array((1/radius,scale_time/radius)),"time_label":"τ","state_labels":["r = R/R₀","dr/dτ"],"frequency_label":"frequency (cycles / τ)"}
        return {"time_scale":1e6,"state_scale":np.array((1e6,1.)),"time_label":"t (µs)","state_labels":["R (µm)","dR/dt (m/s)"],"frequency_label":"frequency (MHz)"}
    return {"time_scale":1.,"state_scale":None,"time_label":"t","state_labels":[],"frequency_label":"frequency"}

def _display_series(series: TimeSeries, metadata: dict[str,Any]) -> tuple[np.ndarray,np.ndarray,dict[str,Any]]:
    return _sources.display_series(series,metadata)

def _display_analysis_series(series: TimeSeries, source_metadata: dict[str,Any], analysis_metadata: dict[str,Any]) -> tuple[np.ndarray,np.ndarray,dict[str,Any]]:
    """Display derived states once when they are already in integration coordinates."""
    preprocessing=analysis_metadata.get("preprocessing") or {}
    if preprocessing.get("coordinate_system")!="dimensionless": return _display_series(series,source_metadata)
    context=_axis_context(source_metadata)
    return series.time*float(context["time_scale"]),series.values,context
