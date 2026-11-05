"""Read-only source loading and display transforms for :mod:`viewer`."""
from __future__ import annotations
from pathlib import Path
from typing import Any
import json, numpy as np
from archiver.gateway import VIEWER_ARCHIVE
from time_series import TimeSeries

def manifest(path:Path)->dict[str,Any]:
    with VIEWER_ARCHIVE.open_read(path) as data:
        if "metadata_json" not in data.files: raise ValueError(f"{path} is not a dynsyslab archive.")
        return json.loads(str(data["metadata_json"]))
def solution_record(path:Path)->tuple[TimeSeries,dict[str,Any]]:
    with VIEWER_ARCHIVE.open_read(path) as data: metadata=json.loads(str(data["metadata_json"]))
    metadata["source_archive"]=str(path); return VIEWER_ARCHIVE.read_time_series(path),metadata
def control_target_series(metadata:dict[str,Any])->TimeSeries|None:
    """Read an archived target only when the solve recorded one."""
    history=metadata.get("control_target_history") or {}
    if not history.get("present") or metadata.get("controller") is None: return None
    source=metadata.get("source_archive")
    if source is None: return None
    return VIEWER_ARCHIVE.read_control_target_time_series(source)
def control_series(metadata:dict[str,Any])->TimeSeries|None:
    """Read the actual archived controller output on the solution grid."""
    history=metadata.get("control_history") or {}
    if not history.get("present") or metadata.get("controller") is None: return None
    source=metadata.get("source_archive")
    if source is None: return None
    return VIEWER_ARCHIVE.read_control_time_series(source)
def trajectory_paths(manifest_data:dict[str,Any],system_index:int,cluster_index:int|None=None,trajectory_index:int|None=None)->list[Path]:
    systems=manifest_data.get("systems")
    if systems is None: raise ValueError("The requested level requires a tiered hierarchy manifest.")
    system=systems[system_index]; clusters=system["clusters"] if cluster_index is None else [system["clusters"][cluster_index]]; paths=[]
    for cluster in clusters:
        records=cluster["trajectories"] if trajectory_index is None else [cluster["trajectories"][trajectory_index]]; paths.extend(Path(item["solution_archive"]) for item in records)
    if not paths or any(not item.is_file() for item in paths): raise FileNotFoundError("A referenced child trajectory archive is unavailable.")
    return paths
def load_archive_source(source:Any)->list[tuple[str,TimeSeries,dict[str,Any]]]:
    path=Path(source.path).expanduser().resolve()
    if not path.is_file(): raise FileNotFoundError(path)
    try: manifest_data=manifest(path)
    except ValueError:
        if source.level!="trajectory": raise
        return [(source.label or path.stem,*solution_record(path))]
    if manifest_data.get("format")!="dynsyslab.hierarchy/v1":
        if source.level!="trajectory": raise ValueError("Non-hierarchy archives can only be viewed at trajectory level.")
        return [(source.label or path.stem,*solution_record(path))]
    selected=range(len(manifest_data["systems"])) if source.level=="batch" else [source.system_index or 0]; output=[]
    for index in selected:
        paths=trajectory_paths(manifest_data,index,source.cluster_index,source.trajectory_index) if source.level=="trajectory" else trajectory_paths(manifest_data,index,source.cluster_index if source.level=="cluster" else None)
        for item in paths:
            series,metadata=solution_record(item); metadata["system_index"]=index; output.append((f"system-{index}: {item.stem}",series,metadata))
    return output
def axis_context(metadata:dict[str,Any])->dict[str,Any]:
    system=metadata.get("system",""); form=metadata.get("integration_form","dimensional"); parameters=metadata.get("system_parameters",{})
    if system=="rayleigh-plesset":
        if form=="dimensionless":
            preset=parameters.get("preset",{}); radius=float(parameters["equilibrium_radius_m"])
            scaling=metadata.get("nondimensionalization") or {}
            if scaling.get("time_scale_s") is not None: scale_time=float(scaling["time_scale_s"])
            else:
                # Legacy dimensionless archives used the inertial-pressure
                # time scale and predate explicit scaling provenance.
                pressure=float(parameters["ambient_pressure_pa"]); density=float(preset["rho_kg_m3"]); scale_time=radius*np.sqrt(density/pressure)
            return {"time_scale":1/scale_time,"state_scale":np.array((1/radius,scale_time/radius)),"time_label":"τ","state_labels":["r = R/R₀","dr/dτ"],"frequency_label":"frequency (cycles / τ)"}
        return {"time_scale":1e6,"state_scale":np.array((1e6,1.)),"time_label":"t (µs)","state_labels":["R (µm)","dR/dt (m/s)"],"frequency_label":"frequency (MHz)"}
    if system=="marmottant":
        if form=="dimensionless": return {"time_scale":1.,"state_scale":None,"time_label":"τ","state_labels":["r = R/R₀","dr/dτ"],"frequency_label":"frequency (cycles / τ)"}
        return {"time_scale":1e6,"state_scale":np.array((1e6,1.)),"time_label":"t (µs)","state_labels":["R (µm)","dR/dt (m/s)"],"frequency_label":"frequency (MHz)"}
    if system=="reddy-szeri":
        return {"time_scale":1.,"state_scale":None,"time_label":"τ","state_labels":["r = R/R₀","dr/dτ","x/R₀","dx/dτ","a₂","da₂/dτ","a₃","da₃/dτ","a₄","da₄/dτ","a₅","da₅/dτ"],"frequency_label":"frequency (cycles / τ)"}
    if system=="duffing":
        return {"time_scale":1.,"state_scale":None,"time_label":"t","state_labels":["x","dx/dt"],"frequency_label":"frequency"}
    return {"time_scale":1.,"state_scale":None,"time_label":"t","state_labels":[],"frequency_label":"frequency"}
def display_series(series:TimeSeries,metadata:dict[str,Any]):
    context=axis_context(metadata); scale=context["state_scale"]; return series.time*context["time_scale"],series.values if scale is None else series.values*np.asarray(scale),context

def derived_array(path: str|Path,name: str,*,columns=None) -> np.ndarray:
    """Resolve a stored, aliased, or virtual derived array without writes."""
    from postprocessing.core import load_derived_array
    return load_derived_array(path,name,columns=columns)
