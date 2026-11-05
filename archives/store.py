"""Atomic per-system solution archives and optional live history checkpoints."""
from __future__ import annotations
import hashlib, io, json, os, time, uuid, zipfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import numpy as np
from provenance import serializable_system
from time_series import TimeSeries

# This package owns numerical data independently of both solver implementation
# and dynamical-system definitions. Override the location for any installed
# environment with DYNSYSLAB_ARCHIVE_DIR. Editable/source installs retain the
# project-local archive tree; wheel installs default to a user-writable folder.
_package_root = Path(__file__).resolve().parent
_source_root = _package_root.parent
_default_root = _package_root if (_source_root / "pyproject.toml").exists() else Path.home() / "dynsyslab_archives"
ARCHIVE_ROOT = Path(os.environ.get("DYNSYSLAB_ARCHIVE_DIR", _default_root)).expanduser().resolve()

def _stamp(): return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
def _safe(value: str): return "".join(c if c.isalnum() or c in "-_" else "-" for c in value)
def _metadata_json(metadata: dict[str,Any]) -> str: return json.dumps(metadata,default=str,separators=(",",":"),sort_keys=True)
def _uniform_time_steps(time: np.ndarray | None) -> bool:
    """Return whether a saved trajectory uses a single constant time increment."""
    if time is None:
        return False
    samples = np.asarray(time, dtype=float).reshape(-1)
    if samples.size < 2:
        return True
    increments = np.diff(samples)
    return bool(np.allclose(increments, increments[0], rtol=1e-12, atol=1e-15))
@contextmanager
def _path_lock(path: Path, *, timeout_seconds: float=60.0, stale_after_seconds: float=24*3600):
    """Small cross-process lock for a single archive record.

    Readers never need the lock because writers publish with ``os.replace``.
    The lock only serializes read-modify-write operations, avoiding lost
    metadata updates when the solver, JobManager, and postprocessors finish at
    nearly the same time.
    """
    lock_path=path.with_suffix(path.suffix+".lock")
    deadline=time.monotonic()+timeout_seconds
    handle=None
    while handle is None:
        try:
            handle=lock_path.open("x",encoding="utf-8")
            handle.write(f"pid={os.getpid()}\nutc={datetime.now(timezone.utc).isoformat()}\n"); handle.flush(); os.fsync(handle.fileno())
        except FileExistsError:
            # The lock is intentionally a tiny sidecar, not a process-bound
            # descriptor, so an OS crash can leave one behind. A writer holds
            # it only for a read-modify-write transaction; a day-old sidecar
            # is therefore safe to recover and should not permanently block
            # an otherwise healthy campaign.
            try:
                if time.time()-lock_path.stat().st_mtime>stale_after_seconds:
                    lock_path.unlink(); continue
            except FileNotFoundError:
                continue
            if time.monotonic()>=deadline:
                raise TimeoutError(f"Timed out waiting for archive lock: {lock_path}")
            time.sleep(0.05)
    try:
        yield
    finally:
        handle.close()
        try: lock_path.unlink()
        except FileNotFoundError: pass

def _atomic_npz(path: Path, **data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary = path.with_suffix(path.suffix+f".{uuid.uuid4().hex}.partial")
    # LZMA is lossless and normally denser than numpy's deflate default for
    # repeated smooth trajectory data. NPZ remains directly readable by
    # numpy.load; no custom reader is required.
    codec={"lzma":zipfile.ZIP_LZMA,"deflate":zipfile.ZIP_DEFLATED}.get(os.environ.get("DYNSYSLAB_COMPRESSION","lzma").lower())
    if codec is None: raise ValueError("DYNSYSLAB_COMPRESSION must be 'lzma' or 'deflate'.")
    with zipfile.ZipFile(temporary,"w",compression=codec) as archive:
        for name,value in data.items():
            buffer=io.BytesIO(); np.save(buffer,np.asarray(value),allow_pickle=False)
            archive.writestr(f"{name}.npy",buffer.getvalue())
    # Make the file durable before publishing it as a completed archive.
    # Windows requires a writable descriptor for fsync, unlike many POSIX
    # filesystems where a read-only descriptor is accepted.
    with temporary.open("rb+") as handle: os.fsync(handle.fileno())
    os.replace(temporary, path)

def recompress_archive(path: str | Path) -> Path:
    """Losslessly rewrite one dynsyslab archive with the standard LZMA codec."""
    target=Path(path)
    with _path_lock(target):
        with np.load(target,allow_pickle=False) as source: payload={name:source[name] for name in source.files}
        _atomic_npz(target,**payload)
    return target

def update_archive_metadata(path: str | Path, patch: dict[str,Any]) -> Path:
    """Atomically persist post-solve provenance without changing numeric arrays."""
    target=Path(path)
    with _path_lock(target):
        with np.load(target,allow_pickle=False) as source: payload={name:source[name] for name in source.files}
        metadata=json.loads(str(payload["metadata_json"])); metadata.update(patch); payload["metadata_json"]=_metadata_json(metadata)
        _atomic_npz(target,**payload)
    return target

def rewrite_archive_metadata(path: str | Path, transform: Any) -> Path:
    """Atomically replace metadata while preserving every numerical array."""
    target=Path(path)
    with _path_lock(target):
        with np.load(target,allow_pickle=False) as source: payload={name:source[name] for name in source.files}
        metadata=json.loads(str(payload["metadata_json"])); updated=transform(metadata)
        if not isinstance(updated,dict): raise TypeError("metadata transform must return a dictionary.")
        payload["metadata_json"]=_metadata_json(updated); _atomic_npz(target,**payload)
    return target

def _time_grid_key(time: np.ndarray) -> str:
    samples=np.ascontiguousarray(np.asarray(time,dtype=float).reshape(-1))
    digest=hashlib.sha256(); digest.update(str(samples.shape).encode()); digest.update(samples.dtype.str.encode()); digest.update(samples.tobytes())
    return digest.hexdigest()

def _time_grid_path(reference: str, metadata: dict[str,Any] | None=None) -> Path:
    """Resolve a deduplicated grid relative to the solution's own archive root.

    Older records did not carry ``time_grid_root`` and retain the primary-root
    fallback, so compacting an old archive remains backwards compatible.
    """
    root=Path((metadata or {}).get("time_grid_root") or (metadata or {}).get("archive_root") or ARCHIVE_ROOT).expanduser().resolve()
    target=(root/reference).resolve()
    if not target.is_file(): raise FileNotFoundError(f"Referenced archive time grid is unavailable: {target}")
    return target

def _time_from_archive(data: Any, metadata: dict[str,Any], *, raw: bool=False) -> np.ndarray:
    direct="raw_time" if raw else "time"
    if direct in data.files: return np.asarray(data[direct])
    reference=metadata.get("raw_time_grid_ref" if raw else "time_grid_ref")
    if reference is None: raise KeyError(f"Archive has neither {direct!r} nor a time-grid reference.")
    with np.load(_time_grid_path(reference,metadata),allow_pickle=False) as grid: return np.asarray(grid["time"])

def load_time_series(path: str | Path) -> TimeSeries:
    """Load the primary archived solution as the shared TimeSeries object."""
    with np.load(Path(path),allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        time=_time_from_archive(data,metadata) if "time" not in data.files else data["time"]
        return TimeSeries(time,data["values"] if "values" in data.files else data["state"],int(metadata.get("pooled_members",0)))

def _extend_solution_unlocked(path: str | Path, continuation: TimeSeries, *, raw_continuation: TimeSeries | None=None,
                    control_continuation: TimeSeries | None=None, control_target_continuation: TimeSeries | None=None,
                    forcing_continuation: TimeSeries | None=None) -> Path:
    """Append a continuation directly to its original solution archive.

    This is the sole archive-level extension operation: it never creates a
    second solution record.  Aligned control/forcing/raw histories must be
    supplied when the original archive contains them, preventing silent
    metadata/data misalignment.
    """
    target=Path(path).resolve()
    with np.load(target,allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"])); time=_time_from_archive(data,metadata); values=np.asarray(data["values"] if "values" in data.files else data["state"])
        payload={name:np.asarray(data[name]) for name in data.files if name not in {"metadata_json","time","values","state"}}
    def join(old_time,old_values,new: TimeSeries, name: str):
        new_time=np.asarray(new.time); new_values=np.asarray(new.values)
        if old_values.shape[1:]!=new_values.shape[1:]: raise ValueError(f"{name} state shape does not match the original solution.")
        if new_time.size==0 or new_time[0]<old_time[-1]: raise ValueError(f"{name} must start at or after the original final time.")
        start=1 if np.isclose(new_time[0],old_time[-1],rtol=1e-12,atol=1e-15) else 0
        return np.concatenate((old_time,new_time[start:])),np.concatenate((old_values,new_values[start:]),axis=0)
    combined_time,combined_values=join(time,values,continuation,"continuation")
    for key, supplied in (("control_values",control_continuation),("control_target_values",control_target_continuation),("forcing_values",forcing_continuation)):
        if key in payload:
            if supplied is None: raise ValueError(f"Original archive contains {key}; provide its aligned continuation.")
            _,payload[key]=join(time,payload[key],supplied,key)
    if "raw_values" in payload:
        if raw_continuation is None: raise ValueError("Original archive contains raw data; provide raw_continuation when extending it.")
        raw_time=np.asarray(payload.pop("raw_time")); raw_values=np.asarray(payload.pop("raw_values")); payload["raw_time"],payload["raw_values"]=join(raw_time,raw_values,raw_continuation,"raw continuation")
    metadata.pop("time_grid_ref",None); metadata.pop("time_grid_root",None); metadata["time_grid_shared"]=False; metadata["uniform_time_steps"]=TimeSeries(combined_time,combined_values).uniform_time_steps; metadata["delta_t"]=TimeSeries(combined_time,combined_values).delta_t
    metadata.setdefault("extensions",[]).append({"extended_utc":_stamp(),"prior_samples":int(time.size),"appended_samples":int(combined_time.size-time.size),"final_samples":int(combined_time.size)})
    metadata["time_span"]=[float(combined_time[0]),float(combined_time[-1])]; metadata["metadata_revision"]=int(metadata.get("metadata_revision",0))+1
    _atomic_npz(target,time=combined_time,values=combined_values,metadata_json=_metadata_json(metadata),**payload)
    from .index import index_solution
    index_solution(target)
    from infrastructure.source_events import notify_source_update
    notify_source_update(target,source_type="solution",reason="trajectory_extended")
    return target

def extend_solution(path: str | Path, continuation: TimeSeries, *, raw_continuation: TimeSeries | None=None,
                    control_continuation: TimeSeries | None=None, control_target_continuation: TimeSeries | None=None,
                    forcing_continuation: TimeSeries | None=None) -> Path:
    """Append under an exclusive lock so no concurrent metadata update is lost."""
    target=Path(path).resolve()
    with _path_lock(target):
        return _extend_solution_unlocked(target,continuation,raw_continuation=raw_continuation,
                                         control_continuation=control_continuation,
                                         control_target_continuation=control_target_continuation,
                                         forcing_continuation=forcing_continuation)

def extend_solution_from_archive(path: str | Path, continuation_archive: str | Path) -> Path:
    """Append a completed continuation archive, preserving forcing/control history."""
    source=Path(continuation_archive)
    with np.load(source,allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"])); time=_time_from_archive(data,metadata); values=np.asarray(data["values"] if "values" in data.files else data["state"])
        continuation=TimeSeries(time,values)
        raw=TimeSeries(np.asarray(data["raw_time"]),np.asarray(data["raw_values"])) if {"raw_time","raw_values"}.issubset(data.files) else None
        control=TimeSeries(time,np.asarray(data["control_values"])) if "control_values" in data.files else None
        control_target=TimeSeries(time,np.asarray(data["control_target_values"])) if "control_target_values" in data.files else None
        forcing=TimeSeries(time,np.asarray(data["forcing_values"])) if "forcing_values" in data.files else None
    return extend_solution(path,continuation,raw_continuation=raw,control_continuation=control,control_target_continuation=control_target,forcing_continuation=forcing)

def load_control_time_series(path: str | Path, *, raw: bool=False) -> TimeSeries:
    """Load solve-time controller values aligned with the archived solution."""
    with np.load(Path(path),allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        value_key="raw_control_values" if raw else "control_values"
        if value_key not in data.files: raise KeyError("This archive does not contain a controller history.")
        return TimeSeries(_time_from_archive(data,metadata,raw=raw),data[value_key])

def load_control_target_time_series(path: str | Path, *, raw: bool=False) -> TimeSeries:
    """Load the controller reference aligned with the archived solution."""
    with np.load(Path(path),allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        value_key="raw_control_target_values" if raw else "control_target_values"
        if value_key not in data.files: raise KeyError("This archive does not contain a controller-target history.")
        return TimeSeries(_time_from_archive(data,metadata,raw=raw),data[value_key])

def load_hierarchy_time_series(path: str | Path, *, system_index: int=0, cluster_index: int | None=None) -> TimeSeries:
    """Reconstruct a pooled System or Cluster TimeSeries from canonical children."""
    with np.load(Path(path),allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        # Legacy hierarchy files may contain redundant pooled arrays. Prefer
        # canonical child archives whenever their references are present.
        try: system=metadata["systems"][system_index]
        except IndexError as error: raise IndexError(f"System index {system_index} is not present in this hierarchy.") from error
        clusters=system["clusters"] if cluster_index is None else [system["clusters"][cluster_index]]
        paths=[trajectory.get("solution_archive") for cluster in clusters for trajectory in cluster["trajectories"]]
        if not paths or any(item is None for item in paths):
            raise KeyError("This hierarchy has no complete child solution references to reconstruct its pooled TimeSeries.")
    children=[load_time_series(item) for item in paths]
    reference=children[0].time
    if any(item.time.shape!=reference.shape or not np.allclose(item.time,reference,rtol=1e-12,atol=1e-15) for item in children[1:]):
        raise ValueError("Cannot reconstruct pooled hierarchy data: child trajectory time grids differ.")
    shapes={item.values.shape[1:] for item in children}
    if len(shapes)!=1: raise ValueError("Cannot reconstruct pooled hierarchy data: child solution shapes differ.")
    return TimeSeries(reference,np.stack([item.values for item in children],axis=1),pooled_members=len(children))

def save_hierarchy(target: Any) -> Path:
    """Persist a completed Batch or System hierarchy and its pooled arrays."""
    from hierarchy import hierarchy_payload
    metadata,arrays=hierarchy_payload(target)
    sources=[trajectory.get("solution_archive") for system in metadata["systems"] for cluster in system["clusters"] for trajectory in cluster["trajectories"]]
    deduplicate_solution_time_grids([Path(item) for item in sources if item is not None])
    systems=target.systems if hasattr(target,"systems") else [target]
    configured_roots={Path(getattr(system,"archive_root",None) or ARCHIVE_ROOT).expanduser().resolve() for system in systems}
    # A single-System hierarchy, or a Batch whose Systems share storage, stays
    # inside that registered archive root.  A cross-root Batch manifest falls
    # back to the primary archive because it references multiple stores.
    archive_root=configured_roots.pop() if len(configured_roots)==1 else Path(ARCHIVE_ROOT).resolve()
    root=archive_root/"hierarchies"; root.mkdir(parents=True,exist_ok=True)
    label=_safe(str(metadata["label"])); path=root/f"{label}__{_stamp()}.npz"
    _atomic_npz(path,metadata_json=_metadata_json(metadata),**arrays)
    return path

def compact_hierarchy(path: str | Path) -> Path:
    """Losslessly remove legacy duplicated pooled arrays from a hierarchy file."""
    target=Path(path)
    with np.load(target,allow_pickle=False) as data:
        if "metadata_json" not in data.files: raise ValueError("Not a dynsyslab hierarchy archive.")
        metadata=json.loads(str(data["metadata_json"]))
    if metadata.get("format")!="dynsyslab.hierarchy/v1": raise ValueError("Not a dynsyslab hierarchy archive.")
    # Verify references before removing legacy duplicated arrays.
    for system in metadata.get("systems",[]):
        for cluster in system.get("clusters",[]):
            for trajectory in cluster.get("trajectories",[]):
                source=trajectory.get("solution_archive")
                if source is None or not Path(source).is_file():
                    raise FileNotFoundError("Cannot compact hierarchy: a canonical child solution archive is missing.")
    _atomic_npz(target,metadata_json=_metadata_json(metadata))
    return target

def deduplicate_solution_time_grids(paths: list[Path]) -> int:
    """Losslessly share only repeated *final* time meshes among solutions.

    Raw adaptive meshes stay with their own solve when interpolation occurred;
    they are distinct numerical evidence. A shared grid is written only for a
    group of two or more byte-identical final meshes, so singleton archives do
    not gain an extra sidecar file.
    """
    groups: dict[tuple[str,str],list[tuple[Path,np.ndarray,dict[str,Any]]]]={}
    for path in paths:
        if not path.is_file(): continue
        with np.load(path,allow_pickle=False) as data:
            metadata=json.loads(str(data["metadata_json"]))
            if "time" not in data.files: continue
            time=np.asarray(data["time"])
        root=Path(metadata.get("archive_root") or ARCHIVE_ROOT).expanduser().resolve()
        groups.setdefault((str(root),_time_grid_key(time)),[]).append((path,time,metadata))
    rewritten=0
    for (root_text,key),items in groups.items():
        if len(items)<2: continue
        archive_root=Path(root_text)
        root=archive_root/"time_grids"; root.mkdir(parents=True,exist_ok=True)
        grid=root/f"{key}.npz"; reference=str(grid.relative_to(archive_root)).replace("\\","/")
        if not grid.exists(): _atomic_npz(grid,time=items[0][1])
        for path,_,metadata in items:
            metadata["time_grid_ref"]=reference
            metadata["time_grid_root"]=str(archive_root)
            metadata["time_grid_shared"]=True
            with _path_lock(path):
                # Reload under the lock: a solver may have just extended this
                # trajectory after the grouping pass, in which case it is no
                # longer eligible for this shared grid.
                with np.load(path,allow_pickle=False) as source:
                    current={name:source[name] for name in source.files}
                if "time" not in current or _time_grid_key(current["time"])!=key:
                    continue
                current_metadata=json.loads(str(current["metadata_json"])); current_metadata.update({"time_grid_ref":reference,"time_grid_root":str(archive_root),"time_grid_shared":True})
                payload={name:value for name,value in current.items() if name not in {"time","metadata_json"}}
                payload["metadata_json"]=_metadata_json(current_metadata)
                _atomic_npz(path,**payload); rewritten+=1
    return rewritten

def _resample_series(time: np.ndarray, values: np.ndarray, target: np.ndarray, method: str) -> np.ndarray:
    from scipy.interpolate import PchipInterpolator, interp1d

    unique,indices=np.unique(np.asarray(time,dtype=float),return_index=True)
    source=np.asarray(values)[indices]
    if method=="pchip": return np.asarray(PchipInterpolator(unique,source,axis=0)(target))
    if method in {"linear","nearest"}: return np.asarray(interp1d(unique,source,axis=0,kind=method)(target))
    raise ValueError(f"Unsupported forcing interpolation method: {method}.")

class ArchiveSession:
    def __init__(self, request: Any):
        self.request=request; self.system_name=_safe(getattr(request.system,"name",type(request.system).__name__))
        self.archive_root=Path(getattr(request,"archive_root",None) or ARCHIVE_ROOT).resolve()
        self.root=self.archive_root/"solutions"/self.system_name; self.root.mkdir(parents=True,exist_ok=True)
        self.run_id=f"{_safe(request.label)}__{_stamp()}"; self.history_dir=self.archive_root/"history"/self.system_name/self.run_id
        self.history_count=0
    def telemetry(self,event: str, **details: Any) -> None:
        """Append tiny lifecycle telemetry; never serializes state arrays."""
        if not getattr(self.request.history,"telemetry_enabled",True): return
        directory=self.archive_root/"solver_logs"; directory.mkdir(parents=True,exist_ok=True)
        path=directory/f"{self.run_id}.telemetry.jsonl"
        with path.open("a",encoding="utf-8") as handle:
            handle.write(_metadata_json({"time_utc":datetime.now(timezone.utc).isoformat(),"event":event,"label":self.request.label,"kind":self.request.kind,"job_id":getattr(self.request,"job_id",None),**details})+"\n")
    def metadata(self, status: str, message: str=""):
        metadata={"system":self.system_name,"kind":self.request.kind,"label":self.request.label,"status":status,"message":message,
                "options":self.request.options.__dict__,"history":self.request.history.__dict__,"time_span":self.request.time_span,
                "initial_state_override":None if self.request.initial_state is None else np.asarray(self.request.initial_state,dtype=float).tolist(),
                "interpolation_request":self.request.interpolation.__dict__,"spatial_domain":self.request.spatial_domain,"mesh":str(self.request.mesh),"postprocessors_requested":len(self.request.postprocessors),"system_parameters":serializable_system(self.request.system),
                "controller":None if self.request.controller is None else serializable_system(self.request.controller),
                "estimator":None if self.request.estimator is None else serializable_system(self.request.estimator),"archive_root":str(self.archive_root)}
        metadata.update(getattr(self.request,"metadata",{}) or {})
        return metadata
    def checkpoint(self, time: np.ndarray, values: np.ndarray, note: str="running"):
        if not self.request.history.track_during_solve and self.request.history.mode not in {"sampled","checkpoint"}: return None
        self.history_dir.mkdir(parents=True,exist_ok=True); self.history_count+=1
        path=self.history_dir/f"checkpoint-{self.history_count:06d}.npz"
        # Checkpoint metadata belongs once in the run manifest, not copied into
        # every incrementally written numerical checkpoint.
        _atomic_npz(path,time=np.asarray(time),values=np.asarray(values))
        (self.history_dir/"manifest.json").write_text(_metadata_json(self.metadata(note)),encoding="utf-8")
        return path
    def _forcing(self, time: np.ndarray, values: np.ndarray):
        """Sample built-in drivers and external additive controls on solution points."""
        system,control=self.request.system,self.request.control; columns=[]; names=[]; details=[]
        driver=getattr(system,"driver",None)
        driver_active=False
        if driver is not None and callable(getattr(driver,"value",None)):
            driver_values=np.asarray([driver.value(float(t)) for t in time],dtype=float).reshape(-1,1)
            # A zero-amplitude built-in driver is recorded as a system parameter,
            # but is correctly classified as free evolution.
            if np.any(np.abs(driver_values)>np.finfo(float).eps):
                # Driver samples are exactly regenerable from its compact
                # parameter metadata, unlike arbitrary external controls.
                driver_active=True; details.append({"source":"system_driver","parameters":serializable_system(driver),"samples_archived":False})
        if control is not None:
            samples=[]
            for t,state in zip(time,values):
                sample=control.value(float(t),state) if hasattr(control,"value") else control(float(t),state)
                samples.append(np.asarray(sample,dtype=float).reshape(-1))
            external=np.vstack(samples); columns.append(external)
            names.extend([f"external_control_{index}" for index in range(external.shape[1])])
            details.append({"source":"external_control","description":repr(control),"columns":int(external.shape[1])})
        # Controllers may have internal integrator/observer state.  Their
        # actual solve-time values are saved separately by ``complete``;
        # never reevaluate a stateful controller while archiving.
        controller=self.request.controller
        if controller is not None:
            details.append({"source":"controller","parameters":serializable_system(controller),"samples_archived":"control_values"})
        forced=bool(columns) or controller is not None or driver_active
        info={"forcing_type":"forced" if forced else "free","forcing_present":forced,"sources":details,"columns":names,
              "time_aligned_with_solution":True}
        return info,(np.hstack(columns) if columns else None)
    def complete(self, result: Any):
        path=self.root/f"{self.run_id}.npz"; metadata=self.metadata(result.status,result.message)|result.metadata
        if not isinstance(result.time_series,TimeSeries):
            raise TypeError("ArchiveSession.complete requires NumericalResult.time_series to be a TimeSeries.")
        series=result.time_series; raw=result.raw_time_series
        # Binary, archive-level property of the actual time array written below.
        metadata["uniform_time_steps"]=series.uniform_time_steps
        metadata["delta_t"]=series.delta_t
        metadata["time_series_format"]="dynsyslab.TimeSeries/v1"
        metadata["pooled_members"]=series.pooled_members
        metadata["is_pooled"]=series.is_pooled
        metadata.setdefault("interpolated",False)
        if raw is not None:
            metadata["raw_uniform_time_steps"]=raw.uniform_time_steps; metadata["raw_delta_t"]=raw.delta_t
        payload={"values":series.values,"time":series.time,"metadata_json":_metadata_json(metadata)}
        if raw is not None:
            payload["raw_time"]=raw.time; payload["raw_values"]=raw.values
        control=result.control_time_series; raw_control=result.raw_control_time_series
        if control is not None:
            if control.time.shape!=series.time.shape or not np.allclose(control.time,series.time,rtol=1e-12,atol=1e-15):
                raise ValueError("Control history must align with the archived solution time grid.")
            payload["control_values"]=control.values
            metadata["control_history"]={"present":True,"time_aligned_with_solution":True,"columns":int(control.values.shape[1]),"coordinate_system":result.metadata.get("control_input_coordinate_system"),"units":result.metadata.get("control_input_units"),"scale_to_dimensional":result.metadata.get("control_input_scale_to_dimensional"),"uniform_time_steps":control.uniform_time_steps,"delta_t":control.delta_t}
            if raw_control is not None:
                if raw is None or raw_control.time.shape!=raw.time.shape or not np.allclose(raw_control.time,raw.time,rtol=1e-12,atol=1e-15):
                    raise ValueError("Raw control history must align with the raw archived solution time grid.")
                payload["raw_control_values"]=raw_control.values
                metadata["control_history"]["raw_samples"]=raw_control.sample_count
        else: metadata["control_history"]={"present":False}
        if result.control_call_time is not None or result.control_call_values is not None:
            call_time=np.asarray(result.control_call_time,dtype=float).reshape(-1); call_values=np.asarray(result.control_call_values)
            if call_values.ndim!=2 or call_values.shape[0]!=call_time.size:
                raise ValueError("Controller call history time and values must have matching sample counts.")
            payload["control_call_time"]=call_time; payload["control_call_values"]=call_values
            metadata["control_history"]["controller_calls_recorded"]=int(call_time.size)
            metadata["control_history"]["evaluation_policy"]=result.metadata.get("control_evaluation")
            metadata["control_history"]["call_log_semantics"]="actual_controller_evaluations"
            metadata["control_history"]["call_log_contains_repeated_times"]=bool(np.unique(call_time).size!=call_time.size)
        control_target=result.control_target_time_series; raw_control_target=result.raw_control_target_time_series
        if control_target is not None:
            if control_target.time.shape!=series.time.shape or not np.allclose(control_target.time,series.time,rtol=1e-12,atol=1e-15):
                raise ValueError("Control-target history must align with the archived solution time grid.")
            payload["control_target_values"]=control_target.values
            metadata["control_target_history"]={"present":True,"time_aligned_with_solution":True,"columns":int(control_target.values.shape[1]),"coordinate_system":"canonical_archive_state","uniform_time_steps":control_target.uniform_time_steps,"delta_t":control_target.delta_t}
            if raw_control_target is not None:
                if raw is None or raw_control_target.time.shape!=raw.time.shape or not np.allclose(raw_control_target.time,raw.time,rtol=1e-12,atol=1e-15):
                    raise ValueError("Raw control-target history must align with the raw archived solution time grid.")
                payload["raw_control_target_values"]=raw_control_target.values
                metadata["control_target_history"]["raw_samples"]=raw_control_target.sample_count
        else: metadata["control_target_history"]={"present":False}
        if series is not None:
            force_info,forcing=self._forcing(raw.time if raw is not None else series.time,raw.values if raw is not None else series.values)
            metadata["forcing"]=force_info; metadata["forcing_type"]=force_info["forcing_type"]
            if forcing is not None:
                if raw is not None:
                    payload["raw_forcing_values"]=forcing
                    interpolation=result.metadata.get("interpolation",{}); method=interpolation.get("method","pchip")
                    payload["forcing_values"]=_resample_series(raw.time,forcing,series.time,method)
                    metadata["forcing"]["interpolation"]={"method":method,"raw_samples":int(forcing.shape[0]),"output_samples":series.sample_count}
                else: payload["forcing_values"]=forcing
        payload["metadata_json"]=_metadata_json(metadata)
        _atomic_npz(path,**payload)
        from .index import index_solution
        index_solution(path)
        return path
