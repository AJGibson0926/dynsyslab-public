"""Archive-backed trajectory and derived-analysis data contracts."""
from __future__ import annotations
import hashlib, io, json, os, uuid, zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from archiver.gateway import POSTPROCESSING_ARCHIVE, VIEWER_ARCHIVE
from time_series import TimeSeries

DERIVED_DEFAULT_COMPRESSION="deflate"
_COMPRESSION_TYPES={"deflate":zipfile.ZIP_DEFLATED,"lzma":zipfile.ZIP_LZMA,"stored":zipfile.ZIP_STORED}

def _archive_compression(path: Path) -> str:
    with zipfile.ZipFile(path,"r") as archive:
        members=archive.infolist()
        if not members: raise ValueError(f"Derived archive has no members: {path}")
        codec=members[0].compress_type
    return next((name for name,value in _COMPRESSION_TYPES.items() if value==codec),"stored")

def _write_derived_npz(path: Path, payload: dict[str,object], *, compression: str=DERIVED_DEFAULT_COMPRESSION) -> Path:
    if compression not in _COMPRESSION_TYPES: raise ValueError("compression must be 'deflate', 'lzma', or 'stored'.")
    temporary=path.with_suffix(path.suffix+f".{uuid.uuid4().hex}.partial")
    options={"compression":_COMPRESSION_TYPES[compression]}
    if compression=="deflate": options["compresslevel"]=1
    with zipfile.ZipFile(temporary,"w",**options) as archive:
        for name,value in payload.items():
            buffer=io.BytesIO(); np.save(buffer,np.asarray(value),allow_pickle=False); archive.writestr(f"{name}.npy",buffer.getvalue())
    with temporary.open("rb+") as handle: os.fsync(handle.fileno())
    os.replace(temporary,path); return path

@dataclass(frozen=True)
class Trajectory:
    time_series: TimeSeries
    metadata: dict
    source: Path | None

    @property
    def time(self): return self.time_series.time
    @property
    def values(self): return self.time_series.values

@dataclass(frozen=True)
class DerivedRecord:
    path: Path
    analysis: str
    source_archive: str
    created_utc: str | None
    publication_id: str | None
    analysis_signature_sha256: str | None

def load_solution(path: str | Path) -> Trajectory:
    """Read a completed archive without importing solver or system-library code."""
    path=Path(path).resolve()
    with VIEWER_ARCHIVE.open_read(path) as data:
        metadata=json.loads(str(data["metadata_json"]))
    return Trajectory(VIEWER_ARCHIVE.read_time_series(path),metadata,path)

def integration_coordinates(source: Trajectory) -> Trajectory:
    """Recover the coordinate system used by a dimensionless solver run.

    Solver archives intentionally store canonical dimensional time/state values.
    Analyses that must reproduce the generating equations can opt back into the
    archived run's integration coordinates without changing the source archive.
    """
    form=str(source.metadata.get("integration_form","dimensional"))
    metadata=dict(source.metadata)
    preprocessing=dict(metadata.get("postprocessing_preprocessing") or {})
    if form!="dimensionless":
        preprocessing["coordinate_system"]="dimensional"
        preprocessing["coordinate_transform"]="identity"
        metadata["postprocessing_preprocessing"]=preprocessing
        return Trajectory(source.time_series,metadata,source.source)
    try:
        from dynsyslib import resolve_archived_system
        system=resolve_archived_system(metadata.get("system"),metadata.get("system_parameters"))
    except Exception as error:
        raise ValueError("Dimensionless integration coordinates require a resolvable archived dynamical system.") from error
    if not hasattr(system,"to_dimensionless_time") or not hasattr(system,"to_dimensionless_state"):
        raise ValueError(f"Archived system {metadata.get('system')!r} does not expose its dimensionless coordinate transform.")
    time=np.asarray(system.to_dimensionless_time(source.time),dtype=float)
    values=np.asarray([system.to_dimensionless_state(row) for row in source.values],dtype=float)
    preprocessing.update({"coordinate_system":"dimensionless","coordinate_transform":"archive-dimensional-to-integration",
                          "source_integration_form":form})
    metadata["postprocessing_preprocessing"]=preprocessing
    metadata["delta_t"]=float(time[1]-time[0])
    return Trajectory(TimeSeries(time,values),metadata,source.source)

def _column_selection(values: np.ndarray, columns) -> np.ndarray:
    if columns is None: return np.asarray(values)
    indices=np.asarray(columns,dtype=int).reshape(-1)
    if values.ndim!=2: raise ValueError("Column selection requires a two-dimensional derived array.")
    if np.any(indices<0) or np.any(indices>=values.shape[1]): raise IndexError("Derived-array column index is out of range.")
    return np.asarray(values[:,indices])

def _virtual_derived_array(data, specification: dict, columns) -> np.ndarray:
    if specification.get("encoding")!="polynomial-dictionary-evaluation/v1":
        raise ValueError(f"Unsupported virtual derived-array encoding {specification.get('encoding')!r}.")
    dictionary=specification.get("dictionary") or {}
    if dictionary.get("encoding")!="polynomial-dictionary/v1":
        raise ValueError("Virtual eigenfunction reconstruction requires a polynomial dictionary specification.")
    from spectral.evaluation import polynomial_eigenfunction_values
    try:
        return polynomial_eigenfunction_values(data[specification["states"]],data[specification["coefficients"]],
            dictionary=dictionary,state_preprocessing=specification.get("state_preprocessing") or {},columns=columns)
    except IndexError as error:
        raise IndexError("Virtual derived-array column index is out of range.") from error

def load_derived_array(path: str | Path, name: str, *, columns=None) -> np.ndarray:
    """Read, alias, or lazily reconstruct a canonical derived array."""
    target=Path(path)
    with np.load(target,allow_pickle=False) as data:
        if name in data.files: return _column_selection(np.asarray(data[name]),columns)
        if "metadata_json" not in data.files: raise KeyError(name)
        metadata=json.loads(str(data["metadata_json"])); analysis=metadata.get("analysis_metadata") or {}
        aliases=(analysis.get("storage") or {}).get("array_aliases") or {}; canonical=aliases.get(name)
        if canonical is not None:
            if canonical not in data.files: raise KeyError(name)
            return _column_selection(np.asarray(data[canonical]),columns)
        virtual=((analysis.get("storage") or {}).get("virtual_arrays") or {}).get(name)
        if virtual is None: raise KeyError(name)
        return _virtual_derived_array(data,virtual,columns)

def save_derived(source: Trajectory, analysis: str, analysis_metadata: dict | None = None, *, compression: str=DERIVED_DEFAULT_COMPRESSION, publication: str="immutable", **arrays) -> Path:
    """Publish a derived record without replacing an earlier analysis by default."""
    if source.source is None: raise ValueError("Derived archive publication requires a source archive path.")
    if publication not in {"immutable","replace"}: raise ValueError("publication must be 'immutable' or 'replace'.")
    root=POSTPROCESSING_ARCHIVE.archive_root/"derived"/analysis; root.mkdir(parents=True,exist_ok=True)
    payload={name:np.asarray(value) for name,value in arrays.items()}
    source_metadata=json.dumps(source.metadata,default=str,separators=(",",":"),sort_keys=True)
    source_metadata_hash=hashlib.sha256(source_metadata.encode()).hexdigest()
    signature_payload=json.dumps({"analysis":analysis,"source_archive":str(Path(source.source).resolve()),"source_metadata_sha256":source_metadata_hash,"analysis_metadata":analysis_metadata or {}},default=str,separators=(",",":"),sort_keys=True)
    signature=hashlib.sha256(signature_payload.encode()).hexdigest(); publication_id=str(uuid.uuid4()); created=datetime.now(timezone.utc)
    if publication=="replace": path=root/f"{source.source.stem}__{analysis}.npz"
    else:
        stamp=created.strftime("%Y%m%dT%H%M%S.%fZ"); path=root/f"{source.source.stem}__{analysis}__{stamp}__{publication_id[:8]}.npz"
        if path.exists(): raise FileExistsError(f"Refusing to replace immutable derived record: {path}")
    payload["metadata_json"]=json.dumps({"format":"dynsyslab.derived/v1","analysis":analysis,"source_archive":str(Path(source.source).resolve()),"source_metadata_sha256":source_metadata_hash,
        "analysis_signature_sha256":signature,"analysis_metadata":analysis_metadata or {},"storage_compression":f"npz-{compression}",
        "publication":{"mode":publication,"id":publication_id,"created_utc":created.isoformat(),"immutable":publication=="immutable"}},default=str,separators=(",",":"),sort_keys=True)
    return _write_derived_npz(path,payload,compression=compression)

def discover_derived(*, analysis: str|None=None, source_archive: str|Path|None=None) -> list[DerivedRecord]:
    """Discover published derived records newest-first through compact metadata."""
    root=POSTPROCESSING_ARCHIVE.archive_root/"derived"
    if not root.exists(): paths=[]
    elif analysis: paths=list((root/analysis).glob("*.npz"))
    else: paths=list(root.rglob("*.npz"))
    expected=None if source_archive is None else str(Path(source_archive).resolve()); records=[]
    for path in paths:
        try:
            with np.load(path,allow_pickle=False) as data: metadata=json.loads(str(data["metadata_json"]))
        except (KeyError,OSError,ValueError,zipfile.BadZipFile): continue
        source=str(Path(metadata.get("source_archive","")).resolve())
        if expected is not None and source!=expected: continue
        publication=metadata.get("publication") or {}
        records.append(DerivedRecord(path.resolve(),str(metadata.get("analysis") or path.parent.name),metadata.get("source_archive",""),
            publication.get("created_utc"),publication.get("id"),metadata.get("analysis_signature_sha256")))
    return sorted(records,key=lambda item:(item.created_utc or "",item.path.name),reverse=True)

def compact_derived_provenance(path: str | Path) -> Path:
    """Remove a legacy duplicated source-metadata snapshot, losslessly."""
    target=Path(path)
    with np.load(target,allow_pickle=False) as source: payload={name:source[name] for name in source.files}
    metadata=json.loads(str(payload["metadata_json"])); snapshot=metadata.pop("source_metadata",None)
    if snapshot is not None:
        encoded=json.dumps(snapshot,default=str,separators=(",",":"),sort_keys=True)
        metadata["source_metadata_sha256"]=__import__("hashlib").sha256(encoded.encode()).hexdigest()
        payload["metadata_json"]=json.dumps(metadata,default=str,separators=(",",":"),sort_keys=True)
        _write_derived_npz(target,payload,compression=_archive_compression(target))
    return target

def recompress_derived(path: str | Path, *, algorithm: str="lzma") -> Path:
    """Losslessly rewrite an existing NPZ using the selected ZIP codec.

    This changes only the container encoding, never values or archive metadata.
    """
    path=Path(path)
    if algorithm not in _COMPRESSION_TYPES: raise ValueError("algorithm must be 'lzma', 'deflate', or 'stored'.")
    with np.load(path,allow_pickle=False) as source: arrays={name:source[name] for name in source.files}
    if "metadata_json" in arrays:
        metadata=json.loads(str(arrays["metadata_json"])); metadata["storage_compression"]=f"npz-{algorithm}"
        arrays["metadata_json"]=json.dumps(metadata,default=str,separators=(",",":"),sort_keys=True)
    return _write_derived_npz(path,arrays,compression=algorithm)


__all__=["DERIVED_DEFAULT_COMPRESSION","DerivedRecord","Trajectory",
    "compact_derived_provenance","discover_derived","integration_coordinates",
    "load_derived_array","load_solution","recompress_derived","save_derived"]
