"""KROM-aware phase-space errors for approximate Koopman input maps."""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Sequence
from uuid import uuid4

import numpy as np


@dataclass(frozen=True)
class ControlApproximationErrorField:
    x_values: np.ndarray
    y_values: np.ndarray
    exact_norm: np.ndarray
    absolute_error: np.ndarray
    relative_error: np.ndarray
    approximation: str
    metadata: dict[str,Any]

    def __post_init__(self):
        x=np.array(self.x_values,dtype=float,copy=True).reshape(-1); y=np.array(self.y_values,dtype=float,copy=True).reshape(-1)
        expected=(y.size,x.size)
        exact=np.array(self.exact_norm,dtype=float,copy=True); absolute=np.array(self.absolute_error,dtype=float,copy=True)
        relative=np.array(self.relative_error,dtype=float,copy=True)
        if x.size<2 or y.size<2 or np.any(np.diff(x)<=0.) or np.any(np.diff(y)<=0.):
            raise ValueError("Error-field coordinates must be increasing and contain at least two points per axis.")
        if exact.shape!=expected or absolute.shape!=expected or relative.shape!=expected:
            raise ValueError("Error-field arrays must align with the y-by-x grid.")
        if any(np.any(~np.isfinite(item)) for item in (x,y,exact,absolute,relative)):
            raise ValueError("Control-error fields must be finite.")
        for item in (x,y,exact,absolute,relative): item.setflags(write=False)
        object.__setattr__(self,"x_values",x); object.__setattr__(self,"y_values",y)
        object.__setattr__(self,"exact_norm",exact); object.__setattr__(self,"absolute_error",absolute)
        object.__setattr__(self,"relative_error",relative); object.__setattr__(self,"metadata",dict(self.metadata))

    def summary(self) -> dict[str,Any]:
        def distribution(values):
            flat=np.asarray(values,dtype=float).reshape(-1)
            return {"minimum":float(np.min(flat)),"median":float(np.median(flat)),
                "mean":float(np.mean(flat)),"p95":float(np.quantile(flat,.95)),"maximum":float(np.max(flat))}
        return {"encoding":"koopman-control-approximation-error-field/v1","approximation":self.approximation,
            "grid_shape":list(self.absolute_error.shape),"rectangle":{"x":[float(self.x_values[0]),float(self.x_values[-1])],
            "y":[float(self.y_values[0]),float(self.y_values[-1])]},"absolute":distribution(self.absolute_error),
            "relative":distribution(self.relative_error),"metadata":dict(self.metadata)}


@dataclass(frozen=True)
class ControlApproximationEvaluation:
    """One fitted input map, its controlled-data lineage, and its error field."""

    evaluation_id: str
    label: str
    approximation: str
    representation: Any
    controlled_dataset: dict[str,Any]
    error_field: ControlApproximationErrorField
    fit_diagnostics: dict[str,Any] = field(default_factory=dict)
    metadata: dict[str,Any] = field(default_factory=dict)

    def __post_init__(self):
        from .control_inputs import BilinearControlMatrix, ConstantControlMatrix
        expected="constant" if isinstance(self.representation,ConstantControlMatrix) else (
            "bilinear" if isinstance(self.representation,BilinearControlMatrix) else None)
        if expected is None or expected!=self.approximation:
            raise TypeError("Control evaluation approximation and representation type do not agree.")
        if self.error_field.approximation!=self.approximation:
            raise ValueError("Control evaluation and heat-map approximation names do not agree.")
        if not self.evaluation_id or not self.label.strip():
            raise ValueError("Control approximation evaluations require an id and label.")
        object.__setattr__(self,"controlled_dataset",dict(self.controlled_dataset))
        object.__setattr__(self,"fit_diagnostics",dict(self.fit_diagnostics))
        object.__setattr__(self,"metadata",dict(self.metadata))

    def summary(self) -> dict[str,Any]:
        return {"encoding":"koopman-control-approximation-evaluation/v1",
            "evaluation_id":self.evaluation_id,"label":self.label,
            "approximation":self.approximation,"representation":self.representation.summary(),
            "controlled_dataset":dict(self.controlled_dataset),
            "fit_diagnostics":dict(self.fit_diagnostics),"error_field":self.error_field.summary(),
            "metadata":dict(self.metadata)}


def _dataset_identity(dataset) -> dict[str,Any]:
    if isinstance(dataset,dict):
        summary=dict(dataset)
        encoded=json.dumps(summary,sort_keys=True,default=str,separators=(",",":")).encode()
        summary.setdefault("fingerprint_sha256",sha256(encoded).hexdigest())
        return summary
    trajectories=getattr(dataset,"trajectories",None)
    if trajectories is None or not hasattr(dataset,"summary"):
        raise TypeError("controlled_dataset must be a ControlledDataset or a provenance dictionary.")
    digest=sha256()
    for trajectory in trajectories:
        for value in (trajectory.time,trajectory.states,trajectory.inputs):
            array=np.ascontiguousarray(value)
            digest.update(str(array.dtype).encode()); digest.update(str(array.shape).encode())
            digest.update(array.view(np.uint8))
    summary=dict(dataset.summary()); summary["fingerprint_sha256"]=digest.hexdigest()
    return summary


def control_approximation_error_fields(
    krom,
    bounds: tuple[Sequence[float],Sequence[float]],
    *,
    grid_shape: tuple[int,int]=(81,81),
    representations: Sequence[str]=( "constant","bilinear"),
    gradient_backend: str="auto",
    relative_floor: float|None=None,
) -> dict[str,ControlApproximationErrorField]:
    """Compare any stored input approximation to exact ``B_psi`` on a grid."""
    lower=np.asarray(bounds[0],dtype=float).reshape(-1); upper=np.asarray(bounds[1],dtype=float).reshape(-1)
    if lower.shape!=(2,) or upper.shape!=(2,) or np.any(upper<=lower):
        raise ValueError("A two-dimensional increasing state rectangle is required.")
    ny,nx=(int(value) for value in grid_shape)
    if nx<2 or ny<2: raise ValueError("grid_shape must contain at least two points per axis.")
    x_values=np.linspace(lower[0],upper[0],nx); y_values=np.linspace(lower[1],upper[1],ny)
    states=np.column_stack((np.tile(x_values,ny),np.repeat(y_values,nx)))
    exact=np.asarray([krom.B_psi(state,output="complex",representation="exact",
        gradient_backend=gradient_backend) for state in states]).reshape(ny,nx,krom.lifted_dimension,krom.input_dimension)
    exact_norm=np.linalg.norm(exact,axis=(-2,-1)); scale=float(np.max(exact_norm))
    floor=np.sqrt(np.finfo(float).eps)*max(scale,1.) if relative_floor is None else float(relative_floor)
    if floor<=0. or not np.isfinite(floor): raise ValueError("relative_floor must be finite and positive.")
    fields={}
    available=krom.available_control_representations()
    for name in representations:
        if name not in {"constant","bilinear"}: raise ValueError("Approximation names must be 'constant' or 'bilinear'.")
        if name not in available: raise ValueError(f"KROM has no stored {name!r} control representation.")
        approximation=np.asarray([krom.B_psi(state,output="complex",representation=name) for state in states])
        approximation=approximation.reshape(exact.shape)
        absolute=np.linalg.norm(exact-approximation,axis=(-2,-1)); relative=absolute/np.maximum(exact_norm,floor)
        metadata={"norm":"frobenius","gradient_backend":gradient_backend,
            "exact_eigenfunction_representation_kinds":[item.representation_kind for item in krom.eigenfunctions],
            "relative_denominator":"max(norm(exact_Bpsi), relative_floor)","relative_floor":floor}
        fields[name]=ControlApproximationErrorField(x_values,y_values,exact_norm,absolute,relative,name,metadata)
    return fields


def evaluate_control_approximation(
    krom,
    representation,
    controlled_dataset,
    bounds: tuple[Sequence[float],Sequence[float]],
    *,
    label: str,
    grid_shape: tuple[int,int]=(81,81),
    gradient_backend: str="auto",
    relative_floor: float|None=None,
    fit_diagnostics: dict[str,Any]|None=None,
    metadata: dict[str,Any]|None=None,
) -> ControlApproximationEvaluation:
    """Evaluate a candidate fit without replacing the KROM's active input map."""
    from .control_inputs import BilinearControlMatrix, ConstantControlMatrix
    if not isinstance(label,str) or not label.strip():
        raise ValueError("A control-approximation evaluation label is required.")
    if isinstance(representation,ConstantControlMatrix): name="constant"
    elif isinstance(representation,BilinearControlMatrix): name="bilinear"
    else: raise TypeError("representation must be ConstantControlMatrix or BilinearControlMatrix.")
    candidate=krom.with_control_representation(name,representation)
    field_value=control_approximation_error_fields(candidate,bounds,grid_shape=grid_shape,
        representations=(name,),gradient_backend=gradient_backend,
        relative_floor=relative_floor)[name]
    source=_dataset_identity(controlled_dataset)
    field_value=ControlApproximationErrorField(field_value.x_values,field_value.y_values,
        field_value.exact_norm,field_value.absolute_error,field_value.relative_error,name,
        {**field_value.metadata,"controlled_dataset_fingerprint_sha256":source["fingerprint_sha256"]})
    diagnostics=dict(fit_diagnostics or getattr(representation,"provenance",{}) or {})
    return ControlApproximationEvaluation(str(uuid4()),label.strip(),name,representation,source,
        field_value,diagnostics,{"comparison_truth":"analytic chain-rule B_psi of fitted eigenfunctions",
        **dict(metadata or {})})


def compare_control_approximation_evaluations(
    evaluations: Sequence[ControlApproximationEvaluation],
    *, metric: str="relative.mean",
) -> dict[str,Any]:
    """Rank retained fits, including fits of the same family from different data."""
    items=tuple(evaluations)
    if len(items)<2: raise ValueError("Comparison requires at least two control evaluations.")
    reference=items[0].error_field
    if any(not np.array_equal(item.error_field.x_values,reference.x_values)
           or not np.array_equal(item.error_field.y_values,reference.y_values)
           or not np.allclose(item.error_field.exact_norm,reference.exact_norm,rtol=1e-12,atol=1e-14)
           for item in items[1:]):
        raise ValueError("Control evaluations must share one grid and exact-B_psi truth field before ranking.")
    allowed={"relative.mean","relative.p95","relative.maximum",
        "absolute.mean","absolute.p95","absolute.maximum"}
    if metric not in allowed: raise ValueError(f"metric must be one of {sorted(allowed)}.")
    category,key=metric.split(".")
    rows=[]
    for item in items:
        summary=item.error_field.summary()[category]
        rows.append({"evaluation_id":item.evaluation_id,"label":item.label,
            "approximation":item.approximation,"controlled_dataset_fingerprint_sha256":
            item.controlled_dataset.get("fingerprint_sha256"),"metric":metric,"value":float(summary[key])})
    rows.sort(key=lambda row:row["value"])
    best=rows[0]["value"]
    for rank,row in enumerate(rows,1):
        row["rank"]=rank; row["ratio_to_best"]=(float(row["value"]/best) if best>0. else (1. if row["value"]==0. else float("inf")))
    return {"encoding":"koopman-control-approximation-comparison/v1","metric":metric,
        "evaluation_count":len(rows),"ranking":rows}


def save_control_approximation_evaluations(
    path: str|Path,
    evaluations: Sequence[ControlApproximationEvaluation],
    *, metadata: dict[str,Any]|None=None,
) -> Path:
    """Persist multiple candidate fits and heat maps without selecting a winner."""
    from .control_inputs import BilinearControlMatrix, ConstantControlMatrix
    items=tuple(evaluations)
    if not items: raise ValueError("At least one control approximation evaluation is required.")
    destination=Path(path).expanduser().resolve()
    if destination.suffix.lower()!=".npz": destination=destination.with_suffix(".npz")
    destination.parent.mkdir(parents=True,exist_ok=True)
    payload={}; descriptions=[]
    for index,item in enumerate(items):
        prefix=f"evaluation_{index}"; field_value=item.error_field
        payload[f"{prefix}_x_values"]=field_value.x_values; payload[f"{prefix}_y_values"]=field_value.y_values
        payload[f"{prefix}_exact_norm"]=field_value.exact_norm
        payload[f"{prefix}_absolute_error"]=field_value.absolute_error
        payload[f"{prefix}_relative_error"]=field_value.relative_error
        if isinstance(item.representation,ConstantControlMatrix):
            payload[f"{prefix}_matrix"]=item.representation.matrix
            if item.representation.discrete_matrix is not None:
                payload[f"{prefix}_discrete_matrix"]=item.representation.discrete_matrix
        elif isinstance(item.representation,BilinearControlMatrix):
            payload[f"{prefix}_operators"]=item.representation.operators
            payload[f"{prefix}_offsets"]=item.representation.offsets
            if item.representation.conjugate_operators is not None:
                payload[f"{prefix}_conjugate_operators"]=item.representation.conjugate_operators
        descriptions.append({"prefix":prefix,"summary":item.summary(),
            "field_metadata":dict(field_value.metadata)})
    payload["metadata_json"]=np.asarray(json.dumps({
        "format":"dynsyslab.control-approximation-evaluations/v1",
        "evaluations":descriptions,"metadata":dict(metadata or {})},sort_keys=True,default=str))
    np.savez_compressed(destination,**payload); return destination


def load_control_approximation_evaluations(path: str|Path) -> tuple[list[ControlApproximationEvaluation],dict[str,Any]]:
    """Restore independently executable fits and their full heat-map fields."""
    from .control_inputs import BilinearControlMatrix, ConstantControlMatrix
    source=Path(path).expanduser().resolve(); result=[]
    with np.load(source,allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        if metadata.get("format")!="dynsyslab.control-approximation-evaluations/v1":
            raise ValueError(f"{source} is not a control-approximation evaluation archive.")
        for description in metadata["evaluations"]:
            prefix=description["prefix"]; summary=description["summary"]; rep=summary["representation"]
            if summary["approximation"]=="constant":
                discrete=(np.asarray(data[f"{prefix}_discrete_matrix"]) if f"{prefix}_discrete_matrix" in data.files else None)
                representation=ConstantControlMatrix(np.asarray(data[f"{prefix}_matrix"]),
                    rep.get("time_domain","continuous"),rep.get("delta_t"),
                    dict(rep.get("provenance") or {}),discrete)
            else:
                conjugate=(np.asarray(data[f"{prefix}_conjugate_operators"]) if f"{prefix}_conjugate_operators" in data.files else None)
                representation=BilinearControlMatrix(np.asarray(data[f"{prefix}_operators"]),
                    np.asarray(data[f"{prefix}_offsets"]),rep.get("time_domain","continuous"),
                    rep.get("delta_t"),dict(rep.get("provenance") or {}),conjugate)
            field_value=ControlApproximationErrorField(np.asarray(data[f"{prefix}_x_values"]),
                np.asarray(data[f"{prefix}_y_values"]),np.asarray(data[f"{prefix}_exact_norm"]),
                np.asarray(data[f"{prefix}_absolute_error"]),np.asarray(data[f"{prefix}_relative_error"]),
                summary["approximation"],dict(description.get("field_metadata") or {}))
            result.append(ControlApproximationEvaluation(summary["evaluation_id"],summary["label"],
                summary["approximation"],representation,dict(summary["controlled_dataset"]),field_value,
                dict(summary.get("fit_diagnostics") or {}),dict(summary.get("metadata") or {})))
    metadata["source_archive"]=str(source); return result,metadata


def save_control_approximation_error_fields(path: str|Path, fields: dict[str,ControlApproximationErrorField], *,
                                            metadata: dict[str,Any]|None=None) -> Path:
    """Archive typed heat-map fields for read-only native Viewer display."""
    if not fields: raise ValueError("At least one control-approximation error field is required.")
    destination=Path(path).expanduser().resolve(); destination.parent.mkdir(parents=True,exist_ok=True)
    names=list(fields); first=fields[names[0]]
    if any(not np.array_equal(field.x_values,first.x_values) or not np.array_equal(field.y_values,first.y_values)
           for field in fields.values()):
        raise ValueError("Archived error fields must share one phase-space grid.")
    payload={"x_values":first.x_values,"y_values":first.y_values}
    descriptions=[]
    for index,name in enumerate(names):
        field=fields[name]; prefix=f"field_{index}"
        payload[f"{prefix}_exact_norm"]=field.exact_norm; payload[f"{prefix}_absolute_error"]=field.absolute_error
        payload[f"{prefix}_relative_error"]=field.relative_error
        descriptions.append({"name":name,"approximation":field.approximation,"prefix":prefix,"metadata":field.metadata})
    payload["metadata_json"]=np.asarray(json.dumps({"format":"dynsyslab.control-approximation-error-fields/v1",
        "fields":descriptions,"metadata":dict(metadata or {})},sort_keys=True))
    np.savez_compressed(destination,**payload); return destination


def load_control_approximation_error_fields(path: str|Path) -> tuple[dict[str,ControlApproximationErrorField],dict[str,Any]]:
    """Load a control-error archive without importing Viewer."""
    source=Path(path).expanduser().resolve()
    with np.load(source,allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        if metadata.get("format")!="dynsyslab.control-approximation-error-fields/v1":
            raise ValueError(f"{source} is not a DynSysLab control-error archive.")
        x=np.asarray(data["x_values"]); y=np.asarray(data["y_values"]); fields={}
        for description in metadata["fields"]:
            prefix=description["prefix"]
            fields[description["name"]]=ControlApproximationErrorField(x,y,
                np.asarray(data[f"{prefix}_exact_norm"]),np.asarray(data[f"{prefix}_absolute_error"]),
                np.asarray(data[f"{prefix}_relative_error"]),description["approximation"],description.get("metadata") or {})
    metadata["source_archive"]=str(source); return fields,metadata


__all__=["ControlApproximationErrorField","ControlApproximationEvaluation",
    "compare_control_approximation_evaluations","control_approximation_error_fields",
    "evaluate_control_approximation","load_control_approximation_error_fields",
    "load_control_approximation_evaluations","save_control_approximation_error_fields",
    "save_control_approximation_evaluations"]
