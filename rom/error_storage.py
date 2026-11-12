"""Persistence for attached eigenfunction and KROM control-error records."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .control_error import (load_control_approximation_evaluations,
                            save_control_approximation_evaluations)


def _complex(value: complex) -> dict[str,float]:
    return {"real":float(np.real(value)),"imag":float(np.imag(value))}


def _signature(krom) -> dict[str,Any]:
    return {"encoding":"krom-error-evaluation-signature/v1","method":krom.method,
        "coordinate_system":krom.coordinate_system,"lifted_dimension":krom.lifted_dimension,
        "continuous_eigenvalues":[_complex(value.continuous_time) for value in krom.eigenvalues]}


def save_krom_error_evaluations(path: str|Path, krom, *, metadata: dict[str,Any]|None=None) -> Path:
    """Save all attached free-data metrics, fitted maps, and heat-map arrays."""
    free=[]
    for index,function in enumerate(krom.eigenfunctions):
        free.append({"retained_index":index,"continuous_eigenvalue":_complex(function.eigenvalue.continuous_time),
            "evaluations":[item.summary() for item in function.error_evaluations.values()]})
    bundle={"encoding":"dynsyslab.krom-error-evaluations/v1","krom":_signature(krom),
        "free_trajectory_evaluations":free,"metadata":dict(metadata or {})}
    controls=list(krom.control_error_evaluations.values())
    destination=Path(path).expanduser().resolve()
    if destination.suffix.lower()!=".npz": destination=destination.with_suffix(".npz")
    destination.parent.mkdir(parents=True,exist_ok=True)
    if controls:
        return save_control_approximation_evaluations(destination,controls,metadata=bundle)
    np.savez_compressed(destination,metadata_json=np.asarray(json.dumps({
        "format":"dynsyslab.control-approximation-evaluations/v1","evaluations":[],
        "metadata":bundle},sort_keys=True,default=str)))
    return destination


def load_krom_error_evaluations(path: str|Path, krom, *, attach: bool=True,
                                replace: bool=False) -> dict[str,Any]:
    """Load and validate an error bundle against one exact KROM spectrum."""
    from spectral.errors import EigenfunctionErrorEvaluation
    controls,outer=load_control_approximation_evaluations(path)
    bundle=dict(outer.get("metadata") or {})
    if bundle.get("encoding")!="dynsyslab.krom-error-evaluations/v1":
        raise ValueError("Archive does not contain a KROM error-evaluation bundle.")
    expected=_signature(krom); actual=dict(bundle.get("krom") or {})
    if actual.get("lifted_dimension")!=expected["lifted_dimension"]:
        raise ValueError("Error bundle lifted dimension does not match the KROM.")
    archived=np.asarray([complex(item["real"],item["imag"])
        for item in actual.get("continuous_eigenvalues",[])])
    current=np.asarray([item.continuous_time for item in krom.eigenvalues])
    if archived.shape!=current.shape or not np.allclose(archived,current,rtol=1e-10,atol=1e-12):
        raise ValueError("Error bundle eigenvalues do not match the KROM.")
    free={}
    for group in bundle.get("free_trajectory_evaluations",[]):
        index=int(group["retained_index"])
        if index<0 or index>=krom.lifted_dimension: raise IndexError("Archived eigenfunction index is invalid.")
        evaluations=[EigenfunctionErrorEvaluation.from_summary(item) for item in group.get("evaluations",[])]
        free[index]=evaluations
        if attach:
            for evaluation in evaluations:
                krom.eigenfunctions[index].attach_error_evaluation(evaluation,replace=replace)
    if attach:
        for evaluation in controls: krom.attach_control_error_evaluation(evaluation,replace=replace)
    return {"encoding":bundle["encoding"],"source_archive":outer["source_archive"],
        "free_trajectory_evaluations":free,"control_approximation_evaluations":controls,
        "metadata":dict(bundle.get("metadata") or {})}


__all__=["load_krom_error_evaluations","save_krom_error_evaluations"]
