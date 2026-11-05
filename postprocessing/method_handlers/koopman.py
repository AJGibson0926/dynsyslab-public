"""Koopman-family processor handlers and publication contracts."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from ..core import Trajectory, load_solution, save_derived
from ..koopman.dictionary import PolynomialDictionary
from ..koopman.dmd import dmd, hankel_dmd
from ..koopman.edmd import edmd
from ..koopman.koopman_sindy import koopman_with_sindy
from ..koopman.pipeline import _save
from ..koopman.spectrum import DynamicalSystemReference, KoopmanEigenvalue, KoopmanSpectrum
from ..koopman.storage import sindy_archive_payload
from .base import MethodHandler, WorkflowContext, source_archive


def _archive_spectrum(path) -> KoopmanSpectrum:
    with np.load(Path(path),allow_pickle=False) as data:
        metadata=json.loads(str(data["metadata_json"]))
        values=np.asarray(data["eigenvalues"],complex) if "eigenvalues" in data.files else None
    if values is None:
        raise KeyError("This spectrum archive stores eigenvalues by reference; provide its source spectrum archive.")
    source=load_solution(metadata["source_archive"])
    reference=DynamicalSystemReference.from_archive_metadata(source.metadata,str(source.source))
    delta_t=float(np.diff(source.time)[0])
    eigenvalues=[KoopmanEigenvalue.from_discrete(value,delta_t,reference,
        {"method":"archive-reference","index":index}) for index,value in enumerate(values)]
    return KoopmanSpectrum(reference,eigenvalues,[],metadata={"method":"archive-reference"})


def _run_dmd(analysis: Trajectory, source: Trajectory, settings: dict[str,Any], context: WorkflowContext):
    return dmd(analysis.values,time=analysis.time,metadata=analysis.metadata,
               source_archive=source_archive(source),**settings)


def _run_hdmd(analysis: Trajectory, source: Trajectory, settings: dict[str,Any], context: WorkflowContext):
    return hankel_dmd(analysis.values,time=analysis.time,metadata=analysis.metadata,
                      source_archive=source_archive(source),**settings)


def _run_edmd(analysis: Trajectory, source: Trajectory, settings: dict[str,Any], context: WorkflowContext):
    dictionary=settings.pop("dictionary",None) or PolynomialDictionary()
    return edmd(analysis.values,dictionary,time=analysis.time,metadata=analysis.metadata,
                source_archive=source_archive(source),**settings)


def _run_sindy(analysis: Trajectory, source: Trajectory, settings: dict[str,Any], context: WorkflowContext):
    # Preserve the generic/manual processor's historical materialized-array
    # behavior.  The immutable standard H-DMD -> SINDy policy opts into the
    # coefficient-first representation explicitly.
    archive_storage=settings.pop("archive_storage","full")
    settings.setdefault("compute_training_values",archive_storage=="full")
    supplied=settings.pop("spectrum","previous")
    spectrum=context.spectrum if supplied is None or supplied=="previous" else (
        _archive_spectrum(supplied) if isinstance(supplied,(str,Path)) else supplied)
    if spectrum is None:
        raise ValueError("Koopman-with-SINDy needs a preceding or explicit spectrum.")
    result=koopman_with_sindy(spectrum=spectrum,trajectories=[analysis],**settings)
    result.metadata["archive_storage_policy"]=archive_storage
    if supplied is None or supplied=="previous":
        archive=context.spectrum_archive
        result.metadata["eigenvalue_source_archive"]=None if archive is None else str(archive)
    elif isinstance(supplied,(str,Path)):
        result.metadata["eigenvalue_source_archive"]=str(Path(supplied).expanduser().resolve())
    return result


def _publish_spectrum(source: Trajectory, result):
    return _save(source,result)


def _publish_sindy(source: Trajectory, result):
    arrays,metadata=sindy_archive_payload(result,
        storage_policy=result.metadata.get("archive_storage_policy","full"))
    return save_derived(source,"koopman-with-sindy",analysis_metadata=metadata,**arrays)


KOOPMAN_HANDLERS={
    "dmd":MethodHandler(_run_dmd,_publish_spectrum),
    "hdmd":MethodHandler(_run_hdmd,_publish_spectrum),
    "edmd":MethodHandler(_run_edmd,_publish_spectrum),
    "sindy":MethodHandler(_run_sindy,_publish_sindy),
}

__all__=["KOOPMAN_HANDLERS"]
