"""Composable eigenvalue, eigenfunction, and control-identification stages."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

from .spectrum import KoopmanSpectrum

EigenvalueMethod=Literal["hdmd","optimization","provided"]
EigenfunctionMethod=Literal["hdmd-interpolation","sindy","boundary-functions"]


@dataclass(frozen=True)
class KoopmanPipelineSpecification:
    eigenvalue_method: EigenvalueMethod
    eigenfunction_method: EigenfunctionMethod


@dataclass
class EigenvalueStageResult:
    spectrum: KoopmanSpectrum
    method: str
    artifact: Any | None=None


@dataclass
class EigenfunctionStageResult:
    spectrum: KoopmanSpectrum
    krom: Any
    method: str
    artifact: Any


@dataclass
class ComposableKoopmanResult:
    specification: KoopmanPipelineSpecification
    eigenvalues: EigenvalueStageResult
    eigenfunctions: EigenfunctionStageResult


def eigenvalues_from_hdmd(result) -> EigenvalueStageResult:
    """Expose H-DMD eigenvalues without coupling later stages to its functions."""
    source=result.spectrum
    spectrum=KoopmanSpectrum(source.system,list(source.eigenvalues),[],[],{
        "method":"hdmd-eigenvalues","stage_kind":"eigenvalue-extraction",
        "source_spectrum_id":source.metadata.get("spectrum_id"),
        "source_method":source.metadata.get("method"),
        "preprocessing":dict(source.metadata.get("preprocessing") or {}),
    },dynamical_system=source.dynamical_system,eigenvalue_source=source)
    return EigenvalueStageResult(spectrum,"hdmd",result)


def provided_eigenvalues(spectrum: KoopmanSpectrum) -> EigenvalueStageResult:
    """Declare an existing paired eigenvalue spectrum as a pipeline input."""
    return EigenvalueStageResult(spectrum,"provided",None)


def construct_eigenfunctions(stage: EigenvalueStageResult, method: EigenfunctionMethod, **kwargs) -> EigenfunctionStageResult:
    """Construct eigenfunctions from an independently selected eigenvalue stage."""
    if method=="hdmd-interpolation":
        if stage.method!="hdmd" or stage.artifact is None:
            raise ValueError("hdmd-interpolation requires eigenvalues_from_hdmd from the same H-DMD result.")
        artifact=stage.artifact
        return EigenfunctionStageResult(artifact.spectrum,artifact.krom,method,artifact)
    if method=="sindy":
        from .koopman_sindy import koopman_with_sindy
        artifact=koopman_with_sindy(spectrum=stage.spectrum,**kwargs)
        return EigenfunctionStageResult(artifact.spectrum,artifact.krom,method,artifact)
    if method=="boundary-functions":
        from .boundary_eigenfunctions import fit_boundary_eigenfunctions
        artifact=fit_boundary_eigenfunctions(stage.spectrum,**kwargs)
        return EigenfunctionStageResult(artifact.spectrum,artifact.krom,method,artifact)
    raise ValueError("Unknown eigenfunction construction method.")


def run_composable_koopman_pipeline(
    specification: KoopmanPipelineSpecification,
    *,
    eigenvalue_options: dict[str,Any],
    eigenfunction_options: dict[str,Any],
) -> ComposableKoopmanResult:
    """Run the algorithms selected independently in ``specification``."""
    if specification.eigenvalue_method=="hdmd":
        stage=eigenvalues_from_hdmd(eigenvalue_options["result"])
    elif specification.eigenvalue_method=="optimization":
        from .eigenvalue_optimization import optimize_trajectory_eigenvalues
        artifact=optimize_trajectory_eigenvalues(**eigenvalue_options)
        stage=EigenvalueStageResult(artifact.spectrum,"optimization",artifact)
    elif specification.eigenvalue_method=="provided":
        stage=provided_eigenvalues(eigenvalue_options["spectrum"])
    else: raise ValueError("Unknown eigenvalue extraction method.")
    functions=construct_eigenfunctions(stage,specification.eigenfunction_method,**eigenfunction_options)
    functions.spectrum.metadata["composable_pipeline"]={
        "eigenvalue_method":specification.eigenvalue_method,
        "eigenfunction_method":specification.eigenfunction_method,
        "eigenvalue_spectrum_id":stage.spectrum.metadata.get("spectrum_id"),
        "eigenfunction_spectrum_id":functions.spectrum.metadata.get("spectrum_id"),
    }
    return ComposableKoopmanResult(specification,stage,functions)


__all__=["ComposableKoopmanResult","EigenfunctionStageResult","EigenvalueStageResult",
    "KoopmanPipelineSpecification","construct_eigenfunctions","eigenvalues_from_hdmd",
    "provided_eigenvalues","run_composable_koopman_pipeline"]
