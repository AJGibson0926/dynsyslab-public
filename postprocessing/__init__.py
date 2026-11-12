"""Independent archive-data post-processing package.

Archive primitives stay lightweight. Numerical method families and workflow
orchestration are resolved only when their public objects are requested.
"""
from __future__ import annotations

from importlib import import_module

from .core import (
    DERIVED_DEFAULT_COMPRESSION,
    DerivedRecord,
    Trajectory,
    compact_derived_provenance,
    discover_derived,
    integration_coordinates,
    load_derived_array,
    load_solution,
    recompress_derived,
    save_derived,
)
from time_series import TimeSeries


_LAZY_EXPORTS={
    "run_processor":(".pipeline","run_processor"),
    "plan_trajectory_hdmd_sindy":(".planning","plan_trajectory_hdmd_sindy"),
    "summarize_processing_runs":(".planning","summarize_processing_runs"),
    "SVDResult":(".pca_pod","SVDResult"),
    "HAVOKResult":(".pca_pod","HAVOKResult"),
    "PODMode":(".pca_pod","PODMode"),
    "PROM":(".pca_pod","PROM"),
    "PODSystemReference":(".pca_pod","PODSystemReference"),
    "svd":(".pca_pod","svd"),
    "havok":(".pca_pod","havok"),
    "svd_archive":(".pca_pod","svd_archive"),
    "havok_archive":(".pca_pod","havok_archive"),
    "processor":(".processor","processor"),
    "PostprocessingPipeline":(".processor","PostprocessingPipeline"),
    "ProcessingMethod":(".processor","ProcessingMethod"),
    "ProcessingMethodResult":(".processor","ProcessingMethodResult"),
    "ProcessingRun":(".processor","ProcessingRun"),
    "PreparedAnalysis":(".workflow","PreparedAnalysis"),
    "SourceResolver":(".workflow","SourceResolver"),
    "WorkflowContext":(".workflow","WorkflowContext"),
    "WorkflowPolicy":(".policies","WorkflowPolicy"),
    "TRAJECTORY_HDMD_SINDY_POLICY":(".policies","TRAJECTORY_HDMD_SINDY_POLICY"),
    "RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY":(".policies","RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY"),
    "LinearAlgebraBackend":(".linalg_backends","LinearAlgebraBackend"),
    "available_linear_algebra_backends":(".linalg_backends","available_linear_algebra_backends"),
    "select_linear_algebra_backend":(".linalg_backends","select_linear_algebra_backend"),
    "ResamplingOptions":(".resampling","ResamplingOptions"),
    "ResamplingDecisionRequired":(".resampling","ResamplingDecisionRequired"),
    "require_uniform_time_steps":(".resampling","require_uniform_time_steps"),
    "resample":(".resampling","resample"),
}


def __getattr__(name: str):
    target=_LAZY_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module=import_module(target[0],__name__)
    value=getattr(module,target[1])
    globals()[name]=value
    return value


def __dir__():
    return sorted(set(globals())|set(_LAZY_EXPORTS))


__all__=["Trajectory","DerivedRecord","load_solution","integration_coordinates","load_derived_array","save_derived","discover_derived","recompress_derived","compact_derived_provenance","DERIVED_DEFAULT_COMPRESSION","run_processor","plan_trajectory_hdmd_sindy","summarize_processing_runs","TimeSeries","SVDResult","HAVOKResult","PODMode","PROM","PODSystemReference","svd","havok","svd_archive","havok_archive","PostprocessingPipeline","processor","ProcessingMethod","ProcessingMethodResult","ProcessingRun","PreparedAnalysis","SourceResolver","WorkflowContext","WorkflowPolicy","TRAJECTORY_HDMD_SINDY_POLICY","RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY","LinearAlgebraBackend","available_linear_algebra_backends","select_linear_algebra_backend","ResamplingOptions","ResamplingDecisionRequired","require_uniform_time_steps","resample"]
