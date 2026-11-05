"""PCA/POD/SVD/HAVOK processor handlers and publication contracts."""
from __future__ import annotations

from typing import Any

from ..core import Trajectory, save_derived
from ..pca_pod.algorithms import bagged_pca, hankel_pod, havok, pca, robust_pca, snapshot_pod, svd
from .base import MethodHandler, WorkflowContext, source_archive


def _simple(function):
    def run(analysis: Trajectory, source: Trajectory, settings: dict[str,Any], context: WorkflowContext):
        return function(analysis.time_series,**settings)
    return run


def _with_provenance(function):
    def run(analysis: Trajectory, source: Trajectory, settings: dict[str,Any], context: WorkflowContext):
        return function(analysis.time_series,source_metadata=analysis.metadata,
                        source_archive=source_archive(source),**settings)
    return run


def _publish_components(source: Trajectory, result):
    return save_derived(source,result.method,analysis_metadata=result.metadata,
        mean=result.mean,components=result.components,singular_values=result.singular_values,
        scores=result.scores,explained_variance_ratio=result.explained_variance_ratio)


def _publish_robust(source: Trajectory, result):
    return save_derived(source,"robust-pca",**result)


def _publish_svd(source: Trajectory, result):
    return save_derived(source,result.metadata["method"],analysis_metadata=result.metadata,
        left_singular_vectors=result.left_singular_vectors,singular_values=result.singular_values,
        right_singular_vectors=result.right_singular_vectors,mean=result.mean)


def _publish_havok(source: Trajectory, result):
    return save_derived(source,"havok",analysis_metadata=result.metadata,A=result.A,B=result.B,
        forcing=result.forcing,coordinates=result.coordinates,coordinate_time=result.coordinate_time,
        left_singular_vectors=result.left_singular_vectors,singular_values=result.singular_values,
        right_singular_vectors=result.right_singular_vectors,observable=result.observable)


DECOMPOSITION_HANDLERS={
    "pca":MethodHandler(_simple(pca),_publish_components),
    "pod":MethodHandler(_simple(snapshot_pod),_publish_components),
    "bagged_pca":MethodHandler(_simple(bagged_pca),_publish_components),
    "robust_pca":MethodHandler(_simple(robust_pca),_publish_robust),
    "svd":MethodHandler(_with_provenance(svd),_publish_svd),
    "hankel_pod":MethodHandler(_with_provenance(hankel_pod),_publish_svd),
    "havok":MethodHandler(_with_provenance(havok),_publish_havok),
}

__all__=["DECOMPOSITION_HANDLERS"]
