"""Immutable named policies for reproducible post-processing workflows."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Mapping

from .koopman.dmd import (
    HDMD_DEFAULT_MATRIX_SIZE,
    HDMD_DEFAULT_SAMPLING_WINDOW,
    HDMD_DEFAULT_SVD_TOLERANCE,
)


@dataclass(frozen=True)
class WorkflowPolicy:
    """Named immutable defaults for one ordered workflow."""

    name: str
    hdmd: Mapping[str,Any]
    sindy: Mapping[str,Any]

    def hdmd_settings(self) -> dict[str,Any]:
        return dict(self.hdmd)

    def sindy_settings(self) -> dict[str,Any]:
        return dict(self.sindy)


TRAJECTORY_HDMD_SETTINGS=MappingProxyType({
    "matrix_size":HDMD_DEFAULT_MATRIX_SIZE,
    "sampling_window":HDMD_DEFAULT_SAMPLING_WINDOW,
    "observable":0,
    "matrix_layout":"window-snapshots",
    "mode_sort":"terminal-amplitude-energy",
    "rank":None,
    "svd_tolerance":HDMD_DEFAULT_SVD_TOLERANCE,
    "svd_driver":"gesvd",
    "delay_stride":1,
    "sample_stride":1,
    "svd_threads":None,
    "eigendecomposition_threads":None,
    "linear_algebra_backend":"auto",
    "max_memory_mb":8192.,
    "analysis_coordinates":"integration",
})

TRAJECTORY_SINDY_SETTINGS=MappingProxyType({
    "spectrum":"previous",
    "polynomial_degree":4,
    "threshold":.05,
    "max_iterations":20,
    "normalize":True,
    "gamma_method":"system_rhs",
    "compute_training_errors":False,
    "compute_training_values":False,
    "reuse_conjugate_fits":True,
    "null_space_backend":"gram-eigh",
    "max_memory_mb":8192.,
    "analysis_coordinates":"integration",
    "archive_storage":"coefficient-first",
})

TRAJECTORY_HDMD_SINDY_POLICY=WorkflowPolicy(
    "trajectory-hdmd-sindy",TRAJECTORY_HDMD_SETTINGS,TRAJECTORY_SINDY_SETTINGS)

# Reproducible RPE control-isolation policy.  This is deliberately named and
# separate from the general-purpose trajectory policy: it preserves the exact
# polynomial conventions used by the validated large-target KLQR experiment
# without silently changing defaults for unrelated systems.
RPE_CLEAN_ISOLATION_SINDY_SETTINGS=MappingProxyType({
    **dict(TRAJECTORY_SINDY_SETTINGS),
    "polynomial_degree":4,
    "threshold":1e-6,
    "max_iterations":1,
    "normalize":False,
    "gamma_method":"system_rhs",
    "compute_training_errors":True,
    "compute_training_values":False,
    "reuse_conjugate_fits":True,
    "null_space_backend":"gram-eigh",
    "archive_storage":"coefficient-first",
})
RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY=WorkflowPolicy(
    "rpe-clean-isolation-hdmd-sindy",TRAJECTORY_HDMD_SETTINGS,RPE_CLEAN_ISOLATION_SINDY_SETTINGS)

_POLICY_ALIASES={
    "default":TRAJECTORY_HDMD_SINDY_POLICY,
    TRAJECTORY_HDMD_SINDY_POLICY.name:TRAJECTORY_HDMD_SINDY_POLICY,
    RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY.name:RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY,
}


def workflow_policy(name: str="trajectory-hdmd-sindy") -> WorkflowPolicy:
    try:
        return _POLICY_ALIASES[name]
    except KeyError as error:
        raise ValueError("Available presets: default, trajectory-hdmd-sindy, rpe-clean-isolation-hdmd-sindy.") from error


__all__=["WorkflowPolicy","TRAJECTORY_HDMD_SETTINGS","TRAJECTORY_SINDY_SETTINGS",
         "TRAJECTORY_HDMD_SINDY_POLICY","RPE_CLEAN_ISOLATION_SINDY_SETTINGS",
         "RPE_CLEAN_ISOLATION_HDMD_SINDY_POLICY","workflow_policy"]
