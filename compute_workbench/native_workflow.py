"""Compatibility import for the dependency-free native Python client.

New code should import :mod:`dynsyslab_client`. This module remains so source
workflows do not need an abrupt import-path change. It owns no implementation.
"""

from dynsyslab_client.workflow import (
    NativeProgress,
    NativeWorkflow,
    NativeWorkflowCancelled,
    NativeWorkflowError,
    NativeWorkflowResourceLimit,
    find_library,
    library_candidates,
    plan_job,
    run_job,
)

__all__ = [
    "NativeProgress",
    "NativeWorkflow",
    "NativeWorkflowCancelled",
    "NativeWorkflowError",
    "NativeWorkflowResourceLimit",
    "find_library",
    "library_candidates",
    "plan_job",
    "run_job",
]
