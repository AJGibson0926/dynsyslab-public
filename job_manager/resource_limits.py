"""Numerical thread limits derived from explicit Job resource requests."""
from __future__ import annotations

from contextlib import contextmanager
import os
from typing import Any, Iterator, Mapping


THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def numerical_thread_policy(cores: int) -> dict[str, Any]:
    """Return the portable policy recorded with a Job and worker result."""
    count = int(cores)
    if count < 1:
        raise ValueError("A numerical thread policy requires at least one core.")
    return {
        "schema": "dynsyslab.numerical-thread-policy/v1",
        "allocated_cores": count,
        "threads_per_process": count,
        "environment_variables": list(THREAD_VARIABLES),
        "purpose": "prevent nested numerical-library oversubscription",
    }


def worker_thread_environment(cores: int, base: Mapping[str, str] | None = None) -> dict[str, str]:
    """Build an isolated worker environment without mutating the parent."""
    count = int(cores)
    numerical_thread_policy(count)
    environment = dict(os.environ if base is None else base)
    for name in THREAD_VARIABLES:
        environment[name] = str(count)
    environment["DYNSYSLAB_ALLOCATED_CORES"] = str(count)
    return environment


@contextmanager
def numerical_thread_limit(cores: int) -> Iterator[dict[str, Any]]:
    """Limit already-loaded BLAS/OpenMP pools for in-process serial Jobs."""
    policy = numerical_thread_policy(cores)
    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=int(cores)):
        yield policy


__all__ = [
    "THREAD_VARIABLES",
    "numerical_thread_limit",
    "numerical_thread_policy",
    "worker_thread_environment",
]
