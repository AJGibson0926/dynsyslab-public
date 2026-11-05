"""Solver/Post-processor application boundary for DynSysLab.

The package is deliberately a thin orchestration layer.  Numerical work stays
in :mod:`solver`, :mod:`control`, and :mod:`postprocessing`; this package owns
versioned job documents, discovery, planning, and archive/UI integration.
"""

from .contracts import JOB_SCHEMA, RESULT_SCHEMA, load_job, validate_job
from .execution import execute_job, plan_job
from .registry import build_registry

__all__ = [
    "JOB_SCHEMA", "RESULT_SCHEMA", "build_registry", "execute_job", "load_job",
    "plan_job", "validate_job",
]
