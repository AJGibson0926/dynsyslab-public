"""Shared post-integration lifecycle for every solver execution pathway.

Numerical kernels return :class:`NumericalResult` objects.  This module owns
the application boundary that follows: interpolation, durable archival,
telemetry, Viewer invalidation, postprocessing, and partial-history exposure.
Keeping that sequence in one place lets scalar and batched kernels share
identical externally visible behavior.
"""
from __future__ import annotations

from archiver.gateway import SOLVER_ARCHIVE
from infrastructure.source_events import notify_source_update

from .interpolation import apply_interpolation
from .models import NumericalResult, SolveRequest


def finalize_result(request: SolveRequest, session, result: NumericalResult) -> NumericalResult:
    """Apply the canonical lifecycle to a kernel-produced result."""
    if result.status in {"complete", "paused"}:
        result = apply_interpolation(request, result)
        result.archive_path = str(session.complete(result))
        if result.status == "paused":
            session.telemetry("paused", checkpoint_path=result.metadata.get("checkpoint_path"))
            result.history_path = str(session.history_dir)
            notify_source_update(
                result.archive_path,
                source_type="solution",
                reason="solve_paused_at_checkpoint",
            )
            return result

        session.telemetry("completed", archive_path=result.archive_path)
        notify_source_update(result.archive_path, source_type="solution", reason="solve_completed")
        outcomes = []
        for processor in request.postprocessors:
            name = getattr(processor, "__name__", type(processor).__name__)
            try:
                derived = str(processor(result.archive_path))
                outcomes.append({"processor": name, "result": derived})
                notify_source_update(derived, source_type="derived", reason="postprocessing_completed")
            except Exception as exc:
                outcomes.append({"processor": name, "error": str(exc)})
        if outcomes:
            result.metadata["postprocessing"] = outcomes
            SOLVER_ARCHIVE.update_metadata(result.archive_path, {"postprocessing": outcomes})
    elif request.history.track_during_solve and request.history.keep_partial_history_on_failure:
        result.history_path = str(session.history_dir)

    if result.status == "failed":
        session.telemetry("failed", message=result.message)
    return result


__all__ = ["finalize_result"]
