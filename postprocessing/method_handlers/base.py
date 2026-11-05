"""Method-handler contract for ordered post-processing workflows."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from ..core import Trajectory
from ..workflow import WorkflowContext

RunMethod = Callable[[Trajectory, Trajectory, dict[str, Any], WorkflowContext], Any]
PublishMethod = Callable[[Trajectory, Any], Path | None]


@dataclass(frozen=True)
class MethodHandler:
    """Own one method's execution and derived-archive representation."""
    run: RunMethod
    publish: PublishMethod


def source_archive(source: Trajectory) -> str | None:
    return None if source.source is None else str(source.source)


__all__=["MethodHandler","WorkflowContext","source_archive"]
