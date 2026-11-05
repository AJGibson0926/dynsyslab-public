"""Composable archive post-processing entry points."""
from __future__ import annotations
from typing import Callable, Any
from .core import Trajectory, load_solution

ArchiveProcessor = Callable[[str], Any]

def run_processor(archive_path: str, processor: Callable[[Trajectory], Any]) -> Any:
    """Load exactly one completed solution archive and process it."""
    return processor(load_solution(archive_path))
