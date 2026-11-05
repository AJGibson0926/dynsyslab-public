"""Canonical project paths and exclusions for tools that inspect dynsyslab."""
from __future__ import annotations

from pathlib import Path

PROJECT_ROOT=Path(__file__).resolve().parent
# Build products are intentionally excluded from every source-oriented tool.
GENERATED_DIRECTORIES=frozenset({
    ".pytest_cache", ".test-deps", "__pycache__", "build", "dist",
    "dynsyslab.egg-info", "test_reports", "testing", "tmp",
    "viewer_compute_archive", "viewer_outputs",
})

def is_generated(path: str|Path) -> bool:
    """Whether a path is inside a generated, non-source project directory."""
    target=Path(path).resolve()
    try:
        return any(part in GENERATED_DIRECTORIES for part in target.relative_to(PROJECT_ROOT).parts)
    except ValueError:
        return False

def source_files(pattern: str="*.py") -> list[Path]:
    """Source discovery that never walks generated package output."""
    return [path for path in PROJECT_ROOT.rglob(pattern) if not is_generated(path)]
