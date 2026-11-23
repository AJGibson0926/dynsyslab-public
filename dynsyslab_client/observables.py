"""Thin native observable evaluation over supplied state samples."""
from collections.abc import Mapping
from pathlib import Path
from typing import Any
from .periodic import _query


def evaluate_batch(request: Mapping[str, Any], *, executable: str | Path | None = None,
                   timeout: float = 120) -> dict[str, Any]:
    """Evaluate named real observables on flat sample-major states in native C++.

    Required: system, observables (names), states. No trajectory integration,
    refitting or differentiation is performed. Outputs are sample-major.
    """
    return _query("observable", {"operation": "evaluate_batch", **dict(request)},
                  "dynsyslab.observable-query/v1", "dynsyslab.observable-result/v1",
                  executable, timeout)
