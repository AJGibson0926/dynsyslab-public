"""Read-only extraction of successful paths for conditioned candidate fits."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from analysis_data import load_solution
from .common import ACTIVE_WINDOW


def prior_successful_exact_path(
    output_root: Path, indices: list[int]
) -> tuple[np.ndarray, str]:
    """Return the latest matching exact-control path from a campaign report."""
    report_path = output_root / "campaign-report.json"
    if not report_path.is_file():
        raise FileNotFoundError(
            "Path-conditioned fitting needs a prior campaign report."
        )
    report = json.loads(report_path.read_text(encoding="utf-8"))
    for attempt in reversed(report.get("attempts", [])):
        if attempt.get("candidate", {}).get("retained_indices") != indices:
            continue
        exact = attempt.get("cases", {}).get("exact", {})
        if exact.get("status") != "complete" or not exact.get("archive_path"):
            continue
        trajectory = load_solution(exact["archive_path"])
        time = np.asarray(trajectory.time, dtype=float)
        values = np.asarray(trajectory.values, dtype=float)
        mask = (time >= ACTIVE_WINDOW[0]) & (time <= ACTIVE_WINDOW[1])
        path = values[mask][::5]
        if path.shape[0] >= 3 and np.all(np.isfinite(path)):
            return path, str(Path(exact["archive_path"]).resolve())
    raise ValueError(
        f"No prior successful exact-control path exists for modes {indices}."
    )
