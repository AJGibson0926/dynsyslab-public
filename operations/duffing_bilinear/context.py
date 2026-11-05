"""Validated immutable-by-convention context for one Duffing search run."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from rom import load_krom_archive


class DuffingSearchContext:
    """Load shared sources once and retain bounded-search configuration."""

    def __init__(
        self,
        parent_archive: Path,
        controlled_dataset,
        controlled_states: np.ndarray,
        output_root: Path,
        *,
        output_dt: float,
        grid_size: int,
        case_timeout_seconds: float,
    ) -> None:
        self.parent_archive = parent_archive.resolve()
        self.parent = load_krom_archive(self.parent_archive)
        self.controlled_dataset = controlled_dataset
        self.controlled_states = np.asarray(controlled_states, dtype=float)
        self.output_root = output_root.resolve()
        self.output_dt = float(output_dt)
        self.grid_size = int(grid_size)
        self.case_timeout_seconds = float(case_timeout_seconds)
