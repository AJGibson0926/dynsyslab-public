"""Minimal estimator contracts shared by solver and estimator plug-ins."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np


class Estimator(Protocol):
    """Structural contract for a solver-compatible state estimator."""

    def estimate(self,t: float,true_state: np.ndarray) -> np.ndarray: ...


@dataclass(frozen=True)
class IdentityEstimator:
    """Default estimator that returns an isolated copy of the true state."""

    def estimate(self,t,true_state):
        return np.asarray(true_state,dtype=float).copy()


__all__=["Estimator","IdentityEstimator"]
