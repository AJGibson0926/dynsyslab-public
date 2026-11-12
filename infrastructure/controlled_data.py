"""Neutral contracts for aligned controlled-trajectory data."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np


@dataclass(frozen=True)
class ControlledTrajectory:
    """One state trajectory with interval-held control samples."""

    time: np.ndarray
    states: np.ndarray
    inputs: np.ndarray
    metadata: dict[str, Any] | None = None

    def __post_init__(self):
        time = np.array(self.time, dtype=float, copy=True).reshape(-1)
        states = np.array(self.states, dtype=float, copy=True)
        inputs = np.array(self.inputs, dtype=float, copy=True)
        if inputs.ndim == 1:
            inputs = inputs[:, None]
        if time.size < 3 or np.any(~np.isfinite(time)) or np.any(np.diff(time) <= 0):
            raise ValueError(
                "Controlled trajectory time must be finite and strictly increasing."
            )
        if states.ndim != 2 or states.shape[0] != time.size or np.any(~np.isfinite(states)):
            raise ValueError(
                "Controlled trajectory states must align with time and be finite."
            )
        if inputs.ndim != 2 or inputs.shape[0] != time.size - 1 or np.any(~np.isfinite(inputs)):
            raise ValueError(
                "Controlled inputs must have one row per state-transition interval."
            )
        time.setflags(write=False)
        states.setflags(write=False)
        inputs.setflags(write=False)
        object.__setattr__(self, "time", time)
        object.__setattr__(self, "states", states)
        object.__setattr__(self, "inputs", inputs)
        object.__setattr__(self, "metadata", dict(self.metadata or {}))

    @property
    def delta_t(self) -> float:
        """Return the required uniform sampling interval."""
        differences = np.diff(self.time)
        delta = float(differences[0])
        if not np.allclose(differences, delta, rtol=1e-9, atol=1e-12):
            raise ValueError(
                "Controlled-data regression requires a uniform sample interval."
            )
        return delta

    @property
    def input_dimension(self) -> int:
        """Number of independently sampled input channels."""
        return int(self.inputs.shape[1])


@dataclass(frozen=True)
class ControlledDataset:
    """Validated ensemble of compatible controlled trajectories."""

    trajectories: tuple[ControlledTrajectory, ...]
    metadata: dict[str, Any] | None = None

    def __init__(
        self,
        trajectories: Sequence[ControlledTrajectory],
        metadata: dict[str, Any] | None = None,
    ):
        items = tuple(trajectories)
        if not items:
            raise ValueError("ControlledDataset requires at least one trajectory.")
        state_dimension = items[0].states.shape[1]
        inputs = items[0].input_dimension
        delta_t = items[0].delta_t
        if any(
            item.states.shape[1] != state_dimension or item.input_dimension != inputs
            for item in items
        ):
            raise ValueError(
                "Controlled trajectories must share state and input dimensions."
            )
        if any(
            not np.isclose(item.delta_t, delta_t, rtol=1e-9, atol=1e-12)
            for item in items
        ):
            raise ValueError("Controlled trajectories must share one sample interval.")
        object.__setattr__(self, "trajectories", items)
        object.__setattr__(self, "metadata", dict(metadata or {}))

    @property
    def delta_t(self) -> float:
        """Shared trajectory sampling interval."""
        return self.trajectories[0].delta_t

    @property
    def input_dimension(self) -> int:
        """Shared input dimension."""
        return self.trajectories[0].input_dimension

    @property
    def transition_count(self) -> int:
        """Total state-transition samples in the ensemble."""
        return int(sum(item.inputs.shape[0] for item in self.trajectories))

    def summary(self) -> dict[str, Any]:
        """Return compact, JSON-compatible dataset provenance."""
        return {
            "encoding": "controlled-trajectory-ensemble/v1",
            "trajectory_count": len(self.trajectories),
            "transition_count": self.transition_count,
            "delta_t": self.delta_t,
            "input_dimension": self.input_dimension,
            "metadata": dict(self.metadata),
        }


__all__ = ["ControlledDataset", "ControlledTrajectory"]
