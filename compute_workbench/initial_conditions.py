"""Deterministic initial-condition generation and cluster packaging."""
from __future__ import annotations

from typing import Any
import numpy as np

from hierarchy import Cluster, Trajectory


def generate_initial_conditions(system: Any, spec: dict, *, controlled: bool) -> np.ndarray:
    """Generate deterministic physical-coordinate initial states for a job."""
    mode = spec.get("mode", "default")
    count = int(spec.get("count", 1))
    dimension = int(np.asarray(system.initial_state).size)
    if mode == "explicit":
        values = np.asarray(spec.get("values"), dtype=float)
        if values.ndim == 1:
            values = values[None, :]
        if values.shape != (count, dimension):
            raise ValueError(f"Explicit initial conditions must have shape ({count}, {dimension}).")
        return values
    default = np.asarray(system.initial_state, dtype=float)
    # The solver gate independently enforces the controlled-RPE equilibrium
    # default.  Keeping physical ICs here makes all hierarchy entries truthful.
    if mode == "default":
        return np.repeat(default[None, :], count, axis=0)
    ranges = np.asarray(spec.get("ranges"), dtype=float)
    if ranges.shape != (dimension, 2) or np.any(ranges[:, 0] > ranges[:, 1]):
        raise ValueError(f"Generated initial conditions require {dimension} ordered [low, high] ranges.")
    generator = np.random.default_rng(int(spec.get("seed", 0)))
    if mode == "uniform":
        return generator.uniform(ranges[:, 0], ranges[:, 1], size=(count, dimension))
    if mode == "normal":
        center = ranges.mean(axis=1)
        sigma = (ranges[:, 1] - ranges[:, 0]) / 6.0
        return np.clip(generator.normal(center, sigma, size=(count, dimension)),
                       ranges[:, 0], ranges[:, 1])
    raise ValueError(f"Unsupported initial-condition mode {mode!r}.")


def make_clusters(values: np.ndarray, count: int) -> list[Cluster]:
    """Distribute initial states round-robin into hierarchy clusters."""
    groups: list[list[Trajectory]] = [[] for _ in range(count)]
    for index, row in enumerate(np.asarray(values, dtype=float)):
        groups[index % count].append(Trajectory(row, label=f"trajectory-{index:03d}"))
    return [Cluster(items, label=f"cluster-{index:03d}") for index, items in enumerate(groups)]
