"""Classical and finite-horizon observability analysis for Koopman ROMs.

The routines in this module deliberately distinguish the raw classical
observability matrix from numerically safer finite-horizon Gramians.  They
operate on real/imaginary split Koopman coordinates, matching ``KROM.A`` and
``ObservableReconstruction.split_matrix`` exactly.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Literal

import numpy as np


ObservabilityMethod = Literal[
    "matrix", "discrete-gramian", "continuous-gramian", "continuous-infinite-gramian"
]


def _square_matrix(value: np.ndarray, name: str) -> np.ndarray:
    result = np.array(value, dtype=float, copy=True)
    if result.ndim != 2 or result.shape[0] != result.shape[1] or result.shape[0] < 1:
        raise ValueError(f"{name} must be a nonempty square matrix.")
    if np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values.")
    return result


def _measurement_matrix(value: np.ndarray, state_dimension: int) -> np.ndarray:
    result = np.array(value, dtype=float, copy=True)
    if result.ndim != 2 or result.shape[0] < 1 or result.shape[1] != state_dimension:
        raise ValueError("C must have shape (measurement_dimension, state_dimension).")
    if np.any(~np.isfinite(result)):
        raise ValueError("C must contain only finite values.")
    return result


def classical_observability_matrix(
    A: np.ndarray,
    C: np.ndarray,
    *,
    horizon: int | None = None,
    maximum_absolute_entry: float = 1e150,
) -> np.ndarray:
    """Return ``[C; C A; ...; C A**(horizon-1)]`` with overflow gates.

    The default horizon is the state dimension.  Large or unstable systems can
    make the power matrix meaningless in floating point; in that case this
    function fails explicitly and callers should use a finite-horizon Gramian.
    """
    state = _square_matrix(A, "A")
    measurement = _measurement_matrix(C, state.shape[0])
    steps = state.shape[0] if horizon is None else int(horizon)
    if steps < 1 or (horizon is not None and steps != horizon):
        raise ValueError("horizon must be a positive integer or None.")
    if not np.isfinite(maximum_absolute_entry) or maximum_absolute_entry <= 0.0:
        raise ValueError("maximum_absolute_entry must be finite and positive.")
    blocks = []
    block = measurement.copy()
    for _ in range(steps):
        if np.any(~np.isfinite(block)) or np.max(np.abs(block)) > maximum_absolute_entry:
            raise FloatingPointError(
                "The classical observability matrix exceeded its numerical safety bound; "
                "use a finite-horizon Gramian or rescale the lifted coordinates."
            )
        blocks.append(block.copy())
        block = block @ state
    return np.vstack(blocks)


def discrete_observability_gramian(
    A: np.ndarray, C: np.ndarray, *, horizon: int
) -> np.ndarray:
    """Return the discrete finite-horizon observability Gramian."""
    state = _square_matrix(A, "A")
    measurement = _measurement_matrix(C, state.shape[0])
    steps = int(horizon)
    if steps < 1 or steps != horizon:
        raise ValueError("horizon must be a positive integer.")
    gramian = np.zeros_like(state)
    block = measurement.copy()
    for _ in range(steps):
        gramian += block.T @ block
        block = block @ state
        if np.any(~np.isfinite(block)) or np.any(~np.isfinite(gramian)):
            raise FloatingPointError(
                "The discrete observability Gramian overflowed; shorten the horizon or rescale coordinates."
            )
    return 0.5 * (gramian + gramian.T)


def continuous_observability_gramian(
    A: np.ndarray, C: np.ndarray, *, horizon: float
) -> np.ndarray:
    """Return ``integral_0^T exp(A.T t) C.T C exp(A t) dt``.

    A Van Loan block exponential avoids time-grid quadrature and therefore
    provides a deterministic result for a requested physical/model horizon.
    """
    from scipy.linalg import expm

    state = _square_matrix(A, "A")
    measurement = _measurement_matrix(C, state.shape[0])
    duration = float(horizon)
    if not np.isfinite(duration) or duration <= 0.0:
        raise ValueError("horizon must be finite and positive.")
    forcing = measurement.T @ measurement
    dimension = state.shape[0]
    van_loan = np.block(
        [[-state.T, forcing], [np.zeros_like(state), state]]
    )
    transition = expm(van_loan * duration)
    gramian = expm(state.T * duration) @ transition[:dimension, dimension:]
    gramian = np.real_if_close(0.5 * (gramian + gramian.T), tol=1000)
    result = np.asarray(gramian, dtype=float)
    if np.any(~np.isfinite(result)):
        raise FloatingPointError(
            "The continuous observability Gramian overflowed; shorten the horizon or rescale coordinates."
        )
    return result


def continuous_infinite_observability_gramian(A: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Return the infinite-horizon Gramian for a strictly stable generator."""
    from scipy.linalg import solve_continuous_lyapunov

    state = _square_matrix(A, "A")
    measurement = _measurement_matrix(C, state.shape[0])
    poles = np.linalg.eigvals(state)
    if np.any(np.real(poles) >= 0.0):
        raise ValueError("An infinite-horizon observability Gramian requires a strictly stable A.")
    gramian = solve_continuous_lyapunov(state.T, -(measurement.T @ measurement))
    return np.asarray(np.real_if_close(0.5 * (gramian + gramian.T), tol=1000), dtype=float)


def _effective_rank(strengths: np.ndarray) -> float:
    energy = np.square(np.maximum(np.asarray(strengths, dtype=float), 0.0))
    total = float(np.sum(energy))
    if total <= np.finfo(float).tiny:
        return 0.0
    probabilities = energy / total
    nonzero = probabilities > 0.0
    return float(np.exp(-np.sum(probabilities[nonzero] * np.log(probabilities[nonzero]))))


@dataclass(frozen=True)
class ObservabilityAnalysis:
    """Immutable, Viewer-ready observability result in split KROM coordinates."""

    method: ObservabilityMethod
    A: np.ndarray
    C: np.ndarray
    gramian: np.ndarray
    singular_values: np.ndarray
    directions: np.ndarray
    numerical_rank: int
    tolerance: float
    horizon: float | int | None
    matrix: np.ndarray | None = None
    mode_participation: np.ndarray | None = None
    metadata: dict[str, Any] | None = None

    def __post_init__(self):
        state = _square_matrix(self.A, "A")
        measurement = _measurement_matrix(self.C, state.shape[0])
        gramian = _square_matrix(self.gramian, "gramian")
        strengths = np.asarray(self.singular_values, dtype=float).reshape(-1)
        directions = np.asarray(self.directions, dtype=float)
        if gramian.shape != state.shape or strengths.size != state.shape[0]:
            raise ValueError("Observability arrays do not match the state dimension.")
        if directions.shape != state.shape:
            raise ValueError("directions must contain observability directions as columns.")
        matrix = None if self.matrix is None else np.asarray(self.matrix, dtype=float)
        if matrix is not None and (matrix.ndim != 2 or matrix.shape[1] != state.shape[0]):
            raise ValueError("The classical observability matrix has the wrong width.")
        participation = None if self.mode_participation is None else np.asarray(
            self.mode_participation, dtype=float
        )
        if participation is not None and participation.ndim != 2:
            raise ValueError("mode_participation must be a two-dimensional array.")
        for item in (state, measurement, gramian, strengths, directions, matrix, participation):
            if item is not None:
                item.setflags(write=False)
        object.__setattr__(self, "A", state)
        object.__setattr__(self, "C", measurement)
        object.__setattr__(self, "gramian", gramian)
        object.__setattr__(self, "singular_values", strengths)
        object.__setattr__(self, "directions", directions)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "mode_participation", participation)
        object.__setattr__(self, "metadata", dict(self.metadata or {}))

    @property
    def state_dimension(self) -> int:
        return int(self.A.shape[0])

    @property
    def measurement_dimension(self) -> int:
        return int(self.C.shape[0])

    @property
    def condition_number(self) -> float:
        positive = self.singular_values[self.singular_values > self.tolerance]
        return float(np.inf if positive.size == 0 else positive[0] / positive[-1])

    @property
    def effective_rank(self) -> float:
        return _effective_rank(self.singular_values)

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "koopman-observability-analysis/v1",
            "method": self.method,
            "state_dimension": self.state_dimension,
            "measurement_dimension": self.measurement_dimension,
            "horizon": self.horizon,
            "numerical_rank": int(self.numerical_rank),
            "tolerance": float(self.tolerance),
            "condition_number": self.condition_number,
            "effective_rank": self.effective_rank,
            "singular_values": self.singular_values.tolist(),
            "directions_are_columns": True,
            "metadata": dict(self.metadata),
        }


def analyze_observability(
    A: np.ndarray,
    C: np.ndarray,
    *,
    method: ObservabilityMethod = "continuous-gramian",
    horizon: float | int | None = None,
    sample_time: float | None = None,
    coordinate_groups: tuple[tuple[int, ...], ...] | None = None,
    metadata: dict[str, Any] | None = None,
) -> ObservabilityAnalysis:
    """Analyze an explicit real state-space pair without hidden conversions."""
    from scipy.linalg import expm

    state = _square_matrix(A, "A")
    measurement = _measurement_matrix(C, state.shape[0])
    matrix = None
    if method == "matrix":
        steps = state.shape[0] if horizon is None else int(horizon)
        matrix = classical_observability_matrix(state, measurement, horizon=steps)
        _, strengths, right = np.linalg.svd(matrix, full_matrices=False)
        directions = right.T
        gramian = matrix.T @ matrix
        used_horizon: float | int | None = steps
    elif method == "discrete-gramian":
        steps = state.shape[0] if horizon is None else int(horizon)
        if sample_time is not None:
            dt = float(sample_time)
            if not np.isfinite(dt) or dt <= 0.0:
                raise ValueError("sample_time must be finite and positive.")
            operator = expm(state * dt)
        else:
            operator = state
        gramian = discrete_observability_gramian(operator, measurement, horizon=steps)
        used_horizon = steps
        eigenvalues, directions = np.linalg.eigh(gramian)
        order = np.argsort(eigenvalues)[::-1]
        strengths = np.sqrt(np.maximum(eigenvalues[order], 0.0))
        directions = directions[:, order]
    elif method == "continuous-gramian":
        duration = 1.0 if horizon is None else float(horizon)
        gramian = continuous_observability_gramian(state, measurement, horizon=duration)
        used_horizon = duration
        eigenvalues, directions = np.linalg.eigh(gramian)
        order = np.argsort(eigenvalues)[::-1]
        strengths = np.sqrt(np.maximum(eigenvalues[order], 0.0))
        directions = directions[:, order]
    elif method == "continuous-infinite-gramian":
        if horizon is not None:
            raise ValueError("horizon is not used for an infinite-horizon Gramian.")
        gramian = continuous_infinite_observability_gramian(state, measurement)
        used_horizon = None
        eigenvalues, directions = np.linalg.eigh(gramian)
        order = np.argsort(eigenvalues)[::-1]
        strengths = np.sqrt(np.maximum(eigenvalues[order], 0.0))
        directions = directions[:, order]
    else:
        raise ValueError(f"Unsupported observability method {method!r}.")
    largest = float(strengths[0]) if strengths.size else 0.0
    tolerance = float(max(measurement.shape[0], state.shape[0]) * np.finfo(float).eps * largest)
    rank = int(np.sum(strengths > tolerance))
    participation = None
    if coordinate_groups is not None:
        rows = []
        used = set()
        for group in coordinate_groups:
            indices = tuple(int(item) for item in group)
            if not indices or any(item < 0 or item >= state.shape[0] for item in indices):
                raise ValueError("coordinate_groups contains an invalid split-coordinate index.")
            if used.intersection(indices):
                raise ValueError("coordinate_groups must not overlap.")
            used.update(indices)
            rows.append(np.sum(np.square(directions[list(indices), :]), axis=0))
        participation = np.vstack(rows)
    return ObservabilityAnalysis(
        method, state, measurement, gramian, strengths, directions, rank,
        tolerance, used_horizon, matrix, participation, metadata
    )


def _krom_coordinate_groups(krom) -> tuple[tuple[int, ...], ...]:
    groups = []
    offset = 0
    for function in krom.eigenfunctions:
        width = 1 if function.is_real else 2
        groups.append(tuple(range(offset, offset + width)))
        offset += width
    return tuple(groups)


def analyze_krom_observability(
    krom,
    reconstruction=None,
    *,
    C: np.ndarray | None = None,
    method: ObservabilityMethod = "continuous-gramian",
    horizon: float | int | None = None,
    sample_time: float | None = None,
) -> ObservabilityAnalysis:
    """Analyze a KROM using a declared observable reconstruction or explicit C."""
    if reconstruction is not None and C is not None:
        raise ValueError("Pass reconstruction or C, not both.")
    selected = reconstruction if reconstruction is not None else getattr(
        krom, "observable_reconstruction", None
    )
    if isinstance(selected, str):
        selected = krom.observable_map(selected)
    if C is None:
        if selected is None:
            raise ValueError("KROM observability requires an ObservableReconstruction or explicit C.")
        C = selected.split_matrix(krom.eigenfunctions)
    observable = None if selected is None else selected.observable.summary()
    return analyze_observability(
        krom.A(output="split"), C, method=method, horizon=horizon,
        sample_time=sample_time, coordinate_groups=_krom_coordinate_groups(krom),
        metadata={
            "coordinate_system": getattr(krom, "coordinate_system", "model"),
            "observable": observable,
            "split_coordinate_groups": [list(group) for group in _krom_coordinate_groups(krom)],
        },
    )


def save_observability_analysis(path: str | Path, analysis: ObservabilityAnalysis) -> Path:
    """Save a self-contained, non-executable numerical analysis artifact."""
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        "A": analysis.A, "C": analysis.C, "gramian": analysis.gramian,
        "singular_values": analysis.singular_values, "directions": analysis.directions,
        "metadata_json": np.asarray(json.dumps(analysis.summary(), sort_keys=True)),
    }
    if analysis.matrix is not None:
        arrays["matrix"] = analysis.matrix
    if analysis.mode_participation is not None:
        arrays["mode_participation"] = analysis.mode_participation
    np.savez_compressed(target, **arrays)
    return target


def load_observability_analysis(path: str | Path) -> ObservabilityAnalysis:
    source = Path(path).expanduser().resolve()
    with np.load(source, allow_pickle=False) as data:
        summary = json.loads(str(data["metadata_json"]))
        return ObservabilityAnalysis(
            summary["method"], data["A"], data["C"], data["gramian"],
            data["singular_values"], data["directions"], summary["numerical_rank"],
            summary["tolerance"], summary.get("horizon"),
            None if "matrix" not in data.files else data["matrix"],
            None if "mode_participation" not in data.files else data["mode_participation"],
            summary.get("metadata") or {},
        )


__all__ = [
    "ObservabilityAnalysis", "ObservabilityMethod", "analyze_observability",
    "analyze_krom_observability", "classical_observability_matrix",
    "continuous_infinite_observability_gramian", "continuous_observability_gramian",
    "discrete_observability_gramian", "load_observability_analysis",
    "save_observability_analysis",
]
