"""Matrix-norm error fields for opt-in Koopman-control approximations."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def _readonly(value: np.ndarray, *, dtype) -> np.ndarray:
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class ControlApproximationErrorField:
    """Absolute and relative Frobenius errors on a rectangular phase grid."""

    x_values: np.ndarray
    y_values: np.ndarray
    exact_norm: np.ndarray
    absolute_error: np.ndarray
    relative_error: np.ndarray
    approximation_kind: str
    inside_training_hull: np.ndarray | None
    metadata: dict[str, Any]

    def __post_init__(self):
        x_values = _readonly(self.x_values, dtype=float).reshape(-1)
        y_values = _readonly(self.y_values, dtype=float).reshape(-1)
        expected = (y_values.size, x_values.size)
        exact_norm = _readonly(self.exact_norm, dtype=float)
        absolute = _readonly(self.absolute_error, dtype=float)
        relative = _readonly(self.relative_error, dtype=float)
        if exact_norm.shape != expected or absolute.shape != expected or relative.shape != expected:
            raise ValueError("Error-field arrays must align with the y-by-x rectangular grid.")
        if np.any(~np.isfinite(exact_norm)) or np.any(~np.isfinite(absolute)) or np.any(~np.isfinite(relative)):
            raise ValueError("Error-field arrays must be finite.")
        mask = None
        if self.inside_training_hull is not None:
            mask = _readonly(self.inside_training_hull, dtype=bool)
            if mask.shape != expected:
                raise ValueError("Training-hull mask must align with the error grid.")
        object.__setattr__(self, "x_values", x_values)
        object.__setattr__(self, "y_values", y_values)
        object.__setattr__(self, "exact_norm", exact_norm)
        object.__setattr__(self, "absolute_error", absolute)
        object.__setattr__(self, "relative_error", relative)
        object.__setattr__(self, "inside_training_hull", mask)
        object.__setattr__(self, "metadata", dict(self.metadata))

    def summary(self) -> dict[str, Any]:
        mask = np.ones_like(self.absolute_error, dtype=bool)
        if self.inside_training_hull is not None:
            mask = self.inside_training_hull
        return {
            "encoding": "candidate-control-approximation-error-field/v1",
            "candidate_only": True,
            "approximation_kind": self.approximation_kind,
            "grid_shape": list(self.absolute_error.shape),
            "rectangle": {
                "x": [float(self.x_values[0]), float(self.x_values[-1])],
                "y": [float(self.y_values[0]), float(self.y_values[-1])],
            },
            "absolute": _distribution_summary(self.absolute_error),
            "relative": _distribution_summary(self.relative_error),
            "inside_training_hull": {
                "point_count": int(np.sum(mask)),
                "absolute": _distribution_summary(self.absolute_error[mask]),
                "relative": _distribution_summary(self.relative_error[mask]),
            },
            "metadata": dict(self.metadata),
        }


def _distribution_summary(values: np.ndarray) -> dict[str, float | int]:
    flattened = np.asarray(values, dtype=float).reshape(-1)
    return {
        "count": int(flattened.size),
        "minimum": float(np.min(flattened)),
        "median": float(np.median(flattened)),
        "mean": float(np.mean(flattened)),
        "p90": float(np.quantile(flattened, 0.90)),
        "p95": float(np.quantile(flattened, 0.95)),
        "maximum": float(np.max(flattened)),
    }


def compare_control_approximations_on_grid(
    x_values: np.ndarray,
    y_values: np.ndarray,
    exact_Bpsi: np.ndarray,
    *,
    lifted_values: np.ndarray | None = None,
    constant_matrix: np.ndarray | None = None,
    bilinear_operators: np.ndarray | None = None,
    bilinear_offsets: np.ndarray | None = None,
    inside_training_hull: np.ndarray | None = None,
    relative_floor: float | None = None,
    metadata: dict[str, Any] | None = None,
) -> dict[str, ControlApproximationErrorField]:
    """Compare constant and/or bilinear matrices to exact ``B_psi(x)``.

    ``exact_Bpsi`` is ``(ny,nx,modes,inputs)``.  Bilinear lifting is
    ``(ny,nx,modes)`` and operators are ``(inputs,modes,modes)``.
    Errors use the Frobenius matrix norm at every phase-space point.
    """
    x = np.asarray(x_values, dtype=float).reshape(-1)
    y = np.asarray(y_values, dtype=float).reshape(-1)
    exact = np.asarray(exact_Bpsi, dtype=complex)
    expected_prefix = (y.size, x.size)
    if exact.ndim != 4 or exact.shape[:2] != expected_prefix:
        raise ValueError("exact_Bpsi must have shape (len(y), len(x), modes, inputs).")
    if np.any(~np.isfinite(exact)):
        raise ValueError("exact_Bpsi must be finite.")
    modes, inputs = exact.shape[-2:]
    exact_norm = np.linalg.norm(exact, axis=(-2, -1))
    scale = float(np.max(exact_norm))
    floor = (
        np.sqrt(np.finfo(float).eps) * max(scale, 1.0)
        if relative_floor is None
        else float(relative_floor)
    )
    if floor <= 0.0 or not np.isfinite(floor):
        raise ValueError("relative_floor must be finite and positive.")
    mask = None if inside_training_hull is None else np.asarray(inside_training_hull, dtype=bool)
    common = {
        **dict(metadata or {}),
        "norm": "frobenius",
        "relative_denominator": "max(norm(exact_Bpsi), relative_floor)",
        "relative_floor": floor,
    }
    results: dict[str, ControlApproximationErrorField] = {}

    if constant_matrix is not None:
        constant = np.asarray(constant_matrix, dtype=complex)
        if constant.shape != (modes, inputs):
            raise ValueError("constant_matrix must align with exact B_psi matrix dimensions.")
        approximation = np.broadcast_to(constant, exact.shape)
        results["constant"] = _make_field(
            x, y, exact, approximation, "constant", mask, common
        )

    if (bilinear_operators is None) != (bilinear_offsets is None):
        raise ValueError("Bilinear operators and offsets must be provided together.")
    if bilinear_operators is not None:
        if lifted_values is None:
            raise ValueError("Bilinear comparison requires lifted_values on the phase grid.")
        lifted = np.asarray(lifted_values, dtype=complex)
        operators = np.asarray(bilinear_operators, dtype=complex)
        offsets = np.asarray(bilinear_offsets, dtype=complex)
        if lifted.shape != expected_prefix + (modes,):
            raise ValueError("lifted_values must have shape (len(y), len(x), modes).")
        if operators.shape != (inputs, modes, modes) or offsets.shape != (modes, inputs):
            raise ValueError("Bilinear arrays do not align with exact B_psi dimensions.")
        approximation = np.einsum("imn,yxn->yxmi", operators, lifted) + offsets
        results["bilinear"] = _make_field(
            x, y, exact, approximation, "bilinear", mask, common
        )

    if not results:
        raise ValueError("At least one constant or bilinear approximation is required.")
    return results


def _make_field(
    x: np.ndarray,
    y: np.ndarray,
    exact: np.ndarray,
    approximation: np.ndarray,
    kind: str,
    mask: np.ndarray | None,
    metadata: dict[str, Any],
) -> ControlApproximationErrorField:
    exact_norm = np.linalg.norm(exact, axis=(-2, -1))
    absolute = np.linalg.norm(exact - approximation, axis=(-2, -1))
    floor = float(metadata["relative_floor"])
    relative = absolute / np.maximum(exact_norm, floor)
    return ControlApproximationErrorField(
        x,
        y,
        exact_norm,
        absolute,
        relative,
        kind,
        mask,
        metadata,
    )


__all__ = ["ControlApproximationErrorField", "compare_control_approximations_on_grid"]
