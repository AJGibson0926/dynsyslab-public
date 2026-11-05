"""Named Koopman control-input representations retained by a KROM."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np


TimeDomain = Literal["continuous", "discrete"]


def complex_to_split(matrix: np.ndarray, eigenfunctions) -> np.ndarray:
    """Convert rows aligned to nonredundant complex coordinates to real rows."""
    values = np.asarray(matrix, dtype=complex)
    if values.ndim != 2 or values.shape[0] != len(eigenfunctions):
        raise ValueError("Control matrix rows must align with the KROM eigenfunctions.")
    rows = []
    for index, function in enumerate(eigenfunctions):
        row = values[index]
        rows.append(row.real[None, :] if function.is_real else np.vstack((row.real, row.imag)))
    return np.vstack(rows)


@dataclass(frozen=True)
class ConstantControlMatrix:
    """State-independent lifted input matrix identified from controlled data."""

    matrix: np.ndarray
    time_domain: TimeDomain = "continuous"
    delta_t: float | None = None
    provenance: dict[str, Any] | None = None
    discrete_matrix: np.ndarray | None = None

    def __post_init__(self):
        matrix = np.array(self.matrix, dtype=complex, copy=True)
        if matrix.ndim != 2 or not matrix.size or np.any(~np.isfinite(matrix)):
            raise ValueError("Constant Koopman B must be a finite nonempty matrix.")
        if self.time_domain not in {"continuous", "discrete"}:
            raise ValueError("time_domain must be 'continuous' or 'discrete'.")
        if self.time_domain == "discrete" and (
            self.delta_t is None or self.delta_t <= 0 or not np.isfinite(self.delta_t)
        ):
            raise ValueError("A discrete constant B requires a positive finite delta_t.")
        discrete = None
        if self.discrete_matrix is not None:
            discrete = np.array(self.discrete_matrix, dtype=complex, copy=True)
            if discrete.shape != matrix.shape or np.any(~np.isfinite(discrete)):
                raise ValueError("discrete_matrix must match matrix shape and be finite.")
            discrete.setflags(write=False)
        matrix.setflags(write=False)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "discrete_matrix", discrete)
        object.__setattr__(self, "provenance", dict(self.provenance or {}))

    @property
    def lifted_dimension(self) -> int:
        return int(self.matrix.shape[0])

    @property
    def input_dimension(self) -> int:
        return int(self.matrix.shape[1])

    def evaluate(self, lifted_state: np.ndarray | None = None) -> np.ndarray:
        return np.asarray(self.matrix)

    def subset(self, positions) -> "ConstantControlMatrix":
        selected = np.asarray(positions, dtype=int)
        discrete = None if self.discrete_matrix is None else self.discrete_matrix[selected]
        return ConstantControlMatrix(
            self.matrix[selected], self.time_domain, self.delta_t,
            {**self.provenance, "subset_positions": selected.tolist()}, discrete,
        )

    def summary(self) -> dict[str, Any]:
        return {
            "kind": "constant",
            "encoding": "constant-koopman-input/v1",
            "shape": list(self.matrix.shape),
            "time_domain": self.time_domain,
            "delta_t": self.delta_t,
            "has_discrete_matrix": self.discrete_matrix is not None,
            "provenance": dict(self.provenance),
        }


@dataclass(frozen=True)
class BilinearControlMatrix:
    """Closed affine input map on a nonredundant complex lift.

    The standard form is ``M[j] psi + b[:,j]``.  ``conjugate_operators``
    optionally adds ``N[j] conjugate(psi)``; this is exactly a general real-
    linear bilinear map in real/imaginary split coordinates and corresponds to
    retaining both members of every conjugate eigenfunction pair explicitly.
    """

    operators: np.ndarray
    offsets: np.ndarray
    time_domain: TimeDomain = "continuous"
    delta_t: float | None = None
    provenance: dict[str, Any] | None = None
    conjugate_operators: np.ndarray | None = None

    def __post_init__(self):
        operators = np.array(self.operators, dtype=complex, copy=True)
        offsets = np.array(self.offsets, dtype=complex, copy=True)
        conjugate = None if self.conjugate_operators is None else np.array(
            self.conjugate_operators, dtype=complex, copy=True
        )
        if operators.ndim != 3:
            raise ValueError("Bilinear operators must have shape (inputs, modes, modes).")
        inputs, rows, columns = operators.shape
        if rows != columns or offsets.shape != (rows, inputs):
            raise ValueError("Bilinear offsets/operators have incompatible dimensions.")
        if np.any(~np.isfinite(operators)) or np.any(~np.isfinite(offsets)):
            raise ValueError("Bilinear Koopman input arrays must be finite.")
        if conjugate is not None and (
            conjugate.shape != operators.shape or np.any(~np.isfinite(conjugate))
        ):
            raise ValueError("Conjugate-bilinear operators must match standard operators and be finite.")
        if self.time_domain not in {"continuous", "discrete"}:
            raise ValueError("time_domain must be 'continuous' or 'discrete'.")
        if self.time_domain == "discrete" and (
            self.delta_t is None or self.delta_t <= 0 or not np.isfinite(self.delta_t)
        ):
            raise ValueError("A discrete bilinear B requires a positive finite delta_t.")
        operators.setflags(write=False); offsets.setflags(write=False)
        if conjugate is not None: conjugate.setflags(write=False)
        object.__setattr__(self, "operators", operators)
        object.__setattr__(self, "offsets", offsets)
        object.__setattr__(self, "conjugate_operators", conjugate)
        object.__setattr__(self, "provenance", dict(self.provenance or {}))

    @property
    def lifted_dimension(self) -> int:
        return int(self.offsets.shape[0])

    @property
    def input_dimension(self) -> int:
        return int(self.offsets.shape[1])

    def evaluate(self, lifted_state: np.ndarray) -> np.ndarray:
        psi = np.asarray(lifted_state, dtype=complex).reshape(-1)
        if psi.size != self.lifted_dimension:
            raise ValueError("Bilinear input evaluation has the wrong lifted dimension.")
        result = self.offsets.copy()
        for input_index in range(self.input_dimension):
            result[:, input_index] += self.operators[input_index] @ psi
            if self.conjugate_operators is not None:
                result[:, input_index] += self.conjugate_operators[input_index] @ np.conjugate(psi)
        return result

    def subset(self, positions) -> "BilinearControlMatrix":
        """Return an independent closed principal restriction.

        Both controlled output rows and lifted predictor columns are restricted
        to the selected coordinates.  The result has no dependency on the
        source KROM.  Re-regression, when desired, is a separate explicit step.
        """
        selected = np.asarray(positions, dtype=int)
        if selected.ndim != 1 or not selected.size or len(np.unique(selected)) != selected.size:
            raise ValueError("Bilinear subset positions must be a nonempty unique vector.")
        if np.any(selected < 0) or np.any(selected >= self.lifted_dimension):
            raise IndexError("Bilinear subset position is outside the source representation.")
        operators = self.operators[:, selected][:, :, selected]
        conjugate = (None if self.conjugate_operators is None else
                     self.conjugate_operators[:, selected][:, :, selected])
        offsets = self.offsets[selected]
        return BilinearControlMatrix(
            operators, offsets, self.time_domain, self.delta_t,
            {**self.provenance, "subset_positions": selected.tolist(),
             "inheritance_mode": "independent-closed-principal-restriction"},
            conjugate,
        )

    def summary(self) -> dict[str, Any]:
        return {
            "kind": "bilinear",
            "encoding": "widely-linear-bilinear-koopman-input/v1" if self.conjugate_operators is not None else "bilinear-koopman-input/v1",
            "operator_shape": list(self.operators.shape),
            "offset_shape": list(self.offsets.shape),
            "has_conjugate_operators": self.conjugate_operators is not None,
            "time_domain": self.time_domain,
            "delta_t": self.delta_t,
            "provenance": dict(self.provenance),
        }


__all__ = [
    "BilinearControlMatrix", "ConstantControlMatrix", "TimeDomain", "complex_to_split"
]
