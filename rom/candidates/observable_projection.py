"""Candidate full-spectrum observable projection for parity experiments.

This module deliberately does not alter :class:`rom.KROM`.  The accepted
nonredundant KROM path remains :func:`rom.fit_observable_reconstruction`.
These routines cover MATLAB-style full-complex arrays where ``Y = C Psi``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def _readonly(value: np.ndarray, *, dtype) -> np.ndarray:
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _metrics(actual: np.ndarray, predicted: np.ndarray) -> dict[str, Any]:
    residual = predicted - actual
    real_residual = predicted.real - actual.real
    centered = actual.real - np.mean(actual.real, axis=0, keepdims=True)
    denominator = np.sum(centered * centered, axis=0)
    numerator = np.sum(real_residual * real_residual, axis=0)
    r2 = np.where(
        denominator > np.finfo(float).eps,
        1.0 - numerator / denominator,
        np.where(numerator <= np.finfo(float).eps, 1.0, -np.inf),
    )
    return {
        "sample_count": int(actual.shape[0]),
        "complex_rmse_by_output": np.sqrt(np.mean(np.abs(residual) ** 2, axis=0)).tolist(),
        "real_rmse_by_output": np.sqrt(np.mean(real_residual**2, axis=0)).tolist(),
        "r2_by_output": r2.tolist(),
        "maximum_imaginary_leakage_by_output": np.max(np.abs(predicted.imag), axis=0).tolist(),
    }


@dataclass(frozen=True)
class ObservableProjectionCandidate:
    """Unpromoted affine ``Y = Psi C.T + intercept`` candidate."""

    coefficients: np.ndarray
    intercept: np.ndarray
    diagnostics: dict[str, Any]

    def __post_init__(self):
        coefficients = _readonly(self.coefficients, dtype=complex)
        intercept = _readonly(self.intercept, dtype=complex).reshape(-1)
        if coefficients.ndim != 2 or coefficients.shape[0] != intercept.size:
            raise ValueError("Projection C must have shape (outputs, lifted modes).")
        if not coefficients.size or np.any(~np.isfinite(coefficients)) or np.any(~np.isfinite(intercept)):
            raise ValueError("Projection C and intercept must be finite and nonempty.")
        object.__setattr__(self, "coefficients", coefficients)
        object.__setattr__(self, "intercept", intercept)
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    @property
    def output_dimension(self) -> int:
        return int(self.coefficients.shape[0])

    @property
    def lifted_dimension(self) -> int:
        return int(self.coefficients.shape[1])

    def reconstruct(self, lifted_samples: np.ndarray) -> np.ndarray:
        values = np.asarray(lifted_samples, dtype=complex)
        one = values.ndim == 1
        if one:
            values = values.reshape(1, -1)
        if values.ndim != 2 or values.shape[1] != self.lifted_dimension:
            raise ValueError("Lifted samples must be sample-major and align with projection C.")
        result = values @ self.coefficients.T + self.intercept
        return result[0] if one else result

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "candidate-full-complex-observable-projection/v1",
            "candidate_only": True,
            "output_dimension": self.output_dimension,
            "lifted_dimension": self.lifted_dimension,
            "diagnostics": dict(self.diagnostics),
        }


def fit_observable_projection_candidate(
    lifted_samples: np.ndarray,
    observable_samples: np.ndarray,
    *,
    validation_lifted_samples: np.ndarray | None = None,
    validation_observable_samples: np.ndarray | None = None,
    fit_intercept: bool = False,
    ridge: float = 0.0,
) -> ObservableProjectionCandidate:
    """Fit the MATLAB-style projection ``Y = C Psi`` without promotion.

    Arrays are sample-major.  If no independent validation arrays are given,
    every fifth sample is held out deterministically.  ``fit_intercept=False``
    reproduces ``C = Y * pinv(Psi)`` (up to the explicitly chosen ridge).
    """
    lifted = np.asarray(lifted_samples, dtype=complex)
    observed = np.asarray(observable_samples, dtype=complex)
    if observed.ndim == 1:
        observed = observed[:, None]
    if lifted.ndim != 2 or observed.ndim != 2 or lifted.shape[0] != observed.shape[0]:
        raise ValueError("Lifted and observable samples must be aligned sample-major matrices.")
    if lifted.shape[0] < 3 or lifted.shape[1] < 1 or observed.shape[1] < 1:
        raise ValueError("Observable projection requires at least three samples and nonempty dimensions.")
    if np.any(~np.isfinite(lifted)) or np.any(~np.isfinite(observed)):
        raise ValueError("Observable projection samples must be finite.")
    if ridge < 0.0 or not np.isfinite(ridge):
        raise ValueError("ridge must be finite and nonnegative.")
    if (validation_lifted_samples is None) != (validation_observable_samples is None):
        raise ValueError("Both independent validation arrays must be provided together.")

    if validation_lifted_samples is None:
        indices = np.arange(lifted.shape[0])
        held_mask = indices % 5 == 0
        if np.sum(~held_mask) < 2:
            raise ValueError("Observable projection leaves too few fit samples.")
        fit_lifted = lifted[~held_mask]
        fit_observed = observed[~held_mask]
        validation_lifted = lifted[held_mask]
        validation_observed = observed[held_mask]
        validation_source = "deterministic-interleaved-20-percent"
    else:
        fit_lifted = lifted
        fit_observed = observed
        validation_lifted = np.asarray(validation_lifted_samples, dtype=complex)
        validation_observed = np.asarray(validation_observable_samples, dtype=complex)
        if validation_observed.ndim == 1:
            validation_observed = validation_observed[:, None]
        if (
            validation_lifted.ndim != 2
            or validation_observed.ndim != 2
            or validation_lifted.shape[0] != validation_observed.shape[0]
            or validation_lifted.shape[1] != lifted.shape[1]
            or validation_observed.shape[1] != observed.shape[1]
        ):
            raise ValueError("Independent validation arrays do not align with the fit arrays.")
        validation_source = "independent"

    if fit_intercept:
        design = np.column_stack((np.ones(fit_lifted.shape[0], dtype=complex), fit_lifted))
        penalty = np.diag(np.r_[0.0, np.ones(fit_lifted.shape[1])])
    else:
        design = fit_lifted
        penalty = np.eye(fit_lifted.shape[1])
    gram = design.conj().T @ design + ridge * penalty
    if ridge == 0.0:
        solution = np.linalg.pinv(design, rcond=1e-12) @ fit_observed
        solve_method = "svd-pseudoinverse"
    else:
        solution = np.linalg.pinv(gram, rcond=1e-12) @ design.conj().T @ fit_observed
        solve_method = "ridge-normal-equations"
    if fit_intercept:
        intercept = solution[0]
        coefficients = solution[1:].T
    else:
        intercept = np.zeros(observed.shape[1], dtype=complex)
        coefficients = solution.T

    provisional = ObservableProjectionCandidate(coefficients, intercept, {})
    diagnostics = {
        "method": "candidate-full-complex-observable-projection/v1",
        "candidate_only": True,
        "fit_intercept": bool(fit_intercept),
        "solve_method": solve_method,
        "ridge": float(ridge),
        "design_shape": list(design.shape),
        "design_rank": int(np.linalg.matrix_rank(design)),
        "condition_number": float(np.linalg.cond(gram)),
        "validation_source": validation_source,
        "fit": _metrics(fit_observed, provisional.reconstruct(fit_lifted)),
        "validation": _metrics(
            validation_observed, provisional.reconstruct(validation_lifted)
        ),
        "promotion": "none; caller must explicitly adapt this full-spectrum candidate",
    }
    return ObservableProjectionCandidate(coefficients, intercept, diagnostics)


__all__ = ["ObservableProjectionCandidate", "fit_observable_projection_candidate"]
