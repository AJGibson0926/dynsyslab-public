"""Exact evaluators shared by in-memory and archived spectral representations."""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np


def polynomial_eigenfunction_values(states: np.ndarray, coefficients: np.ndarray, *,
                                    dictionary: dict[str,Any],
                                    state_preprocessing: dict[str,Any],
                                    columns: Sequence[int] | np.ndarray | None=None) -> np.ndarray:
    """Evaluate serialized polynomial eigenfunctions, optionally by column.

    The returned columns correspond to coefficient rows.  This is the sole
    evaluator for the versioned ``polynomial-dictionary-evaluation/v1`` archive
    encoding, preventing publication and Viewer from drifting numerically.
    """
    if dictionary.get("encoding")!="polynomial-dictionary/v1":
        raise ValueError("Polynomial eigenfunction evaluation requires a polynomial-dictionary/v1 specification.")
    x=np.asarray(states,dtype=float); coefficient_matrix=np.asarray(coefficients,dtype=complex)
    if x.ndim!=2 or coefficient_matrix.ndim!=2:
        raise ValueError("Polynomial eigenfunction states and coefficients must be two-dimensional.")
    indices=np.arange(coefficient_matrix.shape[0],dtype=int) if columns is None else np.asarray(columns,dtype=int).reshape(-1)
    if np.any(indices<0) or np.any(indices>=coefficient_matrix.shape[0]):
        raise IndexError("Polynomial eigenfunction column index is out of range.")
    required={"center","scale","multiplier"}
    if not required.issubset(state_preprocessing):
        raise ValueError("Polynomial eigenfunction evaluation requires complete state-preprocessing metadata.")
    center=np.asarray(state_preprocessing["center"],dtype=float)
    scale=np.asarray(state_preprocessing["scale"],dtype=float)
    multiplier=np.asarray(state_preprocessing["multiplier"],dtype=float)
    if center.shape!=scale.shape or center.shape!=multiplier.shape or center.shape!=(x.shape[1],):
        raise ValueError("Polynomial eigenfunction state-preprocessing vectors do not match the state dimension.")
    if np.any(scale<=0) or np.any(multiplier<=0) or not np.all(np.isfinite(center)) or not np.all(np.isfinite(scale)) or not np.all(np.isfinite(multiplier)):
        raise ValueError("Polynomial eigenfunction state preprocessing must be finite with positive scales and multipliers.")
    transformed=(x-center)/scale*multiplier
    from .dictionary import PolynomialDictionary
    features=PolynomialDictionary(degree=int(dictionary["degree"]),
        include_constant=bool(dictionary.get("include_constant",True))).transform(transformed)
    if features.shape[1]!=coefficient_matrix.shape[1]:
        raise ValueError("Serialized polynomial dictionary and coefficient feature counts are inconsistent.")
    return np.asarray(features,dtype=complex)@coefficient_matrix[indices].T


__all__=["polynomial_eigenfunction_values"]
