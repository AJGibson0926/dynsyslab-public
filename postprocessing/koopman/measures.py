"""Carefully labelled finite-data spectral-weight estimates."""
from __future__ import annotations
import numpy as np
from .result import KoopmanResult

def empirical_spectral_weight(result: KoopmanResult, eigenvalue_index: int) -> np.ndarray:
    """Return a normalized empirical weight, not a rigorously identified eigenmeasure.

    It weights each observed sample by squared magnitude of an estimated Koopman
    eigenfunction.  This is useful for exploratory spectral localization only.
    """
    if result.eigenfunction_values is None: raise ValueError("Eigenfunction samples are unavailable.")
    values=np.abs(result.eigenfunction_values[:,eigenvalue_index])**2; total=values.sum()
    return values/total if total else np.full(values.size,1/values.size)
