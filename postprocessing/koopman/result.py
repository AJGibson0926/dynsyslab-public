from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
import numpy as np
from .spectrum import KoopmanSpectrum

@dataclass
class KoopmanResult:
    method: str
    eigenvalues: np.ndarray
    modes: np.ndarray | None
    eigenfunction_values: np.ndarray | None
    operator: np.ndarray
    metadata: dict = field(default_factory=dict)
    spectrum: KoopmanSpectrum | None = None
    krom: Any | None = None
    modal_amplitudes: np.ndarray | None = None
