from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

@dataclass
class PCAResult:
    method: str
    mean: np.ndarray
    components: np.ndarray
    singular_values: np.ndarray
    scores: np.ndarray
    explained_variance_ratio: np.ndarray
    metadata: dict=field(default_factory=dict)

@dataclass
class SVDResult:
    """Economy SVD of a samples-by-features data matrix."""
    left_singular_vectors: np.ndarray
    singular_values: np.ndarray
    right_singular_vectors: np.ndarray
    mean: np.ndarray
    metadata: dict=field(default_factory=dict)
    pod_modes: object | None=None
    @property
    def rank(self) -> int: return int(self.singular_values.size)
    @property
    def prom(self): return self.pod_modes

@dataclass
class HAVOKResult:
    """A delay-coordinate linear model with an identified forcing signal."""
    A: np.ndarray
    B: np.ndarray
    forcing: np.ndarray
    coordinates: np.ndarray
    coordinate_time: np.ndarray
    left_singular_vectors: np.ndarray
    singular_values: np.ndarray
    right_singular_vectors: np.ndarray
    observable: np.ndarray
    metadata: dict=field(default_factory=dict)
    prom: object | None=None
    @property
    def state_coordinates(self) -> np.ndarray: return self.coordinates[:, :-1]
