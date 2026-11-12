"""POD--Galerkin reduced-order models, independent of Koopman machinery."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass
class PODGalerkinROM:
    """A modal Galerkin ROM with explicit POD modes and reference state.

    This class is deliberately separate from :class:`KROM`: POD modal
    coordinates and their reduced operators need not be Koopman observables.
    """
    modes: np.ndarray
    reference_state: np.ndarray
    A: np.ndarray
    B: np.ndarray | None=None
    method: str="pod-galerkin"

    def __post_init__(self):
        self.modes=np.asarray(self.modes,dtype=float)
        self.reference_state=np.asarray(self.reference_state,dtype=float).reshape(-1)
        self.A=np.asarray(self.A,dtype=float)
        self.B=None if self.B is None else np.asarray(self.B,dtype=float)
        if self.modes.ndim!=2 or self.modes.shape[0]!=self.reference_state.size:
            raise ValueError("POD modes must have one row per physical-state component.")
        if self.A.shape!=(self.modes.shape[1],self.modes.shape[1]):
            raise ValueError("PODGalerkin A must be square in the retained modal dimension.")
        if self.B is not None and self.B.shape[0]!=self.modes.shape[1]:
            raise ValueError("PODGalerkin B rows must match the retained modal dimension.")

    def project(self,state: np.ndarray) -> np.ndarray:
        return self.modes.T@(np.asarray(state,dtype=float).reshape(-1)-self.reference_state)

    def reconstruct(self,coordinates: np.ndarray) -> np.ndarray:
        return self.reference_state+self.modes@np.asarray(coordinates,dtype=float).reshape(-1)
