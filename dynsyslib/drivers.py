"""Forcing terms that can be attached to library systems."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class SinusoidalDriver:
    frequency_hz: float
    amplitude: float
    phase_rad: float = 0.0

    @property
    def period(self) -> float:
        return 1.0 / self.frequency_hz

    def value(self, t: float) -> float:
        return self.amplitude * np.cos(2.0 * np.pi * self.frequency_hz * t + self.phase_rad)

@dataclass(frozen=True)
class ZeroDriver:
    """Explicit free-evolution driver compatible with every driven library model."""
    @property
    def period(self) -> float: return float("inf")
    def value(self,t: float) -> float: return 0.

@dataclass(frozen=True)
class TabulatedDriver:
    """Piecewise-linear, archive-serializable scalar input time series."""
    time: np.ndarray
    values: np.ndarray
    def __post_init__(self):
        object.__setattr__(self,"time",np.asarray(self.time,dtype=float).reshape(-1)); object.__setattr__(self,"values",np.asarray(self.values,dtype=float).reshape(-1))
        if self.time.size<2 or self.time.size!=self.values.size or np.any(np.diff(self.time)<=0): raise ValueError("TabulatedDriver requires matching strictly increasing time and value arrays.")
    @property
    def period(self) -> float: return float("inf")
    def value(self,t: float) -> float: return float(np.interp(t,self.time,self.values,left=self.values[0],right=self.values[-1]))
