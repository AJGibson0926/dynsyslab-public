"""Estimator interfaces independent of solvers, controllers, and system libraries."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable
import numpy as np

from infrastructure.estimation import Estimator, IdentityEstimator

@dataclass(frozen=True)
class StaticMapEstimator:
    model: Callable[[np.ndarray],np.ndarray]
    def estimate(self,t,true_state): return np.asarray(self.model(np.asarray(true_state,dtype=float)),dtype=float)

@dataclass
class LinearLuenbergerEstimator:
    """Continuous observer updated on increasing solver evaluation times."""
    A: np.ndarray; C: np.ndarray; L: np.ndarray; initial_estimate: np.ndarray
    _state: np.ndarray=field(init=False); _time: float | None=field(init=False,default=None)
    requires_sample_clock: bool=field(init=False,default=True,repr=False)
    def __post_init__(self): self._state=np.asarray(self.initial_estimate,dtype=float).copy()
    def estimate(self,t,true_state):
        observation=np.asarray(true_state,dtype=float); current=float(t)
        if self._time is not None and current>self._time:
            dt=current-self._time; self._state=self._state+dt*(np.asarray(self.A)@self._state+np.asarray(self.L)@(observation-np.asarray(self.C)@self._state))
        self._time=current if self._time is None or current>self._time else self._time
        return self._state.copy()

@dataclass
class LQEEstimator:
    """Continuous-time linear-quadratic estimator with solver-supplied measurements.

    The solver currently supplies the plant state as its measurement source;
    ``C`` selects or maps that state to the measurement vector used by the
    observer.  The steady-state Kalman/LQE gain is obtained from the continuous
    algebraic Riccati equation defined by ``A, C, Q, R``.
    """
    A: np.ndarray; C: np.ndarray; Q: np.ndarray; R: np.ndarray; initial_estimate: np.ndarray
    _state: np.ndarray=field(init=False); _time: float | None=field(init=False,default=None); gain: np.ndarray=field(init=False)
    requires_sample_clock: bool=field(init=False,default=True,repr=False)
    def __post_init__(self):
        from scipy.linalg import solve_continuous_are
        A=np.asarray(self.A,dtype=float); C=np.asarray(self.C,dtype=float); Q=np.asarray(self.Q,dtype=float); R=np.asarray(self.R,dtype=float)
        covariance=solve_continuous_are(A.T,C.T,Q,R)
        self.gain=np.linalg.solve(R,C@covariance).T
        self._state=np.asarray(self.initial_estimate,dtype=float).copy()
    def estimate(self,t,true_state):
        current=float(t); observation=np.asarray(self.C,dtype=float)@np.asarray(true_state,dtype=float)
        if self._time is not None and current>self._time:
            dt=current-self._time; A=np.asarray(self.A,dtype=float); C=np.asarray(self.C,dtype=float)
            self._state=self._state+dt*(A@self._state+self.gain@(observation-C@self._state))
        self._time=current if self._time is None or current>self._time else self._time
        return self._state.copy()
