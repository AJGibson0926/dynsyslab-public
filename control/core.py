"""Classical targets and controllers; independent of solver and system libraries."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, ClassVar, Protocol
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_continuous_are
from time_series import TimeSeries

class Target(Protocol):
    """State-reference protocol shared by feedback controllers."""
    def evaluate(self,t: float,state_dimension: int) -> np.ndarray: ...

@dataclass(frozen=True)
class ConstantTarget:
    """Time-invariant state target."""
    state: np.ndarray
    def evaluate(self,t,state_dimension):
        value=np.asarray(self.state,dtype=float).reshape(-1)
        if value.size==1: value=np.full(state_dimension,value.item())
        if value.size!=state_dimension: raise ValueError("Target state dimension does not match the controlled state.")
        return value

@dataclass(frozen=True)
class CallableTarget:
    """Target backed by a caller-supplied time function."""
    function: Callable[[float],np.ndarray]
    def evaluate(self,t,state_dimension): return ConstantTarget(self.function(t)).evaluate(t,state_dimension)

@dataclass(frozen=True)
class TimeSeriesTarget:
    """Target interpolated from a supplied state time series."""
    series: TimeSeries
    interpolation: str="pchip"
    def evaluate(self,t,state_dimension):
        if self.series.values.ndim!=2: raise ValueError("A target TimeSeries must have values shaped (time, state).")
        if self.series.values.shape[1]!=state_dimension: raise ValueError("Target state dimension does not match the controlled state.")
        if t<self.series.time[0] or t>self.series.time[-1]: raise ValueError("Controller target was requested outside its time range.")
        if self.interpolation=="pchip": return np.asarray(PchipInterpolator(self.series.time,self.series.values,axis=0)(t))
        if self.interpolation=="linear": return np.asarray([np.interp(t,self.series.time,self.series.values[:,i]) for i in range(state_dimension)])
        raise ValueError("Target interpolation must be 'pchip' or 'linear'.")

@dataclass
class WaveformController:
    """Multi-tone waveform.  Default actuation replaces a system driver input."""
    frequencies_hz: np.ndarray
    amplitudes: np.ndarray
    phases_rad: np.ndarray | None=None
    offset: float=0.
    actuation: str="system_input"  # system_input or additive_state
    output_matrix: np.ndarray | None=None
    active_time_window: tuple[float,float] | None=None
    archive_replay_safe: ClassVar[bool]=True
    coordinate_system: str=field(init=False,default="dimensional")
    time_units: str=field(init=False,default="dimensional")
    state_units: str=field(init=False,default="dimensional")
    control_units: str=field(init=False,default="system input")
    _bound_system: object | None=field(init=False,default=None,repr=False)
    _physical_active_time_window: tuple[float,float] | None=field(init=False,default=None,repr=False)
    def __post_init__(self):
        self.frequencies_hz=np.asarray(self.frequencies_hz,dtype=float).reshape(-1); self.amplitudes=np.asarray(self.amplitudes,dtype=float).reshape(-1)
        if self.frequencies_hz.size!=self.amplitudes.size or not self.frequencies_hz.size: raise ValueError("Waveform frequencies and amplitudes must be nonempty and equally sized.")
        self.phases_rad=np.zeros_like(self.frequencies_hz) if self.phases_rad is None else np.asarray(self.phases_rad,dtype=float).reshape(-1)
        if self.phases_rad.size!=self.frequencies_hz.size: raise ValueError("Waveform phases must match the number of frequencies.")
        if self.active_time_window is not None and self.active_time_window[1] < self.active_time_window[0]: raise ValueError("active_time_window must satisfy start <= end.")
        self._physical_active_time_window=self.active_time_window
    def bind_coordinate_system(self,coordinate_system: str,*,system=None) -> None:
        if coordinate_system not in {"dimensionless","dimensional"}: raise ValueError("coordinate_system must be 'dimensionless' or 'dimensional'.")
        if coordinate_system=="dimensionless":
            if system is None or not callable(getattr(system,"to_dimensional_time",None)):
                raise TypeError("A dimensionless waveform controller requires an explicit system time conversion.")
            if self.actuation=="system_input" and not callable(getattr(system,"to_dimensionless_input",None)):
                raise TypeError("A dimensionless system-input waveform requires system.to_dimensionless_input().")
            if self.actuation!="system_input":
                raise TypeError("Generic additive-state waveform forcing cannot be nondimensionalized without an explicit state-input scaling.")
            self.active_time_window=None if self._physical_active_time_window is None else tuple(float(np.asarray(system.to_dimensionless_time(item))) for item in self._physical_active_time_window)
        else: self.active_time_window=self._physical_active_time_window
        self._bound_system=system; self.coordinate_system=coordinate_system
        self.time_units=coordinate_system; self.state_units=coordinate_system
        self.control_units="dimensionless system input" if coordinate_system=="dimensionless" else "system input"
    def value(self,t,state):
        physical_time=float(np.asarray(self._bound_system.to_dimensional_time(t))) if self.coordinate_system=="dimensionless" else float(t)
        raw=np.array([self.offset+np.sum(self.amplitudes*np.sin(2*np.pi*self.frequencies_hz*physical_time+self.phases_rad))])
        if self.coordinate_system=="dimensionless": raw=np.asarray(self._bound_system.to_dimensionless_input(raw),dtype=float)
        return raw if self.output_matrix is None else np.asarray(self.output_matrix,dtype=float)@raw

def _matrix(value,state): return np.asarray(value(state),dtype=float) if callable(value) else np.asarray(value,dtype=float)

@dataclass
class LQRController:
    """Continuous LQR or SDRE-style LQR when A/B are state-dependent callables."""
    A: np.ndarray | Callable[[np.ndarray],np.ndarray]
    B: np.ndarray | Callable[[np.ndarray],np.ndarray]
    Q: np.ndarray
    R: np.ndarray
    target: Target | None=None
    actuation: str="additive_state"
    active_time_window: tuple[float,float] | None=None
    last_gain: np.ndarray | None=field(init=False,default=None)
    _fixed_gain: np.ndarray | None=field(init=False,default=None,repr=False)
    # Repeated evaluations at adaptive solver trial states do not advance
    # controller dynamics. ``last_gain`` is diagnostic only.
    continuous_rhs_safe: ClassVar[bool]=True
    @property
    def requires_sample_clock(self) -> bool:
        """State-dependent SDRE gains must be updated on a defined clock."""
        return callable(self.A) or callable(self.B)
    def value(self,t,state):
        x=np.asarray(state,dtype=float).reshape(-1); reference=np.zeros_like(x) if self.target is None else self.target.evaluate(t,x.size)
        A,B=_matrix(self.A,x),_matrix(self.B,x); Q,R=np.asarray(self.Q,dtype=float),np.asarray(self.R,dtype=float)
        if self._fixed_gain is None or self.requires_sample_clock:
            P=solve_continuous_are(A,B,Q,R); K=np.linalg.solve(R,B.T@P)
            if not self.requires_sample_clock: self._fixed_gain=K
        else: K=self._fixed_gain
        self.last_gain=K
        return B@(-K@(x-reference))

@dataclass
class MPCController(LQRController):
    """Unconstrained finite-horizon linear MPC with optional input limits.

    The local continuous model is discretized with a forward-Euler hold over
    ``sample_time`` and a backward Riccati recursion supplies the first move.
    It is deliberately small and dependency-free; constrained optimization
    backends can be attached later without changing the controller contract.
    """
    horizon: int=20
    sample_time: float=1e-3
    control_limit: np.ndarray | None=None
    continuous_rhs_safe: ClassVar[bool]=False
    @property
    def requires_sample_clock(self) -> bool: return True
    def value(self,t,state):
        x=np.asarray(state,dtype=float).reshape(-1); reference=np.zeros_like(x) if self.target is None else self.target.evaluate(t,x.size)
        A,B=_matrix(self.A,x),_matrix(self.B,x); Q,R=np.asarray(self.Q,dtype=float),np.asarray(self.R,dtype=float)
        if self.horizon<1 or self.sample_time<=0: raise ValueError("MPC horizon and sample_time must be positive.")
        Ad=np.eye(A.shape[0])+self.sample_time*A; Bd=self.sample_time*B; P=Q.copy()
        K=None
        for _ in range(self.horizon):
            K=np.linalg.solve(R+Bd.T@P@Bd,Bd.T@P@Ad)
            P=Q+Ad.T@P@(Ad-Bd@K)
        self.last_gain=K
        input_value=-K@(x-reference)
        if self.control_limit is None: return B@input_value
        limits=np.asarray(self.control_limit,dtype=float).reshape(-1)
        if limits.size==1: limits=np.full(input_value.size,limits.item())
        if limits.size!=input_value.size: raise ValueError("MPC control_limit must have one value per control input.")
        input_value=np.clip(input_value,-limits,limits)
        return B@input_value

@dataclass
class PIDController:
    """Sample-clock proportional-integral-derivative controller."""
    kp: np.ndarray | float
    ki: np.ndarray | float=0.
    kd: np.ndarray | float=0.
    target: Target | None=None
    actuation: str="additive_state"
    active_time_window: tuple[float,float] | None=None
    _integral: np.ndarray | None=field(init=False,default=None)
    _previous_error: np.ndarray | None=field(init=False,default=None)
    _previous_time: float | None=field(init=False,default=None)
    requires_sample_clock: ClassVar[bool]=True
    def value(self,t,state):
        x=np.asarray(state,dtype=float).reshape(-1); reference=np.zeros_like(x) if self.target is None else self.target.evaluate(t,x.size); error=reference-x
        if self._integral is None: self._integral=np.zeros_like(error); derivative=np.zeros_like(error)
        else:
            # ``0.0`` is a valid previous sample time, not a missing value.
            dt=float(t-(self._previous_time if self._previous_time is not None else t));
            if dt>0: self._integral+=error*dt; derivative=(error-self._previous_error)/dt
            else: derivative=np.zeros_like(error)
        self._previous_error=error.copy(); self._previous_time=float(t)
        return np.asarray(self.kp)*error+np.asarray(self.ki)*self._integral+np.asarray(self.kd)*derivative
