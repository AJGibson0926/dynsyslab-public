"""Shared numerical time-series value object for dynsyslab.

This module deliberately contains solution arrays only.  System identity,
solver settings, and analysis provenance remain separate metadata.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

@dataclass(frozen=True)
class FourierTransform:
    """One-sided or full FFT data derived from a uniformly sampled TimeSeries."""
    frequencies: np.ndarray
    coefficients: np.ndarray
    one_sided: bool
    delta_t: float

@dataclass(frozen=True)
class TimeSeries:
    """Numerical values indexed by strictly increasing time values."""
    time: np.ndarray
    values: np.ndarray
    pooled_members: int = 0
    uniform_time_steps: bool = field(init=False)
    delta_t: float | None = field(init=False)
    is_pooled: bool = field(init=False)
    fourier_transform: FourierTransform | None = field(init=False,default=None,compare=False)

    def __post_init__(self):
        time=np.asarray(self.time,dtype=float).reshape(-1)
        values=np.asarray(self.values)
        if time.size == 0: raise ValueError("A TimeSeries requires at least one time value.")
        if not np.all(np.isfinite(time)) or np.any(np.diff(time)<=0):
            raise ValueError("TimeSeries time values must be finite and strictly increasing.")
        if values.ndim < 1 or values.shape[0] != time.size:
            raise ValueError("TimeSeries values must have one leading sample axis matching time.")
        members=int(self.pooled_members)
        if members<0 or (members and (values.ndim<2 or values.shape[1]!=members)):
            raise ValueError("pooled_members must be zero, or match values.shape[1] for pooled data.")
        delta=float(time[1]-time[0]) if time.size>1 else None
        uniform=bool(time.size<3 or np.allclose(np.diff(time),delta,rtol=1e-12,atol=1e-15))
        object.__setattr__(self,"time",time)
        object.__setattr__(self,"values",values)
        object.__setattr__(self,"pooled_members",members)
        object.__setattr__(self,"is_pooled",bool(members))
        object.__setattr__(self,"uniform_time_steps",uniform)
        object.__setattr__(self,"delta_t",delta if uniform else None)

    @property
    def sample_count(self) -> int: return int(self.time.size)

    @property
    def sample_shape(self) -> tuple[int,...]: return tuple(self.values.shape[1:])

    @property
    def ensemble_size(self) -> int: return self.pooled_members if self.is_pooled else 1

    def members(self) -> tuple["TimeSeries",...]:
        """Return one member per pooled trajectory, or this series if individual."""
        if not self.is_pooled: return (self,)
        return tuple(TimeSeries(self.time,self.values[:,index]) for index in range(self.pooled_members))

    def fft(self, *, one_sided: bool=True) -> FourierTransform:
        """Compute and cache the FFT along time, requiring a uniform time grid.

        The returned coefficients retain every non-time value axis.  Real-valued
        data use ``rfft`` by default; complex data automatically use the full
        transform because a one-sided spectrum would discard information.
        """
        if not self.uniform_time_steps or self.delta_t is None:
            raise ValueError("FFT requires uniform_time_steps=True; resample first.")
        if self.sample_count < 2:
            raise ValueError("FFT requires at least two time samples.")
        values=np.asarray(self.values)
        use_one_sided=bool(one_sided and not np.iscomplexobj(values))
        if use_one_sided:
            frequencies=np.fft.rfftfreq(self.sample_count,d=self.delta_t)
            coefficients=np.fft.rfft(values,axis=0)
        else:
            frequencies=np.fft.fftfreq(self.sample_count,d=self.delta_t)
            coefficients=np.fft.fft(values,axis=0)
        result=FourierTransform(frequencies,coefficients,use_one_sided,self.delta_t)
        object.__setattr__(self,"fourier_transform",result)
        return result
