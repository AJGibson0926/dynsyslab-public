"""User-approved, in-memory resampling shared by all post-processing methods."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d
from .core import Trajectory
from time_series import TimeSeries

class ResamplingDecisionRequired(ValueError):
    prompt=("Postprocessing is blocked because the trajectory has nonuniform time steps. Would you like to resample it in memory for this analysis? Choose a uniform time_step or explicit uniform times and pchip (recommended), linear, or nearest. The original archive will not be altered.")
    def __init__(self,reason): super().__init__(f"{reason} {self.prompt}"); self.reason=reason; self.interpolation_methods=("pchip","linear","nearest")
@dataclass(frozen=True)
class ResamplingOptions:
    time_step: float|None=None
    times: np.ndarray|None=None
    method: str="pchip"
def require_uniform_time_steps(time:np.ndarray,metadata:dict[str,Any]|None=None):
    samples=np.asarray(time,float).reshape(-1)
    if samples.size<3: raise ValueError("Fixed-lag postprocessing requires at least three time samples.")
    increments=np.diff(samples); verified=bool(np.all(np.isfinite(increments)) and np.all(increments>0) and np.allclose(increments,increments[0],rtol=1e-10,atol=1e-13)); declared=None if metadata is None else metadata.get("uniform_time_steps")
    if not verified: raise ResamplingDecisionRequired("Postprocessing blocked: timestamps are nonuniform or metadata does not verify uniformity.")
    return {"metadata_uniform_time_steps":declared,"verified_uniform_time_steps":verified,"delta_t":float(increments[0])}
def resample(source:Trajectory,options:ResamplingOptions)->Trajectory:
    if options.time_step is not None and options.times is not None: raise ValueError("Choose either time_step or explicit times, not both.")
    raw=np.asarray(source.time,float); unique,indices=np.unique(raw,return_index=True)
    target=np.asarray(options.times,float).reshape(-1) if options.times is not None else np.arange(raw[0],raw[-1]+options.time_step*1e-12,options.time_step,float) if options.time_step is not None and options.time_step>0 else None
    tolerance=max(1.,abs(raw[-1]))*1e-12
    if target is None or target.size<3 or np.any(np.diff(target)<=0) or target[0]<raw[0]-tolerance or target[-1]>raw[-1]+tolerance: raise ValueError("Resampling requires increasing target times within the original span.")
    target=np.clip(target,raw[0],raw[-1])
    method=options.method.lower(); values=np.asarray(source.values)[indices]
    interpolator=PchipInterpolator(unique,values,axis=0) if method=="pchip" else interp1d(unique,values,axis=0,kind=method) if method in {"linear","nearest"} else None
    if interpolator is None: raise ValueError("Resampling method must be pchip, linear, or nearest.")
    metadata=dict(source.metadata); metadata["uniform_time_steps"]=True; metadata["postprocessing_resampling"]={"in_memory_only":True,"method":method,"raw_samples":int(raw.size),"output_samples":int(target.size),"time_step":float(target[1]-target[0])}
    return Trajectory(TimeSeries(target,np.asarray(interpolator(target))),metadata,source.source)
