"""Post-solve resampling for continuous dynsyslab trajectories."""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator, interp1d
from .models import InterpolationOptions, NumericalResult, SolveRequest
from time_series import TimeSeries

def target_times(request: SolveRequest) -> np.ndarray | None:
    """Resolve an explicit target grid or create a uniform grid when requested."""
    options=request.interpolation
    supplied=[item for item in (request.sample_times, options.times) if item is not None]
    if len(supplied)>1:
        raise ValueError("Specify target times through either sample_times or interpolation.times, not both.")
    if supplied:
        target=np.asarray(supplied[0],dtype=float).reshape(-1)
    elif options.enabled:
        if options.time_step is None or request.time_span is None:
            raise ValueError("InterpolationOptions(enabled=True) requires time_step and time_span when no target times are supplied.")
        if options.time_step<=0: raise ValueError("interpolation.time_step must be positive.")
        t0,tf=request.time_span
        # Use half a step of numerical slack so a mathematically aligned final
        # endpoint is retained despite binary floating-point roundoff.
        target=np.arange(t0,tf+options.time_step*.5,options.time_step,dtype=float)
        if np.isclose(target[-1],tf,rtol=1e-12,atol=1e-15): target[-1]=tf
    else:
        return None
    if target.size==0 or not np.all(np.isfinite(target)) or np.any(np.diff(target)<=0):
        raise ValueError("Interpolation target times must be finite, nonempty, and strictly increasing.")
    return target

def apply_interpolation(request: SolveRequest, result: NumericalResult) -> NumericalResult:
    """Keep raw integration output and replace result output only when requested."""
    target=target_times(request)
    if target is None or result.status!="complete": return result
    if request.kind=="lattice":
        raise ValueError("Interpolation is not supported for discrete lattice trajectories.")
    if result.time is None or len(result.time)<2:
        raise ValueError("At least two raw time samples are required for interpolation.")
    raw_time=np.asarray(result.time,dtype=float)
    raw_values=np.asarray(result.values)
    # Checkpointed adaptive integrations can repeat a boundary sample. Retain
    # the raw archive untouched while using the first sample per unique time for interpolation.
    unique_time,indices=np.unique(raw_time,return_index=True)
    if unique_time.size<2 or np.any(np.diff(unique_time)<=0):
        raise ValueError("Raw integration times must contain at least two increasing samples.")
    endpoint_tolerance=32*np.finfo(float).eps*max(1.,abs(float(unique_time[0])),abs(float(unique_time[-1])))
    if target[0] < unique_time[0]-endpoint_tolerance or target[-1] > unique_time[-1]+endpoint_tolerance:
        raise ValueError("Interpolation target times must lie within the solved time span.")
    # Preserve an explicitly requested grid while snapping only roundoff-level
    # endpoint differences to the actual adaptive-solver interval.
    target=np.clip(target,unique_time[0],unique_time[-1])
    samples=raw_values[indices]
    method=request.interpolation.method.lower()
    if method=="pchip": interpolator=PchipInterpolator(unique_time,samples,axis=0)
    elif method=="linear": interpolator=interp1d(unique_time,samples,axis=0,kind="linear")
    elif method=="nearest": interpolator=interp1d(unique_time,samples,axis=0,kind="nearest")
    else: raise ValueError("interpolation.method must be one of: pchip, linear, nearest.")
    result.raw_time_series=TimeSeries(raw_time,raw_values)
    result.time_series=TimeSeries(target,np.asarray(interpolator(target)))
    if result.control_time_series is not None:
        control=result.control_time_series
        if control.time.shape!=raw_time.shape or not np.allclose(control.time,raw_time,rtol=1e-12,atol=1e-15):
            raise ValueError("Control history must be aligned with the raw solution before interpolation.")
        result.raw_control_time_series=control
        if method=="pchip": control_interpolator=PchipInterpolator(unique_time,control.values[indices],axis=0)
        else: control_interpolator=interp1d(unique_time,control.values[indices],axis=0,kind=method)
        result.control_time_series=TimeSeries(target,np.asarray(control_interpolator(target)))
    if result.control_target_time_series is not None:
        control_target=result.control_target_time_series
        if control_target.time.shape!=raw_time.shape or not np.allclose(control_target.time,raw_time,rtol=1e-12,atol=1e-15):
            raise ValueError("Control-target history must be aligned with the raw solution before interpolation.")
        result.raw_control_target_time_series=control_target
        if method=="pchip": target_interpolator=PchipInterpolator(unique_time,control_target.values[indices],axis=0)
        else: target_interpolator=interp1d(unique_time,control_target.values[indices],axis=0,kind=method)
        result.control_target_time_series=TimeSeries(target,np.asarray(target_interpolator(target)))
    result.metadata["interpolated"]=True
    result.metadata["interpolation"]={"method":method,"raw_samples":int(raw_time.size),"output_samples":int(target.size)}
    return result
