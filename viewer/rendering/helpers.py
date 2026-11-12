"""Pure display-shaping helpers shared by Viewer renderer families."""
from __future__ import annotations

from typing import Any, Sequence
import numpy as np

from ..contracts import PlotRequest, PlotStyle, _DEFAULT_TRACE_COLORS, _KROM_MODE_COLORS
from ..settings import get_interpolation_sampling

def _sequence(value: Any, count: int) -> list[Any]:
    if isinstance(value,str) or not isinstance(value,Sequence): return [value]*count
    if not value: raise ValueError("Style sequence cannot be empty.")
    return [value[index % len(value)] for index in range(count)]

def _control_display(metadata: dict[str,Any], component: int) -> tuple[str,float]:
    history=metadata.get("control_history") or {}
    if str(metadata.get("system"))=="rayleigh-plesset" and component==0:
        if history.get("coordinate_system")=="dimensionless":
            return "dimensionless pressure control u* = p_control/p0",1.
        pressure_scale=float((metadata.get("nondimensionalization") or {}).get("pressure_scale_pa") or 0.)
        if np.isfinite(pressure_scale) and pressure_scale>0:
            return "dimensionless pressure control p_c* = p_c/P0",1./pressure_scale
        return "pressure control p_c (Pa)",1.
    units=history.get("units")
    if units:
        prefix="dimensionless control" if history.get("coordinate_system")=="dimensionless" else "control"
        return f"{prefix} u{component} ({units})",1.
    return f"control u{component}",1.

def _krom_display_indices(request: PlotRequest, retained: tuple[int,...]) -> tuple[int,...]:
    """Select dashboard panels without altering the retained model spectrum."""
    if request.krom_mode_indices is not None:
        selected=tuple(int(index) for index in request.krom_mode_indices)
        if not selected: raise ValueError("krom_mode_indices cannot be empty; use a Koopman spectrum plot for spectrum-only display.")
        if len(set(selected))!=len(selected): raise ValueError("krom_mode_indices must be unique.")
        invalid=[index for index in selected if index not in retained]
        if invalid: raise ValueError(f"KROM dashboard indices are not retained model modes: {invalid}.")
        return selected
    limit=request.krom_mode_limit
    if limit is None: return retained
    if int(limit)!=limit or limit<1: raise ValueError("krom_mode_limit must be a positive integer or None.")
    return retained[:int(limit)]

def _krom_colors(style: PlotStyle, count: int) -> list[str]:
    if isinstance(style.colors,str): return [style.colors]*count
    if tuple(style.colors)==_DEFAULT_TRACE_COLORS:
        if count<=len(_KROM_MODE_COLORS): return list(_KROM_MODE_COLORS[:count])
        return [f"hsl({(index*137.507764)%360:.3f},72%,46%)" for index in range(count)]
    return [str(value) for value in _sequence(style.colors,count)]

def _interpolation_kind(metadata: dict[str,Any] | None, fallback: str="pod") -> str:
    method=str((metadata or {}).get("method","")).lower()
    if "sindy" in method: return "sindy"
    if "havok" in method: return "havok"
    if "pod" in method: return "pod"
    if "hankel" in method or "hdmd" in method: return "hdmd"
    return fallback

def _display_samples(samples: np.ndarray, kind: str, requested: int | None) -> tuple[np.ndarray,dict[str,Any]]:
    settings=get_interpolation_sampling(kind); limit=settings.max_state_samples if requested is None else int(requested)
    if limit is None: limit=samples.shape[0]
    if limit<1: raise ValueError("max_display_samples must be positive when specified.")
    count=samples.shape[0]
    indices=np.arange(count) if count<=limit else np.linspace(0,count-1,limit,dtype=int)
    return samples[indices],{"kind":settings.kind,"recommended_max_state_samples":settings.max_state_samples,
                              "recommended_domain_density":settings.recommended_domain_density,"displayed_state_samples":int(indices.size),
                              "source_state_samples":int(count),"notes":settings.notes}

def _trace_display_indices(arrays: Sequence[np.ndarray], limit: int) -> np.ndarray:
    """Select a bounded display set while retaining endpoints and sharp events."""
    count=int(arrays[0].shape[0])
    if limit<1: raise ValueError("max_display_samples must be positive when specified.")
    if count<=limit: return np.arange(count,dtype=int)
    if limit<=2: return np.linspace(0,count-1,limit,dtype=int)
    selected={0,count-1}
    numeric=[]
    for values in arrays:
        try: vector=np.asarray(values,dtype=float).reshape(-1)
        except (TypeError,ValueError): continue
        if vector.size!=count or not np.any(np.isfinite(vector)): continue
        numeric.append(vector)
        finite=np.flatnonzero(np.isfinite(vector))
        for position in (finite[np.argmin(vector[finite])],finite[np.argmax(vector[finite])]):
            if len(selected)<limit: selected.add(int(position))
    # Preserve the strongest discontinuities or narrow transients before
    # filling the remainder uniformly. This keeps control-engagement spikes
    # visible even when a dense solver archive is aggressively thinned.
    if numeric and len(selected)<limit:
        change_score=np.zeros(count-1,dtype=float)
        for vector in numeric:
            finite=vector[np.isfinite(vector)]
            scale=float(np.ptp(finite)) if finite.size else 0.
            if not np.isfinite(scale) or scale<=np.finfo(float).eps: continue
            differences=np.abs(np.diff(vector))/scale
            change_score=np.maximum(change_score,np.nan_to_num(differences,nan=0.,posinf=0.,neginf=0.))
        transition_count=min(16,max(0,(limit-len(selected))//2),change_score.size)
        if transition_count:
            strongest=np.argpartition(change_score,-transition_count)[-transition_count:]
            strongest=strongest[np.argsort(change_score[strongest])[::-1]]
            for position in strongest:
                for candidate in (int(position),int(position)+1):
                    if len(selected)<limit: selected.add(candidate)
    for position in np.linspace(0,count-1,limit,dtype=int):
        if len(selected)>=limit: break
        selected.add(int(position))
    if len(selected)<limit:
        for position in range(count):
            if len(selected)>=limit: break
            selected.add(position)
    return np.asarray(sorted(selected),dtype=int)

def _limit_display_traces(traces: Sequence[dict[str,Any]], requested: int | None) -> tuple[list[dict[str,Any]],dict[str,Any]]:
    """Downsample plot traces only; archived numerical data remain untouched."""
    if requested is None: return list(traces),{}
    limit=int(requested)
    if limit<1: raise ValueError("max_display_samples must be positive when specified.")
    limited=[]; source_counts=[]; display_counts=[]
    for trace in traces:
        arrays={}
        for key in ("x","y","z","value"):
            if key not in trace: continue
            array=np.asarray(trace[key])
            if array.ndim==1: arrays[key]=array
        counts=[array.shape[0] for array in arrays.values()]
        count=max(counts,default=0)
        aligned=[array for array in arrays.values() if array.shape[0]==count]
        if not count or count<=limit or not aligned:
            limited.append(trace); source_counts.append(count); display_counts.append(count); continue
        indices=_trace_display_indices(aligned,limit); result=dict(trace)
        for key,array in arrays.items():
            if array.shape[0]==count: result[key]=array[indices]
        limited.append(result); source_counts.append(count); display_counts.append(int(indices.size))
    metadata={"display_downsampling":{"requested_max_samples":limit,
        "maximum_source_samples":max(source_counts,default=0),
        "maximum_displayed_samples":max(display_counts,default=0),
        "method":"uniform-plus-extrema-and-sharp-events",
        "archive_data_unchanged":True}}
    return limited,metadata

def _focused_spectrum(frequencies: np.ndarray, power: np.ndarray, scale: str) -> tuple[np.ndarray,np.ndarray,tuple[float,float] | None,str]:
    positive=frequencies>0; frequency=frequencies[positive]; values=power[positive]
    if not frequency.size: return frequency,values,None,"spectral power"
    importance=np.maximum(values,0.); total=importance.sum()
    if total>0:
        cutoff=frequency[np.searchsorted(np.cumsum(importance),.995*total)]
        peak=frequency[int(np.argmax(importance))]; xmax=max(cutoff,1.25*peak)
    else: xmax=frequency[-1]
    if scale in {"auto","linear"}: return frequency,values,(0.,float(xmax)),"spectral power"
    if scale!="log_power": raise ValueError("spectrum_scale must be 'linear' or 'log_power'.")
    floor=max(float(np.nanmax(values))*1e-12,np.finfo(float).tiny)
    return frequency,np.log10(np.maximum(values,floor)),(0.,float(xmax)),"log₁₀ spectral power"

def _automatic_time_detail_inset(traces: Sequence[dict[str,Any]], *, disparity_threshold: float=20.) -> dict[str,Any] | None:
    """Describe a detail inset when the main time-series scale hides useful data.

    The decision is deliberately robust to a few extreme samples. Viewer owns
    this display policy; solution and control archives remain untouched.
    """
    records=[]
    for index,trace in enumerate(traces):
        if "x" not in trace or "y" not in trace or trace.get("marker"):
            continue
        x=np.asarray(trace["x"],dtype=float).reshape(-1); y=np.asarray(trace["y"],dtype=float).reshape(-1)
        finite=np.isfinite(x)&np.isfinite(y)
        if np.count_nonzero(finite)<20:
            continue
        x=x[finite]; y=y[finite]
        lower,upper=np.quantile(y,(.01,.99)); robust_span=float(upper-lower)
        full_span=float(np.max(y)-np.min(y))
        records.append({"index":index,"x":x,"y":y,"lower":float(lower),"upper":float(upper),
                        "robust_span":robust_span,"full_span":full_span})
    if not records:
        return None
    candidates=[record for record in records if record["robust_span"]>np.finfo(float).eps and
                record["full_span"]/record["robust_span"]>=disparity_threshold]
    if candidates:
        focus=max(candidates,key=lambda record:record["full_span"]/record["robust_span"])
        selected=[record["index"] for record in records]
        lower=min(record["lower"] for record in records); upper=max(record["upper"] for record in records)
        inlier=(focus["y"]>=focus["lower"])&(focus["y"]<=focus["upper"])
        x_values=focus["x"][inlier]
        xmin,xmax=(float(np.min(x_values)),float(np.max(x_values))) if x_values.size else (float(np.min(focus["x"])),float(np.max(focus["x"])))
        reason="transient-scale-disparity"
        ratio=float(focus["full_span"]/focus["robust_span"])
    elif len(records)>1:
        positive=[record for record in records if record["robust_span"]>np.finfo(float).eps]
        if len(positive)<2:
            return None
        largest=max(record["robust_span"] for record in positive); smallest=min(record["robust_span"] for record in positive)
        if largest/smallest<disparity_threshold:
            return None
        selected=[record["index"] for record in positive if record["robust_span"]<=largest/disparity_threshold]
        detail=[record for record in positive if record["index"] in selected]
        lower=min(record["lower"] for record in detail); upper=max(record["upper"] for record in detail)
        xmin=min(float(np.min(record["x"])) for record in detail); xmax=max(float(np.max(record["x"])) for record in detail)
        reason="trace-scale-disparity"; ratio=float(largest/smallest)
    else:
        return None
    if xmax<=xmin:
        return None
    yspan=upper-lower
    if yspan<=np.finfo(float).eps: yspan=max(abs(lower),1.)*.1
    return {"automatic":True,"reason":reason,"scale_ratio":ratio,"trace_indices":selected,
            "xlim":[xmin,xmax],"ylim":[float(lower-.08*yspan),float(upper+.08*yspan)]}
