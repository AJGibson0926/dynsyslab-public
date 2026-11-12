"""Sparse Identification of Nonlinear Dynamics (SINDy) reduced models."""
from __future__ import annotations
from dataclasses import dataclass, field
from itertools import combinations_with_replacement
from pathlib import Path
from typing import Any
import numpy as np

from time_series import TimeSeries
from analysis_data import load_solution, save_derived
from spectral.dictionary import ObservableDictionary, PolynomialDictionary


@dataclass
class SINDyResult:
    """Sparse discovered vector field in a chosen observable-library basis."""
    coefficients: np.ndarray  # (features, state variables)
    feature_values: np.ndarray
    derivatives: np.ndarray
    feature_names: list[str]
    metadata: dict = field(default_factory=dict)


def _feature_names(library: ObservableDictionary, state_dimension: int, feature_count: int) -> list[str]:
    if isinstance(library, PolynomialDictionary):
        names=[]
        if library.include_constant: names.append("1")
        variables=[f"x{i}" for i in range(state_dimension)]
        for order in range(1,library.degree+1):
            for indices in combinations_with_replacement(range(state_dimension),order):
                names.append("*".join(variables[index] for index in indices))
        return names
    return [f"psi_{index}" for index in range(feature_count)]


def _stlsq(theta: np.ndarray, derivative: np.ndarray, threshold: float, max_iterations: int, normalize_columns: bool) -> np.ndarray:
    if threshold<0: raise ValueError("SINDy threshold must be nonnegative.")
    scales=np.linalg.norm(theta,axis=0)
    scales=np.where(scales>np.finfo(float).eps,scales,1.0) if normalize_columns else np.ones(theta.shape[1])
    scaled=theta/scales
    coefficients=np.linalg.lstsq(scaled,derivative,rcond=None)[0]
    for _ in range(max_iterations):
        previous=coefficients.copy(); coefficients[np.abs(coefficients)<threshold]=0.0
        for state in range(derivative.shape[1]):
            active=coefficients[:,state]!=0
            if np.any(active): coefficients[active,state]=np.linalg.lstsq(scaled[:,active],derivative[:,state],rcond=None)[0]
        if np.array_equal(coefficients==0,previous==0): break
    return coefficients/scales[:,None]


def sindy(states: np.ndarray | TimeSeries, time: np.ndarray|None=None, library: ObservableDictionary|None=None, *, derivatives: np.ndarray|None=None,
          threshold: float=0.1, max_iterations: int=20, normalize_columns: bool=True) -> SINDyResult:
    """Identify the ROM ``dx/dt = Theta(x) Xi`` using STLSQ.

    If derivatives are omitted, second-order finite differences are estimated
    with ``numpy.gradient`` using the supplied timestamps. Uneven timestamps
    are supported; their exact values, not an assumed constant step, are used.
    """
    if isinstance(states,TimeSeries):
        if time is not None: raise ValueError("Pass either a TimeSeries or separate states and time, not both.")
        if states.is_pooled: raise ValueError("This SINDy call received a pooled TimeSeries. Use series.members() to preserve trajectory boundaries.")
        time=states.time; states=states.values
    if time is None: raise ValueError("Array SINDy input requires an explicit time array.")
    x=np.asarray(states,dtype=float)
    if x.ndim==1: x=x[:,None]
    t=np.asarray(time,dtype=float).reshape(-1)
    if x.ndim!=2 or x.shape[0]<3: raise ValueError("SINDy requires at least three time-ordered state samples.")
    if t.size!=x.shape[0] or np.any(~np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError("SINDy requires finite, strictly increasing timestamps matching the state samples.")
    dxdt=np.gradient(x,t,axis=0,edge_order=2) if derivatives is None else np.asarray(derivatives,dtype=float)
    if dxdt.shape!=x.shape: raise ValueError("Supplied SINDy derivatives must match the state array shape.")
    dictionary=library or PolynomialDictionary()
    theta=np.asarray(dictionary.transform(x),dtype=float)
    if theta.ndim!=2 or theta.shape[0]!=x.shape[0]: raise ValueError("SINDy observable library must return one feature row per sample.")
    coefficients=_stlsq(theta,dxdt,threshold,max_iterations,normalize_columns)
    increments=np.diff(t)
    uniform=bool(np.allclose(increments,increments[0],rtol=1e-12,atol=1e-15))
    return SINDyResult(coefficients,theta,dxdt,_feature_names(dictionary,x.shape[1],theta.shape[1]),{
        "threshold":threshold,"max_iterations":max_iterations,"normalize_columns":normalize_columns,
        "derivatives_supplied":derivatives is not None,"uniform_time_steps":uniform,
        "library_class":type(dictionary).__name__,"features":int(theta.shape[1]),"state_dimension":int(x.shape[1]),
        "preprocessing":{"sample_count":int(x.shape[0]),"state_dimension":int(x.shape[1]),"analysis_delta_t":float(increments[0]) if uniform else None,"derivative_method":"supplied" if derivatives is not None else "numpy.gradient"},
        "delay_embedding":{"type":"none","delays":1,"window_samples":1,"window_duration":0.},
    })


def sindy_archive(path: str | Path, library: ObservableDictionary|None=None, *, derivatives: np.ndarray|None=None, threshold: float=.1, max_iterations: int=20, normalize_columns: bool=True):
    """Fit and archive a plain SINDy ROM from one completed solution."""
    source=load_solution(path)
    result=sindy(source.values,source.time,library,derivatives=derivatives,threshold=threshold,max_iterations=max_iterations,normalize_columns=normalize_columns)
    return save_derived(source,"sindy",coefficients=result.coefficients,feature_values=result.feature_values,derivatives=result.derivatives,
                        feature_names=np.asarray(result.feature_names),analysis_metadata=result.metadata)


def sindy_processor(library: ObservableDictionary|None=None, *, threshold: float=.1, max_iterations: int=20, normalize_columns: bool=True):
    """Return a solver-compatible callback that archives a plain SINDy ROM."""
    def processor(path: str): return sindy_archive(path,library,threshold=threshold,max_iterations=max_iterations,normalize_columns=normalize_columns)
    processor.__name__="sindy_rom"; return processor
