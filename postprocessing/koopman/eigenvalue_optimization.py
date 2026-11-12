"""Trajectory-regression optimization of Koopman eigenvalues as one stage."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from scipy.optimize import minimize

from time_series import TimeSeries
from .spectrum import DynamicalSystemReference, KoopmanEigenvalue, KoopmanSpectrum


@dataclass(frozen=True)
class EigenvalueOptimizationResult:
    """Eigenvalue-only output, deliberately independent of eigenfunctions."""

    spectrum: KoopmanSpectrum
    groups: tuple[tuple[int,...],...]
    diagnostics: dict[str,Any]


def _powers_and_fit(eigenvalues: np.ndarray, measurements: np.ndarray):
    steps=np.arange(measurements.shape[0],dtype=float)[:,None]
    powers=np.power(eigenvalues[None,:],steps)
    boundary=np.linalg.lstsq(powers,measurements,rcond=1e-12)[0]
    residual=measurements-powers@boundary
    return powers,boundary,residual


def _objective(parameters: np.ndarray, measurements: np.ndarray) -> tuple[float,np.ndarray]:
    count=parameters.size//2
    eigenvalues=parameters[:count]+1j*parameters[count:]
    powers,boundary,residual=_powers_and_fit(eigenvalues,measurements)
    normalizer=max(float(np.sum(np.abs(measurements)**2)),np.finfo(float).tiny)
    gradient_real=np.empty(count); gradient_imag=np.empty(count)
    steps=np.arange(measurements.shape[0],dtype=float)
    for index,value in enumerate(eigenvalues):
        derivative=np.zeros(measurements.shape[0],dtype=complex)
        if measurements.shape[0]>1:
            derivative[1:]=steps[1:]*np.power(value,steps[1:]-1.)
        direction=derivative[:,None]*boundary[index][None,:]
        gradient_real[index]=-2.*np.real(np.vdot(residual,direction))/normalizer
        gradient_imag[index]=-2.*np.real(np.vdot(residual,1j*direction))/normalizer
    objective=float(np.sum(np.abs(residual)**2)/normalizer)
    return objective,np.concatenate((gradient_real,gradient_imag))


def _symmetrize_from_initial(values: np.ndarray, initial: np.ndarray) -> np.ndarray:
    """Preserve conjugate-dual structure identified in the initialization."""
    result=np.array(values,dtype=complex,copy=True); used=set()
    for index,value in enumerate(initial):
        if index in used or np.isclose(value.imag,0.,atol=1e-12): continue
        partner=next((candidate for candidate in range(index+1,initial.size) if candidate not in used
                      and np.isclose(initial[candidate],np.conjugate(value),rtol=1e-7,atol=1e-9)),None)
        if partner is None: continue
        representative=.5*(result[index]+np.conjugate(result[partner]))
        result[index]=representative; result[partner]=np.conjugate(representative)
        used.update((index,partner))
    return result


def optimize_trajectory_eigenvalues(
    series: TimeSeries,
    *,
    observable,
    initial_eigenvalue_groups: Sequence[Sequence[complex]],
    system: Any,
    metadata: dict[str,Any] | None=None,
    bounds: tuple[float,float]=(-1.05,1.05),
    max_iterations: int=500,
    tolerance: float=1e-10,
    preserve_conjugate_pairs: bool=True,
) -> EigenvalueOptimizationResult:
    """Optimize discrete eigenvalues without constructing eigenfunctions.

    Each observable output is fitted independently by minimizing the residual
    after analytically eliminating trajectory-specific boundary values.  This
    is the separated form of the Korda--Mezic ``eigOptim_grad`` calculation.
    """
    if not isinstance(series,TimeSeries) or not series.is_pooled:
        raise TypeError("Eigenvalue optimization requires a pooled TimeSeries ensemble.")
    from dynsyslib.observables import coerce_observable
    obs=coerce_observable(observable,int(series.values.shape[2]))
    groups=tuple(np.asarray(group,dtype=complex).reshape(-1) for group in initial_eigenvalue_groups)
    if len(groups)!=obs.output_dimension or any(group.size<1 for group in groups):
        raise ValueError("Provide one nonempty initial eigenvalue group per observable output.")
    if max_iterations<1 or tolerance<=0. or not np.isfinite(tolerance):
        raise ValueError("Optimization iteration and tolerance settings are invalid.")
    lower,upper=map(float,bounds)
    if not lower<upper: raise ValueError("Eigenvalue bounds must be increasing.")
    states=np.asarray(series.values,dtype=float)
    samples,trajectories,state_dimension=states.shape
    observed=obs.evaluate(states.transpose(1,0,2).reshape(-1,state_dimension))
    observed=observed.reshape(trajectories,samples,obs.output_dimension).transpose(1,0,2)
    optimized=[]; diagnostics=[]; group_indices=[]; offset=0
    for output_index,initial in enumerate(groups):
        parameters=np.concatenate((initial.real,initial.imag))
        result=minimize(lambda p:_objective(p,observed[:,:,output_index]),parameters,jac=True,
            method="L-BFGS-B",bounds=[(lower,upper)]*parameters.size,
            options={"maxiter":int(max_iterations),"ftol":float(tolerance),"gtol":float(tolerance),"maxls":50})
        values=result.x[:initial.size]+1j*result.x[initial.size:]
        if preserve_conjugate_pairs: values=_symmetrize_from_initial(values,initial)
        _,boundary,residual=_powers_and_fit(values,observed[:,:,output_index])
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(observed[:,:,output_index]),np.finfo(float).tiny))
        optimized.extend(values.tolist()); indices=tuple(range(offset,offset+values.size)); group_indices.append(indices); offset+=values.size
        diagnostics.append({"observable_output_index":output_index,"mode_count":int(values.size),
            "success":bool(result.success),"status":int(result.status),"message":str(result.message),
            "iterations":int(result.nit),"function_evaluations":int(result.nfev),
            "normalized_squared_residual":float(relative**2),"relative_residual":relative,
            "boundary_value_shape":list(boundary.shape)})
    dt=float(series.time[1]-series.time[0])
    source=DynamicalSystemReference.from_archive_metadata(metadata,None)
    values=[KoopmanEigenvalue.from_discrete(value,dt,source,{"method":"trajectory-boundary eigenvalue optimization",
        "index":index}) for index,value in enumerate(optimized)]
    stage_metadata={"method":"trajectory-boundary-eigenvalue-optimization","stage_kind":"eigenvalue-extraction",
        "eigenfunctions_constructed":False,"observable":obs.summary(),"mode_groups":[list(group) for group in group_indices],
        "bounds":[lower,upper],"preserve_conjugate_pairs":bool(preserve_conjugate_pairs),
        "optimization_diagnostics":diagnostics,"preprocessing":{"coordinate_system":str((metadata or {}).get("integration_form","dimensional")),
        "sample_count":int(samples*trajectories),"trajectory_count":int(trajectories),"analysis_delta_t":dt}}
    spectrum=KoopmanSpectrum(source,values,[],[],stage_metadata,dynamical_system=system)
    return EigenvalueOptimizationResult(spectrum,tuple(group_indices),stage_metadata)


__all__=["EigenvalueOptimizationResult","optimize_trajectory_eigenvalues"]
