"""Fitted-boundary eigenfunctions for externally supplied eigenvalues."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

from time_series import TimeSeries
from .spectrum import KoopmanEigenfunction, KoopmanSpectrum, ScatteredEigenfunctionInterpolator


@dataclass
class BoundaryEigenfunctionResult:
    spectrum: KoopmanSpectrum
    training_states: np.ndarray
    training_eigenfunction_values: np.ndarray
    boundary_values: tuple[np.ndarray,...]
    metadata: dict[str,Any]
    krom: Any


def fit_boundary_eigenfunctions(
    eigenvalue_spectrum: KoopmanSpectrum,
    series: TimeSeries,
    *,
    observable,
    mode_groups: Sequence[Sequence[int]] | None=None,
) -> BoundaryEigenfunctionResult:
    """Fit trajectory boundary values and construct scattered eigenfunctions.

    The eigenvalues are immutable inputs to this stage.  For every observable
    output, only the initial boundary values of its assigned eigenfunctions are
    regressed; temporal values then follow exactly as ``mu**k * g0``.
    """
    if not isinstance(series,TimeSeries) or not series.is_pooled:
        raise TypeError("Boundary-eigenfunction fitting requires a pooled TimeSeries ensemble.")
    from dynsyslib.observables import coerce_observable
    from rom import KROM
    from rom.observation import ObservableReconstruction
    obs=coerce_observable(observable,int(series.values.shape[2]))
    groups=mode_groups if mode_groups is not None else eigenvalue_spectrum.metadata.get("mode_groups")
    if groups is None:
        if obs.output_dimension!=1: raise ValueError("Multi-output boundary fitting requires one mode group per output.")
        groups=(tuple(range(len(eigenvalue_spectrum.eigenvalues))),)
    groups=tuple(tuple(int(index) for index in group) for group in groups)
    if len(groups)!=obs.output_dimension or any(not group for group in groups):
        raise ValueError("mode_groups must contain one nonempty group per observable output.")
    flattened=[index for group in groups for index in group]
    if sorted(flattened)!=list(range(len(eigenvalue_spectrum.eigenvalues))) or len(set(flattened))!=len(flattened):
        raise ValueError("mode_groups must partition every supplied eigenvalue exactly once.")
    states=np.asarray(series.values,dtype=float)
    samples,trajectories,state_dimension=states.shape
    training_states=states.transpose(1,0,2).reshape(-1,state_dimension)
    observed=obs.evaluate(training_states).reshape(trajectories,samples,obs.output_dimension).transpose(1,0,2)
    function_values=np.empty((samples*trajectories,len(flattened)),dtype=complex)
    boundaries=[]; diagnostics=[]
    for output_index,group in enumerate(groups):
        eigenvalues=np.asarray([eigenvalue_spectrum.eigenvalues[index].discrete_time for index in group])
        powers=np.power(eigenvalues[None,:],np.arange(samples,dtype=float)[:,None])
        boundary=np.linalg.lstsq(powers,observed[:,:,output_index],rcond=1e-12)[0]
        predicted=powers@boundary; residual=observed[:,:,output_index]-predicted
        for local,index in enumerate(group):
            modal_values=powers[:,local,None]*boundary[local,None,:]
            function_values[:,index]=modal_values.T.reshape(-1)
        boundaries.append(boundary)
        diagnostics.append({"observable_output_index":output_index,"mode_indices":list(group),
            "boundary_value_shape":list(boundary.shape),
            "relative_residual":float(np.linalg.norm(residual)/max(np.linalg.norm(observed[:,:,output_index]),np.finfo(float).tiny))})
    functions=[]
    for index,(value,sampled) in enumerate(zip(eigenvalue_spectrum.eigenvalues,function_values.T)):
        interpolator=ScatteredEigenfunctionInterpolator(training_states,sampled)
        functions.append(KoopmanEigenfunction(value,"fitted boundary values with scattered state interpolation",
            values=sampled,state_callable=interpolator,representation_kind="scattered-interpolation",
            metadata={"method":"fitted-boundary-eigenfunctions","index":index,
                "interpolator":interpolator.summary(),"gradient_backend":"finite_difference"}))
    metadata={"method":"fitted-boundary-eigenfunctions","stage_kind":"eigenfunction-construction",
        "eigenvalue_algorithm":eigenvalue_spectrum.metadata.get("method"),"mode_groups":[list(group) for group in groups],
        "observable":obs.summary(),"fit_diagnostics":diagnostics,"preprocessing":{"coordinate_system":
        str((eigenvalue_spectrum.metadata.get("preprocessing") or {}).get("coordinate_system","dimensional")),
        "sample_count":int(training_states.shape[0]),"trajectory_count":int(trajectories),
        "analysis_delta_t":float(series.time[1]-series.time[0])},
        "training_domain":{"encoding":"axis-aligned-training-domain/v1","minimum":np.min(training_states,axis=0).tolist(),
        "maximum":np.max(training_states,axis=0).tolist(),"sample_count":int(training_states.shape[0]),
        "coordinate_system":str((eigenvalue_spectrum.metadata.get("preprocessing") or {}).get("coordinate_system","dimensional")),
        "source":"fitted-boundary training states"}}
    spectrum=KoopmanSpectrum(eigenvalue_spectrum.system,list(eigenvalue_spectrum.eigenvalues),functions,[],metadata,
        dynamical_system=eigenvalue_spectrum.dynamical_system,eigenvalue_source=eigenvalue_spectrum)
    krom=KROM(spectrum,eigenvalue_spectrum.dynamical_system)
    full_C=np.zeros((obs.output_dimension,len(functions)),dtype=complex)
    for output_index,group in enumerate(groups): full_C[output_index,list(group)]=1.
    C=full_C[:,list(krom._retained_indices)]
    reconstruction=ObservableReconstruction(obs,C,np.zeros(obs.output_dimension),
        {"passed":True,"reason":"boundary functions were jointly regressed to this observable",
         "training_relative_residual_by_output":[item["relative_residual"] for item in diagnostics]},
        {"method":"fitted-boundary modal coefficients","definition":"unit nu per assigned boundary eigenfunction"})
    krom.observable_reconstruction=reconstruction; krom._validate_observable_reconstruction()
    spectrum.metadata["observable_reconstruction"]=reconstruction.summary()
    return BoundaryEigenfunctionResult(spectrum,training_states,function_values,tuple(boundaries),metadata,krom)


__all__=["BoundaryEigenfunctionResult","fit_boundary_eigenfunctions"]
