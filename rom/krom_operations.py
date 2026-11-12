"""KROM extraction, control-system rebinding, and autonomous rollout."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any
import numpy as np

from time_series import TimeSeries
from .backends import ROMBackendPreference, select_krom_backend

@dataclass(frozen=True)
class KROMRollout:
    """Autonomous lifted and reconstructed trajectories with provenance."""
    reconstructed: TimeSeries
    lifted: TimeSeries
    metadata: dict[str, Any]


class _KROMOperationsMixin:
    """Independent model operations."""

    def extract(self, selection, *, include_conjugate_duals: bool=True) -> "KROM":
        """Create an independent KROM from selected retained spectral objects.

        Integer selectors index this KROM's nonredundant ``eigenvalues`` /
        ``eigenfunctions`` view. Exact eigenvalue or eigenfunction objects from
        the full source spectrum are also accepted. Conjugate partners are
        included by default when they exist in the source spectrum.
        """
        return extract_krom(self,selection,include_conjugate_duals=include_conjugate_duals)

    def with_analytic_control_system(self, system: Any, *, reason: str) -> "KROM":
        """Rebind analytic ``B_psi`` evaluation without altering source lineage.

        Autonomous spectral objects remain tied to their training archive. This
        explicit derivative is for a control-equivalent system revision, such
        as a corrected input sign, where the autonomous state/time coordinates
        are unchanged. Fitted linear/bilinear control matrices cannot be safely
        rebound this way.
        """
        if self.default_control_representation!="exact":
            raise ValueError("Only an analytic-control KROM can rebind its live control system.")
        if not isinstance(reason,str) or not reason.strip():
            raise ValueError("A nonempty control-system rebinding reason is required.")
        if type(system) is not type(self.system):
            raise TypeError("A rebound KROM control system must have the same concrete system type.")
        source_dimension=int(np.asarray(self.system.initial_state).size)
        target_dimension=int(np.asarray(system.initial_state).size)
        if source_dimension!=target_dimension:
            raise ValueError("A rebound KROM control system must preserve state dimension.")
        source_model=self._model_state(np.asarray(self.system.initial_state,dtype=float),"physical")
        if self.coordinate_system=="dimensionless":
            transform=getattr(system,"to_dimensionless_state",None)
            if not callable(transform):
                raise TypeError("A rebound dimensionless KROM system must provide to_dimensionless_state.")
            target_model=np.asarray(transform(np.asarray(system.initial_state,dtype=float)),dtype=float)
        else:
            target_model=np.asarray(system.initial_state,dtype=float)
        if not np.allclose(source_model,target_model,rtol=1e-12,atol=1e-14):
            raise ValueError("A rebound KROM control system must preserve model-state coordinates.")
        derived=self.extract(range(len(self.eigenvalues)))
        if self.observation_map is not None:
            derived.observation_map=self.observation_map
            derived._validate_observation_map()
        derived.system=system
        derived.spectrum.dynamical_system=system
        from provenance import serializable
        derived.spectrum.metadata["control_system_binding"]={
            "reason":reason.strip(),
            "source_system":serializable(self.system),
            "bound_system":serializable(system),
            "spectral_source_unchanged":True,
        }
        derived._validate_control_representation()
        return derived

    def rollout(self, initial_state: np.ndarray, time: np.ndarray, *,
                execution_backend: ROMBackendPreference="auto",
                physical_state_only: bool=True,state_coordinates: str="model",time_coordinates: str="model",
                output_coordinates: str="model") -> KROMRollout:
        """Propagate an autonomous KROM on arbitrary increasing sample times.

        Controlled rollout remains a separate future contract because analytic
        ``B_psi(x)`` is state dependent.  This method intentionally implements
        only the exact diagonal autonomous lifted flow.
        """
        external_time=np.asarray(time,dtype=float).reshape(-1)
        samples=self._model_time(external_time,time_coordinates)
        if samples.size < 1 or np.any(~np.isfinite(samples)) or np.any(np.diff(samples)<=0):
            raise ValueError("KROM rollout time must contain finite, strictly increasing samples.")
        initial=self._model_state(initial_state,state_coordinates)
        expected=int(np.asarray(self.system.initial_state).size)
        if initial.size != expected or np.any(~np.isfinite(initial)):
            raise ValueError(f"KROM initial_state must be a finite vector of length {expected}.")
        coordinates=self.psi(initial,output="complex",state_coordinates="model")
        eigenvalues=np.asarray([value.continuous_time for value in self.eigenvalues],dtype=complex)
        backend=select_krom_backend(execution_backend)
        if self.observation_map is not None:
            if execution_backend=="native":
                raise ValueError("Native learned-observation-map rollout is not implemented; use execution_backend='python' or 'auto'.")
            backend=select_krom_backend("python")
            lifted=np.exp(np.multiply.outer(samples-samples[0],eigenvalues))*coordinates
            reconstructed=self.reconstruct(lifted,physical_state_only=physical_state_only)
        elif backend.compiled:
            modes=self.mode_matrix(physical_state_only=physical_state_only)
            from importlib import import_module
            lifted,reconstructed=import_module("_dynsyslab_native").krom_rollout(eigenvalues,modes,coordinates,samples)
            lifted=np.asarray(lifted,dtype=complex); reconstructed=np.asarray(reconstructed,dtype=complex)
        else:
            modes=self.mode_matrix(physical_state_only=physical_state_only)
            lifted=np.exp(np.multiply.outer(samples-samples[0],eigenvalues))*coordinates
            reconstructed=lifted@modes.T
        reconstructed=np.real_if_close(reconstructed,tol=1000)
        if output_coordinates=="physical":
            if self.coordinate_system=="dimensionless":
                transform=getattr(self.system,"to_dimensional_state",None)
                if not callable(transform): raise NotImplementedError("This KROM cannot transform reconstructed model states to physical coordinates.")
                reconstructed=np.asarray([transform(row) for row in np.asarray(reconstructed,dtype=float)])
        elif output_coordinates!="model": raise ValueError("output_coordinates must be 'model' or 'physical'.")
        metadata={**self.metadata,"rollout":"autonomous-diagonal","sample_count":int(samples.size),
                  "physical_state_only":bool(physical_state_only),"state_coordinates":state_coordinates,"time_coordinates":time_coordinates,
                  "output_coordinates":output_coordinates,"execution_backend":backend.provenance(execution_backend)}
        return KROMRollout(TimeSeries(external_time,reconstructed),TimeSeries(external_time,lifted),metadata)

def extract_krom(source: KROM, selection, *, include_conjugate_duals: bool=True) -> KROM:
    """Extract an immutable spectral subset into a new, independent KROM."""
    from .koopman import KROM
    if not isinstance(source,KROM): raise TypeError("source must be a KROM.")
    try: requested=list(selection)
    except TypeError as error: raise TypeError("selection must be an iterable of KROM indices, eigenvalues, or eigenfunctions.") from error
    if not requested: raise ValueError("KROM extraction selection cannot be empty.")
    full_values=list(source.spectrum.eigenvalues); full_functions=list(source.spectrum.eigenfunctions)
    selected=[]
    for item in requested:
        if isinstance(item,(int,np.integer)) and not isinstance(item,(bool,np.bool_)):
            position=int(item)
            if position<0 or position>=len(source._retained_indices):
                raise IndexError(f"KROM selection index {position} is outside [0, {len(source._retained_indices)}).")
            full_index=source._retained_indices[position]
        else:
            full_index=next((index for index,value in enumerate(full_values) if value is item),None)
            if full_index is None: full_index=next((index for index,function in enumerate(full_functions) if function is item),None)
            if full_index is None: raise ValueError("A selected eigenvalue/eigenfunction does not belong to the source KROM spectrum.")
        if full_index not in selected: selected.append(full_index)
    expanded=set(selected)
    if include_conjugate_duals:
        for index in tuple(selected):
            dual=full_values[index].dual
            if dual is not None: expanded.add(next(other for other,value in enumerate(full_values) if value is dual))
    indices=tuple(index for index in range(len(full_values)) if index in expanded)
    if source.spectrum.modes and len(source.spectrum.modes)!=len(full_values):
        raise ValueError("KROM extraction cannot preserve reconstruction modes unless one mode is aligned with every full-spectrum eigenvalue.")

    from spectral.contracts import KoopmanEigenfunction, KoopmanEigenvalue, KoopmanSpectrum
    cloned_values=[]; cloned_functions=[]
    for index in indices:
        original=full_values[index]; value_metadata=deepcopy(original.metadata)
        value_metadata.pop("dual_index",None); value_metadata["extraction_source_index"]=index
        cloned_values.append(KoopmanEigenvalue(original.discrete_time,original.continuous_time,original.delta_t,
                                               original.original_domain,original.system,value_metadata))
    for new_index,source_index in enumerate(indices):
        original=full_functions[source_index]; function_metadata=deepcopy(original.metadata)
        function_metadata.pop("dual_index",None); function_metadata["extraction_source_index"]=source_index
        cloned_functions.append(KoopmanEigenfunction(cloned_values[new_index],original.representation,
            values=None if original.values is None else np.array(original.values,copy=True),
            coefficients=None if original.coefficients is None else np.array(original.coefficients,copy=True),
            state_callable=original.state_callable,gradient_callable=original.gradient_callable,
            residual_operator=original.residual_operator,metadata=function_metadata,
            representation_kind=original.representation_kind,
            error_evaluations=deepcopy(original.error_evaluations)))
    modes=[] if not source.spectrum.modes else [np.array(source.spectrum.modes[index],copy=True) for index in indices]
    spectrum_metadata=deepcopy(source.spectrum.metadata); source_spectrum_id=spectrum_metadata.pop("spectrum_id",None)
    spectrum_metadata["extraction"]={"source_spectrum_id":source_spectrum_id,"source_method":source.method,
        "requested_source_indices":selected,"selected_source_indices":list(indices),
        "include_conjugate_duals":bool(include_conjugate_duals),"source_krom_retained_indices":list(source._retained_indices)}
    spectrum=KoopmanSpectrum(source.spectrum.system,cloned_values,cloned_functions,modes,spectrum_metadata,
                             dynamical_system=source.system)
    derived=KROM(spectrum,source.system,coordinate_observable=source.coordinate_observable)

    positions={full_index:position for position,full_index in enumerate(source._retained_indices)}
    derived_source_indices=[value.metadata["extraction_source_index"] for value in derived.eigenvalues]
    retained_positions=[]; conjugated=[]
    for index in derived_source_indices:
        if index in positions:
            retained_positions.append(positions[index]); conjugated.append(False); continue
        dual=full_values[index].dual
        dual_index=None if dual is None else next((other for other,value in enumerate(full_values) if value is dual),None)
        if dual_index not in positions:
            raise ValueError("KROM extraction selected a coordinate absent from its nonredundant representation.")
        retained_positions.append(positions[dual_index]); conjugated.append(True)
    if any(conjugated) and source.control_representations:
        raise ValueError("Conjugate-only extraction of fitted control representations is not supported; include conjugate duals.")
    restricted={name:representation.subset(retained_positions)
        for name,representation in (source.control_representations or {}).items()}
    derived.control_representations=restricted
    if source.default_control_representation=="exact" or source.default_control_representation in derived.control_representations:
        derived.default_control_representation=source.default_control_representation
    derived._validate_control_representation()
    if source.observable_reconstructions:
        from .observation import ObservableReconstruction
        restricted_reconstructions={}
        for name,source_reconstruction in source.observable_reconstructions.items():
            reconstruction=source_reconstruction.subset(retained_positions)
            if any(conjugated):
                coefficients=np.array(reconstruction.coefficients,copy=True)
                for position,needs_conjugate in enumerate(conjugated):
                    if needs_conjugate: coefficients[:,position]=np.conjugate(coefficients[:,position])
                reconstruction=ObservableReconstruction(reconstruction.observable,coefficients,
                    reconstruction.intercept,reconstruction.quality,
                    {**reconstruction.provenance,"conjugate_only_columns":conjugated})
            restricted_reconstructions[name]=reconstruction
        derived.observable_reconstructions=restricted_reconstructions
        derived.default_observable_reconstruction=source.default_observable_reconstruction
        derived.observable_reconstruction=(None if derived.default_observable_reconstruction is None else
                                           restricted_reconstructions[derived.default_observable_reconstruction])
        derived._validate_observable_reconstruction()
    # A validated state decoder is a promoted whole-model capability.  Mode
    # extraction invalidates it; fit a new decoder for the derived KROM.
    derived.observation_map=None
    derived.spectrum.metadata["extraction"]["derived_retained_source_indices"]=[
        value.metadata["extraction_source_index"] for value in derived.eigenvalues]
    return derived
