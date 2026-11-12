"""State and declared-observable reconstruction capabilities."""
from __future__ import annotations

from typing import Any
import numpy as np

class _KROMObservationMixin:
    """State and declared-observable reconstruction capabilities."""

    @property
    def state_representation(self) -> str:
        if self.observation_map is not None: return "learned-state-observation-map"
        embedding=self.spectrum.metadata.get("delay_embedding") or {}
        if embedding.get("type")=="hankel": return "delay-observable"
        if not self.spectrum.modes: return "eigenfunction-only"
        dimension=int(np.asarray(self.system.initial_state).size)
        if len(self.spectrum.modes)==len(self.spectrum.eigenvalues) and all(np.asarray(mode).size==dimension for mode in self.spectrum.modes):
            return "state-modes"
        return "observable-modes"

    def _validate_observation_map(self) -> None:
        if self.observation_map is None: return
        from .observation import StateObservationMap
        if not isinstance(self.observation_map,StateObservationMap):
            raise TypeError("observation_map must be a StateObservationMap.")
        split_dimension=sum(1 if function.is_real else 2 for function in self.eigenfunctions)
        state_dimension=int(getattr(self.system,"state_dimension",0) or np.asarray(self.system.initial_state).size)
        if self.observation_map.lifted_feature_dimension!=split_dimension:
            raise ValueError("Observation map lifted dimension does not match the KROM.")
        if self.observation_map.state_dimension!=state_dimension:
            raise ValueError("Observation map state dimension does not match the KROM system.")

    def _validate_observable_reconstruction(self) -> None:
        from .observation import ObservableReconstruction
        reconstructions=dict(self.observable_reconstructions or {})
        if self.observable_reconstruction is not None:
            if not isinstance(self.observable_reconstruction,ObservableReconstruction):
                raise TypeError("observable_reconstruction must be an ObservableReconstruction.")
            name=self.observable_reconstruction.observable.name
            reconstructions[name]=self.observable_reconstruction
            if self.default_observable_reconstruction is None:
                self.default_observable_reconstruction=name
        for name,reconstruction in reconstructions.items():
            if not str(name).strip():
                raise ValueError("Observable reconstruction names must be nonempty.")
            if not isinstance(reconstruction,ObservableReconstruction):
                raise TypeError("Every observable_reconstructions value must be an ObservableReconstruction.")
            if reconstruction.lifted_dimension!=self.lifted_dimension:
                raise ValueError(f"Observable reconstruction {name!r} does not match the KROM dimension.")
        if self.default_observable_reconstruction is not None:
            if self.default_observable_reconstruction not in reconstructions:
                raise KeyError("default_observable_reconstruction is not present in observable_reconstructions.")
            self.observable_reconstruction=reconstructions[self.default_observable_reconstruction]
        elif reconstructions:
            self.default_observable_reconstruction=next(iter(reconstructions))
            self.observable_reconstruction=reconstructions[self.default_observable_reconstruction]
        self.observable_reconstructions=reconstructions

    def with_observation_map(self, observation_map) -> "KROM":
        """Return an independent KROM carrying a validated state decoder."""
        derived=self.extract(range(len(self.eigenvalues)))
        derived.observation_map=observation_map
        derived._validate_observation_map()
        derived.spectrum.metadata["observation_map"]=observation_map.summary()
        return derived

    def with_observable_reconstruction(self, reconstruction, *, name: str|None=None,
                                       make_default: bool=True) -> "KROM":
        """Return an independent KROM carrying a declared-observable C map."""
        derived=self.extract(range(len(self.eigenvalues)))
        key=str(name or reconstruction.observable.name)
        if not key.strip(): raise ValueError("Observable reconstruction name must be nonempty.")
        derived.observable_reconstructions={**derived.observable_reconstructions,key:reconstruction}
        if make_default:
            derived.default_observable_reconstruction=key
            derived.observable_reconstruction=reconstruction
        derived._validate_observable_reconstruction()
        derived.spectrum.metadata["observable_reconstructions"]={
            item:value.summary() for item,value in derived.observable_reconstructions.items()}
        derived.spectrum.metadata["default_observable_reconstruction"]=derived.default_observable_reconstruction
        derived.spectrum.metadata["observable_reconstruction"]=derived.observable_reconstruction.summary()
        return derived

    def observable_map(self, name: str):
        """Return one named observable projection without changing KROM state."""
        key=str(name)
        if key not in self.observable_reconstructions:
            raise KeyError(f"Observable reconstruction {key!r} is unavailable; available={sorted(self.observable_reconstructions)}.")
        return self.observable_reconstructions[key]

    def mode_matrix(self, *, physical_state_only: bool=True) -> np.ndarray:
        """Return reconstruction modes aligned with retained eigenfunctions."""
        if self.observation_map is not None:
            raise ValueError("A learned state observation map is affine and has no complex Koopman mode matrix; use reconstruct().")
        if physical_state_only and self.state_representation!="state-modes":
            raise ValueError(f"KROM state reconstruction is unavailable for {self.state_representation} representation.")
        if len(self.spectrum.modes) != len(self.spectrum.eigenvalues):
            raise ValueError("KROM reconstruction requires one mode per full-spectrum eigenvalue.")
        matrix=np.column_stack([np.asarray(self.spectrum.modes[index],dtype=complex).reshape(-1) for index in self._retained_indices])
        state_dimension=int(getattr(self.system,"state_dimension",0) or np.asarray(self.system.initial_state).size)
        if physical_state_only:
            if matrix.shape[0] < state_dimension:
                raise ValueError("Koopman modes have fewer rows than the physical state dimension.")
            matrix=matrix[:state_dimension]
        return matrix

    def reconstruct(self, lifted_state: np.ndarray, *, physical_state_only: bool=True) -> np.ndarray:
        """Map one or many retained lifted-coordinate rows to state space."""
        coordinates=np.asarray(lifted_state,dtype=complex)
        if coordinates.shape[-1] != self.lifted_dimension:
            raise ValueError("The final lifted-state axis must match KROM.lifted_dimension.")
        if self.observation_map is not None:
            if not physical_state_only:
                raise ValueError("A learned state observation map reconstructs physical state coordinates only.")
            return self.observation_map.reconstruct(coordinates,self.eigenfunctions)
        values=coordinates@self.mode_matrix(physical_state_only=physical_state_only).T
        return np.real_if_close(values,tol=1000)

    def C(self, observable: str|None=None, *, output: str="split") -> np.ndarray:
        """Return the retained observable reconstruction modes."""
        reconstruction=(self.observable_reconstruction if observable is None
                        else self.observable_map(observable))
        if reconstruction is None:
            raise ValueError("This KROM has no declared observable reconstruction C.")
        if output=="complex": return np.asarray(reconstruction.coefficients)
        if output=="split": return reconstruction.split_matrix(self.eigenfunctions)
        raise ValueError("output must be 'complex' or 'split'.")

    def reconstruct_observable(self,lifted_state: np.ndarray, observable: str|None=None) -> np.ndarray:
        reconstruction=(self.observable_reconstruction if observable is None
                        else self.observable_map(observable))
        if reconstruction is None:
            raise ValueError("This KROM has no declared observable reconstruction C.")
        return reconstruction.reconstruct(lifted_state,self.eigenfunctions)
