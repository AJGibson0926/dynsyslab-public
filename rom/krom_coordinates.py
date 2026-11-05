"""Coordinate transforms, eigenfunction evaluation, and lifted generators."""
from __future__ import annotations

from typing import Any
import numpy as np

from .control_inputs import complex_to_split

class _KROMCoordinateMixin:
    """Coordinate transforms, eigenfunction evaluation, and lifted generators."""

    @property
    def coordinate_system(self) -> str:
        preprocessing=self.spectrum.metadata.get("preprocessing") or {}
        return str(preprocessing.get("coordinate_system","dimensional"))

    @property
    def training_domain(self) -> dict[str,Any] | None:
        domain=self.spectrum.metadata.get("training_domain")
        return None if not isinstance(domain,dict) else dict(domain)

    def domain_report(self, state: np.ndarray, *, state_coordinates: str="model") -> dict[str,Any]:
        """Report axis-aligned training-domain membership without mutation."""
        domain=self.training_domain
        if domain is None:
            return {"available":False,"inside":None,"state_coordinates":"model","violating_axes":[]}
        if domain.get("encoding")!="axis-aligned-training-domain/v1":
            raise ValueError("Unsupported KROM training-domain encoding.")
        model_state=self._model_state(state,state_coordinates)
        minimum=np.asarray(domain.get("minimum"),dtype=float)
        maximum=np.asarray(domain.get("maximum"),dtype=float)
        if minimum.shape!=model_state.shape or maximum.shape!=model_state.shape or np.any(maximum<minimum):
            raise ValueError("KROM training-domain bounds do not match model state coordinates.")
        scale=np.maximum(maximum-minimum,np.finfo(float).eps)
        tolerance=64.*np.finfo(float).eps*np.maximum(1.,np.maximum(np.abs(minimum),np.abs(maximum)))
        below=model_state<minimum-tolerance; above=model_state>maximum+tolerance
        violation=np.maximum((minimum-model_state)/scale,(model_state-maximum)/scale)
        violation=np.maximum(violation,0.)
        axes=np.flatnonzero(below|above)
        return {"available":True,"inside":not bool(axes.size),"state_coordinates":"model",
            "state":model_state.tolist(),"minimum":minimum.tolist(),"maximum":maximum.tolist(),
            "violating_axes":axes.astype(int).tolist(),"relative_excursion":violation.tolist(),
            "maximum_relative_excursion":float(np.max(violation))}

    def _time_scale(self) -> float:
        if self.coordinate_system!="dimensionless": return 1.
        parameters=getattr(self.system,"dimensionless_parameters",None)
        scale=getattr(parameters,"time_scale_s",None)
        if scale is None or not np.isfinite(scale) or scale<=0: raise ValueError("Dimensionless KROM needs a positive system time_scale_s.")
        return float(scale)

    def _model_state(self,state: np.ndarray, coordinates: str) -> np.ndarray:
        value=np.asarray(state,dtype=float).reshape(-1)
        if coordinates=="model" or self.coordinate_system=="dimensional" and coordinates=="physical": return value
        if coordinates!="physical": raise ValueError("state_coordinates must be 'model' or 'physical'.")
        transform=getattr(self.system,"to_dimensionless_state",None)
        if not callable(transform): raise NotImplementedError("This KROM has dimensionless eigenfunctions but its system cannot transform physical states.")
        return np.asarray(transform(value),dtype=float).reshape(-1)

    def _model_time(self,time: np.ndarray, coordinates: str) -> np.ndarray:
        value=np.asarray(time,dtype=float)
        if coordinates=="model" or self.coordinate_system=="dimensional" and coordinates=="physical": return value
        if coordinates!="physical": raise ValueError("time_coordinates must be 'model' or 'physical'.")
        return value/self._time_scale()

    def _input_matrix(self,model_state: np.ndarray,time_coordinates: str) -> np.ndarray | None:
        if self.coordinate_system=="dimensionless":
            function=getattr(self.system,"dimensionless_input_matrix",None)
            if not callable(function): raise NotImplementedError("Dimensionless analytic control requires system.dimensionless_input_matrix(x).")
            matrix=np.asarray(function(model_state),dtype=float)
            if time_coordinates=="physical": matrix=matrix/self._time_scale()
            elif time_coordinates!="model": raise ValueError("time_coordinates must be 'model' or 'physical'.")
            return matrix
        if time_coordinates not in {"model","physical"}: raise ValueError("time_coordinates must be 'model' or 'physical'.")
        function=getattr(self.system,"input_matrix",None)
        return None if not callable(function) else function(model_state)

    @property
    def input_dimension(self) -> int:
        probe=self._model_state(np.asarray(self.system.initial_state,dtype=float),"physical")
        matrix=self._input_matrix(probe,"model")
        if matrix is None: raise ValueError("This KROM's system has no compatible input matrix.")
        return int(np.asarray(matrix,dtype=float).shape[1])

    @property
    def eigenvalues(self):
        """Nonredundant spectral modes retained by this KROM."""
        return [self.spectrum.eigenvalues[index] for index in self._retained_indices]

    @property
    def eigenfunctions(self):
        """Eigenfunctions paired with :attr:`eigenvalues` in KROM order."""
        return [self.spectrum.eigenfunctions[index] for index in self._retained_indices]

    def psi(self,state: np.ndarray, *, output: str="complex", state_coordinates: str="model") -> np.ndarray:
        x=self._model_state(state,state_coordinates)
        if output=="complex": return np.asarray([item.psi_from_state(x) for item in self.eigenfunctions],dtype=complex)
        if output=="split": return np.concatenate([np.array((item.psi_from_state(x).real,)) if item.is_real else item.psi_from_state(x,output="split") for item in self.eigenfunctions]).astype(float,copy=False)
        raise ValueError("output must be 'complex' or 'split'.")

    def nabla_psi(self,state: np.ndarray, *, output: str="complex", state_coordinates: str="model",
                  gradient_backend: str="auto", observable_backend: str="auto") -> np.ndarray:
        """Return ``nabla_x Psi = nabla_y Psi J_h`` for ``y=h(x)``."""
        x=self._model_state(state,state_coordinates)
        complex_result=np.vstack([item.nabla_psi_from_state(x,backend=gradient_backend,
            observable_backend=observable_backend) for item in self.eigenfunctions])
        if output=="complex": return complex_result
        if output=="split": return complex_to_split(complex_result,self.eigenfunctions)
        raise ValueError("output must be 'complex' or 'split'.")

    def A(self, *, output: str="complex", time_coordinates: str="model") -> np.ndarray:
        """Return the continuous lifted generator.

        A complex eigenvalue ``a + i b`` becomes the real block
        ``[[a, -b], [b, a]]`` in split coordinates.  Its imaginary-frequency
        contribution is skew-symmetric; the repeated ``a`` retains physical
        exponential growth or decay.  All entries outside these 2-by-2 blocks
        are zero.
        """
        if time_coordinates not in {"model","physical"}: raise ValueError("time_coordinates must be 'model' or 'physical'.")
        values=np.asarray([item.continuous_time for item in self.eigenvalues],dtype=complex)
        if time_coordinates=="physical" and self.coordinate_system=="dimensionless": values=values/self._time_scale()
        if output=="complex": return np.diag(values)
        if output=="split":
            dimension=sum(1 if value.is_real else 2 for value in self.eigenvalues); result=np.zeros((dimension,dimension),dtype=float); offset=0
            for value in values:
                if np.isclose(value.imag,0.): result[offset,offset]=value.real; offset+=1
                else: result[offset:offset+2,offset:offset+2]=((value.real,-value.imag),(value.imag,value.real)); offset+=2
            return result
        raise ValueError("output must be 'complex' or 'split'.")
