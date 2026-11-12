"""Koopman reduced-order model facade with focused capabilities."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal
import numpy as np

from .krom_control import _KROMControlMixin
from .krom_coordinates import _KROMCoordinateMixin
from .krom_observation import _KROMObservationMixin
from .krom_operations import (KROMRollout, _KROMOperationsMixin, extract_krom)


@dataclass
class KROM(_KROMCoordinateMixin, _KROMObservationMixin, _KROMControlMixin, _KROMOperationsMixin):
    """Spectrum-based Koopman ROM with independently gated capabilities.

    For ``xdot=f(x)+B(x)u``, the lifted input map is assembled as
    ``B_psi(x) = nabla(Psi)(x) @ B(x)``.  This follows the control-affine
    eigenfunction-coordinate construction used by Kaiser--Kutz--Brunton.
    """
    spectrum: Any
    system: Any | None = None
    control_mode: Literal["analytic","bilinear","linear"] = "analytic"
    B1: np.ndarray | None = None
    B2: np.ndarray | None = None
    observation_map: Any | None = None
    coordinate_observable: Any | None = None
    observable_reconstruction: Any | None = None
    control_representations: dict[str, Any] | None = None
    default_control_representation: str | None = None
    control_error_evaluations: dict[str, Any] = field(default_factory=dict, repr=False)
    observable_reconstructions: dict[str, Any] | None = None
    default_observable_reconstruction: str | None = None

    def __post_init__(self):
        if self.system is None:
            self.system=getattr(self.spectrum,"dynamical_system",None)
        if self.system is None:
            raise ValueError("KROM requires the specific dynamical-system object associated with its spectrum.")
        from dynsyslib.observables import Observable
        state_dimension=int(getattr(self.system,"state_dimension",0) or np.asarray(self.system.initial_state).size)
        if self.coordinate_observable is None:
            self.coordinate_observable=Observable.identity(state_dimension)
        elif not isinstance(self.coordinate_observable,Observable):
            raise TypeError("coordinate_observable must be a dynsyslib Observable.")
        elif self.coordinate_observable.state_dimension!=state_dimension:
            raise ValueError("coordinate_observable state dimension does not match the KROM system.")
        for function in self.spectrum.eigenfunctions:
            if function.coordinate_observable is None:
                function.coordinate_observable=self.coordinate_observable
            elif function.coordinate_observable.summary()!=self.coordinate_observable.summary():
                raise ValueError("KROM eigenfunctions do not share the declared coordinate observable.")
        self.spectrum.pair_conjugates()
        functions=list(self.spectrum.eigenfunctions)
        if not functions or any(not item.has_state_callable for item in functions):
            raise ValueError("KROM requires callable eigenfunctions for every retained mode.")
        dt={float(item.delta_t) for item in self.spectrum.eigenvalues}
        if len(dt)!=1: raise ValueError("KROM requires a common delta_t across spectrum eigenvalues.")
        retained=[]
        for index,value in enumerate(self.spectrum.eigenvalues):
            if value.is_real or value.dual is None or value.continuous_time.imag > 0.: retained.append(index)
        self._retained_indices=tuple(sorted(retained,key=lambda index:(0 if self.spectrum.eigenvalues[index].is_real else 1,
            self.spectrum.eigenvalues[index].growth_decay_rate,self.spectrum.eigenvalues[index].angular_frequency,index)))
        self._validate_control_representation()
        self._validate_observation_map()
        self._validate_observable_reconstruction()
        pending=dict(self.control_error_evaluations); self.control_error_evaluations={}
        for key,evaluation in pending.items():
            if str(key)!=evaluation.evaluation_id:
                raise ValueError("KROM control-error evaluation keys must equal their evaluation IDs.")
            self.attach_control_error_evaluation(evaluation)

    @property
    def lifted_dimension(self) -> int: return len(self.eigenfunctions)

    def capabilities(self) -> dict[str,Any]:
        """Declare independently gated KROM operations."""
        gradients=all(item.has_gradient for item in self.eigenfunctions)
        analytic_gradients=all(item.gradient_callable is not None for item in self.eigenfunctions)
        try: control=self._input_matrix(self._model_state(np.asarray(self.system.initial_state,dtype=float),"physical"),"model") is not None
        except (AttributeError,NotImplementedError,TypeError,ValueError): control=False
        reconstruction=self.state_representation in {"state-modes","learned-state-observation-map"}
        validated=bool(self.observation_map is not None and self.observation_map.validated)
        return {"lifted_evolution":True,"state_reconstruction":reconstruction,
                "validated_state_reconstruction":validated,
                "exact_control_coupling":bool(gradients and control),
                "analytic_control_coupling":bool(analytic_gradients and control and self.coordinate_observable.has_analytic_jacobian),
                "eigenfunction_gradient_backends":[item.gradient_backend for item in self.eigenfunctions],
                "control_representations":self.available_control_representations(),"coordinate_system":self.coordinate_system,
                "state_representation":self.state_representation}

    def require_capability(self, capability: str) -> None:
        """Raise when an explicitly named KROM capability is unavailable."""
        available=self.capabilities()
        if capability not in available: raise KeyError(f"Unknown KROM capability {capability!r}.")
        if not bool(available[capability]):
            raise ValueError(f"KROM capability {capability!r} is unavailable; capabilities={available}.")

    @property
    def method(self) -> str: return str(self.spectrum.metadata.get("method") or "unspecified")

    @property
    def system_reference(self): return self.spectrum.system

    @property
    def delta_t(self) -> float: return float(self.spectrum.eigenvalues[0].delta_t)

    @property
    def metadata(self) -> dict[str,Any]:
        return {"method":self.method,"system":self.system_reference,"delta_t":self.delta_t,
                "eigenvalue_domain":"continuous","state_dimension":int(getattr(self.system,"state_dimension",0)),
                "eigenfunction_count":len(self.eigenfunctions),"full_spectrum_eigenfunction_count":len(self.spectrum.eigenfunctions),
                "retained_indices":list(self._retained_indices),"control_mode":self.control_mode,
                "default_control_representation":self.default_control_representation,
                "available_control_representations":self.available_control_representations(),
                "control_representations":{key:value.summary() for key,value in (self.control_representations or {}).items()},
                "control_error_evaluations":[value.summary() for value in self.control_error_evaluations.values()],
                "B1_shape":None if self.B1 is None else list(self.B1.shape),"B2_shape":None if self.B2 is None else list(self.B2.shape),
                "coordinate_system":self.coordinate_system,"state_representation":self.state_representation,
                "coordinate_observable":self.coordinate_observable.summary(),
                "observable_reconstruction":None if self.observable_reconstruction is None else self.observable_reconstruction.summary(),
                "default_observable_reconstruction":self.default_observable_reconstruction,
                "observable_reconstructions":{name:value.summary() for name,value in self.observable_reconstructions.items()},
                "observation_map":None if self.observation_map is None else self.observation_map.summary(),
                "capabilities":self.capabilities(),"spectrum_metadata":self.spectrum.metadata}
