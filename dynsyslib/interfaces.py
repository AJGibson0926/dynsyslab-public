"""Uniform plug-in contracts for dynsyslab system-library entries."""
from __future__ import annotations
from typing import Protocol
import numpy as np

class LibrarySystem:
    """Common metadata exposed by every plug-in, irrespective of numerical kind."""
    kind: str = "unknown"
    name: str
    @property
    def archive_name(self) -> str: return self.name
    def describe(self) -> dict: return {"name": self.name, "kind": self.kind}

class ODESystem(LibrarySystem):
    kind = "ode"
    @property
    def state_dimension(self) -> int: return int(np.asarray(self.initial_state).size)
    def state_component_names(self, coordinate_system: str | None=None) -> tuple[str,...]:
        """Stable names used by the guaranteed state-observable catalog."""
        return tuple(f"x{index}" for index in range(self.state_dimension))
    def state_component_labels(self, coordinate_system: str | None=None) -> tuple[str,...]:
        return self.state_component_names(coordinate_system)
    def model_observables(self, coordinate_system: str | None=None) -> dict:
        """System-specific scalar observables; state entries are added centrally."""
        return {}
    def observable_aliases(self, coordinate_system: str | None=None) -> dict[str,str]:
        return {}
    def observable_coordinate_presets(self, coordinate_system: str | None=None) -> dict[str,tuple[str,...]]:
        """Named ordered observable coordinates available to Koopman methods."""
        return {}
    def observable_catalog(self, coordinate_system: str | None=None):
        """Return state components, full state, and model-specific observables."""
        from .observables import state_observable_catalog
        selected=coordinate_system or "native"
        return state_observable_catalog(
            self.state_dimension,
            component_names=self.state_component_names(selected),
            component_labels=self.state_component_labels(selected),
            coordinate_system=selected,
            extras=self.model_observables(selected),
            aliases=self.observable_aliases(selected),
        )
    def observable(self, name: str, *, coordinate_system: str | None=None):
        return self.observable_catalog(coordinate_system)[name]
    def observable_coordinates(self, selection, *, coordinate_system: str | None=None,
                               name: str | None=None):
        """Resolve an ordered observable coordinate map ``y=h(x)``."""
        from .observables import coerce_coordinate_observable
        return coerce_coordinate_observable(selection,self.state_dimension,system=self,
            coordinate_system=coordinate_system,name=name)
    def input_matrix(self, state: np.ndarray) -> np.ndarray | None:
        """Return constant ``B`` or evaluate a state-dependent ``B(x)``.

        Library systems expose ``control_matrix`` as either a constant matrix
        property or a callable ``control_matrix(x)``.  This normalizer keeps
        solver and equation implementations agnostic to that distinction.
        """
        matrix=getattr(self,"control_matrix",None)
        if matrix is None: return None
        return np.asarray(matrix(np.asarray(state,dtype=float)) if callable(matrix) else matrix,dtype=float)

class DelaySystem(LibrarySystem):
    kind = "dde"

class LatticeSystem(LibrarySystem):
    kind = "lattice"

class ForcedSystem(Protocol):
    name: str
    initial_state: np.ndarray
    period: float
    def rhs(self, t: float, x: np.ndarray) -> np.ndarray: ...
    def jacobian(self, t: float, x: np.ndarray) -> np.ndarray: ...
    def invalid(self, t: float, x: np.ndarray) -> bool: ...
