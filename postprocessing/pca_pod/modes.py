"""Typed POD modes and their trajectory-based state-space representations."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Callable
import numpy as np

@dataclass(frozen=True)
class PODSystemReference:
    """Compact provenance pointer for modes derived from one trajectory."""
    source_archive: str | None
    system_name: str | None
    physical_parameters: dict[str, Any]

    @classmethod
    def from_source(cls, metadata: dict[str, Any] | None, source_archive: str | None=None) -> "PODSystemReference":
        metadata=metadata or {}
        return cls(source_archive, metadata.get("system"), dict(metadata.get("system_parameters") or {}))

def _trajectory_interpolant(states: np.ndarray, values: np.ndarray) -> Callable[[np.ndarray], float]:
    """Local inverse-distance interpolation of a scalar POD amplitude."""
    samples=np.asarray(states,dtype=float)
    amplitudes=np.asarray(values,dtype=float).reshape(-1)
    if samples.ndim==1: samples=samples[:,None]
    if samples.shape[0]!=amplitudes.size: raise ValueError("POD interpolation states and amplitudes must have equal sample counts.")
    def evaluate(state: np.ndarray) -> float:
        point=np.asarray(state,dtype=float).reshape(-1)
        if point.size!=samples.shape[1]: raise ValueError("State dimension does not match the POD interpolation data.")
        distance=np.linalg.norm(samples-point,axis=1); nearest=int(np.argmin(distance))
        if distance[nearest]<=1e-13: return float(amplitudes[nearest])
        count=min(8,amplitudes.size); indices=np.argpartition(distance,count-1)[:count]
        weights=1./np.maximum(distance[indices],1e-13)**2
        return float(np.sum(weights*amplitudes[indices])/np.sum(weights))
    return evaluate

@dataclass
class PODMode:
    """One SVD/POD mode with a state-space-interpolated modal amplitude.

    ``spatial_vector`` is the corresponding right singular vector.  ``values``
    are its signed temporal amplitudes ``sigma_i U[:, i]`` at source states;
    their interpolant gives the mode amplitude at an arbitrary state.  The
    reconstructed feature contribution is ``evaluate(x) * spatial_vector``.
    """
    index: int
    singular_value: float
    spatial_vector: np.ndarray
    values: np.ndarray
    system: PODSystemReference
    state_callable: Callable[[np.ndarray], float] | None=None
    metadata: dict[str, Any]=field(default_factory=dict)

    @property
    def has_state_callable(self) -> bool: return self.state_callable is not None

    def phi(self, x: np.ndarray) -> float:
        """Evaluate the state-interpolated POD amplitude at state ``x``."""
        if self.state_callable is None: raise NotImplementedError("This POD mode has no state-space interpolation callable.")
        return float(self.state_callable(np.asarray(x,dtype=float)))

    evaluate=phi

    def reconstruction(self, x: np.ndarray) -> np.ndarray:
        """Return this mode's reconstructed feature-space contribution at ``x``."""
        return self.phi(x)*np.asarray(self.spatial_vector,dtype=float)

    def summary(self) -> dict[str, Any]:
        return {"index":self.index,"singular_value":self.singular_value,"feature_dimension":int(np.asarray(self.spatial_vector).size),
                "value_samples":int(np.asarray(self.values).size),"has_phi":self.has_state_callable,"metadata":self.metadata}

@dataclass
class PROM:
    """Proper Orthogonal Decomposition Reduced-Order Model from one SVD."""
    system: PODSystemReference
    modes: list[PODMode]
    mean: np.ndarray
    metadata: dict[str, Any]=field(default_factory=dict)

    def __post_init__(self):
        if any(mode.system is not self.system for mode in self.modes):
            raise ValueError("Every POD mode must share the collection's system reference instance.")

    def summary(self) -> dict[str, Any]:
        return {"system":{"source_archive":self.system.source_archive,"system_name":self.system.system_name,
                           "physical_parameters":self.system.physical_parameters},"mode_count":len(self.modes),
                "modes":[mode.summary() for mode in self.modes],"metadata":self.metadata}

def pod_modes_from_svd(left: np.ndarray, singular: np.ndarray, right: np.ndarray, states: np.ndarray, *,
                       mean: np.ndarray, source_metadata: dict[str, Any] | None=None,
                       source_archive: str | None=None, metadata: dict[str, Any] | None=None) -> PROM:
    """Build a nonredundant POD-mode collection from direct SVD factors."""
    left=np.asarray(left,dtype=float); singular=np.asarray(singular,dtype=float).reshape(-1); right=np.asarray(right,dtype=float)
    samples=np.asarray(states,dtype=float)
    if samples.ndim==1: samples=samples[:,None]
    rank=singular.size
    if left.shape[0]!=samples.shape[0] or left.shape[1]<rank or right.shape[0]<rank:
        raise ValueError("SVD factors and source state samples have incompatible dimensions.")
    system=PODSystemReference.from_source(source_metadata,source_archive)
    shared=dict(metadata or {})
    modes=[]
    for index in range(rank):
        values=left[:,index]*singular[index]
        modes.append(PODMode(index,float(singular[index]),right[index].copy(),values,system,
                             _trajectory_interpolant(samples,values),
                             {**shared,"index":index,"representation":"trajectory-amplitude interpolation"}))
    return PROM(system,modes,np.asarray(mean,dtype=float).copy(),
                             {**shared,"representation":"direct SVD right-singular-vector spatial modes"})

def pod_modes_from_hankel_svd(left: np.ndarray, singular: np.ndarray, right: np.ndarray, states: np.ndarray, *,
                              mean: np.ndarray, source_metadata: dict[str, Any] | None=None,
                              source_archive: str | None=None, metadata: dict[str, Any] | None=None) -> PROM:
    """Build POD modes when a Hankel SVD stores temporal coordinates in ``V``."""
    left=np.asarray(left,dtype=float); singular=np.asarray(singular,dtype=float).reshape(-1); right=np.asarray(right,dtype=float)
    samples=np.asarray(states,dtype=float)
    if samples.ndim==1: samples=samples[:,None]
    rank=singular.size
    if left.shape[1]<rank or right.shape[0]<rank or right.shape[1]!=samples.shape[0]:
        raise ValueError("Hankel SVD factors and source state samples have incompatible dimensions.")
    system=PODSystemReference.from_source(source_metadata,source_archive); shared=dict(metadata or {}); modes=[]
    for index in range(rank):
        values=singular[index]*right[index,:]
        modes.append(PODMode(index,float(singular[index]),left[:,index].copy(),values,system,
                             _trajectory_interpolant(samples,values),
                             {**shared,"index":index,"representation":"Hankel temporal-coordinate interpolation"}))
    return PROM(system,modes,np.asarray(mean,dtype=float).copy(),
                {**shared,"representation":"Hankel SVD delay modes with temporal-coordinate interpolation"})
