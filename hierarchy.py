"""Shared tiered experiment objects: Batch > System > Cluster > Trajectory."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import Any
import numpy as np
from time_series import TimeSeries

def _same_time(reference: np.ndarray, candidate: np.ndarray) -> bool:
    return reference.shape==candidate.shape and bool(np.allclose(reference,candidate,rtol=1e-12,atol=1e-15))

def _serializable(value: Any) -> Any:
    if is_dataclass(value): return _serializable(asdict(value))
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,dict): return {str(key):_serializable(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)): return [_serializable(item) for item in value]
    if isinstance(value,(str,int,float,bool)) or value is None: return value
    return repr(value)

def _pool(children: list[TimeSeries], tier: str) -> TimeSeries:
    if not children:
        raise ValueError(f"A {tier} needs at least one completed child TimeSeries to pool.")
    reference=children[0].time
    if any(not _same_time(reference,item.time) for item in children[1:]):
        raise ValueError(f"Cannot pool {tier}: every child must use identical time values.")
    shapes={item.values.shape[1:] for item in children}
    if len(shapes)!=1:
        raise ValueError(f"Cannot pool {tier}: child state/field shapes differ.")
    return TimeSeries(reference,np.stack([item.values for item in children],axis=1),pooled_members=len(children))

@dataclass
class Trajectory:
    """One initial condition and its eventual solver-produced TimeSeries."""
    initial_condition: np.ndarray
    label: str = "trajectory"
    time_series: TimeSeries | None = None
    archive_path: str | None = None
    metadata: dict[str,Any] = field(default_factory=dict)

    def __post_init__(self):
        self.initial_condition=np.asarray(self.initial_condition,dtype=float).reshape(-1)
        if not self.initial_condition.size or not np.all(np.isfinite(self.initial_condition)):
            raise ValueError("Trajectory initial_condition must be a nonempty finite vector.")

@dataclass
class Cluster:
    """A collection of trajectories sharing one System's dynamics and time mesh."""
    trajectories: list[Trajectory] = field(default_factory=list)
    label: str = "cluster"
    time_series: TimeSeries | None = None
    pooled_initial_conditions: np.ndarray | None = None

    def pool_initial_conditions(self) -> np.ndarray:
        """Stack child initial conditions as ``(trajectory, state)``."""
        if not self.trajectories: raise ValueError("A cluster needs trajectories to pool initial conditions.")
        dimensions={item.initial_condition.shape for item in self.trajectories}
        if len(dimensions)!=1: raise ValueError("Cannot pool cluster initial conditions with different state dimensions.")
        self.pooled_initial_conditions=np.stack([item.initial_condition for item in self.trajectories])
        return self.pooled_initial_conditions

    def pool(self) -> TimeSeries:
        self.pool_initial_conditions()
        missing=[item.label for item in self.trajectories if item.time_series is None]
        if missing: raise ValueError(f"Cannot pool cluster with unsolved trajectories: {missing}.")
        self.time_series=_pool([item.time_series for item in self.trajectories if item.time_series is not None],"cluster")
        return self.time_series

@dataclass
class System:
    """A dynsyslib system plus its shared prospective integration time sequence."""
    dynamical_system: Any
    clusters: list[Cluster] = field(default_factory=list)
    time_span: tuple[float,float] | None = None
    times: np.ndarray | None = None
    uniform_time_steps: bool = False
    delta_t: float | None = None
    kind: str | None = None
    label: str = "system"
    archive_root: str | None = None
    time_series: TimeSeries | None = None
    pooled_initial_conditions: np.ndarray | None = None
    metadata: dict[str,Any] = field(default_factory=dict)

    def set_initial_conditions(self, initial_conditions: np.ndarray, *, cluster_label: str="explicit-initial-conditions") -> Cluster:
        """Create this System's sole trajectory cluster from explicit state rows."""
        if self.clusters: raise ValueError("System already has clusters; do not mix stored and explicit initial conditions.")
        values=np.asarray(initial_conditions,dtype=float)
        if values.ndim==1: values=values[None,:]
        dimension=int(np.asarray(self.dynamical_system.initial_state).size)
        if values.ndim!=2 or values.shape[0]<1 or values.shape[1]!=dimension or not np.all(np.isfinite(values)):
            raise ValueError(f"initial_conditions must have shape (trajectory_count, {dimension}).")
        cluster=Cluster([Trajectory(row,label=f"trajectory-{index:03d}") for index,row in enumerate(values)],label=cluster_label)
        self.clusters=[cluster]; self.pool_initial_conditions()
        return cluster

    def pool_initial_conditions(self) -> np.ndarray:
        """Concatenate all child cluster initial-condition arrays by trajectory."""
        if not self.clusters: raise ValueError("A system needs clusters to pool initial conditions.")
        arrays=[item.pool_initial_conditions() for item in self.clusters]
        dimensions={item.shape[1] for item in arrays}
        if len(dimensions)!=1: raise ValueError("Cannot pool system initial conditions with different state dimensions.")
        self.pooled_initial_conditions=np.concatenate(arrays,axis=0)
        return self.pooled_initial_conditions

    def integration_times(self) -> np.ndarray:
        if self.times is not None:
            times=np.asarray(self.times,dtype=float).reshape(-1)
            if times.size<2 or np.any(~np.isfinite(times)) or np.any(np.diff(times)<=0):
                raise ValueError("System.times must contain at least two finite, strictly increasing values.")
            return times
        if not self.uniform_time_steps:
            raise ValueError("Specify System.times, or set uniform_time_steps=True with time_span and delta_t.")
        if self.time_span is None or self.delta_t is None or self.delta_t<=0:
            raise ValueError("A uniform System needs time_span and positive delta_t.")
        start,stop=map(float,self.time_span); steps=(stop-start)/self.delta_t
        if start>=stop or not np.isclose(steps,round(steps),rtol=1e-11,atol=1e-11):
            raise ValueError("time_span must be an increasing integer multiple of delta_t.")
        return np.linspace(start,stop,int(round(steps))+1)

    @property
    def solver_kind(self) -> str:
        return self.kind or getattr(self.dynamical_system,"kind","ode")

    def pool(self) -> TimeSeries:
        self.pool_initial_conditions()
        missing=[item.label for item in self.clusters if item.time_series is None]
        if missing: raise ValueError(f"Cannot pool system with unpooled clusters: {missing}.")
        children=[item.time_series for item in self.clusters if item.time_series is not None]
        if not children: raise ValueError("A system needs at least one pooled cluster to pool.")
        reference=children[0].time
        if any(not _same_time(reference,item.time) for item in children[1:]):
            raise ValueError("Cannot pool system: every cluster must use identical time values.")
        shapes={item.values.shape[2:] for item in children}
        if len(shapes)!=1: raise ValueError("Cannot pool system: cluster state/field shapes differ.")
        # Flatten the cluster/trajectory hierarchy only at the data level so a
        # pooled System is uniformly (time, total_trajectory, state_or_field).
        values=np.concatenate([item.values for item in children],axis=1)
        self.time_series=TimeSeries(reference,values,pooled_members=int(values.shape[1]))
        return self.time_series

    @classmethod
    def rpe_uniform_cluster(cls, dynamical_system: Any, *, trajectory_count: int=30,
                            radius_fraction_range: tuple[float,float]=(.5,2.), time_span: tuple[float,float],
                            delta_t: float, seed: int | None=None, label: str="rpe-system") -> "System":
        """Build the requested RPE ensemble with radii uniformly sampled about equilibrium."""
        if trajectory_count<1 or radius_fraction_range[0]<=0 or radius_fraction_range[1]<radius_fraction_range[0]:
            raise ValueError("Provide a positive trajectory_count and an ordered positive radius fraction range.")
        equilibrium=float(dynamical_system.equilibrium_radius_m); generator=np.random.default_rng(seed)
        radii=generator.uniform(*radius_fraction_range,size=trajectory_count)*equilibrium
        trajectories=[Trajectory(np.array((radius,0.)),label=f"trajectory-{index:03d}") for index,radius in enumerate(radii)]
        return cls(dynamical_system,[Cluster(trajectories,label="initial-condition-ensemble")],time_span,None,True,delta_t,label=label)

@dataclass
class Batch:
    """Top-level container for independently configured Systems."""
    systems: list[System] = field(default_factory=list)
    label: str = "batch"
    pooled_initial_conditions: list[np.ndarray] = field(default_factory=list)
    metadata: dict[str,Any] = field(default_factory=dict)

    def pool_initial_conditions(self) -> list[np.ndarray]:
        """Pool each System independently; systems may have different state dimensions."""
        self.pooled_initial_conditions=[system.pool_initial_conditions() for system in self.systems]
        return self.pooled_initial_conditions

    def pool(self) -> list[TimeSeries]:
        """Pool each system independently; systems need not share a time mesh."""
        self.pool_initial_conditions()
        return [system.pool() for system in self.systems]

def hierarchy_payload(target: System | Batch) -> tuple[dict[str,Any],dict[str,np.ndarray]]:
    """Return a compact hierarchy manifest, not a duplicate numerical archive.

    Individual trajectory archives are the sole canonical time-series storage.
    The manifest records their relationship and can reconstruct pooled views on
    demand, so a System/Cluster does not repeat its children's solution arrays.
    """
    systems=target.systems if isinstance(target,Batch) else [target]
    metadata={"format":"dynsyslab.hierarchy/v1","tier":"batch" if isinstance(target,Batch) else "system","label":target.label,"metadata":_serializable(getattr(target,"metadata",{})),"systems":[]}
    arrays={}
    for system_index,system in enumerate(systems):
        entry={"label":system.label,"kind":system.solver_kind,"dynamical_system":_serializable(system.dynamical_system),"metadata":_serializable(system.metadata),"clusters":[]}
        for cluster_index,cluster in enumerate(system.clusters):
            cluster_entry={"label":cluster.label,"pooled_initial_conditions":None if cluster.pooled_initial_conditions is None else cluster.pooled_initial_conditions.tolist(),"trajectories":[]}
            for trajectory in cluster.trajectories:
                cluster_entry["trajectories"].append({"label":trajectory.label,"initial_condition":trajectory.initial_condition.tolist(),"solution_archive":trajectory.archive_path,"metadata":trajectory.metadata})
            entry["clusters"].append(cluster_entry)
        entry["pooled_initial_conditions"]=None if system.pooled_initial_conditions is None else system.pooled_initial_conditions.tolist()
        metadata["systems"].append(entry)
    return metadata,arrays
