"""Serializable mesh-parallel plans.  These describe work; they do not solve it."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from typing import Literal

@dataclass(frozen=True)
class StructuredMesh:
    shape: tuple[int,...]
    spacing: tuple[float,...] | None=None
    origin: tuple[float,...] | None=None
    boundary: tuple[str,...] | None=None
    def __post_init__(self):
        if not self.shape or any(int(value)<2 for value in self.shape): raise ValueError("A structured mesh needs at least two cells on each axis.")
        if self.spacing is not None and len(self.spacing)!=len(self.shape): raise ValueError("spacing must match mesh dimension.")
        if self.boundary is not None and len(self.boundary)!=len(self.shape): raise ValueError("boundary must match mesh dimension.")
    @property
    def dimension(self) -> int: return len(self.shape)
    def payload(self) -> dict: return asdict(self)

@dataclass(frozen=True)
class MeshPartition:
    rank: int
    interior: tuple[tuple[int,int],...]
    neighbors: dict[str,int]=field(default_factory=dict)
    weight: float=1.0
    def __post_init__(self):
        if self.rank<0 or any(end<=start for start,end in self.interior): raise ValueError("A partition needs positive half-open interior ranges.")
    @property
    def shape(self) -> tuple[int,...]: return tuple(end-start for start,end in self.interior)
    def payload(self) -> dict: return asdict(self)

@dataclass(frozen=True)
class AdaptiveMeshPolicy:
    enabled: bool=False
    refinement_threshold: float=0.8
    coarsening_threshold: float=0.1
    rebalance_imbalance: float=1.25
    checkpoint_before_repartition: bool=True

@dataclass(frozen=True)
class ParallelMeshPlan:
    """Complete PDE-internal parallelization contract, stored with a PDE run."""
    mesh: StructuredMesh
    partition_shape: tuple[int,...]
    ghost_width: int=1
    communication: Literal["local","mpi"]="local"
    adaptive: AdaptiveMeshPolicy=field(default_factory=AdaptiveMeshPolicy)
    partitions: tuple[MeshPartition,...]=()
    load_balance: Literal["static","adaptive"]="static"
    local_workers: int | None=None
    def __post_init__(self):
        if len(self.partition_shape)!=self.mesh.dimension or any(value<1 for value in self.partition_shape): raise ValueError("partition_shape must match mesh dimension and be positive.")
        if self.ghost_width<1: raise ValueError("ghost_width must be at least one.")
        if self.local_workers is not None and self.local_workers<1: raise ValueError("local_workers must be positive when supplied.")
        if self.partitions and len(self.partitions)!=self.worker_count: raise ValueError("partition count must equal the cartesian worker count.")
    @property
    def worker_count(self) -> int:
        product=1
        for value in self.partition_shape: product*=value
        return product
    def payload(self) -> dict: return asdict(self)
