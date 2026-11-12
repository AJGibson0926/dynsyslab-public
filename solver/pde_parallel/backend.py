"""PDE-only parallel backend facade, ready to be attached to solve_pde later."""
from __future__ import annotations
from abc import ABC
from dataclasses import dataclass
from typing import Mapping
import numpy as np
from .exchange import GhostExchange, MPIExchange
from .models import ParallelMeshPlan
from .partition import StructuredPartitioner

class PDEParallelBackend(ABC):
    """Owns intra-solve mesh communication, never JobManager scheduling."""
    def __init__(self, plan: ParallelMeshPlan):
        self.plan=StructuredPartitioner.build(plan) if not plan.partitions else plan
    def metadata(self) -> dict: return {"parallel_mesh_plan":self.plan.payload(),"backend":type(self).__name__}

@dataclass
class LocalPDEParallelBackend(PDEParallelBackend):
    plan: ParallelMeshPlan
    def __post_init__(self): PDEParallelBackend.__init__(self,self.plan); self.exchange=GhostExchange(self.plan)
    def exchange_ghosts(self, fields: Mapping[int,np.ndarray]) -> dict[int,np.ndarray]: return self.exchange.padded(fields)

@dataclass
class MPIPDEParallelBackend(PDEParallelBackend):
    plan: ParallelMeshPlan
    communicator: object | None=None
    def __post_init__(self): PDEParallelBackend.__init__(self,self.plan); self.exchange=MPIExchange(self.plan,self.communicator)
    def exchange_rank_field(self, field: np.ndarray) -> np.ndarray: return self.exchange.exchange_rank_field(field)
