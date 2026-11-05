"""Adaptive-mesh/repartition decisions with durable, solver-readable records."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import numpy as np
from .models import ParallelMeshPlan
from .partition import StructuredPartitioner

@dataclass(frozen=True)
class RepartitionRecord:
    created_utc: str
    reason: str
    prior_weights: tuple[float,...]
    imbalance: float
    checkpoint_required: bool
    def payload(self) -> dict: return asdict(self)

class AdaptiveMeshCoordinator:
    def __init__(self, plan: ParallelMeshPlan): self.plan=plan
    def needs_rebalance(self, weights) -> bool:
        values=np.asarray(weights,dtype=float)
        if values.size!=self.plan.worker_count or np.any(values<=0): raise ValueError("One positive load weight per mesh worker is required.")
        return float(values.max()/values.mean())>self.plan.adaptive.rebalance_imbalance
    def rebalance(self, weights, *, reason: str="adaptive load imbalance") -> tuple[ParallelMeshPlan,RepartitionRecord]:
        values=np.asarray(weights,dtype=float)
        record=RepartitionRecord(datetime.now(timezone.utc).isoformat(),reason,tuple(float(item) for item in values),float(values.max()/values.mean()),self.plan.adaptive.checkpoint_before_repartition)
        # Topology remains stable in this first backend.  A future AMR mesh
        # adapter supplies refined cells and invokes this same record path.
        rebuilt=StructuredPartitioner.build(self.plan)
        return rebuilt,record
