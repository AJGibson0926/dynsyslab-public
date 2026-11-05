"""Deterministic structured-domain partitioning and neighbor discovery."""
from __future__ import annotations
from itertools import product
from .models import MeshPartition, ParallelMeshPlan

class StructuredPartitioner:
    @staticmethod
    def build(plan: ParallelMeshPlan) -> ParallelMeshPlan:
        axes=[StructuredPartitioner._cuts(cells,parts) for cells,parts in zip(plan.mesh.shape,plan.partition_shape)]
        coordinates=list(product(*(range(count) for count in plan.partition_shape)))
        rank_of={coordinate:index for index,coordinate in enumerate(coordinates)}; partitions=[]
        for rank,coordinate in enumerate(coordinates):
            interior=tuple(axis[index] for axis,index in zip(axes,coordinate)); neighbors={}
            for axis,index in enumerate(coordinate):
                for offset,name in ((-1,f"axis{axis}-"),(1,f"axis{axis}+")):
                    nearby=list(coordinate); nearby[axis]+=offset
                    if 0<=nearby[axis]<plan.partition_shape[axis]: neighbors[name]=rank_of[tuple(nearby)]
            partitions.append(MeshPartition(rank,interior,neighbors))
        return ParallelMeshPlan(plan.mesh,plan.partition_shape,plan.ghost_width,plan.communication,plan.adaptive,tuple(partitions),plan.load_balance,plan.local_workers)
    @staticmethod
    def _cuts(cells: int, parts: int) -> tuple[tuple[int,int],...]:
        base,remainder=divmod(cells,parts); start=0; output=[]
        for index in range(parts):
            width=base+(1 if index<remainder else 0); end=start+width
            if width<1: raise ValueError("More partitions than cells along an axis.")
            output.append((start,end)); start=end
        return tuple(output)
