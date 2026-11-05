"""Ghost-region exchange for local test backends and optional MPI deployments."""
from __future__ import annotations
from typing import Mapping
import numpy as np
from .models import ParallelMeshPlan

class GhostExchange:
    """Exchange nearest-neighbor halos from rank-indexed interior arrays."""
    def __init__(self, plan: ParallelMeshPlan):
        self.plan=plan
        if not plan.partitions: raise ValueError("Partition the plan before constructing GhostExchange.")
    def padded(self, fields: Mapping[int,np.ndarray]) -> dict[int,np.ndarray]:
        width=self.plan.ghost_width; output={}
        for partition in self.plan.partitions:
            field=np.asarray(fields[partition.rank]);
            if field.shape[:len(partition.shape)]!=partition.shape: raise ValueError(f"Rank {partition.rank} field shape does not match its interior partition.")
            output[partition.rank]=np.pad(field,[(width,width)]*len(partition.shape)+[(0,0)]*(field.ndim-len(partition.shape)),mode="edge")
        for partition in self.plan.partitions:
            target=output[partition.rank]
            for axis in range(len(partition.shape)):
                for sign,name in ((-1,f"axis{axis}-"),(1,f"axis{axis}+")):
                    neighbor=partition.neighbors.get(name)
                    if neighbor is None: continue
                    source=np.asarray(fields[neighbor]); slices=[slice(None)]*source.ndim; slices[axis]=slice(-width,None) if sign<0 else slice(0,width)
                    destination=[slice(None)]*target.ndim; destination[axis]=slice(0,width) if sign<0 else slice(-width,None)
                    for other in range(len(partition.shape)):
                        if other!=axis: destination[other]=slice(width,-width)
                    target[tuple(destination)]=source[tuple(slices)]
        return output

class MPIExchange(GhostExchange):
    """Optional mpi4py halo exchange; imported only for an MPI PDE run."""
    def __init__(self, plan: ParallelMeshPlan, communicator=None):
        super().__init__(plan)
        if communicator is None:
            try:
                from mpi4py import MPI
            except ImportError as error: raise RuntimeError("MPIExchange requires optional dependency mpi4py.") from error
            communicator=MPI.COMM_WORLD
        self.communicator=communicator
        if communicator.Get_size()!=plan.worker_count: raise ValueError("MPI communicator size must match ParallelMeshPlan.worker_count.")
    def exchange_rank_field(self, field: np.ndarray) -> np.ndarray:
        rank=self.communicator.Get_rank(); width=self.plan.ghost_width; partition=self.plan.partitions[rank]
        padded=np.pad(np.asarray(field),[(width,width)]*len(partition.shape)+[(0,0)]*(field.ndim-len(partition.shape)),mode="edge")
        requests=[]; received=[]
        for axis in range(len(partition.shape)):
            for sign,name in ((-1,f"axis{axis}-"),(1,f"axis{axis}+")):
                neighbor=partition.neighbors.get(name)
                if neighbor is None: continue
                send=np.ascontiguousarray(np.take(field,range(-width,0) if sign<0 else range(width),axis=axis)); receive=np.empty_like(send)
                destination=[slice(None)]*padded.ndim; destination[axis]=slice(0,width) if sign<0 else slice(-width,None)
                for other in range(len(partition.shape)):
                    if other!=axis: destination[other]=slice(width,-width)
                requests.extend((self.communicator.Irecv(receive,source=neighbor,tag=axis*2+(1 if sign<0 else 0)),self.communicator.Isend(send,dest=neighbor,tag=axis*2+(0 if sign<0 else 1))))
                received.append((tuple(destination),receive))
        for request in requests: request.Wait()
        for destination,receive in received: padded[destination]=receive
        return padded
