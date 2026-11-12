"""PDE pathway via method-of-lines systems exposing rhs(t, field, mesh)."""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from time_series import TimeSeries
from scipy.integrate import solve_ivp
from typing import Any
from .models import NumericalResult, SolveRequest
from .backends import select_reference_backend


def _fixed_step_pde(rhs, t0: float, tf: float, field0: np.ndarray, step_size: float,
                    method: str, request: SolveRequest, session: Any):
    """Reference explicit integration that preserves arbitrary field shape."""
    times=[float(t0)]; fields=[np.asarray(field0,dtype=float).copy()]
    current=float(t0); field=np.asarray(field0,dtype=float).copy()
    while current < tf:
        step=min(float(step_size),tf-current)
        if method=="Euler":
            field=field+step*rhs(current,field.ravel()).reshape(field.shape)
        else:
            k1=rhs(current,field.ravel()).reshape(field.shape)
            k2=rhs(current+step/2,(field+step*k1/2).ravel()).reshape(field.shape)
            k3=rhs(current+step/2,(field+step*k2/2).ravel()).reshape(field.shape)
            k4=rhs(current+step,(field+step*k3).ravel()).reshape(field.shape)
            field=field+step*(k1+2*k2+2*k3+k4)/6
        current+=step; times.append(current); fields.append(field.copy())
        if request.history.track_during_solve:
            session.checkpoint(np.asarray((current,)),field[None,...])
    return np.asarray(times),np.asarray(fields)

def _partitioned_pde(request: SolveRequest, session: Any, field0: np.ndarray) -> NumericalResult:
    """Synchronous PDE subdomain pathway for systems with ``rhs_subdomain``.

    It is deliberately separate from method-of-lines ``solve_ivp``: each
    explicit stage exchanges ghost regions before independent subdomain RHS
    evaluations.  MPI uses the same plan objects but must be launched through
    an MPI-aware driver; this local pathway never creates distributed ranks.
    """
    from .pde_parallel import LocalPDEParallelBackend, ParallelMeshPlan
    numerical_backend=select_reference_backend(request.options.execution_backend,"pde")
    plan=request.parallel_mesh_plan
    if not isinstance(plan,ParallelMeshPlan): raise TypeError("parallel_mesh_plan must be a ParallelMeshPlan.")
    backend_name=request.pde_parallel_backend or plan.communication
    if backend_name!="local": raise ValueError("solve_pde currently attaches the local mesh backend only. Launch MPI ranks externally with MPIPDEParallelBackend.")
    if tuple(plan.mesh.shape)!=tuple(field0.shape[:plan.mesh.dimension]): raise ValueError("ParallelMeshPlan mesh shape must match the leading spatial shape of initial_field(mesh).")
    if not hasattr(request.system,"rhs_subdomain"): raise TypeError("A parallel PDE system must implement rhs_subdomain(t, ghosted_field, partition, mesh, parallel_plan).")
    if request.options.step_size is None or request.options.step_size<=0: raise ValueError("Partitioned PDE integration requires SolverOptions(step_size=...) for synchronized ghost exchanges.")
    if request.options.method not in {"Euler","RK4"}: raise ValueError("Partitioned PDE integration currently supports explicit Euler or RK4 only.")
    backend=LocalPDEParallelBackend(plan); plan=backend.plan; dt=float(request.options.step_size); t0,tf=request.time_span
    fields={partition.rank:np.array(field0[tuple(slice(start,end) for start,end in partition.interior)],copy=True) for partition in plan.partitions}
    def evaluate(time: float, values: dict[int,np.ndarray]) -> dict[int,np.ndarray]:
        ghosts=backend.exchange_ghosts(values)
        def rhs_for(partition):
            derivative=np.asarray(request.system.rhs_subdomain(time,ghosts[partition.rank],partition,request.mesh,plan),dtype=float)
            if derivative.shape!=values[partition.rank].shape: raise ValueError(f"rhs_subdomain returned {derivative.shape} for rank {partition.rank}; expected {values[partition.rank].shape}.")
            return partition.rank,derivative
        workers=plan.local_workers or 1
        if workers==1: return dict(rhs_for(partition) for partition in plan.partitions)
        with ThreadPoolExecutor(max_workers=min(workers,plan.worker_count)) as pool: return dict(pool.map(rhs_for,plan.partitions))
    def assemble(values: dict[int,np.ndarray]) -> np.ndarray:
        output=np.empty_like(field0)
        for partition in plan.partitions: output[tuple(slice(start,end) for start,end in partition.interior)]=values[partition.rank]
        return output
    times=[float(t0)]; snapshots=[assemble(fields)]; current=float(t0)
    try:
        while current < tf:
            step=min(dt,tf-current)
            if request.options.method=="Euler":
                derivative=evaluate(current,fields); fields={rank:value+step*derivative[rank] for rank,value in fields.items()}
            else:
                k1=evaluate(current,fields); stage={rank:value+0.5*step*k1[rank] for rank,value in fields.items()}
                k2=evaluate(current+0.5*step,stage); stage={rank:value+0.5*step*k2[rank] for rank,value in fields.items()}
                k3=evaluate(current+0.5*step,stage); stage={rank:value+step*k3[rank] for rank,value in fields.items()}
                k4=evaluate(current+step,stage); fields={rank:value+(step/6.)*(k1[rank]+2*k2[rank]+2*k3[rank]+k4[rank]) for rank,value in fields.items()}
            current+=step; times.append(current); frame=assemble(fields); snapshots.append(frame)
            if request.history.track_during_solve: session.checkpoint(np.asarray((current,)),frame[None,...])
        metadata={"mesh_shape":field0.shape,"parallel_mesh":backend.metadata(),"parallel_method":request.options.method,"parallel_step_size":dt,"execution_backend":numerical_backend.provenance(requested=request.options.execution_backend)}
        return NumericalResult(TimeSeries(np.asarray(times),np.asarray(snapshots)),"complete","parallel mesh solve completed",metadata=metadata)
    except Exception as exc:
        return NumericalResult(None,"failed",str(exc),metadata={"mesh_shape":field0.shape,"parallel_mesh":backend.metadata(),"execution_backend":numerical_backend.provenance(requested=request.options.execution_backend)})

def solve_pde(request: SolveRequest, session: Any) -> NumericalResult:
    if request.time_span is None or request.mesh is None: raise ValueError("PDE solve requires time_span and mesh.")
    system,mesh=request.system,request.mesh; field0=np.asarray(system.initial_field(mesh),dtype=float); shape=field0.shape
    if request.parallel_mesh_plan is not None: return _partitioned_pde(request,session,field0)
    numerical_backend=select_reference_backend(request.options.execution_backend,"pde")
    def rhs(t,flat): return np.asarray(system.rhs(t,flat.reshape(shape),mesh),dtype=float).ravel()
    try:
        if request.options.method in {"Euler","RK4"}:
            if request.options.step_size is None or request.options.step_size<=0:
                raise ValueError(f"{request.options.method} PDE integration requires positive options.step_size.")
            time,fields=_fixed_step_pde(rhs,*request.time_span,field0,request.options.step_size,request.options.method,request,session)
        elif request.history.track_during_solve:
            t0,tf=request.time_span; interval=request.history.checkpoint_interval or (tf-t0)/20
            current_t,current= t0,field0.ravel(); time_parts=[np.asarray((t0,))]; field_parts=[field0[None,...]]
            while current_t < tf:
                end=min(tf,current_t+interval); sol=solve_ivp(rhs,(current_t,end),current,method=request.options.method,rtol=request.options.rtol,atol=request.options.atol,max_step=request.options.max_step or np.inf)
                if not sol.success: raise RuntimeError(sol.message)
                block=sol.y.T.reshape((-1,)+shape); session.checkpoint(sol.t,block)
                # Each solve_ivp segment repeats its initial boundary. Retain
                # it once so archived TimeSeries time remains strictly increasing.
                time_parts.append(sol.t[1:]); field_parts.append(block[1:]); current_t,current=end,sol.y[:,-1]
            time=np.concatenate(time_parts); fields=np.concatenate(field_parts)
        else:
            sol=solve_ivp(rhs,request.time_span,field0.ravel(),method=request.options.method,rtol=request.options.rtol,atol=request.options.atol,max_step=request.options.max_step or np.inf)
            if not sol.success: raise RuntimeError(sol.message)
            time,fields=sol.t,sol.y.T.reshape((-1,)+shape)
        return NumericalResult(TimeSeries(time,fields),"complete","solve completed",metadata={"mesh_shape":shape,"execution_backend":numerical_backend.provenance(requested=request.options.execution_backend)})
    except Exception as exc: return NumericalResult(None,"failed",str(exc),metadata={"mesh_shape":shape,"execution_backend":numerical_backend.provenance(requested=request.options.execution_backend)})
