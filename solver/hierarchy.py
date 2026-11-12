"""Solver orchestration for shared Batch/System/Cluster/Trajectory objects."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from datetime import datetime, timezone
from time import perf_counter
import json
import numpy as np
from archiver.gateway import SOLVER_ARCHIVE
from hierarchy import Batch, System
from time_series import TimeSeries
from .backends import select_ode_backend
from .lifecycle import finalize_result
from .models import HistoryOptions, InterpolationOptions, NumericalResult, SolveRequest, SolverOptions
from .native_ode import integrate_fixed_native, native_fixed_eligibility

@dataclass
class HierarchyResult:
    """Aggregate outcome for one solved System or Batch hierarchy."""
    target: System | Batch
    trajectory_results: list[Any]
    archive_path: str
    status: str="complete"


def _native_cluster_results(requests: list[SolveRequest]) -> list[NumericalResult] | None:
    """Integrate one eligible cluster in one native call.

    ``None`` means the ordinary scalar route should be used.  Each successful
    native trajectory still receives its own archive session and passes through
    the same lifecycle as a scalar solve.
    """
    if len(requests) < 2 or any(request.kind != "ode" for request in requests):
        return None
    first = requests[0]
    method = first.options.method
    if method not in {"Euler", "RK4"} or not first.options.step_size:
        return None
    requested_form = first.options.use_dimensionless_form
    has_dimensionless = hasattr(first.system, "dimensionless_rhs") and hasattr(first.system, "dimensionless_initial_state")
    dimensionless = (has_dimensionless and bool(getattr(first.system,"prefer_dimensionless_form",True))) if requested_form is None else requested_form
    eligible, reason = native_fixed_eligibility(first, dimensionless=dimensionless)
    backend = select_ode_backend(
        first.options.execution_backend,
        method,
        native_eligible=eligible,
        ineligibility_reason=reason,
    )
    if not backend.compiled or not backend.batched_ode:
        return None

    t0, tf = first.time_span
    initial_states = np.stack([np.asarray(request.initial_state, dtype=float) for request in requests])
    sessions = [SOLVER_ARCHIVE.new_session(request) for request in requests]
    for request, session in zip(requests, sessions):
        session.telemetry("started", time_span=request.time_span, method=method)
        if request.on_event is not None:
            request.on_event({"event": "started", "kind": "ode", "time_span": request.time_span})
    try:
        time, values = integrate_fixed_native(
            first.system,
            float(t0),
            float(tf),
            initial_states,
            float(first.options.step_size),
            method,
        )
        if values.shape != (len(requests), time.size, initial_states.shape[1]):
            raise RuntimeError("The native batched ODE kernel returned an invalid result shape.")
        output = []
        for index, (request, session) in enumerate(zip(requests, sessions)):
            result = NumericalResult(
                TimeSeries(time, values[index]),
                "complete",
                "solve completed",
                metadata={
                    "state_dimension": int(initial_states.shape[1]),
                    "integration_form": "dimensional",
                    "control_sample_time": request.options.control_sample_time,
                    "control_evaluation": request.options.control_evaluation,
                    "checkpoint_pause": False,
                    "checkpoint_path": None,
                    "execution_backend": backend.provenance(requested=request.options.execution_backend),
                    "batch_execution": {
                        "enabled": True,
                        "batch_size": len(requests),
                        "trajectory_index": index,
                    },
                },
            )
            if request.on_event is not None:
                request.on_event({"event": "completed", "samples": result.time_series.sample_count, "checkpoint_path": None})
            output.append(finalize_result(request, session, result))
        return output
    except Exception as exc:
        output = []
        for request, session in zip(requests, sessions):
            if request.on_event is not None:
                request.on_event({"event": "failed", "message": str(exc)})
            result = NumericalResult(
                None,
                "failed",
                str(exc),
                metadata={
                    "state_dimension": int(initial_states.shape[1]),
                    "execution_backend": backend.provenance(requested=request.options.execution_backend),
                    "batch_execution": {"enabled": True, "batch_size": len(requests)},
                },
            )
            output.append(finalize_result(request, session, result))
        return output

def solve_hierarchy(target: System | Batch, *, options: SolverOptions | None=None,
                    history: HistoryOptions | None=None, interpolation_method: str="pchip",
                    postprocessors: tuple=(), initial_conditions=None, controller=None, estimator=None,
                    controllers: list[Any] | tuple[Any,...] | None=None,
                    on_event=None, cancel_token=None, pause_token=None, job_id: str | None=None) -> HierarchyResult:
    """Solve every trajectory, pool completed children, and archive the hierarchy.

    Each System supplies its single target mesh.  It is passed as an explicit
    interpolation grid, so every trajectory can be safely pooled afterward.
    """
    systems=target.systems if isinstance(target,Batch) else [target]
    if controllers is not None and len(controllers)!=len(systems):
        raise ValueError("controllers must contain exactly one controller per System.")
    if initial_conditions is not None:
        if isinstance(target,Batch):
            raise ValueError("For a Batch, configure initial conditions on each System; explicit initial_conditions is for a single System.")
        target.set_initial_conditions(initial_conditions)
    results=[]
    from . import solve  # delayed to avoid the package's import cycle
    total_trajectories=sum(len(cluster.trajectories) for system in systems for cluster in system.clusters)
    completed_trajectories=0
    started=perf_counter()
    label=getattr(target,"label","hierarchy-solve"); safe="".join(item if item.isalnum() or item in "-_" else "-" for item in label)
    log_dir=SOLVER_ARCHIVE.archive_root/"solver_logs"; log_dir.mkdir(parents=True,exist_ok=True); log_path=log_dir/f"{safe}.jsonl"; emitted=set()
    aggression=(options or SolverOptions()).estimate_aggression
    interval={"low":25,"medium":10,"high":5}.get(aggression)
    if interval is None: raise ValueError("estimate_aggression must be 'low', 'medium', or 'high'.")
    def log_progress(*,force: bool=False):
        percent=int(100*completed_trajectories/total_trajectories) if total_trajectories else 100; threshold=(percent//10)*10
        threshold=(percent//interval)*interval
        if not ((threshold and threshold not in emitted) or (force and 100 not in emitted)): return
        emitted.add(threshold if threshold else 100); paths=[Path(item.archive_path) for item in results if getattr(item,"archive_path",None) and Path(item.archive_path).is_file()]
        elapsed=perf_counter()-started; rate=None if not completed_trajectories else elapsed/completed_trajectories
        with log_path.open("a",encoding="utf-8") as handle: handle.write(json.dumps({"event":"progress","time_utc":datetime.now(timezone.utc).isoformat(),"completed_trajectories":completed_trajectories,"planned_trajectories":total_trajectories,"percent":percent,"threshold_percent":threshold,"estimate_aggression":aggression,"elapsed_seconds":elapsed,"seconds_per_trajectory":rate,"estimated_remaining_seconds":None if rate is None else rate*(total_trajectories-completed_trajectories),"solution_file_count":len(paths),"solution_bytes":sum(path.stat().st_size for path in paths)})+"\n")
    selected_options=options or SolverOptions()
    selected_history=history or HistoryOptions()
    for system_index,system in enumerate(systems):
        system_controller=controllers[system_index] if controllers is not None else controller
        times=system.integration_times()
        for cluster_index,cluster in enumerate(system.clusters):
            requests=[]
            for trajectory_index,trajectory in enumerate(cluster.trajectories):
                request_kwargs={
                    "system":system.dynamical_system,"kind":system.solver_kind,
                    "time_span":(float(times[0]),float(times[-1])),
                    "initial_state":trajectory.initial_condition,"options":selected_options,"history":selected_history,
                    "interpolation":InterpolationOptions(times=times,method=interpolation_method),"postprocessors":postprocessors,
                    "controller":system_controller,"on_event":on_event,"cancel_token":cancel_token,"pause_token":pause_token,
                    "archive_root":system.archive_root,"job_id":job_id,
                    "metadata":dict(system.metadata)|dict(trajectory.metadata),
                    "label":f"{system.label}-c{cluster_index:03d}-t{trajectory_index:03d}",
                }
                # Omission intentionally preserves SolveRequest's identity-estimator
                # default; only an explicitly supplied estimator replaces it.
                if estimator is not None: request_kwargs["estimator"]=estimator
                requests.append(SolveRequest(**request_kwargs))
            cluster_results=_native_cluster_results(requests)
            if cluster_results is None:
                cluster_results=[solve(request) for request in requests]
            for trajectory_index,(trajectory,result) in enumerate(zip(cluster.trajectories,cluster_results)):
                results.append(result)
                if result.status=="paused":
                    trajectory.time_series=result.time_series; trajectory.archive_path=result.archive_path
                    return HierarchyResult(target,results,str(result.archive_path),status="paused")
                if result.status!="complete" or result.time_series is None:
                    paths=[Path(item.archive_path) for item in results if getattr(item,"archive_path",None) and Path(item.archive_path).is_file()]
                    with log_path.open("a",encoding="utf-8") as handle: handle.write(json.dumps({"event":"failed","time_utc":datetime.now(timezone.utc).isoformat(),"completed_trajectories":completed_trajectories,"planned_trajectories":total_trajectories,"failed_trajectory":trajectory.label,"reason":result.message,"solution_file_count":len(paths),"solution_bytes":sum(path.stat().st_size for path in paths)})+"\n")
                    raise RuntimeError(f"Trajectory {trajectory.label} failed: {result.message}")
                trajectory.time_series=result.time_series; trajectory.archive_path=result.archive_path
                trajectory.metadata={"status":result.status,"message":result.message,"solver_metadata":result.metadata}
                completed_trajectories+=1
                log_progress()
                if on_event is not None:
                    on_event({"event":"hierarchy_progress","completed_trajectories":completed_trajectories,"total_trajectories":total_trajectories,
                              "system_index":system_index,"cluster_index":cluster_index,"trajectory_index":trajectory_index,"archive_path":result.archive_path})
            cluster.pool()
        system.pool()
    path=SOLVER_ARCHIVE.save_hierarchy(target)
    log_progress(force=True)
    return HierarchyResult(target,results,str(path))
