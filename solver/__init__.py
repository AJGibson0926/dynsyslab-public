"""Standalone numerical solver entry point for dynsyslab."""
from __future__ import annotations

from typing import Any

from .models import HistoryOptions, InterpolationOptions, NumericalResult, SolveRequest, SolverOptions
from .backends import SolverBackendCapabilities, available_solver_backends


class _InMemorySession:
    """Minimal non-persistent session accepted by the ODE numerical kernel."""
    def __init__(self,request: SolveRequest): self.request=request
    def telemetry(self,*args,**kwargs): return None
    def checkpoint(self,*args,**kwargs):
        raise RuntimeError("In-memory ODE solves do not support checkpoints.")


def solve(request: Any, **kwargs) -> NumericalResult | Any:
    """Solve once and archive every completed result under archives/<system>/.

    With history tracking enabled, partial checkpoints are written as the solver
    advances; otherwise no data are written until successful completion.
    """
    from hierarchy import Batch, System

    if isinstance(request,(System,Batch)):
        from .hierarchy import solve_hierarchy as _solve_hierarchy
        return _solve_hierarchy(request,**kwargs)
    if kwargs: raise TypeError("Keyword arguments are only accepted when solving a hierarchy.")
    if not isinstance(request,SolveRequest):
        raise TypeError("solve expects a SolveRequest, hierarchy.System, or hierarchy.Batch.")
    from archiver.gateway import SOLVER_ARCHIVE
    from .dde import solve_dde
    from .lattice import solve_lattice
    from .lifecycle import finalize_result
    from .ode import solve_ode
    from .pde import solve_pde

    session=SOLVER_ARCHIVE.new_session(request)
    session.telemetry("started",time_span=request.time_span,method=request.options.method)
    handlers={"ode":solve_ode,"dde":solve_dde,"pde":solve_pde,"lattice":solve_lattice}
    result=handlers[request.kind](request,session)
    return finalize_result(request,session,result)


def solve_hierarchy(*args,**kwargs):
    """Lazy public adapter for :func:`solver.hierarchy.solve_hierarchy`."""
    from .hierarchy import solve_hierarchy as implementation
    return implementation(*args,**kwargs)


def solve_ode_in_memory(request: SolveRequest) -> NumericalResult:
    """Run the canonical ODE kernel without archive or Viewer side effects.

    This boundary is intended for bounded diagnostics and campaign screening,
    not production results. Durable histories and postprocessors are rejected.
    """
    if not isinstance(request,SolveRequest) or request.kind!="ode":
        raise TypeError("solve_ode_in_memory requires an ODE SolveRequest.")
    if request.history.track_during_solve or request.history.mode!="telemetry":
        raise ValueError("In-memory ODE solves require non-tracking telemetry history.")
    if request.postprocessors:
        raise ValueError("In-memory ODE solves cannot run postprocessors.")
    from .ode import solve_ode
    return solve_ode(request,_InMemorySession(request))


def resume_ode_from_pause(*args,**kwargs):
    """Lazy public adapter for durable ODE checkpoint resumption."""
    from .checkpoints import resume_ode_from_pause as implementation
    return implementation(*args,**kwargs)


def __getattr__(name: str):
    """Load public value types without importing numerical engines eagerly."""
    if name=="HierarchyResult":
        from .hierarchy import HierarchyResult
        return HierarchyResult
    if name=="TimeSeries":
        from time_series import TimeSeries
        return TimeSeries
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__=["solve","solve_hierarchy","solve_ode_in_memory","resume_ode_from_pause","HierarchyResult","SolveRequest","SolverOptions","HistoryOptions","InterpolationOptions","NumericalResult","TimeSeries","SolverBackendCapabilities","available_solver_backends"]
