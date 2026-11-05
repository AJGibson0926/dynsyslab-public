"""Delay-differential pathway; systems own history-aware stepping."""
from __future__ import annotations
import numpy as np
from time_series import TimeSeries
from typing import Any
from .models import NumericalResult, SolveRequest
from .backends import select_reference_backend

def solve_dde(request: SolveRequest, session: Any) -> NumericalResult:
    """Execute one history-aware delay-system request through its model API."""
    if request.time_span is None: raise ValueError("DDE solve requires time_span.")
    backend=select_reference_backend(request.options.execution_backend,"dde")
    duration=request.time_span[1]-request.time_span[0]; dt=request.options.step_size or .05
    try:
        interval=request.history.checkpoint_interval or max(dt, duration/20)
        checkpoint_steps=max(1,int(round(interval/dt)))
        callback=(lambda time,values: session.checkpoint(time+request.time_span[0],values)) if request.history.track_during_solve else None
        try:
            time,values=request.system.integrate(duration,dt,checkpoint_steps=checkpoint_steps,on_checkpoint=callback)
        except TypeError:  # Third-party DDE systems may expose only the minimal API.
            time,values=request.system.integrate(duration,dt)
        time=np.asarray(time)+request.time_span[0]; values=np.asarray(values)
        if request.history.track_during_solve: session.checkpoint(time,values)
        return NumericalResult(TimeSeries(time,values),"complete","solve completed",metadata={"execution_backend":backend.provenance(requested=request.options.execution_backend)})
    except Exception as exc: return NumericalResult(None,"failed",str(exc),metadata={"execution_backend":backend.provenance(requested=request.options.execution_backend)})
