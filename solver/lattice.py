"""Stochastic/lattice pathway for models exposing run(steps)."""
from __future__ import annotations
import numpy as np
from time_series import TimeSeries
from typing import Any
from .models import NumericalResult, SolveRequest
from .backends import select_lattice_backend
from .native_lattice import native_lattice_eligibility, run_ising2d_native

def solve_lattice(request: SolveRequest, session: Any) -> NumericalResult:
    """Execute one stochastic lattice request with optional checkpoints."""
    if request.time_span is None: raise ValueError("Lattice solve uses time_span=(start_sweep,end_sweep).")
    eligible,reason=native_lattice_eligibility(request)
    backend=select_lattice_backend(request.options.execution_backend,native_eligible=eligible,ineligibility_reason=reason)
    try:
        steps=int(request.time_span[1]-request.time_span[0]); start=request.time_span[0]
        if request.history.track_during_solve and hasattr(request.system,"sweep"):
            interval=max(1,int(request.history.checkpoint_interval or max(1,steps//20))); pieces=[]
            for offset in range(0,steps,interval):
                count=min(interval,steps-offset); chunk=request.system.run(count); pieces.append(chunk)
                session.checkpoint(np.arange(offset,offset+count)+start,chunk)
            values=np.vstack(pieces)
        elif backend.compiled:
            values,lattice_metadata=run_ising2d_native(request.system,steps)
        else: values=np.asarray(request.system.run(steps))
        time=np.arange(steps)+start
        return NumericalResult(TimeSeries(time,values),"complete","solve completed",metadata={"lattice_backend":locals().get("lattice_metadata"),"execution_backend":backend.provenance(requested=request.options.execution_backend)})
    except Exception as exc: return NumericalResult(None,"failed",str(exc),metadata={"execution_backend":backend.provenance(requested=request.options.execution_backend)})
