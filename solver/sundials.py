"""Optional adaptive ODE adapter for the official SUNDIALS Python wheel."""
from __future__ import annotations

from importlib import import_module
from typing import Callable

import numpy as np


def adaptive_native_eligibility(request) -> tuple[bool, str | None]:
    """Restrict the first adapter to semantics it can preserve exactly."""
    if request.history.track_during_solve or request.history.mode in {"sampled", "checkpoint"}:
        return False, "SUNDIALS checkpoint/history integration is not implemented"
    if request.cancel_token is not None or request.pause_token is not None:
        return False, "SUNDIALS cooperative cancellation and pause are not implemented"
    return True, None


def integrate_adaptive_sundials(
    rhs: Callable[[float, np.ndarray], np.ndarray],
    t0: float,
    tf: float,
    initial_state: np.ndarray,
    *,
    method: str,
    rtol: float,
    atol: float,
    max_step: float | None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Return accepted CVODE steps using zero-copy N_Vector array views."""
    core = import_module("sundials4py.core")
    cvodes = import_module("sundials4py.cvodes")
    status, context = core.SUNContext_Create(core.SUN_COMM_NULL)
    if status != core.SUN_SUCCESS:
        raise RuntimeError(f"SUNContext_Create failed with status {status}.")
    initial = np.asarray(initial_state, dtype=float).reshape(-1)
    vector = core.N_VNew_Serial(initial.size, context)
    if vector is None:
        raise RuntimeError("N_VNew_Serial could not allocate the solution vector.")
    state_view = core.N_VGetArrayPointer(vector)
    state_view[:] = initial
    callback_error: list[Exception] = []

    def callback(time, state_vector, derivative_vector, _unused):
        try:
            state = core.N_VGetArrayPointer(state_vector)
            derivative = np.asarray(rhs(float(time), state), dtype=float).reshape(-1)
            if derivative.shape != state.shape or np.any(~np.isfinite(derivative)):
                raise ValueError("ODE RHS returned an invalid derivative for SUNDIALS.")
            core.N_VGetArrayPointer(derivative_vector)[:] = derivative
            return 0
        except Exception as error:
            callback_error.append(error)
            return -1

    family = cvodes.CV_BDF if method in {"Radau", "BDF", "LSODA"} else cvodes.CV_ADAMS
    memory = cvodes.CVodeCreate(family, context)
    if memory is None:
        raise RuntimeError("CVodeCreate failed.")
    if cvodes.CVodeInit(memory.get(), callback, float(t0), vector) != cvodes.CV_SUCCESS:
        raise RuntimeError("CVodeInit failed.")
    if cvodes.CVodeSStolerances(memory.get(), float(rtol), float(atol)) != cvodes.CV_SUCCESS:
        raise RuntimeError("CVodeSStolerances failed.")
    matrix = core.SUNDenseMatrix(initial.size, initial.size, context)
    linear_solver = core.SUNLinSol_Dense(vector, matrix, context)
    if matrix is None or linear_solver is None:
        raise RuntimeError("SUNDIALS dense linear solver allocation failed.")
    if cvodes.CVodeSetLinearSolver(memory.get(), linear_solver, matrix) != cvodes.CV_SUCCESS:
        raise RuntimeError("CVodeSetLinearSolver failed.")
    if max_step is not None and cvodes.CVodeSetMaxStep(memory.get(), float(max_step)) != cvodes.CV_SUCCESS:
        raise RuntimeError("CVodeSetMaxStep failed.")
    if cvodes.CVodeSetStopTime(memory.get(), float(tf)) != cvodes.CV_SUCCESS:
        raise RuntimeError("CVodeSetStopTime failed.")

    times = [float(t0)]
    values = [initial.copy()]
    while times[-1] < tf:
        step_status, reached = cvodes.CVode(memory.get(), float(tf), vector, cvodes.CV_ONE_STEP)
        if step_status < 0:
            if callback_error:
                raise RuntimeError(f"SUNDIALS RHS callback failed: {callback_error[-1]}") from callback_error[-1]
            raise RuntimeError(f"CVode failed with status {step_status}.")
        reached = float(reached)
        if reached <= times[-1]:
            raise RuntimeError("CVode returned a non-increasing accepted-step time.")
        times.append(reached)
        values.append(np.asarray(state_view, dtype=float).copy())

    _, steps = cvodes.CVodeGetNumSteps(memory.get())
    _, evaluations = cvodes.CVodeGetNumRhsEvals(memory.get())
    return np.asarray(times), np.asarray(values), {
        "library": "SUNDIALS",
        "interface": "sundials4py",
        "algorithm": "CVODE-BDF" if family == cvodes.CV_BDF else "CVODE-Adams",
        "accepted_steps": int(steps),
        "rhs_evaluations": int(evaluations),
    }


__all__ = ["adaptive_native_eligibility", "integrate_adaptive_sundials"]
