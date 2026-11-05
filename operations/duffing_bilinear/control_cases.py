"""Bounded isolated KLQR case execution for Duffing KROM candidates."""
from __future__ import annotations

import json
import multiprocessing as mp
from pathlib import Path
from time import perf_counter
from typing import Any
import numpy as np

from control import ConstantTarget, KoopmanLQRController
from rom import load_krom_archive
from solver import InterpolationOptions, SolveRequest, SolverOptions, solve
from .common import (ACTIVE_WINDOW, FINAL_TIME, INITIAL_STATE, TARGET_STATE,
                     _control_metrics, _sha, _write_json)

def _solve_case_local(
    attempt_number: int,
    key: str,
    label: str,
    model,
    representation: str,
    *,
    output_dt: float,
) -> tuple[dict[str, Any], Any]:
    split_dimension = model.A(output="split").shape[0]
    controller = KoopmanLQRController(
        model,
        Q=12500.0 * np.eye(split_dimension),
        R=np.asarray(((0.01,),)),
        target=ConstantTarget(TARGET_STATE),
        target_space="state",
        active_time_window=ACTIVE_WINDOW,
        target_domain_policy="ignore",
        control_representation=representation,
    )
    started = perf_counter()
    result = solve(
        SolveRequest(
            system=model.system,
            kind="ode",
            time_span=(0.0, FINAL_TIME),
            initial_state=INITIAL_STATE,
            controller=controller,
            options=SolverOptions(
                method="BDF",
                rtol=1e-8,
                atol=1e-10,
                max_step=0.01,
                use_dimensionless_form=True,
                control_evaluation="continuous_rhs",
                execution_backend="python",
            ),
            interpolation=InterpolationOptions(
                enabled=True, time_step=output_dt, method="pchip"
            ),
            label=f"duffing-bilinear-search-a{attempt_number:02d}-{key}",
        )
    )
    record: dict[str, Any] = {
        "key": key,
        "label": label,
        "representation": representation,
        "status": result.status,
        "message": result.message,
        "solve_seconds": perf_counter() - started,
    }
    if result.status == "complete":
        record.update(
            {
                "archive_path": str(Path(result.archive_path).resolve()),
                "archive_sha256": _sha(result.archive_path),
                "metrics": _control_metrics(result),
            }
        )
    return record, result

def _solve_case_process(payload: dict[str, Any], record_path: str) -> None:
    """Spawn-safe worker that leaves a durable record before exiting."""
    destination = Path(record_path)
    try:
        model = load_krom_archive(payload["krom_archive"])
        record, _ = _solve_case_local(
            int(payload["attempt_number"]),
            str(payload["key"]),
            str(payload["label"]),
            model,
            str(payload["representation"]),
            output_dt=float(payload["output_dt"]),
        )
    except BaseException as error:
        record = {
            "key": str(payload["key"]),
            "label": str(payload["label"]),
            "representation": str(payload["representation"]),
            "status": "failed",
            "message": f"worker exception: {type(error).__name__}: {error}",
            "solve_seconds": None,
        }
    _write_json(destination, record)

def _solve_case_bounded(
    payload: dict[str, Any], record_path: Path, *, timeout_seconds: float
) -> dict[str, Any]:
    """Evaluate one unchanged controller case in an isolated bounded worker."""
    if timeout_seconds <= 0.0 or not np.isfinite(timeout_seconds):
        raise ValueError("case timeout must be finite and positive.")
    if record_path.exists():
        record_path.unlink()
    started = perf_counter()
    process = mp.get_context("spawn").Process(
        target=_solve_case_process, args=(payload, str(record_path.resolve()))
    )
    process.start()
    process.join(timeout_seconds)
    if process.is_alive():
        process.terminate()
        process.join(10.0)
        return {
            "key": str(payload["key"]),
            "label": str(payload["label"]),
            "representation": str(payload["representation"]),
            "status": "timeout",
            "message": f"controller solve exceeded campaign wall-time limit of {timeout_seconds:g} s",
            "solve_seconds": perf_counter() - started,
            "timeout_seconds": float(timeout_seconds),
        }
    if not record_path.is_file():
        return {
            "key": str(payload["key"]),
            "label": str(payload["label"]),
            "representation": str(payload["representation"]),
            "status": "failed",
            "message": f"controller worker exited with code {process.exitcode} without a record",
            "solve_seconds": perf_counter() - started,
        }
    return json.loads(record_path.read_text(encoding="utf-8"))
