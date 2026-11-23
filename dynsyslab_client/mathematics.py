"""Thin client for the compiled, stateless mathematics query API.

No arithmetic, graph, Fourier or distribution algorithms live in Python.
"""
from __future__ import annotations

from collections.abc import Mapping
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any

from .workflow import installation_root


def evolve_fourier_channels(channels, time, *, initial_phase=0.0, initial_time=0.0,
                            executable=None, timeout=120):
    """Advance supplied finite Fourier channels through their native clocks.

    Returns channel-major [real, imaginary] samples. Native C++ owns evolution;
    this adapter only batches calls. Each channel may have its own frequency.
    The initial phase and time are explicit and shared. No ODE is integrated.
    """
    times=list(time)
    return [query(dict(operation="fourier",action="evolve",series=s,time=times,
                       initial_phase=initial_phase,initial_time=initial_time),
                  executable=executable,timeout=timeout)["result"]["values"] for s in channels]


def query(request: Mapping[str, Any], *, executable: str | Path | None = None,
          timeout: float = 120) -> dict[str, Any]:
    """Run a versioned math query using the installed native CLI.

    Pass executable to select a source-build CLI explicitly. Exact rationals
    are strings; this adapter never infers rationality from floating-point data.
    Native errors and subprocess timeouts propagate; temporary input is removed.
    """
    cli = (Path(executable).expanduser().resolve() if executable else
           installation_root() / "bin" / ("dynsyslab-cli.exe" if os.name == "nt"
                                          else "dynsyslab-cli"))
    if not cli.is_file():
        raise FileNotFoundError(f"Native mathematics CLI is missing: {cli}")
    payload = {"schema": "dynsyslab.math-query/v1", **dict(request)}
    with tempfile.TemporaryDirectory(prefix="dynsyslab-math-") as temporary:
        path = Path(temporary) / "query.json"
        path.write_text(json.dumps(payload, allow_nan=False), encoding="utf-8")
        result = subprocess.run([str(cli), "mathematics", str(path)],
                                capture_output=True, text=True, encoding="utf-8",
                                timeout=timeout, check=False)
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip() or
                           f"Native mathematics query exited {result.returncode}")
    response = json.loads(result.stdout)
    if response.get("schema") != "dynsyslab.mathematics/v1":
        raise RuntimeError("Unexpected native mathematics response schema")
    return response
