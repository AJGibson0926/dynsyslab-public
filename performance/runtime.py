"""Compact, serializable identity for numerical benchmark environments."""
from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import os
import platform
import sys
from typing import Any

import numpy as np
import scipy
from threadpoolctl import threadpool_info


def runtime_fingerprint() -> dict[str, Any]:
    """Describe the Python, numerical libraries, hardware, and backends in use."""
    details = StringIO()
    with redirect_stdout(details):
        np.show_runtime()
    from dataclasses import asdict
    from solver.backends import available_solver_backends
    from postprocessing.linalg_backends import available_linear_algebra_backends
    return {
        "schema": "dynsyslab.performance-runtime/v1",
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "python_compiler": platform.python_compiler(),
        "executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "numpy_runtime": details.getvalue().strip(),
        "threadpools": threadpool_info(),
        "solver_backends": {name: asdict(backend) for name,backend in available_solver_backends().items()},
        "linear_algebra_backends": {name:asdict(backend) for name,backend in available_linear_algebra_backends().items()},
    }
