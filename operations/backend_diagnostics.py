"""Cross-package CLI/API for optional numerical execution capabilities."""
from __future__ import annotations

from dataclasses import asdict
import json
import os

from rom.backends import native_rom_backend
from postprocessing.linalg_backends import available_linear_algebra_backends

from solver.backends import available_solver_backends


def backend_diagnostics() -> dict:
    return {
        "schema": "dynsyslab.backend-diagnostics/v1",
        "solver": {name: asdict(backend) for name,backend in available_solver_backends().items()},
        "rom": {"python": {"name": "numpy-reference", "available": True},
                "native": asdict(native_rom_backend())},
        "postprocessing": {"linear_algebra": {
            name:asdict(backend) for name,backend in available_linear_algebra_backends().items()}},
        "resource_environment": {
            name: os.environ.get(name) for name in (
                "DYNSYSLAB_ALLOCATED_CORES", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
            )
        },
    }


def main() -> int:
    print(json.dumps(backend_diagnostics(),indent=2,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())


__all__=["backend_diagnostics","main"]
