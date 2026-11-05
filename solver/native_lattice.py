"""Model-specific adapters for explicit native lattice execution."""
from __future__ import annotations

from importlib import import_module

import numpy as np

from dynsyslib.lattice import IsingModel


def native_lattice_eligibility(request) -> tuple[bool, str | None]:
    """Explain whether a lattice request matches the admitted native kernel."""
    if not isinstance(request.system, IsingModel):
        return False, f"{type(request.system).__name__} has no native lattice adapter"
    if request.history.track_during_solve:
        return False, "native lattice checkpoint/history integration is not implemented"
    if request.cancel_token is not None or request.pause_token is not None:
        return False, "native lattice cancellation and pause are not implemented"
    return True, None


def run_ising2d_native(system: IsingModel, sweeps: int) -> tuple[np.ndarray, dict]:
    """Run a reproducible native stream and update the live model state."""
    module = import_module("_dynsyslab_native")
    # Draw one seed from the model-owned generator so repeated calls advance
    # deterministically while the native engine remains fully documented.
    seed = int(system.rng.integers(0, np.iinfo(np.uint64).max, dtype=np.uint64))
    values, final_spins = module.ising2d_run(
        np.ascontiguousarray(system.spins, dtype=np.int64),
        float(system.coupling),
        float(system.field),
        float(system.temperature),
        seed,
        int(sweeps),
    )
    system.spins = np.asarray(final_spins, dtype=int)
    return np.asarray(values, dtype=float), {
        "algorithm": "Metropolis single-site",
        "random_engine": "std::mt19937_64",
        "native_seed": seed,
        "trajectory_parity": "statistical-not-bitwise",
    }


__all__ = ["native_lattice_eligibility", "run_ising2d_native"]
