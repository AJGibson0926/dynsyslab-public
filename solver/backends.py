"""Numerical execution-backend selection for solver-owned kernels.

This module describes capabilities and policy.  It deliberately owns no job
scheduling, archive I/O, or visualization behavior.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import Any, Literal

BackendPreference = Literal["auto", "python", "native"]


@dataclass(frozen=True)
class SolverBackendCapabilities:
    """Serializable capability declaration for one numerical backend."""
    name: str
    available: bool
    compiled: bool
    fixed_step_ode: bool
    adaptive_ode: bool = False
    batched_ode: bool = False
    fixed_step_models: tuple[str, ...] = ()
    pde: bool = False
    lattice: bool = False
    version: str | None = None
    reason: str | None = None
    auto_admitted: bool = True

    def provenance(self, *, requested: BackendPreference) -> dict[str, Any]:
        return {"requested": requested, **asdict(self)}


PYTHON_BACKEND = SolverBackendCapabilities(
    name="python-reference",
    available=True,
    compiled=False,
    fixed_step_ode=True,
    adaptive_ode=True,
    pde=True,
    lattice=True,
)


def native_backend() -> SolverBackendCapabilities:
    """Describe the optional extension without making it an installation requirement."""
    try:
        module = import_module("_dynsyslab_native")
    except ImportError as error:
        return SolverBackendCapabilities(
            name="dynsyslab-native",
            available=False,
            compiled=True,
            fixed_step_ode=False,
            reason=f"optional extension is not installed: {error}",
        )
    payload = module.capabilities()
    if not isinstance(payload, dict):
        raise TypeError("_dynsyslab_native.capabilities() must return a dictionary.")
    return SolverBackendCapabilities(
        name="dynsyslab-native",
        available=True,
        compiled=True,
        fixed_step_ode=bool(payload.get("fixed_step_ode", False)),
        adaptive_ode=bool(payload.get("adaptive_ode", False)),
        batched_ode=bool(payload.get("batched_ode", False)),
        fixed_step_models=tuple(str(item) for item in payload.get("fixed_step_models", ())),
        pde=bool(payload.get("pde", False)),
        lattice=bool(payload.get("lattice", False)),
        version=str(payload.get("version", "unknown")),
    )


def sundials_backend() -> SolverBackendCapabilities:
    """Describe the maintained official adaptive backend when installed."""
    try:
        import_module("sundials4py.cvodes")
        installed = version("sundials4py")
    except (ImportError, PackageNotFoundError) as error:
        return SolverBackendCapabilities(
            name="sundials-cvode",
            available=False,
            compiled=True,
            fixed_step_ode=False,
            reason=f"optional sundials4py backend is not installed: {error}",
            auto_admitted=False,
        )
    return SolverBackendCapabilities(
        name="sundials-cvode",
        available=True,
        compiled=True,
        fixed_step_ode=False,
        adaptive_ode=True,
        version=installed,
        auto_admitted=False,
    )


def available_solver_backends() -> dict[str, SolverBackendCapabilities]:
    """Inspect all maintained solver backends without selecting one."""
    return {"python": PYTHON_BACKEND, "native": native_backend(), "sundials": sundials_backend()}


def select_ode_backend(
    preference: BackendPreference,
    method: str,
    *,
    native_eligible: bool = True,
    ineligibility_reason: str | None = None,
) -> SolverBackendCapabilities:
    """Resolve one backend without silently weakening an explicit request."""
    if preference not in {"auto", "python", "native"}:
        raise ValueError("execution_backend must be 'auto', 'python', or 'native'.")
    fixed = method in {"Euler", "RK4"}
    if preference == "python":
        return PYTHON_BACKEND
    candidate = native_backend() if fixed else sundials_backend()
    supported = candidate.fixed_step_ode if fixed else candidate.adaptive_ode
    if preference == "native":
        if not candidate.available:
            raise RuntimeError(candidate.reason or "The native solver backend is unavailable.")
        if not supported:
            family = "fixed-step" if fixed else "adaptive"
            raise RuntimeError(f"The native solver backend does not support {family} ODE method {method!r}.")
        if not native_eligible:
            raise RuntimeError(ineligibility_reason or "This solve request is not eligible for the native backend.")
        return candidate
    return candidate if candidate.available and supported and native_eligible and candidate.auto_admitted else PYTHON_BACKEND


def select_reference_backend(preference: BackendPreference, kind: str) -> SolverBackendCapabilities:
    """Resolve kinds that do not yet have an admitted native implementation."""
    if preference not in {"auto", "python", "native"}:
        raise ValueError("execution_backend must be 'auto', 'python', or 'native'.")
    if preference == "native":
        candidate = native_backend()
        if not candidate.available:
            raise RuntimeError(candidate.reason or "The native solver backend is unavailable.")
        raise RuntimeError(f"The native solver backend does not support {kind!r} solves yet.")
    return PYTHON_BACKEND


def select_lattice_backend(
    preference: BackendPreference,
    *,
    native_eligible: bool,
    ineligibility_reason: str | None = None,
) -> SolverBackendCapabilities:
    """Select explicit native lattice execution without changing RNG semantics in auto mode."""
    if preference not in {"auto", "python", "native"}:
        raise ValueError("execution_backend must be 'auto', 'python', or 'native'.")
    if preference in {"auto", "python"}:
        return PYTHON_BACKEND
    candidate = native_backend()
    if not candidate.available:
        raise RuntimeError(candidate.reason or "The native solver backend is unavailable.")
    if not candidate.lattice:
        raise RuntimeError("The native solver backend does not support lattice solves.")
    if not native_eligible:
        raise RuntimeError(ineligibility_reason or "This lattice request is not eligible for the native backend.")
    return candidate


__all__ = [
    "BackendPreference",
    "SolverBackendCapabilities",
    "PYTHON_BACKEND",
    "available_solver_backends",
    "native_backend",
    "select_ode_backend",
    "select_lattice_backend",
    "select_reference_backend",
    "sundials_backend",
]
