"""Thin adapters from built-in Python system definitions to native arrays."""
from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from typing import Literal

import numpy as np

from dynsyslib.drivers import SinusoidalDriver, ZeroDriver
from dynsyslib.lib.duffing import DuffingSystem
from dynsyslib.lib.rpe import RayleighPlessetSystem
from dynsyslib.lib.van_der_pol import VanDerPolSystem

CoordinateSystem = Literal["dimensional", "dimensionless"]


@dataclass(frozen=True)
class NativeODEModel:
    identifier: str
    parameters: np.ndarray
    coordinate_system: CoordinateSystem


def _driver_parameters(driver) -> tuple[float, float, float] | None:
    if isinstance(driver, ZeroDriver):
        return 0.0, 0.0, 0.0
    if isinstance(driver, SinusoidalDriver):
        return float(driver.amplitude), float(driver.frequency_hz), float(driver.phase_rad)
    return None


def native_model(system, *, coordinate_system: CoordinateSystem = "dimensional") -> NativeODEModel | None:
    if coordinate_system not in {"dimensional", "dimensionless"}:
        raise ValueError("coordinate_system must be 'dimensional' or 'dimensionless'.")
    forcing = _driver_parameters(getattr(system, "driver", None))
    if forcing is None:
        return None
    if coordinate_system == "dimensional" and isinstance(system, DuffingSystem):
        amplitude,frequency,phase=forcing
        return NativeODEModel(
            "duffing",
            np.asarray((system.delta, system.alpha, system.beta,
                        system.input_gain*amplitude,frequency,phase), dtype=np.float64),
            coordinate_system,
        )
    if coordinate_system == "dimensional" and isinstance(system, VanDerPolSystem):
        return NativeODEModel(
            "van_der_pol", np.asarray((system.mu, *forcing), dtype=np.float64), coordinate_system
        )
    if coordinate_system == "dimensionless" and isinstance(system, RayleighPlessetSystem):
        groups = system.dimensionless_parameters
        # Native 0.2 adds its scalar drive. Encode the public external-pressure
        # sign in the adapter without changing the stable native ABI.
        rpe_forcing=(system.pressure_input_sign*forcing[0],forcing[1],forcing[2])
        return NativeODEModel(
            "rpe_dimensionless",
            np.asarray(
                (
                    system.preset.kappa,
                    groups.gas_pressure,
                    groups.vapor_pressure,
                    groups.ambient_pressure,
                    groups.capillarity,
                    groups.viscosity,
                    groups.pressure_scale_pa,
                    groups.time_scale_s,
                    *rpe_forcing,
                ),
                dtype=np.float64,
            ),
            coordinate_system,
        )
    return None


def _installed_model_support(specification: NativeODEModel) -> tuple[bool, str | None]:
    """Reject adapter/binary version skew before backend selection.

    Absence is handled by :func:`solver.backends.select_ode_backend`. Native
    0.1 predates explicit model-list provenance and is known to contain only
    the two original dimensional kernels.
    """
    try:
        module = import_module("_dynsyslab_native")
    except ImportError:
        return True, None
    payload = module.capabilities()
    if not isinstance(payload, dict):
        return False, "installed native extension returned invalid capability metadata"
    advertised = payload.get("fixed_step_models")
    models = {"duffing", "van_der_pol"} if advertised is None else {
        str(item) for item in advertised
    }
    if specification.identifier not in models:
        version = payload.get("version", "unknown")
        return False, (
            f"installed native extension {version} does not provide model "
            f"{specification.identifier!r}"
        )
    return True, None


def native_fixed_eligibility(request, *, dimensionless: bool) -> tuple[bool, str | None]:
    if request.control is not None or request.controller is not None:
        return False, "control and estimation remain on the Python reference backend"
    if request.history.track_during_solve:
        return False, "in-solve checkpoint/history callbacks require the Python backend"
    if request.cancel_token is not None or request.pause_token is not None:
        return False, "cooperative cancellation and pause require the Python backend"
    coordinate_system: CoordinateSystem = "dimensionless" if dimensionless else "dimensional"
    specification = native_model(request.system, coordinate_system=coordinate_system)
    if specification is None:
        return False, (
            f"{type(request.system).__name__} has no native model adapter for "
            f"{coordinate_system} coordinates"
        )
    return _installed_model_support(specification)


def integrate_fixed_native(
    system,
    t0: float,
    tf: float,
    initial_state,
    step_size: float,
    method: str,
    *,
    coordinate_system: CoordinateSystem = "dimensional",
):
    specification = native_model(system, coordinate_system=coordinate_system)
    if specification is None:
        raise TypeError(
            f"{type(system).__name__} has no native {coordinate_system} model adapter."
        )
    module = import_module("_dynsyslab_native")
    time, values = module.integrate_fixed(
        specification.identifier,
        specification.parameters,
        np.ascontiguousarray(initial_state, dtype=np.float64),
        float(t0),
        float(tf),
        float(step_size),
        method,
    )
    return np.asarray(time, dtype=float), np.asarray(values, dtype=float)


__all__ = ["NativeODEModel", "integrate_fixed_native", "native_fixed_eligibility", "native_model"]
