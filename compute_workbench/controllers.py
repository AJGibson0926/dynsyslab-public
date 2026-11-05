"""Controller and target factories used by JSON compute jobs."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, ClassVar

import numpy as np

from control import (ConstantTarget, KoopmanLQRController, KoopmanMPCController,
                     WaveformController, rpe_control_target)
from rom.storage import load_krom_archive


@dataclass
class RandomPulseController:
    """Order-independent zero-order-hold random forcing.

    A pulse is a pure function of its time bin and seed, so adaptive solver
    retries cannot advance hidden RNG state or change the forcing history.
    Settings are expressed in physical system-input units and converted at the
    controller/dynamics boundary when a system is solved nondimensionally.
    """
    amplitude: float
    hold_time: float
    seed: int = 0
    offset: float = 0.0
    duty_cycle: float = 1.0
    actuation: str = "system_input"
    active_time_window: tuple[float, float] | None = None
    archive_replay_safe: ClassVar[bool] = True
    continuous_rhs_safe: ClassVar[bool] = True
    coordinate_system: str = field(init=False, default="dimensional")
    time_units: str = field(init=False, default="dimensional")
    state_units: str = field(init=False, default="dimensional")
    control_units: str = field(init=False, default="system input")
    _system: Any = field(init=False, default=None, repr=False)
    _physical_active_time_window: tuple[float, float] | None = field(init=False, default=None, repr=False)

    def __post_init__(self):
        if self.hold_time <= 0 or self.amplitude < 0 or not 0 <= self.duty_cycle <= 1:
            raise ValueError("Random pulses require hold_time>0, amplitude>=0, and duty_cycle in [0,1].")
        self._physical_active_time_window = self.active_time_window

    def bind_coordinate_system(self, coordinate_system: str, *, system=None) -> None:
        if coordinate_system not in {"dimensionless", "dimensional"}:
            raise ValueError("coordinate_system must be dimensionless or dimensional.")
        if coordinate_system == "dimensionless" and (system is None or
                not callable(getattr(system, "to_dimensional_time", None)) or
                not callable(getattr(system, "to_dimensionless_input", None))):
            raise TypeError("Dimensionless random pulses require system time/input conversions.")
        self.coordinate_system = coordinate_system
        self.time_units = self.state_units = coordinate_system
        self.control_units = "dimensionless system input" if coordinate_system == "dimensionless" else "system input"
        self._system = system
        self.active_time_window = (None if self._physical_active_time_window is None else
            tuple(float(np.asarray(system.to_dimensionless_time(value))) for value in self._physical_active_time_window)
            if coordinate_system == "dimensionless" else self._physical_active_time_window)

    def value(self, t: float, state: np.ndarray) -> np.ndarray:
        physical_time = (float(np.asarray(self._system.to_dimensional_time(t)))
                         if self.coordinate_system == "dimensionless" else float(t))
        index = int(np.floor(physical_time / self.hold_time + 1e-12))
        generator = np.random.default_rng(np.random.SeedSequence((int(self.seed), index)))
        active = generator.random() <= self.duty_cycle
        raw = self.offset + (generator.uniform(-self.amplitude, self.amplitude) if active else 0.0)
        value = np.asarray((raw,), dtype=float)
        return (np.asarray(self._system.to_dimensionless_input(value), dtype=float)
                if self.coordinate_system == "dimensionless" else value)

    def archive_metadata(self) -> dict[str, Any]:
        return {"controller": "random-pulse", "amplitude": self.amplitude,
                "hold_time": self.hold_time, "seed": self.seed,
                "offset": self.offset, "duty_cycle": self.duty_cycle,
                "coordinate_system": self.coordinate_system,
                "control_units": self.control_units}


def _matrix(setting: Any, size: int, name: str) -> np.ndarray:
    value = np.asarray(setting, dtype=float)
    if value.ndim == 0:
        return np.eye(size) * float(value)
    if value.ndim == 1:
        if value.size != size:
            raise ValueError(f"{name} diagonal must contain {size} values.")
        return np.diag(value)
    if value.shape != (size, size):
        raise ValueError(f"{name} must be scalar, a {size}-vector, or a {size}x{size} matrix.")
    return value


def build_target(spec: dict[str, Any], state_dimension: int):
    """Construct a control target from one validated UI job specification."""
    kind = spec.get("kind", "constant")
    if kind == "rpe_canonical":
        return rpe_control_target(spec["name"])
    if kind == "constant":
        state = np.asarray(spec.get("state", np.zeros(state_dimension)), dtype=float)
        return ConstantTarget(state)
    raise ValueError(f"Unsupported target kind {kind!r}.")


def build_controller(spec: dict[str, Any], system: Any):
    """Construct a coordinate-aware controller from a compute job section."""
    kind = spec.get("kind", "none")
    settings = dict(spec.get("settings") or {})
    if kind == "none":
        return None
    if kind == "waveform":
        return WaveformController(
            np.asarray(settings.get("frequencies_hz", [1.0]), dtype=float),
            np.asarray(settings.get("amplitudes", [0.0]), dtype=float),
            np.asarray(settings.get("phases_rad", [0.0]), dtype=float),
            offset=float(settings.get("offset", 0.0)),
            active_time_window=tuple(settings["active_time_window"]) if settings.get("active_time_window") else None)
    if kind == "random_pulse":
        return RandomPulseController(
            amplitude=float(settings.get("amplitude", 1.0)),
            hold_time=float(settings.get("hold_time", 0.01)),
            seed=int(settings.get("seed", 0)), offset=float(settings.get("offset", 0.0)),
            duty_cycle=float(settings.get("duty_cycle", 1.0)),
            active_time_window=tuple(settings["active_time_window"]) if settings.get("active_time_window") else None)
    if kind not in {"klqr", "kmpc"}:
        raise ValueError(f"Unsupported controller {kind!r}.")
    source = settings.get("krom")
    if not source:
        raise ValueError(f"{kind.upper()} requires settings.krom.")
    krom = load_krom_archive(source)
    state_dimension = int(np.asarray(system.initial_state).size)
    target = build_target(spec.get("target") or {}, state_dimension)
    lifted = int(krom.A(output="split", time_coordinates="model").shape[0])
    state = np.asarray(krom.system.dimensionless_initial_state if krom.coordinate_system == "dimensionless"
                       else krom.system.initial_state, dtype=float)
    inputs = int(krom.B_psi(state, output="split", state_coordinates="model",
                            time_coordinates="model", representation=settings.get("control_representation", "default")).shape[1])
    Q = _matrix(settings.get("Q", 1.0), lifted, "Q")
    R = _matrix(settings.get("R", 1.0), inputs, "R")
    common = dict(target=target,
                  active_time_window=tuple(settings["active_time_window"]) if settings.get("active_time_window") else None,
                  control_limit=settings.get("control_limit"),
                  target_domain_policy=settings.get("target_domain_policy", "warn"),
                  control_representation=settings.get("control_representation", "default"))
    if kind == "klqr":
        return KoopmanLQRController(krom, Q, R, **common)
    return KoopmanMPCController(krom, Q, R, horizon=int(settings.get("horizon", 20)),
                                sample_time=float(settings.get("sample_time", 0.001)), **common)
