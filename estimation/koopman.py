"""Observable-driven Koopman linear-quadratic estimation.

KLQE is independent of controllers: it can run offline against a measurement
time series, inside a sampled solver workflow, or through immutable design and
result artifacts.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Callable
import hashlib

import numpy as np

from rom import KROM
from time_series import TimeSeries


from .klqe_storage import (load_klqe_design, load_klqe_result,
                           save_klqe_design, save_klqe_result)


def _covariance(value: np.ndarray, dimension: int, name: str, *, definite: bool) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.shape != (dimension, dimension) or np.any(~np.isfinite(result)):
        raise ValueError(f"{name} must be a finite ({dimension}, {dimension}) matrix.")
    if not np.allclose(result, result.T, rtol=1e-11, atol=1e-13):
        raise ValueError(f"{name} must be symmetric.")
    eigenvalues = np.linalg.eigvalsh(0.5 * (result + result.T))
    threshold = 128.0 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(eigenvalues))))
    if definite and np.min(eigenvalues) <= threshold:
        raise ValueError(f"{name} must be positive definite.")
    if not definite and np.min(eigenvalues) < -threshold:
        raise ValueError(f"{name} must be positive semidefinite.")
    return result


def _measurement_matrix(value: np.ndarray, state_dimension: int) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.ndim != 2 or result.shape[0] < 1 or result.shape[1] != state_dimension:
        raise ValueError("C must have shape (measurement_dimension, lifted_dimension).")
    if np.any(~np.isfinite(result)):
        raise ValueError("C must contain only finite values.")
    return result


def _detectability(A: np.ndarray, C: np.ndarray) -> tuple[bool, list[dict[str, Any]]]:
    dimension = A.shape[0]; failures = []
    tolerance = 1e-10 * max(1.0, float(np.linalg.norm(A, ord=2)))
    for pole in np.linalg.eigvals(A):
        if np.real(pole) < -tolerance:
            continue
        pbh = np.vstack((pole * np.eye(dimension) - A, C)).astype(complex)
        rank = int(np.linalg.matrix_rank(pbh))
        if rank < dimension:
            failures.append({"pole_real": float(np.real(pole)), "pole_imag": float(np.imag(pole)), "rank": rank})
    return not failures, failures


def _state_decoder(krom: KROM, reconstruction) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    from rom.observation import ObservableReconstruction, StateObservationMap

    if isinstance(reconstruction, ObservableReconstruction):
        if reconstruction.observable_dimension != np.asarray(krom.system.initial_state).size:
            raise ValueError("The state reconstruction must output the full dynamical state.")
        return reconstruction.split_matrix(krom.eigenfunctions), np.asarray(reconstruction.intercept), reconstruction.summary()
    if isinstance(reconstruction, StateObservationMap):
        return np.asarray(reconstruction.coefficients).T, np.asarray(reconstruction.intercept), reconstruction.summary()
    raise TypeError("state_reconstruction must be ObservableReconstruction or StateObservationMap.")


@dataclass(frozen=True)
class KLQEObservationModel:
    """Explicit measurement and state reconstruction maps for one KROM."""

    rom: KROM
    measurement_reconstruction: Any
    state_reconstruction: Any
    A: np.ndarray = field(init=False)
    C_measurement: np.ndarray = field(init=False)
    measurement_intercept: np.ndarray = field(init=False)
    C_state: np.ndarray = field(init=False)
    state_intercept: np.ndarray = field(init=False)

    def __post_init__(self):
        from rom.observation import ObservableReconstruction

        if not isinstance(self.measurement_reconstruction, ObservableReconstruction):
            raise TypeError("measurement_reconstruction must be an ObservableReconstruction.")
        if self.measurement_reconstruction.lifted_dimension != self.rom.lifted_dimension:
            raise ValueError("Measurement reconstruction does not match the KROM.")
        A = np.asarray(self.rom.A(output="split"), dtype=float)
        C_measurement = self.measurement_reconstruction.split_matrix(self.rom.eigenfunctions)
        C_state, state_intercept, _ = _state_decoder(self.rom, self.state_reconstruction)
        if C_measurement.shape[1] != A.shape[0] or C_state.shape[1] != A.shape[0]:
            raise ValueError("KLQE observation matrices do not match the split lifted dimension.")
        values = {
            "A": A, "C_measurement": C_measurement,
            "measurement_intercept": self.measurement_reconstruction.intercept,
            "C_state": C_state, "state_intercept": state_intercept,
        }
        for name, value in values.items():
            array = np.array(value, dtype=float, copy=True)
            if np.any(~np.isfinite(array)):
                raise ValueError(f"{name} contains non-finite values.")
            array.setflags(write=False); object.__setattr__(self, name, array)

    @property
    def measurement_observable(self):
        return self.measurement_reconstruction.observable

    def measure(self, state: np.ndarray) -> np.ndarray:
        return np.asarray(self.measurement_observable.evaluate(state), dtype=float).reshape(-1)

    def reconstruct_state(self, lifted_state: np.ndarray) -> np.ndarray:
        return self.C_state @ np.asarray(lifted_state, dtype=float) + self.state_intercept

    def summary(self) -> dict[str, Any]:
        _, _, state_summary = _state_decoder(self.rom, self.state_reconstruction)
        return {
            "encoding": "koopman-klqe-observation-model/v1",
            "coordinate_system": self.rom.coordinate_system,
            "lifted_dimension": int(self.A.shape[0]),
            "measurement_dimension": int(self.C_measurement.shape[0]),
            "state_dimension": int(self.C_state.shape[0]),
            "measurement_reconstruction": self.measurement_reconstruction.summary(),
            "state_reconstruction": state_summary,
        }


@dataclass(frozen=True)
class KLQEDesign:
    """Validated steady-state continuous-time KLQE design."""

    observation_model: KLQEObservationModel
    Q: np.ndarray
    R: np.ndarray
    covariance: np.ndarray
    gain: np.ndarray
    observer_poles: np.ndarray
    care_residual_norm: float
    detectable: bool
    detectability_failures: tuple[dict[str, Any], ...] = ()

    def __post_init__(self):
        dimension = self.observation_model.A.shape[0]
        measurements = self.observation_model.C_measurement.shape[0]
        arrays = {
            "Q": _covariance(self.Q, dimension, "Q", definite=False),
            "R": _covariance(self.R, measurements, "R", definite=True),
            "covariance": np.asarray(self.covariance, dtype=float),
            "gain": np.asarray(self.gain, dtype=float),
            "observer_poles": np.asarray(self.observer_poles, dtype=complex).reshape(-1),
        }
        if arrays["covariance"].shape != (dimension, dimension):
            raise ValueError("KLQE covariance has the wrong shape.")
        if arrays["gain"].shape != (dimension, measurements):
            raise ValueError("KLQE gain has the wrong shape.")
        if arrays["observer_poles"].size != dimension:
            raise ValueError("KLQE observer pole count does not match the lifted dimension.")
        for name, value in arrays.items():
            copied = np.array(value, copy=True); copied.setflags(write=False)
            object.__setattr__(self, name, copied)
        object.__setattr__(self, "detectability_failures", tuple(dict(item) for item in self.detectability_failures))

    def estimator(self, initial_state: np.ndarray, *, state_coordinates: str="model") -> "KLQEEstimator":
        model = self.observation_model
        return KLQEEstimator(
            model.rom, model.C_measurement, self.Q, self.R, initial_state,
            model.reconstruct_state, measurement_map=model.measure,
            measurement_offset=model.measurement_intercept,
            gain_override=self.gain, covariance_override=self.covariance,
            initial_state_coordinates=state_coordinates,
            design_metadata=self.summary(),
        )

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "koopman-klqe-design/v1",
            "observation_model": self.observation_model.summary(),
            "Q": self.Q.tolist(), "R": self.R.tolist(),
            "gain": self.gain.tolist(), "covariance": self.covariance.tolist(),
            "observer_poles": [{"real": float(v.real), "imag": float(v.imag)} for v in self.observer_poles],
            "care_residual_norm": float(self.care_residual_norm),
            "detectable": bool(self.detectable),
            "detectability_failures": [dict(item) for item in self.detectability_failures],
            "stable": bool(np.all(np.real(self.observer_poles) < 0.0)),
        }


def design_klqe(
    rom: KROM,
    measurement_reconstruction,
    state_reconstruction,
    Q: np.ndarray,
    R: np.ndarray,
    *,
    require_detectable: bool = True,
    require_stable: bool = True,
    care_residual_tolerance: float = 1e-8,
) -> KLQEDesign:
    """Design a continuous-time steady-state KLQE from declared observables."""
    from scipy.linalg import solve_continuous_are

    model = KLQEObservationModel(rom, measurement_reconstruction, state_reconstruction)
    A = model.A; C = model.C_measurement
    process = _covariance(Q, A.shape[0], "Q", definite=False)
    measurement = _covariance(R, C.shape[0], "R", definite=True)
    detectable, failures = _detectability(A, C)
    if require_detectable and not detectable:
        raise ValueError(f"The KLQE pair (A, C) is not detectable: {failures}")
    covariance = solve_continuous_are(A.T, C.T, process, measurement)
    gain = np.linalg.solve(measurement, C @ covariance).T
    correction = covariance @ C.T @ np.linalg.solve(measurement, C)
    residual = A @ covariance + covariance @ A.T - correction @ covariance + process
    residual_norm = float(np.linalg.norm(residual, ord="fro") / max(1.0, np.linalg.norm(process, ord="fro")))
    if not np.isfinite(care_residual_tolerance) or care_residual_tolerance <= 0.0:
        raise ValueError("care_residual_tolerance must be finite and positive.")
    if not np.isfinite(residual_norm) or residual_norm > care_residual_tolerance:
        raise ValueError(
            f"KLQE CARE residual {residual_norm:.3e} exceeds tolerance {care_residual_tolerance:.3e}."
        )
    poles = np.linalg.eigvals(A - gain @ C)
    if require_stable and np.any(np.real(poles) >= 0.0):
        raise ValueError("The designed KLQE observer has a non-decaying pole.")
    return KLQEDesign(model, process, measurement, covariance, gain, poles,
                      residual_norm, detectable, tuple(failures))


def design_klqe_from_krom(
    rom: KROM,
    measurement: str,
    state_output: str,
    Q: np.ndarray,
    R: np.ndarray,
    **options,
) -> KLQEDesign:
    """Design KLQE by resolving two named observable projections on a KROM."""
    measurement_map=rom.observable_map(measurement)
    if state_output in rom.observable_reconstructions:
        state_map=rom.observable_map(state_output)
    elif state_output=="state" and rom.observation_map is not None:
        state_map=rom.observation_map
    else:
        raise KeyError(
            f"State output {state_output!r} is unavailable; named maps={sorted(rom.observable_reconstructions)}."
        )
    return design_klqe(rom,measurement_map,state_map,Q,R,**options)


@dataclass
class KLQEEstimator:
    """Steady-state KLQE with exact zero-order-held measurement propagation.

    The first seven fields preserve the original lightweight constructor. New
    code should prefer :func:`design_klqe`, which binds explicit observables and
    applies detectability and stability gates.
    """

    rom: KROM
    C: np.ndarray
    Q: np.ndarray
    R: np.ndarray
    initial_state: np.ndarray
    reconstruct: Callable[[np.ndarray], np.ndarray]
    measurement_map: Callable[[np.ndarray], np.ndarray] | None = None
    measurement_offset: np.ndarray | None = None
    gain_override: np.ndarray | None = field(default=None, repr=False)
    covariance_override: np.ndarray | None = field(default=None, repr=False)
    initial_state_coordinates: str = "model"
    design_metadata: dict[str, Any] | None = None
    _state: np.ndarray = field(init=False, repr=False)
    _time: float | None = field(init=False, default=None, repr=False)
    _measurement: np.ndarray | None = field(init=False, default=None, repr=False)
    _transition_cache: dict[float, tuple[np.ndarray, np.ndarray]] = field(init=False, default_factory=dict, repr=False)
    gain: np.ndarray = field(init=False)
    covariance: np.ndarray = field(init=False)
    observer_poles: np.ndarray = field(init=False)
    care_residual_norm: float = field(init=False)
    requires_sample_clock: bool = field(init=False, default=True, repr=False)

    def __post_init__(self):
        from scipy.linalg import solve_continuous_are

        A = np.asarray(self.rom.A(output="split"), dtype=float)
        C = _measurement_matrix(self.C, A.shape[0])
        Q = _covariance(self.Q, A.shape[0], "Q", definite=False)
        R = _covariance(self.R, C.shape[0], "R", definite=True)
        covariance = (solve_continuous_are(A.T, C.T, Q, R) if self.covariance_override is None
                      else np.asarray(self.covariance_override, dtype=float))
        gain = (np.linalg.solve(R, C @ covariance).T if self.gain_override is None
                else np.asarray(self.gain_override, dtype=float))
        if covariance.shape != A.shape or gain.shape != (A.shape[0], C.shape[0]):
            raise ValueError("KLQE covariance or gain override has an incompatible shape.")
        correction = covariance @ C.T @ np.linalg.solve(R, C)
        residual = A @ covariance + covariance @ A.T - correction @ covariance + Q
        self.care_residual_norm = float(np.linalg.norm(residual, ord="fro") / max(1.0, np.linalg.norm(Q, ord="fro")))
        self.observer_poles = np.linalg.eigvals(A - gain @ C)
        self.gain = np.array(gain, copy=True); self.covariance = np.array(covariance, copy=True)
        self.C = C; self.Q = Q; self.R = R
        if self.initial_state_coordinates not in {"model","physical"}:
            raise ValueError("initial_state_coordinates must be 'model' or 'physical'.")
        self._state = np.asarray(self.rom.psi(np.asarray(self.initial_state, dtype=float), output="split",
                                             state_coordinates=self.initial_state_coordinates), dtype=float).reshape(-1)
        if self._state.size != A.shape[0]:
            raise ValueError("KLQE initial lift does not match the generator dimension.")
        if self.measurement_map is None:
            self.measurement_map = lambda x: self.rom.psi(x, output="split")
        offset = np.zeros(C.shape[0]) if self.measurement_offset is None else np.asarray(self.measurement_offset, dtype=float).reshape(-1)
        if offset.size != C.shape[0] or np.any(~np.isfinite(offset)):
            raise ValueError("measurement_offset must match the KLQE measurement dimension.")
        self.measurement_offset = offset

    @property
    def lifted_state(self) -> np.ndarray:
        return self._state.copy()

    @property
    def reconstructed_state(self) -> np.ndarray:
        return np.asarray(self.reconstruct(self._state), dtype=float).reshape(-1)

    @property
    def predicted_measurement(self) -> np.ndarray:
        return self.C @ self._state + self.measurement_offset

    def reset(self, initial_state: np.ndarray | None = None) -> None:
        state = self.initial_state if initial_state is None else initial_state
        self._state = np.asarray(self.rom.psi(np.asarray(state, dtype=float), output="split",
                                             state_coordinates=self.initial_state_coordinates), dtype=float)
        self._time = None; self._transition_cache.clear()
        self._measurement = None

    def _transition(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        from scipy.linalg import expm

        key = float(dt); cached = self._transition_cache.get(key)
        if cached is not None:
            return cached
        A = np.asarray(self.rom.A(output="split"), dtype=float)
        closed = A - self.gain @ self.C
        n = closed.shape[0]; p = self.C.shape[0]
        augmented = np.block([[closed, self.gain], [np.zeros((p, n + p))]])
        propagated = expm(augmented * key)
        result = propagated[:n, :n], propagated[:n, n:]
        self._transition_cache[key] = result
        return result

    def step_measurement(self, t: float, measurement: np.ndarray) -> np.ndarray:
        current = float(t)
        if not np.isfinite(current):
            raise ValueError("KLQE time must be finite.")
        observed = np.asarray(measurement, dtype=float).reshape(-1)
        if observed.size != self.C.shape[0] or np.any(~np.isfinite(observed)):
            raise ValueError("KLQE measurement has the wrong dimension or non-finite values.")
        if self._time is not None and current < self._time:
            raise ValueError("KLQE measurement times must be nondecreasing.")
        if self._time is not None and current > self._time:
            transition, input_map = self._transition(current - self._time)
            held = observed if self._measurement is None else self._measurement
            self._state = transition @ self._state + input_map @ (held - self.measurement_offset)
        self._measurement = observed.copy()
        self._time = current if self._time is None or current > self._time else self._time
        return self.reconstructed_state

    def estimate(self, t: float, true_state: np.ndarray) -> np.ndarray:
        measurement = np.asarray(self.measurement_map(np.asarray(true_state, dtype=float)), dtype=float)
        return self.step_measurement(t, measurement)

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "koopman-klqe-estimator/v2",
            "integration": "exact-zero-order-held-measurement",
            "lifted_dimension": int(self._state.size),
            "measurement_dimension": int(self.C.shape[0]),
            "gain": self.gain.tolist(),
            "observer_poles": [{"real": float(v.real), "imag": float(v.imag)} for v in self.observer_poles],
            "stable": bool(np.all(np.real(self.observer_poles) < 0.0)),
            "care_residual_norm": float(self.care_residual_norm),
            "design":None if self.design_metadata is None else dict(self.design_metadata),
        }


@dataclass(frozen=True)
class KLQEResult:
    """Time-aligned estimates, measurements, innovations, and diagnostics."""

    lifted_estimate: TimeSeries
    state_estimate: TimeSeries
    measurement: TimeSeries
    measurement_estimate: TimeSeries
    innovation: TimeSeries
    metadata: dict[str, Any]
    true_state: TimeSeries | None = None
    state_error: TimeSeries | None = None

    def __post_init__(self):
        time = self.lifted_estimate.time
        series = (self.state_estimate, self.measurement, self.measurement_estimate,
                  self.innovation, self.true_state, self.state_error)
        if any(item is not None and not np.array_equal(item.time, time) for item in series):
            raise ValueError("All KLQE result series must share exactly one time grid.")
        object.__setattr__(self, "metadata", dict(self.metadata))

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "koopman-klqe-result/v1",
            "sample_count": self.lifted_estimate.sample_count,
            "lifted_dimension": int(self.lifted_estimate.values.shape[1]),
            "state_dimension": int(self.state_estimate.values.shape[1]),
            "measurement_dimension": int(self.measurement.values.shape[1]),
            "has_truth": self.true_state is not None,
            "metadata": dict(self.metadata),
        }


def run_klqe(estimator: KLQEEstimator, time: np.ndarray, measurements: np.ndarray, *,
             true_states: np.ndarray | None = None) -> KLQEResult:
    """Run an already designed observer against an explicit measurement series."""
    samples = np.asarray(time, dtype=float).reshape(-1)
    observed = np.asarray(measurements, dtype=float)
    if observed.ndim == 1:
        observed = observed[:, None]
    if samples.size < 1 or observed.shape != (samples.size, estimator.C.shape[0]):
        raise ValueError("KLQE measurements must align with time and C.")
    if samples.size > 1 and np.any(np.diff(samples) <= 0.0):
        raise ValueError("KLQE run times must be strictly increasing.")
    truth = None if true_states is None else np.asarray(true_states, dtype=float)
    if truth is not None and (truth.ndim != 2 or truth.shape[0] != samples.size):
        raise ValueError("true_states must have one row per KLQE sample.")
    lifted = []; reconstructed = []; predicted = []; innovations = []
    for current, measurement in zip(samples, observed):
        reconstructed.append(estimator.step_measurement(float(current), measurement))
        lifted.append(estimator.lifted_state); predicted.append(estimator.predicted_measurement)
        innovations.append(measurement - estimator.predicted_measurement)
    state_values = np.vstack(reconstructed); innovation_values = np.vstack(innovations)
    true_series = None if truth is None else TimeSeries(samples, truth)
    error_series = None if truth is None else TimeSeries(samples, truth - state_values)
    metadata = {
        "estimator": estimator.summary(),
        "state_rmse": None if truth is None else np.sqrt(np.mean((truth - state_values) ** 2, axis=0)).tolist(),
        "innovation_rmse": np.sqrt(np.mean(innovation_values ** 2, axis=0)).tolist(),
    }
    return KLQEResult(
        TimeSeries(samples, np.vstack(lifted)), TimeSeries(samples, state_values),
        TimeSeries(samples, observed), TimeSeries(samples, np.vstack(predicted)),
        TimeSeries(samples, innovation_values), metadata, true_series, error_series,
    )










__all__ = ["KLQEDesign", "KLQEEstimator", "KLQEObservationModel", "KLQEResult",
           "design_klqe", "design_klqe_from_krom", "load_klqe_design", "load_klqe_result",
           "run_klqe", "save_klqe_design", "save_klqe_result"]
