"""Numerical safety screening and publication for generic KROM training plans."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Any
import json
import os

import numpy as np

from campaigner.krom_training import KROMTrainingPlan, StateProbe, _digest
from rom import AxisAlignedStateDomain, load_krom_archive


def _vector(value, dimension: int, label: str) -> np.ndarray | None:
    if value is None:
        return None
    result = np.asarray(value, dtype=float).reshape(-1)
    if result.size == 1:
        result = np.full(dimension, result.item())
    if result.size != dimension or np.any(np.isnan(result)):
        raise ValueError(f"{label} must be non-NaN and scalar or state-dimensional.")
    return result


@dataclass(frozen=True)
class KROMTrainingScreeningPolicy:
    """Integrator settings and optional component-wise safety limits."""

    model_time_end: float = .25
    maximum_model_time_step: float = .002
    state_minimum: tuple[float, ...] | None = None
    state_maximum: tuple[float, ...] | None = None
    maximum_absolute_state: float | tuple[float, ...] | None = None
    method: str = "DOP853"
    rtol: float = 1e-9
    atol: float = 1e-12

    def __post_init__(self):
        numeric = (self.model_time_end, self.maximum_model_time_step, self.rtol, self.atol)
        if any(not np.isfinite(value) or value <= 0. for value in numeric):
            raise ValueError("Screening time and tolerance values must be finite and positive.")
        if self.maximum_model_time_step > self.model_time_end:
            raise ValueError("maximum_model_time_step cannot exceed model_time_end.")

    def limits(self, dimension: int) -> tuple[np.ndarray | None, np.ndarray | None, np.ndarray | None]:
        """Return normalized component-wise lower, upper, and magnitude limits."""
        low = _vector(self.state_minimum, dimension, "state_minimum")
        high = _vector(self.state_maximum, dimension, "state_maximum")
        absolute = _vector(self.maximum_absolute_state, dimension, "maximum_absolute_state")
        if low is not None and high is not None and np.any(high <= low):
            raise ValueError("Every state_maximum must exceed state_minimum.")
        if absolute is not None and np.any(absolute <= 0.):
            raise ValueError("maximum_absolute_state must be positive.")
        return low, high, absolute


@dataclass(frozen=True)
class InitialConditionScreen:
    index: int
    initial_state: tuple[float, ...]
    accepted: bool
    reason: str
    solver_status: str
    sample_count: int
    minimum_state: tuple[float, ...] | None
    maximum_state: tuple[float, ...] | None
    maximum_absolute_state: tuple[float, ...] | None
    final_state: tuple[float, ...] | None
    elapsed_seconds: float

    @property
    def minimum_radius(self) -> float | None:
        return None if self.minimum_state is None else self.minimum_state[0]

    @property
    def maximum_radius(self) -> float | None:
        return None if self.maximum_state is None else self.maximum_state[0]

    @property
    def maximum_absolute_velocity(self) -> float | None:
        return None if self.maximum_absolute_state is None or len(self.maximum_absolute_state) < 2 \
            else self.maximum_absolute_state[1]


@dataclass(frozen=True)
class ProbeScreen:
    probe: StateProbe
    accepted: bool
    reason: str
    minimum_state: tuple[float, ...] | None
    maximum_state: tuple[float, ...] | None
    maximum_absolute_state: tuple[float, ...] | None
    elapsed_seconds: float

    @property
    def target_name(self) -> str:
        return self.probe.requirement

    @property
    def tau(self) -> float:
        return self.probe.parameter

    @property
    def initial_state(self) -> tuple[float, ...]:
        return self.probe.state

    @property
    def minimum_radius(self) -> float | None:
        return None if self.minimum_state is None else self.minimum_state[0]

    @property
    def maximum_radius(self) -> float | None:
        return None if self.maximum_state is None else self.maximum_state[0]

    @property
    def maximum_absolute_velocity(self) -> float | None:
        return None if self.maximum_absolute_state is None or len(self.maximum_absolute_state) < 2 \
            else self.maximum_absolute_state[1]


@dataclass(frozen=True)
class KROMTrainingScreeningReport:
    """Immutable result of screening one fingerprinted generic plan."""

    plan_fingerprint_sha256: str
    policy: KROMTrainingScreeningPolicy
    results: tuple[InitialConditionScreen, ...]
    probe_results: tuple[ProbeScreen, ...]
    publication_ready: bool
    coverage_with_accepted: dict[str, Any]
    elapsed_seconds: float
    format: str = "dynsyslab.krom-training-screen/v1"

    @property
    def accepted_indices(self) -> tuple[int, ...]:
        return tuple(item.index for item in self.results if item.accepted)

    @property
    def rejected_indices(self) -> tuple[int, ...]:
        return tuple(item.index for item in self.results if not item.accepted)

    @property
    def rejected_probe_count(self) -> int:
        return sum(not item.accepted for item in self.probe_results)

    @property
    def target_probe_results(self) -> tuple[ProbeScreen, ...]:
        return self.probe_results

    @property
    def rejected_target_probe_count(self) -> int:
        return self.rejected_probe_count

    @property
    def requirement_summary(self) -> dict[str, Any]:
        output = {}
        for name in sorted({item.probe.requirement for item in self.probe_results}):
            selected = [item for item in self.probe_results if item.probe.requirement == name]
            reasons = {reason: sum(item.reason == reason for item in selected)
                for reason in sorted({item.reason for item in selected})}
            output[name] = {"probe_count": len(selected),
                "accepted_count": sum(item.accepted for item in selected),
                "rejected_count": sum(not item.accepted for item in selected), "reasons": reasons}
        return output

    @property
    def target_summary(self) -> dict[str, Any]:
        return self.requirement_summary

    @property
    def fingerprint(self) -> str:
        return _digest(self.summary(include_fingerprint=False, include_timing=False))

    def summary(self, *, include_results: bool = True, include_fingerprint: bool = True,
                include_timing: bool = True) -> dict[str, Any]:
        payload = {"format": self.format, "plan_fingerprint_sha256": self.plan_fingerprint_sha256,
            "policy": asdict(self.policy), "publication_ready": self.publication_ready,
            "accepted_count": len(self.accepted_indices), "rejected_count": len(self.rejected_indices),
            "accepted_indices": list(self.accepted_indices), "rejected_indices": list(self.rejected_indices),
            "probe_count": len(self.probe_results), "rejected_probe_count": self.rejected_probe_count,
            "requirement_summary": self.requirement_summary,
            "coverage_with_accepted": dict(self.coverage_with_accepted)}
        if include_timing:
            payload["elapsed_seconds"] = self.elapsed_seconds
        if include_results:
            payload["results"] = [asdict(item) for item in self.results]
            payload["probe_results"] = [asdict(item) for item in self.probe_results]
            if not include_timing:
                for item in (*payload["results"], *payload["probe_results"]):
                    item.pop("elapsed_seconds", None)
        if include_fingerprint:
            payload["screen_fingerprint_sha256"] = self.fingerprint
        return payload


def _physical_state(system, state: np.ndarray, coordinates: str) -> np.ndarray:
    if coordinates == "dimensionless":
        transform = getattr(system, "to_dimensional_state", None)
        if not callable(transform):
            raise ValueError("Dimensionless screening requires a physical-state conversion.")
        return np.asarray(transform(state), dtype=float)
    return np.asarray(state, dtype=float)


def screen_krom_training_plan(plan: KROMTrainingPlan, *,
                              policy: KROMTrainingScreeningPolicy | None = None
                              ) -> KROMTrainingScreeningReport:
    """Short-run all planned and probe states through the canonical ODE solver."""
    selected = policy or KROMTrainingScreeningPolicy()
    krom = load_krom_archive(plan.source_krom)
    system = krom.system
    dimension = int(np.asarray(system.initial_state).size)
    low_limit, high_limit, absolute_limit = selected.limits(dimension)
    from solver import SolveRequest, SolverOptions, solve_ode_in_memory
    started = perf_counter()

    def run(initial):
        model_initial = np.asarray(initial, dtype=float)
        physical = _physical_state(system, model_initial, plan.coordinate_system)
        request = SolveRequest(system=system, kind="ode",
            time_span=(0., selected.model_time_end * plan.physical_time_scale),
            initial_state=physical,
            options=SolverOptions(method=selected.method, rtol=selected.rtol, atol=selected.atol,
                max_step=selected.maximum_model_time_step * plan.physical_time_scale,
                use_dimensionless_form=plan.coordinate_system == "dimensionless",
                execution_backend="python"), label="krom-training-screen")
        item_started = perf_counter()
        result = solve_ode_in_memory(request)
        elapsed = perf_counter() - item_started
        accepted, reason, count = False, result.message, 0
        minimum = maximum = max_absolute = final = None
        if result.status == "complete" and result.time_series is not None:
            model = (np.vstack([system.to_dimensionless_state(row) for row in result.values])
                if plan.coordinate_system == "dimensionless" else np.asarray(result.values, dtype=float))
            count = int(model.shape[0])
            minimum, maximum = np.min(model, axis=0), np.max(model, axis=0)
            max_absolute, final = np.max(np.abs(model), axis=0), model[-1]
            if np.any(~np.isfinite(model)):
                reason = "nonfinite_state"
            elif low_limit is not None and np.any(minimum < low_limit):
                reason = "state_below_minimum"
            elif high_limit is not None and np.any(maximum > high_limit):
                reason = "state_above_maximum"
            elif absolute_limit is not None and np.any(max_absolute > absolute_limit):
                reason = "absolute_state_above_maximum"
            else:
                accepted, reason = True, "accepted"
        convert = lambda value: None if value is None else tuple(float(item) for item in value)
        return accepted, reason, result.status, count, convert(minimum), convert(maximum), \
            convert(max_absolute), convert(final), elapsed

    results = []
    for index, initial in enumerate(plan.model_initial_conditions):
        accepted, reason, status, count, minimum, maximum, absolute, final, elapsed = run(initial)
        results.append(InitialConditionScreen(index, tuple(initial), accepted, reason, status,
            count, minimum, maximum, absolute, final, elapsed))
    probe_results = []
    for probe in plan.probes:
        accepted, reason, _, _, minimum, maximum, absolute, _, elapsed = run(probe.state)
        probe_results.append(ProbeScreen(probe, accepted, reason, minimum, maximum, absolute, elapsed))
    accepted_states = np.asarray([item.initial_state for item in results if item.accepted], dtype=float)
    if accepted_states.size:
        union_min = np.minimum(np.asarray(plan.coverage.current.minimum), np.min(accepted_states, axis=0))
        union_max = np.maximum(np.asarray(plan.coverage.current.maximum), np.max(accepted_states, axis=0))
        combined = AxisAlignedStateDomain(tuple(union_min), tuple(union_max), plan.coordinate_system)
        required_covered = bool(np.all(combined.contains(
            (plan.coverage.required.minimum, plan.coverage.required.maximum))))
        anchors_accepted = all(any(np.allclose(item.initial_state, anchor, rtol=0., atol=1e-12)
            and item.accepted for item in results) for anchor in plan.anchors)
        coverage = {"combined_domain": combined.summary(), "required_covered": required_covered,
            "anchors_accepted": anchors_accepted}
    else:
        required_covered = anchors_accepted = False
        coverage = {"combined_domain": None, "required_covered": False, "anchors_accepted": False}
    probes_accepted = all(item.accepted for item in probe_results)
    coverage.update({"probe_states_accepted": probes_accepted,
        "rejected_probe_count": sum(not item.accepted for item in probe_results)})
    ready = bool(accepted_states.size and required_covered and anchors_accepted and probes_accepted)
    return KROMTrainingScreeningReport(plan.fingerprint, selected, tuple(results),
        tuple(probe_results), ready, coverage, perf_counter() - started)


def publish_screened_krom_training_campaign(plan: KROMTrainingPlan,
                                            report: KROMTrainingScreeningReport, *,
                                            label: str | None = None):
    """Publish a free-response campaign only after a matching successful screen."""
    if report.plan_fingerprint_sha256 != plan.fingerprint:
        raise ValueError("Screening report does not belong to this training plan.")
    if not report.publication_ready:
        raise ValueError("KROM training campaign requires a publication-ready screen.")
    krom = load_krom_archive(plan.source_krom)
    accepted = [plan.model_initial_conditions[index] for index in report.accepted_indices]
    physical = [_physical_state(krom.system, np.asarray(state), plan.coordinate_system).tolist()
        for state in accepted]
    from campaigner.designer import Campaigner
    selected_label = label or f"{plan.system_name}-krom-training-domain-expansion"
    campaign, collection = Campaigner().trajectory_ensemble(system_name=plan.system_name,
        preset=plan.preset, parameters=plan.system_parameters, initial_conditions=physical,
        integration={"time_span": [0., plan.model_time_end * plan.physical_time_scale],
            "delta_t": plan.model_time_step * plan.physical_time_scale,
            "solver_options": {"method": "DOP853", "rtol": 1e-9, "atol": 1e-12,
                "use_dimensionless_form": plan.coordinate_system == "dimensionless"}},
        label=selected_label, metadata={"purpose": "krom-training-domain-expansion",
            "source_krom": plan.source_krom, "training_plan": plan.summary(),
            "screening_report": report.summary(), "model_coordinate_initial_conditions": accepted,
            "requirement_names": list(plan.requirement_names), "execution_gate": "screen-passed"})
    return campaign, collection


def save_krom_training_screening_report(plan: KROMTrainingPlan,
                                        report: KROMTrainingScreeningReport, *,
                                        root: str | Path | None = None) -> Path:
    """Atomically persist one immutable, content-addressed screening record."""
    if report.plan_fingerprint_sha256 != plan.fingerprint:
        raise ValueError("Screening report does not belong to this training plan.")
    if root is None:
        from archiver.gateway import MANAGER_ARCHIVE
        directory = Path(MANAGER_ARCHIVE.archive_root) / "campaigns" / "screening"
    else:
        directory = Path(root)
    directory.mkdir(parents=True, exist_ok=True)
    payload = {"format": "dynsyslab.krom-training-screen-record/v1", "plan": plan.summary(),
        "screening": report.summary(include_timing=False)}
    encoded = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    target = directory / f"krom-training-screen__{report.fingerprint}.json"
    if target.exists():
        if target.read_text(encoding="utf-8") != encoded:
            raise RuntimeError("Content-addressed screening report collision.")
        return target.resolve()
    temporary = target.with_suffix(target.suffix + f".{os.getpid()}.partial")
    with temporary.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(encoded)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, target)
    return target.resolve()


__all__ = ["InitialConditionScreen", "KROMTrainingScreeningPolicy",
    "KROMTrainingScreeningReport", "ProbeScreen", "publish_screened_krom_training_campaign",
    "save_krom_training_screening_report", "screen_krom_training_plan"]
