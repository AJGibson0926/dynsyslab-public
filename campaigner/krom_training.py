"""System-agnostic state-domain expansion plans for KROM training data."""
from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable
import json

import numpy as np

from provenance import serializable
from rom import AxisAlignedStateDomain, load_krom_archive, plan_state_domain_expansion


def _digest(payload: dict[str, Any]) -> str:
    return sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode()).hexdigest()


@dataclass(frozen=True)
class StateProbe:
    """One named state that must survive screening before publication."""

    requirement: str
    parameter: float
    state: tuple[float, ...]
    parameter_name: str = "time"

    def summary(self) -> dict[str, Any]:
        return {"requirement": self.requirement, "parameter": self.parameter,
            "parameter_name": self.parameter_name, "state": list(self.state)}


@dataclass(frozen=True)
class StateEnvelopeRequirement:
    """A named state-space envelope with optional anchors and probe states."""

    name: str
    domain: AxisAlignedStateDomain
    anchors: tuple[tuple[float, ...], ...] = ()
    probes: tuple[StateProbe, ...] = ()

    def summary(self) -> dict[str, Any]:
        return {"name": self.name, "domain": self.domain.summary(),
            "anchors": [list(item) for item in self.anchors],
            "probes": [item.summary() for item in self.probes]}


def target_state_requirement(name: str, target, *, sample_parameters: Iterable[float],
                             state_dimension: int, coordinate_system: str = "model",
                             include_initial_anchor: bool = True) -> StateEnvelopeRequirement:
    """Convert any evaluable state target into a generic coverage requirement."""
    parameters = tuple(float(value) for value in sample_parameters)
    if not parameters:
        raise ValueError("Target requirements need at least one sample parameter.")
    states = np.vstack([np.asarray(target.evaluate(value, state_dimension), dtype=float)
        for value in parameters])
    if states.shape != (len(parameters), state_dimension) or np.any(~np.isfinite(states)):
        raise ValueError("Target samples must be finite full-state vectors.")
    envelope = getattr(target, "conservative_state_envelope", None)
    if callable(envelope):
        minimum, maximum = envelope()
    else:
        minimum, maximum = np.min(states, axis=0), np.max(states, axis=0)
    domain = AxisAlignedStateDomain(tuple(minimum), tuple(maximum), coordinate_system)
    probes = tuple(StateProbe(name, parameter, tuple(state))
        for parameter, state in zip(parameters, states))
    anchors = (tuple(states[0]),) if include_initial_anchor else ()
    return StateEnvelopeRequirement(name, domain, anchors, probes)


def _constructor_parameters(system, preset: str | None,
                            overrides: dict[str, Any] | None) -> dict[str, Any]:
    if overrides is not None:
        return dict(overrides)
    if not is_dataclass(system):
        if preset is None:
            raise ValueError("Non-dataclass systems require explicit campaign parameters or a preset.")
        return {}
    excluded = {"driver", "name", "initial_state_override"}
    if preset is not None:
        excluded.add("preset")
    selected_fields = [field for field in fields(system)
        if field.init and field.name not in excluded]
    nested = [field.name for field in selected_fields
        if is_dataclass(getattr(system, field.name))]
    if nested:
        raise ValueError("Supply a registered preset or explicit constructor parameters for "
            f"nested system configuration fields: {nested}.")
    return {field.name: serializable(getattr(system, field.name)) for field in selected_fields}


def _union_requirements(requirements: tuple[StateEnvelopeRequirement, ...],
                        explicit: AxisAlignedStateDomain | None) -> AxisAlignedStateDomain:
    domains = [item.domain for item in requirements]
    if explicit is not None:
        domains.append(explicit)
    if not domains:
        raise ValueError("Supply a required state domain or at least one state-envelope requirement.")
    dimension, coordinates = domains[0].dimension, domains[0].coordinate_system
    if any(item.dimension != dimension or item.coordinate_system != coordinates for item in domains):
        raise ValueError("All requirements must share state dimension and coordinate system.")
    return AxisAlignedStateDomain(tuple(np.min([item.minimum for item in domains], axis=0)),
        tuple(np.max([item.maximum for item in domains], axis=0)), coordinates)


@dataclass(frozen=True)
class KROMTrainingPlan:
    """A non-executing, fingerprinted plan for expanding one KROM's domain."""

    source_krom: str
    system_name: str
    preset: str | None
    system_parameters: dict[str, Any]
    coverage: Any
    model_time_end: float
    model_time_step: float
    physical_time_scale: float
    requirements: tuple[StateEnvelopeRequirement, ...]
    anchors: tuple[tuple[float, ...], ...]
    probes: tuple[StateProbe, ...]
    format: str = "dynsyslab.krom-training-plan/v1"

    @property
    def model_initial_conditions(self) -> tuple[tuple[float, ...], ...]:
        return self.coverage.initial_conditions

    @property
    def coordinate_system(self) -> str:
        return self.coverage.current.coordinate_system

    @property
    def fingerprint(self) -> str:
        return _digest(self.summary(include_fingerprint=False))

    @property
    def requirement_names(self) -> tuple[str, ...]:
        return tuple(item.name for item in self.requirements)

    # Compatibility names for the former RPE-specific contract.
    @property
    def target_names(self) -> tuple[str, ...]:
        return self.requirement_names

    @property
    def target_anchors(self) -> tuple[tuple[float, ...], ...]:
        return self.anchors

    @property
    def target_probes(self) -> tuple[tuple[str, float, tuple[float, ...]], ...]:
        return tuple((item.requirement, item.parameter, item.state) for item in self.probes)

    @property
    def tau_end(self) -> float:
        return self.model_time_end

    @property
    def dtau(self) -> float:
        return self.model_time_step

    @property
    def time_scale_s(self) -> float:
        return self.physical_time_scale

    def summary(self, *, include_fingerprint: bool = True) -> dict[str, Any]:
        payload = {"format": self.format, "source_krom": self.source_krom,
            "system": {"name": self.system_name, "preset": self.preset,
                "parameters": self.system_parameters},
            "coordinate_system": self.coordinate_system,
            "requirements": [item.summary() for item in self.requirements],
            "coverage": self.coverage.summary(), "anchors": [list(item) for item in self.anchors],
            "probes": [item.summary() for item in self.probes],
            "integration": {"model_time_span": [0., self.model_time_end],
                "model_time_step": self.model_time_step,
                "physical_time_span": [0., self.model_time_end * self.physical_time_scale],
                "physical_time_step": self.model_time_step * self.physical_time_scale,
                "physical_time_scale": self.physical_time_scale,
                "use_dimensionless_form": self.coordinate_system == "dimensionless"},
            "risk_gate": {"execution_automatic": False, "requires_screening": True,
                "reason": "New state-domain initial conditions require bounded in-memory screening."}}
        if include_fingerprint:
            payload["plan_fingerprint_sha256"] = _digest(payload)
        return payload


def design_krom_training_plan(*, source_krom: str | Path,
                              required_domain: AxisAlignedStateDomain | None = None,
                              requirements: Iterable[StateEnvelopeRequirement] = (),
                              points_per_axis: Iterable[int], margin_fraction: float = .1,
                              model_time_end: float = 25., model_time_step: float = .001,
                              anchors: Iterable[Iterable[float]] = (),
                              preset: str | None = None,
                              system_parameters: dict[str, Any] | None = None) -> KROMTrainingPlan:
    """Plan domain-expanding free trajectories for any compatible ODE KROM."""
    if model_time_end <= 0. or model_time_step <= 0. or model_time_step > model_time_end:
        raise ValueError("Model time span and step must be positive and ordered.")
    source = Path(source_krom).resolve()
    krom = load_krom_archive(source)
    current = AxisAlignedStateDomain.from_krom(krom)
    selected = tuple(requirements)
    required = _union_requirements(selected, required_domain)
    if required.dimension != current.dimension or required.coordinate_system != current.coordinate_system:
        raise ValueError("Required and KROM training domains must use identical model coordinates.")
    requirement_anchors = tuple(anchor for item in selected for anchor in item.anchors)
    all_anchors = tuple(tuple(float(value) for value in anchor)
        for anchor in (*tuple(anchors), *requirement_anchors))
    probes = tuple(probe for item in selected for probe in item.probes)
    coverage = plan_state_domain_expansion(current, required, points_per_axis=points_per_axis,
        margin_fraction=margin_fraction, outside_current_only=True, anchors=all_anchors)
    system = krom.system
    system_name = str(getattr(system, "name", ""))
    if not system_name:
        raise ValueError("The source KROM system has no registered name.")
    parameters = _constructor_parameters(system, preset, system_parameters)
    from dynsyslib import ZeroDriver, load_system
    campaign_system = load_system(system_name, driver=ZeroDriver(), preset=preset, **parameters)
    if int(np.asarray(campaign_system.initial_state).size) != current.dimension:
        raise ValueError("Campaign system configuration does not match the source KROM dimension.")
    time_scale = 1.
    if current.coordinate_system == "dimensionless":
        converter = getattr(system, "to_dimensional_time", None)
        if not callable(converter):
            raise ValueError("Dimensionless KROM training requires a physical-time conversion.")
        time_scale = float(np.asarray(converter(1.)) - np.asarray(converter(0.)))
    if not np.isfinite(time_scale) or time_scale <= 0.:
        raise ValueError("KROM model-to-physical time scale must be finite and positive.")
    return KROMTrainingPlan(str(source), system_name, preset, parameters, coverage,
        float(model_time_end), float(model_time_step), time_scale, selected,
        all_anchors, probes)


__all__ = ["KROMTrainingPlan", "StateEnvelopeRequirement", "StateProbe",
    "design_krom_training_plan", "target_state_requirement"]
