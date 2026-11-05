"""RPE preset adapter for generic KROM training-domain expansion."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np

from control import RPE_CONTROL_TARGETS, rpe_control_target

from .krom_training import (KROMTrainingPlan, StateEnvelopeRequirement, _digest,
                            design_krom_training_plan, target_state_requirement)


DEFAULT_RPE_KROM = Path(
    "archives/derived/koopman-krom-extract/"
    "rpe-31-free-natural-period-dtau0p001-tau25-c000-t030__20260821T203315.886026Z__"
    "koopman-krom-extract__20260826T202020.762433Z__bdca807e.npz")

RPEKROMTrainingPlan = KROMTrainingPlan


def design_rpe_krom_training_plan(*, source_krom: str | Path = DEFAULT_RPE_KROM,
                                  target_names: Iterable[str] | None = None,
                                  points_per_axis: tuple[int, int] = (9, 11),
                                  margin_fraction: float = .1,
                                  tau_end: float = 25., dtau: float = .001,
                                  target_probes_per_target: int = 33) -> KROMTrainingPlan:
    """Build the canonical RPE preset of the generic expansion workflow."""
    names = tuple(target_names or sorted(RPE_CONTROL_TARGETS))
    if not names:
        raise ValueError("At least one canonical RPE target is required.")
    if target_probes_per_target < 5:
        raise ValueError("target_probes_per_target must be at least five.")
    sample_parameters = np.linspace(0., 2. * np.pi, target_probes_per_target, endpoint=False)
    requirements: tuple[StateEnvelopeRequirement, ...] = tuple(
        target_state_requirement(name, rpe_control_target(name),
            sample_parameters=sample_parameters, state_dimension=2,
            coordinate_system="dimensionless") for name in names)
    plan = design_krom_training_plan(source_krom=source_krom, requirements=requirements,
        points_per_axis=points_per_axis, margin_fraction=margin_fraction,
        model_time_end=tau_end, model_time_step=dtau, preset="air_water_20c")
    if plan.coverage.planned.minimum[0] <= 0.:
        raise ValueError("Planned RPE radius crosses the singular r=0 boundary.")
    return plan


_SCREENING_EXPORTS = {"RPEInitialConditionScreen", "RPEScreeningPolicy",
    "RPETargetProbeScreen", "RPETrainingScreeningReport",
    "publish_screened_rpe_krom_campaign", "save_rpe_krom_screening_report",
    "screen_rpe_krom_training_plan"}


def __getattr__(name: str):
    """Preserve former RPE screening imports through the preset adapter."""
    if name in _SCREENING_EXPORTS:
        from importlib import import_module
        return getattr(import_module("operations.rpe_krom_screening"), name)
    raise AttributeError(name)


__all__ = ["DEFAULT_RPE_KROM", "RPEKROMTrainingPlan", "design_rpe_krom_training_plan",
    *_SCREENING_EXPORTS]
