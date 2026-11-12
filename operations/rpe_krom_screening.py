"""Compatibility preset for RPE-specific KROM training safety limits."""
from __future__ import annotations

from dataclasses import dataclass

from .krom_training_screening import (InitialConditionScreen, KROMTrainingScreeningPolicy,
    KROMTrainingScreeningReport, ProbeScreen, publish_screened_krom_training_campaign,
    save_krom_training_screening_report, screen_krom_training_plan)


@dataclass(frozen=True)
class RPEScreeningPolicy:
    """Legacy names for the canonical two-state RPE screening preset."""

    tau_end: float = .25
    maximum_step_dtau: float = .002
    minimum_radius: float = .05
    maximum_radius: float = 15.
    maximum_absolute_velocity: float = 250.
    method: str = "DOP853"
    rtol: float = 1e-9
    atol: float = 1e-12

    def generic(self) -> KROMTrainingScreeningPolicy:
        unbounded = 1e300
        return KROMTrainingScreeningPolicy(model_time_end=self.tau_end,
            maximum_model_time_step=self.maximum_step_dtau,
            state_minimum=(self.minimum_radius, -unbounded),
            state_maximum=(self.maximum_radius, unbounded),
            maximum_absolute_state=(unbounded, self.maximum_absolute_velocity),
            method=self.method, rtol=self.rtol, atol=self.atol)


RPEInitialConditionScreen = InitialConditionScreen
RPETargetProbeScreen = ProbeScreen
RPETrainingScreeningReport = KROMTrainingScreeningReport


def screen_rpe_krom_training_plan(plan, *, policy: RPEScreeningPolicy | None = None):
    """Apply the RPE safety preset through the generic numerical screener."""
    return screen_krom_training_plan(plan, policy=(policy or RPEScreeningPolicy()).generic())


def publish_screened_rpe_krom_campaign(plan, report, *,
                                       label: str = "rpe-krom-expanded-control-domain"):
    """Publish through the generic campaign operation using the RPE label."""
    return publish_screened_krom_training_campaign(plan, report, label=label)


def save_rpe_krom_screening_report(plan, report, *, root=None):
    """Persist the generic report through the former RPE API name."""
    return save_krom_training_screening_report(plan, report, root=root)


__all__ = ["RPEInitialConditionScreen", "RPEScreeningPolicy", "RPETargetProbeScreen",
    "RPETrainingScreeningReport", "publish_screened_rpe_krom_campaign",
    "save_rpe_krom_screening_report", "screen_rpe_krom_training_plan"]
