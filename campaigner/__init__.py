"""Independent campaign-design package for collections of dynsyslab jobs."""
from .designer import BIFURCATION_PROFILES, CampaignManifest, CampaignRequest, Campaigner
from .krom_training import (KROMTrainingPlan, StateEnvelopeRequirement, StateProbe,
                            design_krom_training_plan, target_state_requirement)
from .rpe_krom import DEFAULT_RPE_KROM, RPEKROMTrainingPlan, design_rpe_krom_training_plan

_GENERIC_SCREENING_EXPORTS = {"KROMTrainingScreeningPolicy", "KROMTrainingScreeningReport",
    "publish_screened_krom_training_campaign", "save_krom_training_screening_report",
    "screen_krom_training_plan"}
_RPE_SCREENING_EXPORTS = {"RPEScreeningPolicy", "RPETrainingScreeningReport",
    "publish_screened_rpe_krom_campaign", "save_rpe_krom_screening_report",
    "screen_rpe_krom_training_plan"}


def __getattr__(name: str):
    """Load numerical screening only when an execution capability is requested."""
    if name in _GENERIC_SCREENING_EXPORTS:
        from importlib import import_module
        return getattr(import_module("operations.krom_training_screening"), name)
    if name in _RPE_SCREENING_EXPORTS:
        from importlib import import_module
        return getattr(import_module("operations.rpe_krom_screening"), name)
    raise AttributeError(name)

__all__ = ["Campaigner", "CampaignManifest", "CampaignRequest", "BIFURCATION_PROFILES",
    "KROMTrainingPlan", "StateEnvelopeRequirement", "StateProbe", "design_krom_training_plan",
    "target_state_requirement", "KROMTrainingScreeningPolicy", "KROMTrainingScreeningReport",
    "publish_screened_krom_training_campaign", "save_krom_training_screening_report",
    "screen_krom_training_plan",
    "DEFAULT_RPE_KROM","RPEKROMTrainingPlan","RPEScreeningPolicy","RPETrainingScreeningReport",
    "design_rpe_krom_training_plan","publish_screened_rpe_krom_campaign","save_rpe_krom_screening_report",
    "screen_rpe_krom_training_plan"]
