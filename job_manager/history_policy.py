"""Campaign-declared history/checkpoint policy resolution for JobManager."""
from __future__ import annotations
from typing import Any

def resolve(policy:dict[str,Any]|None, *, resources:dict[str,Any], options:dict[str,Any], history:dict[str,Any])->tuple[dict[str,Any],dict[str,Any],dict[str,Any]]:
    """Resolve declarative campaign policy into ordinary solver dictionaries.

    Solver never decides whether a trajectory is valuable enough to checkpoint;
    this function makes that decision from manifest metadata before execution.
    """
    policy=dict(policy or {}); resolved={"enabled":False,"reason":"no_policy"}
    if not policy or not policy.get("enabled",False): return options,history,resolved
    estimate=policy.get("estimated_seconds",resources.get("wall_time_s")); threshold=policy.get("if_estimated_seconds_over")
    unknown_policy=str(policy.get("when_estimate_unknown","disable"))
    if unknown_policy not in {"enable","disable"}:
        raise ValueError("when_estimate_unknown must be 'enable' or 'disable'.")
    if threshold is not None and estimate is None and unknown_policy=="disable":
        return options,history,{"enabled":False,"reason":"unknown_estimate","estimated_seconds":estimate,"threshold_seconds":threshold}
    if threshold is not None and estimate is not None and float(estimate)<float(threshold):
        return options,history,{"enabled":False,"reason":"below_estimate_threshold","estimated_seconds":estimate,"threshold_seconds":threshold}
    mode=str(policy.get("mode","checkpoint")); strategy=str(policy.get("strategy","low_level_stepper"))
    if mode not in {"telemetry","sampled","checkpoint"}: raise ValueError(f"Unknown history policy mode {mode!r}.")
    if strategy not in {"segmented","low_level_stepper"}: raise ValueError(f"Unknown checkpoint strategy {strategy!r}.")
    configured={**history,"mode":mode,"track_during_solve":mode!="telemetry","checkpoint_strategy":strategy}
    # Accept the concise names used in campaign planning reports while keeping
    # the solver-facing names identical to HistoryOptions.
    aliases={
        "interval_seconds": "checkpoint_interval",
        "sample_interval_seconds": "sample_interval",
    }
    for source, target in aliases.items():
        if target not in policy and source in policy:
            configured[target]=policy[source]
    for key in ("checkpoint_interval","sample_interval","pause_at_checkpoint","telemetry_interval_seconds","telemetry_enabled","telemetry_heartbeat"):
        if key in policy: configured[key]=policy[key]
    configured_options=dict(options)
    if strategy=="low_level_stepper" and mode in {"sampled","checkpoint"}: configured_options["adaptive_engine"]="low_level_stepper"
    return configured_options,configured,{"enabled":True,"reason":"policy_matched","estimated_seconds":estimate,"estimate_was_unknown":estimate is None,"mode":mode,"strategy":strategy}
