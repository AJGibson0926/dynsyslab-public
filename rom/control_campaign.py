"""Reusable, bounded campaigns for testing KROM control representations."""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence
import numpy as np


MAX_KROM_CONTROL_ATTEMPTS = 25
STANDARD_CONTROL_CASES = (
    "exact", "constant", "bilinear-dynamics", "bilinear-map",
)


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, complex):
        return {"real": float(value.real), "imaginary": float(value.imag)}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)

@dataclass(frozen=True)
class KROMControlCampaignPolicy:
    """Hard bounds and representation-independent control success criteria."""

    max_attempts: int = MAX_KROM_CONTROL_ATTEMPTS
    release_error_tolerance: float = 0.02
    final_equilibrium_error_tolerance: float = 0.10
    required_cases: tuple[str, ...] = STANDARD_CONTROL_CASES
    stop_on_bilinear_success: bool = True

    def __post_init__(self):
        attempts = int(self.max_attempts)
        if attempts < 1 or attempts > MAX_KROM_CONTROL_ATTEMPTS:
            raise ValueError(
                f"max_attempts must be in [1, {MAX_KROM_CONTROL_ATTEMPTS}]."
            )
        if (
            self.release_error_tolerance <= 0.0
            or not np.isfinite(self.release_error_tolerance)
            or self.final_equilibrium_error_tolerance <= 0.0
            or not np.isfinite(self.final_equilibrium_error_tolerance)
        ):
            raise ValueError("Control success tolerances must be finite and positive.")
        cases = tuple(str(item) for item in self.required_cases)
        if not cases or len(set(cases)) != len(cases):
            raise ValueError("required_cases must be nonempty and unique.")
        object.__setattr__(self, "max_attempts", attempts)
        object.__setattr__(self, "required_cases", cases)

    def case_succeeded(self, record: Mapping[str, Any]) -> bool:
        if record.get("status") != "complete":
            return False
        metrics = record.get("metrics") or {}
        release = float(metrics.get("distance_to_target_at_release", np.inf))
        final = float(metrics.get("final_distance_to_left_equilibrium", np.inf))
        return bool(
            np.isfinite(release)
            and np.isfinite(final)
            and release <= self.release_error_tolerance
            and final <= self.final_equilibrium_error_tolerance
        )

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "krom-control-campaign-policy/v1",
            "max_attempts": self.max_attempts,
            "absolute_hard_limit": MAX_KROM_CONTROL_ATTEMPTS,
            "release_error_tolerance": self.release_error_tolerance,
            "final_equilibrium_error_tolerance": self.final_equilibrium_error_tolerance,
            "required_cases": list(self.required_cases),
            "stop_on_bilinear_success": self.stop_on_bilinear_success,
        }

@dataclass(frozen=True)
class KROMCandidate:
    """One pre-screened set of nonredundant parent-KROM coordinates."""

    retained_indices: tuple[int, ...]
    scores: Mapping[str, float]
    rationale: str = ""

    def __post_init__(self):
        indices = tuple(int(item) for item in self.retained_indices)
        if not indices or len(set(indices)) != len(indices) or any(item < 0 for item in indices):
            raise ValueError("A KROM candidate needs unique nonnegative retained indices.")
        normalized = {str(key): float(value) for key, value in self.scores.items()}
        if any(not np.isfinite(value) for value in normalized.values()):
            raise ValueError("Candidate scores must be finite.")
        object.__setattr__(self, "retained_indices", indices)
        object.__setattr__(self, "scores", normalized)

    def summary(self) -> dict[str, Any]:
        return {
            "retained_indices": list(self.retained_indices),
            "mode_count": len(self.retained_indices),
            "scores": dict(self.scores),
            "rationale": self.rationale,
        }

@dataclass(frozen=True)
class KROMControlAttempt:
    """Durable result of one extraction-regression-control sequence."""

    attempt_number: int
    candidate: KROMCandidate
    cases: Mapping[str, Mapping[str, Any]]
    artifacts: Mapping[str, Any]
    diagnostics: Mapping[str, Any]

    def validate(self, policy: KROMControlCampaignPolicy) -> None:
        missing = set(policy.required_cases) - set(self.cases)
        if missing:
            raise ValueError(f"Campaign attempt omitted required control cases: {sorted(missing)}")

    def successful_cases(self, policy: KROMControlCampaignPolicy) -> list[str]:
        return [name for name, record in self.cases.items() if policy.case_succeeded(record)]

    def bilinear_succeeded(self, policy: KROMControlCampaignPolicy) -> bool:
        return any(
            name.startswith("bilinear") and policy.case_succeeded(record)
            for name, record in self.cases.items()
        )

    def summary(self, policy: KROMControlCampaignPolicy) -> dict[str, Any]:
        return {
            "attempt_number": int(self.attempt_number),
            "candidate": self.candidate.summary(),
            "cases": _jsonable(self.cases),
            "successful_cases": self.successful_cases(policy),
            "bilinear_success": self.bilinear_succeeded(policy),
            "artifacts": _jsonable(self.artifacts),
            "diagnostics": _jsonable(self.diagnostics),
        }

@dataclass(frozen=True)
class KROMControlCampaignResult:
    policy: KROMControlCampaignPolicy
    candidates_considered: int
    attempts: tuple[KROMControlAttempt, ...]
    stop_reason: str

    @property
    def successful_attempt(self) -> KROMControlAttempt | None:
        return next(
            (item for item in self.attempts if item.bilinear_succeeded(self.policy)),
            None,
        )

    def summary(self) -> dict[str, Any]:
        winner = self.successful_attempt
        return {
            "encoding": "krom-control-campaign-result/v1",
            "policy": self.policy.summary(),
            "candidates_considered": int(self.candidates_considered),
            "attempt_count": len(self.attempts),
            "stop_reason": self.stop_reason,
            "successful_attempt_number": None if winner is None else winner.attempt_number,
            "attempts": [item.summary(self.policy) for item in self.attempts],
        }

def _write_ledger(path: str | Path | None, result: KROMControlCampaignResult) -> None:
    if path is None:
        return
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(json.dumps(result.summary(), indent=2) + "\n", encoding="utf-8")
    temporary.replace(destination)

def run_krom_control_campaign(
    candidates: Sequence[KROMCandidate],
    evaluator: Callable[[int, KROMCandidate], KROMControlAttempt],
    *,
    policy: KROMControlCampaignPolicy | None = None,
    ledger_path: str | Path | None = None,
    prior_attempts: Sequence[KROMControlAttempt] = (),
) -> KROMControlCampaignResult:
    """Evaluate candidates sequentially with a hard, non-bypassable limit.

    The callback performs the domain-specific extraction, all regressions, all
    required controller solves, and artifact retention for exactly one attempt.
    A validated partial ledger is atomically replaced after every attempt.
    """
    active = policy or KROMControlCampaignPolicy()
    attempts: list[KROMControlAttempt] = list(prior_attempts)
    if len(attempts) > active.max_attempts:
        raise ValueError("Prior attempts exceed the campaign policy limit.")
    for number, attempt in enumerate(attempts, start=1):
        if attempt.attempt_number != number:
            raise ValueError("Prior campaign attempts must be contiguous and one-indexed.")
        attempt.validate(active)
    if any(item.bilinear_succeeded(active) for item in attempts):
        result = KROMControlCampaignResult(
            active, len(candidates), tuple(attempts), "bilinear-success"
        )
        _write_ledger(ledger_path, result)
        return result
    stop_reason = "candidate-list-exhausted"
    remaining = candidates[len(attempts) : active.max_attempts]
    for attempt_number, candidate in enumerate(remaining, start=len(attempts) + 1):
        attempt = evaluator(attempt_number, candidate)
        if attempt.attempt_number != attempt_number or attempt.candidate != candidate:
            raise ValueError("Campaign evaluator returned a misaligned attempt record.")
        attempt.validate(active)
        attempts.append(attempt)
        partial = KROMControlCampaignResult(
            active, len(candidates), tuple(attempts), "in-progress"
        )
        _write_ledger(ledger_path, partial)
        if active.stop_on_bilinear_success and attempt.bilinear_succeeded(active):
            stop_reason = "bilinear-success"
            break
    else:
        if len(candidates) > active.max_attempts:
            stop_reason = "hard-attempt-limit-reached"
    result = KROMControlCampaignResult(
        active, len(candidates), tuple(attempts), stop_reason
    )
    _write_ledger(ledger_path, result)
    return result


from .control_candidate_screening import (
    polynomial_bilinear_closure_quality, rank_krom_candidates,
    resample_trajectory_by_arclength, trajectory_tube_states,
)

__all__ = [
    "KROMCandidate", "KROMControlAttempt", "KROMControlCampaignPolicy",
    "KROMControlCampaignResult", "MAX_KROM_CONTROL_ATTEMPTS",
    "STANDARD_CONTROL_CASES", "rank_krom_candidates",
    "polynomial_bilinear_closure_quality", "run_krom_control_campaign",
    "resample_trajectory_by_arclength", "trajectory_tube_states",
]
