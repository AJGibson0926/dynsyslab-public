"""Top-level persistence and execution of the bounded research campaign."""
from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter
from typing import Any
import numpy as np
from scipy.io import loadmat

from postprocessing.koopman.workflows import (DuffingControlledTrainingSOP,
    generate_duffing_controlled_training)
from rom import (KROMCandidate, KROMControlAttempt, KROMControlCampaignPolicy,
                 load_krom_archive, run_krom_control_campaign)
from .common import (ACTIVE_WINDOW, FINAL_TIME, INITIAL_STATE, LEFT_EQUILIBRIUM,
                     TARGET_STATE, _sha, _write_json)
from .evaluator import DuffingKROMSearch
from .schedule import _candidate_schedule

def _rehydrate_attempts(
    ledger_path: Path, candidates: list[KROMCandidate], policy: KROMControlCampaignPolicy
) -> list[KROMControlAttempt]:
    if not ledger_path.is_file():
        return []
    stored = json.loads(ledger_path.read_text(encoding="utf-8"))
    prior: list[KROMControlAttempt] = []
    for item in stored.get("attempts", []):
        number = int(item["attempt_number"])
        if number > len(candidates):
            raise ValueError("Stored campaign attempt exceeds current candidate schedule.")
        candidate = candidates[number - 1]
        if list(candidate.retained_indices) != item["candidate"]["retained_indices"]:
            raise ValueError("Stored campaign candidate schedule does not match this run.")
        attempt = KROMControlAttempt(
            number,
            candidate,
            item["cases"],
            item.get("artifacts") or {},
            item.get("diagnostics") or {},
        )
        attempt.validate(policy)
        prior.append(attempt)
    return prior

def run(
    manifest_path: Path,
    control_source: Path,
    output_root: Path,
    *,
    max_attempts: int,
    maximum_mode_count: int,
    output_dt: float,
    grid_size: int,
    case_timeout_seconds: float,
    dry_run: bool = False,
    resume: bool = False,
) -> dict[str, Any]:
    started = perf_counter()
    manifest_path = manifest_path.resolve()
    control_source = control_source.resolve()
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    parent_archive = Path(manifest["artifacts"]["sindy_archive"]).resolve()
    parent = load_krom_archive(parent_archive)
    candidates, screening = _candidate_schedule(
        parent, maximum_mode_count=maximum_mode_count
    )
    policy = KROMControlCampaignPolicy(max_attempts=max_attempts)
    _write_json(output_root / "candidate-screening.json", screening)
    if dry_run:
        payload = {
            "dry_run": True,
            "policy": policy.summary(),
            "parent_archive": str(parent_archive),
            "candidate_screening": str((output_root / "candidate-screening.json").resolve()),
            "scheduled_attempts": [
                item.summary() for item in candidates[: policy.max_attempts]
            ],
        }
        _write_json(output_root / "dry-run.json", payload)
        return payload

    source = loadmat(control_source)
    source_controls = np.asarray(source["U_source"], dtype=float)[0]
    controlled, _ = generate_duffing_controlled_training(
        source_controls, sop=DuffingControlledTrainingSOP()
    )
    controlled_states = np.vstack(
        [trajectory.states[:-1] for trajectory in controlled.trajectories]
    )
    evaluator = DuffingKROMSearch(
        parent_archive,
        controlled,
        controlled_states,
        output_root,
        output_dt=output_dt,
        grid_size=grid_size,
        case_timeout_seconds=case_timeout_seconds,
    )
    ledger_path = output_root / "campaign-ledger.json"
    prior = _rehydrate_attempts(ledger_path, candidates, policy) if resume else []
    campaign = run_krom_control_campaign(
        candidates,
        evaluator,
        policy=policy,
        ledger_path=ledger_path,
        prior_attempts=prior,
    )
    payload = {
        **campaign.summary(),
        "schema": "dynsyslab.duffing-bilinear-krom-control-search/v1",
        "contract": {
            "initial_state": INITIAL_STATE,
            "target_state": TARGET_STATE,
            "active_time_window": ACTIVE_WINDOW,
            "release_time": ACTIVE_WINDOW[1],
            "final_time": FINAL_TIME,
            "Q": "12500 I in each extracted split-real lifted space",
            "R": [[0.01]],
            "control_limit": None,
            "control_law": "u=-K(Psi(x)-Psi(target))",
            "solver": "SciPy BDF",
            "rtol": 1e-8,
            "atol": 1e-10,
            "max_step": 0.01,
            "output_dt": output_dt,
            "coordinate_system": "dimensionless",
            "control_evaluation": "continuous_rhs",
            "case_wall_time_limit_seconds": case_timeout_seconds,
        },
        "dataset": controlled.summary(),
        "screening": screening,
        "sources": {
            "manifest": str(manifest_path),
            "manifest_sha256": _sha(manifest_path),
            "parent_krom_archive": str(parent_archive),
            "parent_krom_sha256": _sha(parent_archive),
            "controlled_data": str(control_source),
            "controlled_data_sha256": _sha(control_source),
        },
        "total_seconds": perf_counter() - started,
    }
    _write_json(output_root / "campaign-report.json", payload)
    return payload
