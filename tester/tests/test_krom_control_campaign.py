from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from rom import (
    KROMCandidate,
    KROMControlAttempt,
    KROMControlCampaignPolicy,
    MAX_KROM_CONTROL_ATTEMPTS,
    resample_trajectory_by_arclength,
    run_krom_control_campaign,
    trajectory_tube_states,
)


def _case(release: float, final: float, *, status: str = "complete"):
    return {
        "status": status,
        "metrics": {
            "distance_to_target_at_release": release,
            "final_distance_to_left_equilibrium": final,
        },
    }


class KROMControlCampaignTests(unittest.TestCase):
    def test_policy_rejects_any_attempt_budget_above_absolute_limit(self):
        with self.assertRaisesRegex(ValueError, "max_attempts"):
            KROMControlCampaignPolicy(max_attempts=MAX_KROM_CONTROL_ATTEMPTS + 1)

    def test_campaign_runs_every_case_and_stops_on_bilinear_success(self):
        candidates = [KROMCandidate((index,), {"selection_score": 1.0 / (index + 1)})
                      for index in range(4)]
        calls = []

        def evaluate(number, candidate):
            calls.append((number, candidate.retained_indices))
            success = number == 3
            cases = {
                "exact": _case(0.001, 0.01),
                "constant": _case(0.4, 0.3),
                "bilinear-dynamics": _case(0.01 if success else 0.3, 0.02 if success else 0.4),
                "bilinear-map": _case(0.5, 0.5),
            }
            return KROMControlAttempt(number, candidate, cases, {}, {})

        with tempfile.TemporaryDirectory() as directory:
            ledger = Path(directory) / "ledger.json"
            result = run_krom_control_campaign(candidates, evaluate, ledger_path=ledger)
            stored = json.loads(ledger.read_text(encoding="utf-8"))
        self.assertEqual(calls, [(1, (0,)), (2, (1,)), (3, (2,))])
        self.assertEqual(result.stop_reason, "bilinear-success")
        self.assertEqual(result.successful_attempt.attempt_number, 3)
        self.assertEqual(stored["attempt_count"], 3)
        self.assertEqual(stored["stop_reason"], "bilinear-success")

    def test_campaign_hard_cap_cannot_be_bypassed_by_long_candidate_list(self):
        policy = KROMControlCampaignPolicy(max_attempts=2, stop_on_bilinear_success=False)
        candidates = [KROMCandidate((index,), {"selection_score": 1.0}) for index in range(5)]

        def evaluate(number, candidate):
            cases = {name: _case(1.0, 1.0) for name in policy.required_cases}
            return KROMControlAttempt(number, candidate, cases, {}, {})

        result = run_krom_control_campaign(candidates, evaluate, policy=policy)
        self.assertEqual(len(result.attempts), 2)
        self.assertEqual(result.stop_reason, "hard-attempt-limit-reached")

    def test_attempt_rejects_missing_control_comparison(self):
        policy = KROMControlCampaignPolicy(max_attempts=1)
        candidate = KROMCandidate((0,), {"selection_score": 1.0})

        def evaluate(number, selected):
            return KROMControlAttempt(number, selected, {"exact": _case(0.0, 0.0)}, {}, {})

        with self.assertRaisesRegex(ValueError, "omitted required control cases"):
            run_krom_control_campaign((candidate,), evaluate, policy=policy)

    def test_trajectory_tube_sampling_is_deterministic_bounded_and_nonmutating(self):
        path = np.asarray(((0.0, 0.0), (0.5, -0.5), (1.0, 1.0)))
        original = path.copy()
        sampled = trajectory_tube_states(
            path, np.asarray((-0.25, -0.75)), np.asarray((0.75, 0.75)), radii=(0.1,)
        )
        self.assertEqual(sampled.shape, (15, 2))
        np.testing.assert_array_equal(path, original)
        self.assertTrue(np.all(sampled >= np.asarray((-0.25, -0.75))))
        self.assertTrue(np.all(sampled <= np.asarray((0.75, 0.75))))
        np.testing.assert_allclose(sampled[:3], np.clip(path, (-0.25, -0.75), (0.75, 0.75)))

    def test_arclength_resampling_ignores_dwell_density(self):
        path = np.asarray(((0.0, 0.0), (0.0, 0.0), (1.0, 0.0), (1.0, 1.0)))
        sampled = resample_trajectory_by_arclength(path, 5)
        np.testing.assert_allclose(
            sampled,
            ((0.0, 0.0), (0.5, 0.0), (1.0, 0.0), (1.0, 0.5), (1.0, 1.0)),
        )


if __name__ == "__main__":
    unittest.main()
