from __future__ import annotations
import unittest
from job_manager.history_policy import resolve

class HistoryPolicyTests(unittest.TestCase):
    def test_enabled_checkpoint_policy_configures_low_level_solver(self):
        options,history,result=resolve({"enabled":True,"mode":"checkpoint","strategy":"low_level_stepper","checkpoint_interval":10.,"pause_at_checkpoint":True},resources={"wall_time_s":100.},options={},history={})
        self.assertTrue(result["enabled"]); self.assertEqual(options["adaptive_engine"],"low_level_stepper")
        self.assertEqual(history["mode"],"checkpoint"); self.assertTrue(history["pause_at_checkpoint"])
    def test_threshold_policy_does_not_enable_without_sufficient_estimate(self):
        options,history,result=resolve({"enabled":True,"if_estimated_seconds_over":100.},resources={"wall_time_s":20.},options={},history={})
        self.assertFalse(result["enabled"]); self.assertEqual(options,{}); self.assertEqual(history,{})

    def test_planning_interval_alias_is_converted_to_solver_field(self):
        _, history, result=resolve(
            {"enabled":True,"interval_seconds":15.},
            resources={}, options={}, history={},
        )
        self.assertTrue(result["enabled"])
        self.assertEqual(history["checkpoint_interval"],15.)

    def test_unknown_estimate_can_explicitly_enable_policy(self):
        options, history, result=resolve(
            {"enabled":True,"if_estimated_seconds_over":300.,"when_estimate_unknown":"enable"},
            resources={"wall_time_s":None}, options={}, history={},
        )
        self.assertTrue(result["enabled"])
        self.assertTrue(result["estimate_was_unknown"])
        self.assertEqual(options["adaptive_engine"],"low_level_stepper")
        self.assertEqual(history["mode"],"checkpoint")

    def test_no_policy_preserves_input_dictionaries(self):
        options,history,result=resolve(None,resources={},options={"method":"DOP853"},history={"telemetry_heartbeat":False})
        self.assertEqual(options,{"method":"DOP853"}); self.assertEqual(history,{"telemetry_heartbeat":False}); self.assertEqual(result["reason"],"no_policy")

    def test_invalid_unknown_estimate_choice_is_rejected(self):
        with self.assertRaises(ValueError): resolve({"enabled":True,"when_estimate_unknown":"maybe"},resources={},options={},history={})

    def test_invalid_history_mode_is_rejected(self):
        with self.assertRaises(ValueError): resolve({"enabled":True,"mode":"everything"},resources={},options={},history={})

    def test_invalid_checkpoint_strategy_is_rejected(self):
        with self.assertRaises(ValueError): resolve({"enabled":True,"strategy":"restart_forever"},resources={},options={},history={})

    def test_telemetry_policy_does_not_force_low_level_engine(self):
        options,history,result=resolve({"enabled":True,"mode":"telemetry","strategy":"low_level_stepper"},resources={},options={"method":"RK45"},history={})
        self.assertNotIn("adaptive_engine",options); self.assertFalse(history["track_during_solve"]); self.assertTrue(result["enabled"])
