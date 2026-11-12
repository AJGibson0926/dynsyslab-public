from __future__ import annotations
import tempfile, unittest
import json
from pathlib import Path
import numpy as np
from threading import Event
from solver import HistoryOptions, SolveRequest, SolverOptions, resume_ode_from_pause, solve
from solver import InterpolationOptions

class _Decay:
    name="decay"
    initial_state=np.array([1.])
    def rhs(self,t,x): return -np.asarray(x)

class _SlowDecay(_Decay):
    def rhs(self,t,x):
        from time import sleep
        sleep(.012)
        return super().rhs(t,x)

class SolverCheckpointPauseTests(unittest.TestCase):
    def test_pause_at_checkpoint_then_extend_same_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            token=Event(); token.set()
            request=SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),label="checkpoint-test",archive_root=directory,
                options=SolverOptions(method="RK45"),history=HistoryOptions(track_during_solve=True,checkpoint_interval=.25,pause_at_checkpoint=True),pause_token=token)
            paused=solve(request)
            self.assertEqual(paused.status,"paused"); self.assertTrue(Path(paused.archive_path).is_file())
            resumed=resume_ode_from_pause(paused.archive_path,request)
            self.assertEqual(resumed.status,"complete"); self.assertEqual(Path(resumed.archive_path),Path(paused.archive_path))
            with np.load(paused.archive_path,allow_pickle=False) as archive:
                import json
                metadata=json.loads(str(archive["metadata_json"]))
            self.assertIn("checkpoint_restart",metadata)

    def test_checkpoint_mode_without_pause_token_completes_normally(self):
        with tempfile.TemporaryDirectory() as directory:
            request=SolveRequest(system=_Decay(),kind="ode",time_span=(0.,.5),label="checkpoint-normal",archive_root=directory,
                history=HistoryOptions(track_during_solve=True,checkpoint_interval=.1,pause_at_checkpoint=True))
            result=solve(request)
            self.assertEqual(result.status,"complete")
            self.assertTrue(Path(result.archive_path).is_file())

    def test_resume_rejects_a_nonpaused_solution(self):
        with tempfile.TemporaryDirectory() as directory:
            request=SolveRequest(system=_Decay(),kind="ode",time_span=(0.,.25),label="complete",archive_root=directory)
            complete=solve(request)
            with self.assertRaises(ValueError): resume_ode_from_pause(complete.archive_path,request)

    def test_fixed_step_solver_accepts_checkpoint_history_configuration(self):
        with tempfile.TemporaryDirectory() as directory:
            request=SolveRequest(system=_Decay(),kind="ode",time_span=(0.,.2),label="rk4",archive_root=directory,
                options=SolverOptions(method="RK4",step_size=.05),history=HistoryOptions(track_during_solve=True,checkpoint_interval=.1))
            self.assertEqual(solve(request).status,"complete")

    def test_rk4_decay_is_accurate(self):
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),label="rk4-accuracy",archive_root=directory,options=SolverOptions(method="RK4",step_size=.01)))
            self.assertAlmostEqual(float(result.values[-1,0]),float(np.exp(-1)),places=6)

    def test_euler_decay_has_expected_first_order_accuracy(self):
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),label="euler",archive_root=directory,options=SolverOptions(method="Euler",step_size=.01)))
            self.assertLess(abs(float(result.values[-1,0])-float(np.exp(-1))),.01)

    def test_fixed_step_method_requires_step_size(self):
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),label="bad-rk4",archive_root=directory,options=SolverOptions(method="RK4")))
            self.assertEqual(result.status,"failed")
            self.assertIn("requires",result.message)

    def test_cancellation_token_causes_failed_result_without_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            token=Event(); token.set()
            result=solve(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),label="cancel",archive_root=directory,cancel_token=token))
            self.assertEqual(result.status,"failed")
            self.assertIsNone(result.archive_path)

    def test_explicit_interpolation_grid_is_honored(self):
        with tempfile.TemporaryDirectory() as directory:
            grid=np.linspace(0.,1.,11)
            result=solve(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),label="grid",archive_root=directory,interpolation=InterpolationOptions(times=grid)))
            self.assertTrue(np.array_equal(result.time,grid))

    def test_low_level_stepper_supports_sampled_history(self):
        with tempfile.TemporaryDirectory() as directory:
            request=SolveRequest(system=_Decay(),kind="ode",time_span=(0.,.5),label="sampled",archive_root=directory,
                options=SolverOptions(method="RK45",adaptive_engine="low_level_stepper"),history=HistoryOptions(mode="sampled",sample_interval=.1,checkpoint_strategy="low_level_stepper"))
            result=solve(request)
            self.assertEqual(result.status,"complete")
            self.assertTrue(list((Path(directory)/"history"/"decay").rglob("checkpoint-*.npz")))

    def test_low_level_checkpoint_honors_pause_token(self):
        with tempfile.TemporaryDirectory() as directory:
            token=Event(); token.set()
            request=SolveRequest(system=_Decay(),kind="ode",time_span=(0.,.5),label="low-pause",archive_root=directory,
                options=SolverOptions(method="RK45",adaptive_engine="low_level_stepper"),history=HistoryOptions(mode="checkpoint",checkpoint_interval=.1,pause_at_checkpoint=True,checkpoint_strategy="low_level_stepper"),pause_token=token)
            self.assertEqual(solve(request).status,"paused")

    def test_default_telemetry_writes_lifecycle_log_without_history_arrays(self):
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,.1),label="telemetry",archive_root=directory))
            self.assertEqual(result.status,"complete")
            self.assertTrue(list((Path(directory)/"solver_logs").glob("*.telemetry.jsonl")))
            self.assertFalse((Path(directory)/"history").exists())

    def test_default_heartbeat_options_are_lightweight_and_enabled(self):
        options=HistoryOptions()
        self.assertTrue(options.telemetry_enabled)
        self.assertTrue(options.telemetry_heartbeat)
        self.assertEqual(options.mode,"telemetry")
        self.assertFalse(options.track_during_solve)

    def test_heartbeat_can_be_disabled_without_disabling_solve(self):
        with tempfile.TemporaryDirectory() as directory:
            request=SolveRequest(system=_Decay(),kind="ode",time_span=(0.,.1),label="no-heartbeat",archive_root=directory,
                history=HistoryOptions(telemetry_heartbeat=False))
            self.assertEqual(solve(request).status,"complete")

    def test_heartbeat_reports_scalar_progress_without_history(self):
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_SlowDecay(),kind="ode",time_span=(0.,1.),label="heartbeat",archive_root=directory,
                history=HistoryOptions(telemetry_interval_seconds=.01)))
            telemetry=next((Path(directory)/"solver_logs").glob("*.telemetry.jsonl"))
            events=[json.loads(line) for line in telemetry.read_text(encoding="utf-8").splitlines()]
            heartbeats=[event for event in events if event["event"]=="heartbeat"]
            self.assertEqual(result.status,"complete")
            self.assertTrue(heartbeats)
            self.assertIn("simulation_time",heartbeats[0])
            self.assertIn("progress_fraction",heartbeats[0])
            self.assertFalse((Path(directory)/"history").exists())
