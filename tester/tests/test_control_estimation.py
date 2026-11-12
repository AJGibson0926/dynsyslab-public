from __future__ import annotations
import json
import tempfile
import unittest
import numpy as np

from control import (CallableTarget, ConstantTarget, LQRController, MPCController,
                     PIDController, RPE_CLEAN_ISOLATION_KLQR_SETTINGS, RPE_CONTROL_TARGETS,
                     RPE_OPTIONAL_PRODUCTION_CONTROL_SEQUENCE, TimeSeriesTarget,
                     WaveformController, rpe_control_target)
from estimation import (IdentityEstimator, LinearLuenbergerEstimator, LQEEstimator,
                        StaticMapEstimator)
from dynsyslib import load_system
from solver import SolveRequest, SolverOptions, solve
from solver.ode import _controller_signal
from time_series import TimeSeries


class _Rest:
    name="rest"; initial_state=np.array([0.])
    def rhs(self,t,x): return np.zeros(1)

class _ConstantController:
    actuation="additive_state"
    active_time_window=(.2,.4)
    def value(self,t,state): return np.ones(1)

class _CountingEstimator:
    def __init__(self): self.calls=0
    def estimate(self,t,state): self.calls+=1; return np.asarray(state)+1.


class ControlAndEstimationTests(unittest.TestCase):
    def test_clean_isolation_klqr_policy_is_immutable_and_dimensionless(self):
        policy=RPE_CLEAN_ISOLATION_KLQR_SETTINGS
        self.assertEqual(policy["coordinate_system"],"dimensionless")
        self.assertEqual(policy["input_units"],"p_control/p0")
        self.assertEqual(policy["Q_diagonal"],(5000.,5000.))
        self.assertEqual(policy["R"],1e-4)
        self.assertEqual(policy["control_evaluation"],"continuous_rhs")
        self.assertIsNone(policy["control_limit"])
        self.assertTrue(policy["accept_initial_transient"])
        self.assertFalse(policy["startup_ramp"])
        with self.assertRaises(TypeError): policy["R"]=1e-3

    def test_optional_production_sequence_cannot_change_the_default_policy(self):
        optional=RPE_OPTIONAL_PRODUCTION_CONTROL_SEQUENCE
        self.assertFalse(optional["enabled_by_default"])
        self.assertTrue(optional["opt_in_required"])
        self.assertEqual(optional["reference_engagement"]["tested_duration_tau"],2.)
        self.assertFalse(RPE_CLEAN_ISOLATION_KLQR_SETTINGS["startup_ramp"])
        with self.assertRaises(TypeError): optional["enabled_by_default"]=True
    def test_canonical_rpe_targets_have_analytic_radius_velocity_pairs(self):
        self.assertEqual(len(RPE_CONTROL_TARGETS),7)
        target=rpe_control_target("rpe_multisine_offset_4_periodic")
        value=target.evaluate(0.,2)
        self.assertAlmostEqual(value[0],4.)
        step=1e-7
        derivative=(target.evaluate(step,2)[0]-target.evaluate(-step,2)[0])/(2*step)
        self.assertAlmostEqual(value[1],derivative,places=6)
        with self.assertRaises(ValueError): target.evaluate(0.,1)

    def test_equilibrium_multisine_target_starts_at_equilibrium_and_remains_small(self):
        target=rpe_control_target("rpe_multisine_equilibrium_small")
        np.testing.assert_allclose(target.evaluate(0.,2),np.asarray((1.,0.)),atol=1e-14)
        samples=np.vstack([target.evaluate(t,2) for t in np.linspace(0.,20.,20001)])
        self.assertGreaterEqual(float(samples[:,0].min()),.86-1e-12)
        self.assertLessEqual(float(samples[:,0].max()),1.14+1e-12)

    def test_moderate_identification_target_starts_at_equilibrium(self):
        target=rpe_control_target("rpe_multisine_equilibrium_moderate_identification")
        np.testing.assert_allclose(target.evaluate(0.,2),np.asarray((1.,0.)),atol=1e-14)
        samples=np.vstack([target.evaluate(t,2) for t in np.linspace(0.,20.,20001)])
        self.assertGreaterEqual(float(samples[:,0].min()),.79)
        self.assertLessEqual(float(samples[:,0].max()),1.21)

    def test_constant_target_broadcasts_scalar_and_rejects_wrong_dimension(self):
        self.assertTrue(np.array_equal(ConstantTarget(np.array([2.])).evaluate(0.,3),[2.,2.,2.]))
        with self.assertRaises(ValueError): ConstantTarget(np.array([1.,2.])).evaluate(0.,3)

    def test_callable_target_delegates_and_normalizes(self):
        target=CallableTarget(lambda t: np.array([t]))
        self.assertTrue(np.array_equal(target.evaluate(.25,2),[.25,.25]))

    def test_time_series_target_interpolates_and_checks_domain(self):
        series=TimeSeries(np.array([0.,1.,2.]),np.array([[0.],[1.],[4.]]))
        self.assertAlmostEqual(float(TimeSeriesTarget(series,"linear").evaluate(.5,1)[0]),.5)
        with self.assertRaises(ValueError): TimeSeriesTarget(series).evaluate(3.,1)

    def test_waveform_returns_multitone_value_and_validates_window(self):
        controller=WaveformController(np.array([1.]),np.array([2.]),phases_rad=np.array([np.pi/2]),offset=1.)
        self.assertAlmostEqual(float(controller.value(0.,np.zeros(1))[0]),3.)
        with self.assertRaises(ValueError): WaveformController(np.array([1.]),np.array([1.]),active_time_window=(2.,1.))

    def test_lqr_generates_restoring_control_and_caches_fixed_gain(self):
        controller=LQRController(np.array([[-1.]]),np.array([[1.]]),np.eye(1),np.eye(1))
        first=controller.value(0.,np.array([1.])); gain=controller.last_gain.copy(); second=controller.value(1.,np.array([2.]))
        self.assertLess(first[0],0.); self.assertLess(second[0],first[0]); self.assertTrue(np.array_equal(gain,controller.last_gain))

    def test_explicit_sampled_state_dependent_lqr_requires_sample_time(self):
        controller=LQRController(lambda x:np.array([[-1.-x[0]*0.]]),np.array([[1.]]),np.eye(1),np.eye(1))
        self.assertTrue(controller.requires_sample_clock)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError,"control_sample_time"):
                solve(SolveRequest(system=_Rest(),kind="ode",time_span=(0.,.1),archive_root=directory,
                    options=SolverOptions(control_evaluation="sampled"),controller=controller))

    def test_continuous_rhs_recomputes_stateless_control_and_samples_accepted_states(self):
        class Controller:
            actuation="additive_state"; active_time_window=None; continuous_rhs_safe=True
            def __init__(self): self.calls=0
            def value(self,t,state):
                self.calls+=1
                return np.asarray((1.+float(t)+float(state[0]),))
        controller=Controller()
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Rest(),kind="ode",time_span=(0.,.2),archive_root=directory,
                options=SolverOptions(method="RK4",step_size=.1,control_evaluation="continuous_rhs"),controller=controller))
        self.assertEqual(result.status,"complete")
        self.assertEqual(result.metadata["control_evaluation"],"continuous_rhs")
        self.assertEqual(result.metadata["controller_call_count"],8)
        self.assertEqual(result.control_call_time.size,8)
        self.assertEqual(controller.calls,8+result.time.size)
        expected=1.+result.time+result.values[:,0]
        np.testing.assert_allclose(result.control_time_series.values[:,0],expected)

    def test_continuous_rhs_rejects_stateful_controller(self):
        controller=PIDController(kp=1.,ki=1.,target=ConstantTarget(np.asarray((1.,))))
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError,"stateless RHS-safe"):
                solve(SolveRequest(system=_Rest(),kind="ode",time_span=(0.,.1),archive_root=directory,
                    options=SolverOptions(control_evaluation="continuous_rhs"),controller=controller))

    def test_sampled_clock_records_only_real_controller_updates_without_duplicate_ticks(self):
        class Controller:
            actuation="additive_state"; active_time_window=None
            def __init__(self): self.calls=0
            def value(self,t,state): self.calls+=1; return np.ones(1)
        controller=Controller()
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Rest(),kind="ode",time_span=(0.,.3),archive_root=directory,
                options=SolverOptions(method="RK4",step_size=.1,control_evaluation="sampled",control_sample_time=.1),controller=controller))
        self.assertEqual(result.status,"complete")
        # RK4 evaluates its final stage at the right endpoint, where the
        # right-continuous held command advances to the next integer tick.
        self.assertEqual(controller.calls,4)
        self.assertEqual(result.metadata["controller_call_count"],4)
        np.testing.assert_allclose(result.control_call_time,np.asarray((0.,.1,.2,.3)),atol=1e-14)

    def test_mpc_respects_control_limit(self):
        controller=MPCController(np.array([[-1.]]),np.array([[1.]]),np.eye(1),np.eye(1),control_limit=np.array([.1]),sample_time=.1)
        self.assertLessEqual(abs(float(controller.value(0.,np.array([10.]))[0])),.1000001)

    def test_pid_accumulates_integral_and_derivative_on_increasing_time(self):
        controller=PIDController(kp=1.,ki=1.,kd=1.,target=ConstantTarget(np.array([1.])))
        first=controller.value(0.,np.array([0.])); second=controller.value(1.,np.array([0.])); repeated=controller.value(1.,np.array([0.]))
        self.assertGreater(second[0],first[0]); self.assertAlmostEqual(repeated[0],second[0])

    def test_identity_and_static_estimators_are_separate_from_input(self):
        source=np.array([1.,2.]); identity=IdentityEstimator().estimate(0.,source)
        identity[0]=99.
        self.assertEqual(source[0],1.)
        self.assertTrue(np.array_equal(StaticMapEstimator(lambda x:2*x).estimate(0.,source),[2.,4.]))

    def test_luenberger_observer_only_advances_on_increasing_time(self):
        observer=LinearLuenbergerEstimator(np.array([[0.]]),np.array([[1.]]),np.array([[1.]]),np.array([0.]))
        self.assertEqual(float(observer.estimate(0.,np.array([1.]))[0]),0.)
        self.assertAlmostEqual(float(observer.estimate(1.,np.array([1.]))[0]),1.)
        self.assertAlmostEqual(float(observer.estimate(.5,np.array([0.]))[0]),1.)

    def test_lqe_builds_gain_and_updates_estimate(self):
        observer=LQEEstimator(np.array([[-1.]]),np.array([[1.]]),np.eye(1),np.eye(1),np.array([0.]))
        self.assertEqual(observer.gain.shape,(1,1))
        observer.estimate(0.,np.array([1.])); updated=observer.estimate(.5,np.array([1.]))
        self.assertGreater(updated[0],0.)

    def test_solver_enforces_control_window_and_archives_actual_control_history(self):
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Rest(),kind="ode",time_span=(0.,1.),archive_root=directory,
                options=SolverOptions(method="Euler",step_size=.1),controller=_ConstantController()))
            self.assertEqual(result.status,"complete")
            self.assertAlmostEqual(float(result.values[-1,0]),.3,places=8)
            self.assertIsNotNone(result.control_time_series)
            self.assertEqual(float(result.control_time_series.values[0,0]),0.)

    def test_solver_passes_estimated_state_to_controller(self):
        estimator=_CountingEstimator()
        class Controller:
            actuation="additive_state"; active_time_window=None
            def __init__(self): self.states=[]
            def value(self,t,state): self.states.append(float(state[0])); return np.zeros(1)
        controller=Controller()
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=_Rest(),kind="ode",time_span=(0.,.1),archive_root=directory,controller=controller,estimator=estimator))
        self.assertEqual(result.status,"complete"); self.assertGreater(estimator.calls,0); self.assertTrue(all(value==1. for value in controller.states))

    def test_dimensionless_solver_passes_active_coordinates_to_controller(self):
        class Controller:
            active_time_window=None
            def __init__(self): self.received=None
            def value(self,t,state): self.received=(float(t),np.asarray(state).copy()); return np.zeros(1)
        controller=Controller(); request=SolveRequest(system=_Rest(),kind="ode",time_span=(0.,1.),controller=controller)
        _controller_signal(request,.25,np.asarray((.5,)),dimensionless=True)
        self.assertEqual(controller.received[0],.25)
        np.testing.assert_allclose(controller.received[1],np.asarray((.5,)))

    def test_controlled_rpe_honors_explicit_non_equilibrium_initial_state(self):
        class RPERest(_Rest):
            name="rayleigh-plesset"; equilibrium_radius_m=1e-5; initial_state=np.asarray((1e-5,0.))
            def rhs(self,t,x): return np.zeros(2)
        class Controller:
            active_time_window=None
            def value(self,t,state): return np.zeros(2)
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=RPERest(),kind="ode",time_span=(0.,.1),archive_root=directory,
                initial_state=np.asarray((2e-5,.25)),controller=Controller()))
        self.assertEqual(result.status,"complete")
        np.testing.assert_array_equal(result.values[0],np.asarray((2e-5,.25)))
        self.assertEqual(result.metadata["initial_state_policy"],"controlled-rpe-request-override")

    def test_controlled_rpe_defaults_to_equilibrium_independently_of_free_default(self):
        class RPERest(_Rest):
            name="rayleigh-plesset"; equilibrium_radius_m=1e-5
            initial_state=np.asarray((2.333386598601735e-5,0.))
            def rhs(self,t,x): return np.zeros(2)
        class Controller:
            active_time_window=None
            def value(self,t,state): return np.zeros(2)
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=RPERest(),kind="ode",time_span=(0.,.01),
                archive_root=directory,controller=Controller()))
        self.assertEqual(result.status,"complete")
        np.testing.assert_array_equal(result.values[0],np.asarray((1e-5,0.)))
        self.assertEqual(result.metadata["initial_state_policy"],"controlled-rpe-equilibrium-default")

    def test_controlled_dimensionless_rpe_default_is_one_zero(self):
        class Controller:
            actuation="system_input"; active_time_window=None
            def value(self,t,state): return np.zeros(1)
        system=load_system("rpe")
        scale=system.dimensionless_parameters.time_scale_s
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=system,kind="ode",time_span=(0.,.001*scale),
                archive_root=directory,controller=Controller(),
                options=SolverOptions(method="RK4",step_size=.001*scale,use_dimensionless_form=True)))
            with np.load(result.archive_path,allow_pickle=False) as archive:
                archived=json.loads(str(archive["metadata_json"].item()))
                self.assertEqual(archived["control_history"]["coordinate_system"],"dimensionless")
                self.assertEqual(archived["control_history"]["units"],"p_control/p0")
                self.assertEqual(archived["control_history"]["scale_to_dimensional"],system.ambient_pressure_pa)
        self.assertEqual(result.status,"complete")
        self.assertEqual(result.metadata["controller_coordinate_system"],"dimensionless")
        self.assertEqual(result.metadata["control_input_units"],"p_control/p0")
        np.testing.assert_allclose(system.to_dimensionless_state(result.values[0]),np.asarray((1.,0.)),
                                   rtol=0.,atol=0.)
        self.assertEqual(result.metadata["initial_state_policy"],"controlled-rpe-equilibrium-default")

    def test_free_rpe_uses_its_normal_system_default(self):
        class RPERest(_Rest):
            name="rayleigh-plesset"; equilibrium_radius_m=1e-5
            initial_state=np.asarray((2.333386598601735e-5,-.25))
            def rhs(self,t,x): return np.zeros(2)
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=RPERest(),kind="ode",time_span=(0.,.01),
                archive_root=directory))
        self.assertEqual(result.status,"complete")
        np.testing.assert_array_equal(result.values[0],RPERest.initial_state)
        self.assertEqual(result.metadata["initial_state_policy"],"system-default")

    def test_solver_archives_target_aligned_history_and_exact_controller_call_log(self):
        class Rest2:
            name="rest-2"; initial_state=np.zeros(2)
            def rhs(self,t,x): return np.zeros(2)
        controller=LQRController(-np.eye(2),np.eye(2),np.eye(2),np.eye(2),
                                 target=rpe_control_target("rpe_multisine_offset_4_mixed"))
        with tempfile.TemporaryDirectory() as directory:
            result=solve(SolveRequest(system=Rest2(),kind="ode",time_span=(0.,.05),archive_root=directory,
                options=SolverOptions(method="Euler",step_size=.01),controller=controller))
            self.assertIsNotNone(result.control_target_time_series)
            self.assertEqual(result.control_target_time_series.values.shape,result.values.shape)
            self.assertEqual(result.control_call_time.size,int(result.metadata["controller_call_count"]))
            self.assertGreater(result.control_call_time.size,0)
            with np.load(result.archive_path,allow_pickle=False) as archive:
                self.assertIn("control_values",archive.files)
                self.assertIn("control_target_values",archive.files)
                self.assertIn("control_call_time",archive.files)
                self.assertEqual(archive["control_call_time"].shape[0],archive["control_call_values"].shape[0])
                metadata=json.loads(str(archive["metadata_json"].item()))
                self.assertEqual(metadata["control_history"]["call_log_semantics"],"actual_controller_evaluations")
                self.assertEqual(metadata["control_history"]["evaluation_policy"],"continuous_rhs")
