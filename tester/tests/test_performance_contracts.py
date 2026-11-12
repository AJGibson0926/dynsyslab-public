import unittest
from unittest.mock import patch
import importlib.util
from types import SimpleNamespace

from performance import BenchmarkCase, evaluate_relative_speedups, evaluate_speedups, run_benchmarks, runtime_fingerprint
from performance.cli import _cases, _select_cases
from solver.backends import PYTHON_BACKEND, available_solver_backends, select_ode_backend, select_reference_backend, sundials_backend
from solver.models import HistoryOptions, SolveRequest, SolverOptions
from solver.native_ode import native_fixed_eligibility, native_model
from solver.native_ode import integrate_fixed_native
from solver.ode import _fixed_step, solve_ode
from solver.diagnostics import backend_diagnostics
from dynsyslib import load_system
from dynsyslib.drivers import SinusoidalDriver, TabulatedDriver
import numpy as np


class PerformanceContractTests(unittest.TestCase):
    def test_lazy_sindy_value_reconstruction_is_a_selectable_benchmark(self):
        cases=_cases(1,0); selected=_select_cases(cases,["sindy-lazy-polynomial-values-5000x256"])
        self.assertEqual(len(selected),1)
        self.assertEqual(selected[0].function().shape,(5000,256))
        with self.assertRaisesRegex(ValueError,"Unknown benchmark"):
            _select_cases(cases,["missing-case"])

    def test_runtime_fingerprint_is_versioned_and_serializable(self):
        payload = runtime_fingerprint()
        self.assertEqual(payload["schema"], "dynsyslab.performance-runtime/v1")
        self.assertIn("numpy", payload)
        self.assertIn("scipy", payload)
        self.assertIsInstance(payload["threadpools"], list)
        self.assertIn("scipy",payload["linear_algebra_backends"])

    def test_backend_diagnostics_are_versioned_and_serializable(self):
        payload=backend_diagnostics()
        self.assertEqual(payload["schema"],"dynsyslab.backend-diagnostics/v1")
        self.assertIn("sundials",payload["solver"])
        self.assertTrue(payload["postprocessing"]["linear_algebra"]["scipy"]["available"])
        __import__("json").dumps(payload)

    def test_benchmark_report_records_distribution(self):
        report = run_benchmarks((BenchmarkCase("noop", lambda: None, "test", iterations=3, warmup=0),))
        payload = report.payload()
        self.assertEqual(payload["schema"], "dynsyslab.benchmark-report/v1")
        self.assertEqual(len(payload["samples"][0]["seconds"]), 3)
        self.assertLessEqual(payload["samples"][0]["minimum_seconds"], payload["samples"][0]["maximum_seconds"])

    def test_performance_gate_compares_medians_and_runtime(self):
        runtime = {"machine": "x", "platform": "p", "python": "3", "numpy": "2", "scipy": "1"}
        reference = {"schema": "dynsyslab.benchmark-report/v1", "runtime": runtime, "samples": [{"name": "kernel", "median_seconds": 2.0}]}
        candidate = {"schema": "dynsyslab.benchmark-report/v1", "runtime": runtime, "samples": [{"name": "kernel", "median_seconds": 0.5}]}
        gate = evaluate_speedups(reference, candidate, {"kernel": 3.0})
        self.assertTrue(gate.passed)
        self.assertEqual(gate.results[0].speedup, 4.0)

    def test_performance_gate_rejects_cross_runtime_comparison(self):
        reference = {"schema": "dynsyslab.benchmark-report/v1", "runtime": {"machine": "a"}, "samples": [{"name": "kernel", "median_seconds": 1.0}]}
        candidate = {"schema": "dynsyslab.benchmark-report/v1", "runtime": {"machine": "b"}, "samples": [{"name": "kernel", "median_seconds": 0.1}]}
        gate = evaluate_speedups(reference, candidate, {"kernel": 2.0})
        self.assertFalse(gate.passed)
        self.assertIn("machine", gate.runtime_differences)

    def test_relative_gate_compares_two_cases_from_one_runtime(self):
        report = {"schema": "dynsyslab.benchmark-report/v1", "samples": [
            {"name": "python", "median_seconds": 3.0},
            {"name": "native", "median_seconds": 1.0},
        ]}
        gate = evaluate_relative_speedups(report, {("native", "python"): 2.5})
        self.assertTrue(gate.passed)
        self.assertEqual(gate.results[0].speedup, 3.0)

    def test_invalid_benchmark_counts_are_rejected(self):
        with self.assertRaises(ValueError):
            BenchmarkCase("bad", lambda: None, "test", iterations=0)

    def test_python_backend_is_always_available(self):
        self.assertIs(select_ode_backend("python", "RK4"), PYTHON_BACKEND)
        self.assertTrue(available_solver_backends()["python"].available)

    @unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
    def test_native_backend_reports_coordinate_explicit_rpe_kernel(self):
        backend = available_solver_backends()["native"]
        self.assertIn("rpe_dimensionless", backend.fixed_step_models)

    def test_auto_falls_back_when_native_extension_is_absent(self):
        with patch("solver.backends.import_module", side_effect=ImportError("missing")):
            selected = select_ode_backend("auto", "RK4")
        self.assertIs(selected, PYTHON_BACKEND)

    def test_explicit_native_never_silently_falls_back(self):
        with patch("solver.backends.import_module", side_effect=ImportError("missing")):
            with self.assertRaisesRegex(RuntimeError, "not installed"):
                select_ode_backend("native", "RK4")

    def test_invalid_backend_preference_is_rejected(self):
        with self.assertRaises(ValueError):
            select_ode_backend("gpu", "RK4")

    def test_reference_only_kind_rejects_explicit_native(self):
        with self.assertRaisesRegex(RuntimeError, "does not support 'pde'|not installed"):
            select_reference_backend("native", "pde")

    def test_duffing_has_compact_native_model_specification(self):
        system = load_system("duffing", preset="double_well_reference")
        specification = native_model(system)
        self.assertEqual(specification.identifier, "duffing")
        self.assertEqual(specification.parameters.shape, (6,))
        self.assertEqual(specification.coordinate_system, "dimensional")

    def test_rpe_has_coordinate_explicit_dimensionless_native_specification(self):
        system = load_system(
            "rpe",
            preset="air_water_20c",
            driver=SinusoidalDriver(frequency_hz=25_000.0, amplitude=1_500.0, phase_rad=0.2),
        )
        specification = native_model(system, coordinate_system="dimensionless")
        self.assertEqual(specification.identifier, "rpe_dimensionless")
        self.assertEqual(specification.parameters.shape, (11,))
        self.assertEqual(specification.coordinate_system, "dimensionless")
        self.assertIsNone(native_model(system, coordinate_system="dimensional"))
        groups = system.dimensionless_parameters
        np.testing.assert_allclose(
            specification.parameters,
            np.asarray((
                system.preset.kappa,
                groups.gas_pressure,
                groups.vapor_pressure,
                groups.ambient_pressure,
                groups.capillarity,
                groups.viscosity,
                groups.pressure_scale_pa,
                groups.time_scale_s,
                -1_500.0,
                25_000.0,
                0.2,
            )),
        )

    def test_rpe_native_eligibility_is_dimensionless_only(self):
        system = load_system("rpe", preset="air_water_20c")
        request = SolveRequest(
            system, "ode", time_span=(0.0, 1.0),
            options=SolverOptions(method="RK4", step_size=0.001),
        )
        eligible, reason = native_fixed_eligibility(request, dimensionless=True)
        self.assertTrue(eligible)
        self.assertIsNone(reason)
        eligible, reason = native_fixed_eligibility(request, dimensionless=False)
        self.assertFalse(eligible)
        self.assertIn("dimensional", reason)

    def test_old_native_wheel_rejects_rpe_before_dispatch_but_keeps_legacy_models(self):
        old = SimpleNamespace(capabilities=lambda: {
            "version":"0.1.0","fixed_step_ode":True,"batched_ode":True,
        })
        rpe = load_system("rpe", preset="air_water_20c")
        rpe_request = SolveRequest(rpe,"ode",time_span=(0.,1.),
            options=SolverOptions(method="RK4",step_size=.001))
        duffing = load_system("duffing", preset="double_well_reference")
        duffing_request = SolveRequest(duffing,"ode",time_span=(0.,1.),
            options=SolverOptions(method="RK4",step_size=.001))
        with patch("solver.native_ode.import_module",return_value=old):
            eligible,reason=native_fixed_eligibility(rpe_request,dimensionless=True)
            legacy_eligible,legacy_reason=native_fixed_eligibility(duffing_request,dimensionless=False)
        self.assertFalse(eligible)
        self.assertIn("0.1.0",reason)
        self.assertIn("rpe_dimensionless",reason)
        self.assertTrue(legacy_eligible)
        self.assertIsNone(legacy_reason)

    def test_tabulated_drivers_remain_on_reference_backend(self):
        driver = TabulatedDriver(np.asarray((0.0, 1.0)), np.asarray((0.0, 1.0)))
        system = load_system("duffing", driver=driver)
        request = SolveRequest(system, "ode", time_span=(0.0, 1.0), options=SolverOptions(method="RK4", step_size=0.01))
        eligible, reason = native_fixed_eligibility(request, dimensionless=False)
        self.assertFalse(eligible)
        self.assertIn("no native model adapter", reason)

    def test_adaptive_auto_remains_on_portable_admitted_backend(self):
        selected = select_ode_backend("auto", "RK45")
        self.assertEqual(selected.name, "python-reference")

    @unittest.skipUnless(sundials_backend().available, "optional sundials4py backend is not installed")
    def test_explicit_native_adaptive_solve_uses_sundials(self):
        class Decay:
            name = "adaptive-decay"
            initial_state = np.asarray((1.0,))
            def rhs(self, time, state): return -state
        class Session:
            def telemetry(self, *args, **kwargs): pass
        request = SolveRequest(
            Decay(), "ode", time_span=(0.0, 1.0),
            options=SolverOptions(method="RK45", rtol=1e-9, atol=1e-11, execution_backend="native"),
        )
        result = solve_ode(request, Session())
        self.assertEqual(result.status, "complete")
        self.assertAlmostEqual(result.values[-1, 0], np.exp(-1.0), places=8)
        self.assertEqual(result.metadata["execution_backend"]["name"], "sundials-cvode")
        self.assertEqual(result.metadata["adaptive_backend"]["library"], "SUNDIALS")

    def test_native_path_refuses_in_solve_history_callbacks(self):
        system = load_system("van_der_pol", preset="weakly_nonlinear")
        request = SolveRequest(
            system,
            "ode",
            time_span=(0.0, 1.0),
            options=SolverOptions(method="RK4", step_size=0.01),
            history=HistoryOptions(track_during_solve=True),
        )
        eligible, reason = native_fixed_eligibility(request, dimensionless=False)
        self.assertFalse(eligible)
        self.assertIn("checkpoint", reason)

    @unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
    def test_native_duffing_rk4_matches_python_reference(self):
        class Session:
            request = type("Request", (), {"history": HistoryOptions()})()

        system = load_system("duffing", preset="double_well_reference")
        expected_time, expected_values = _fixed_step(system.rhs, 0.0, 2.0, system.initial_state, 0.002, "RK4", Session())
        actual_time, actual_values = integrate_fixed_native(system, 0.0, 2.0, system.initial_state, 0.002, "RK4")
        np.testing.assert_array_equal(actual_time, expected_time)
        np.testing.assert_allclose(actual_values, expected_values, rtol=2e-13, atol=2e-13)

    @unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
    def test_native_van_der_pol_euler_matches_python_reference(self):
        class Session:
            request = type("Request", (), {"history": HistoryOptions()})()

        system = load_system("van_der_pol", preset="weakly_nonlinear")
        expected_time, expected_values = _fixed_step(system.rhs, 0.0, 0.5, system.initial_state, 0.001, "Euler", Session())
        actual_time, actual_values = integrate_fixed_native(system, 0.0, 0.5, system.initial_state, 0.001, "Euler")
        np.testing.assert_array_equal(actual_time, expected_time)
        np.testing.assert_allclose(actual_values, expected_values, rtol=2e-13, atol=2e-13)

    @unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
    def test_native_dimensionless_rpe_rk4_matches_python_reference_with_physical_driver(self):
        class Session:
            request = type("Request", (), {"history": HistoryOptions()})()

        system = load_system(
            "rpe",
            preset="air_water_20c",
            driver=SinusoidalDriver(frequency_hz=25_000.0, amplitude=1_500.0, phase_rad=0.2),
        )
        initial = np.asarray((1.15, -0.02))
        expected_time, expected_values = _fixed_step(
            system.dimensionless_rhs, 0.0, 0.25, initial, 0.0005, "RK4", Session()
        )
        actual_time, actual_values = integrate_fixed_native(
            system, 0.0, 0.25, initial, 0.0005, "RK4", coordinate_system="dimensionless"
        )
        np.testing.assert_array_equal(actual_time, expected_time)
        np.testing.assert_allclose(actual_values, expected_values, rtol=5e-12, atol=5e-13)

    @unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
    def test_native_kernel_batches_initial_conditions(self):
        system = load_system("duffing", preset="double_well_reference")
        initial = np.asarray(((0.1, 0.0), (0.2, -0.1), (-0.3, 0.2)))
        time, values = integrate_fixed_native(system, 0.0, 0.25, initial, 0.001, "RK4")
        self.assertEqual(values.shape, (3, time.size, 2))
        for index in range(initial.shape[0]):
            _, individual = integrate_fixed_native(system, 0.0, 0.25, initial[index], 0.001, "RK4")
            np.testing.assert_array_equal(values[index], individual)

    @unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
    def test_solve_ode_dispatches_native_and_records_provenance(self):
        class Session:
            def telemetry(self, *args, **kwargs):
                pass

        system = load_system("duffing", preset="double_well_reference")
        request = SolveRequest(
            system,
            "ode",
            time_span=(0.0, 0.1),
            options=SolverOptions(method="RK4", step_size=0.001, execution_backend="native"),
        )
        result = solve_ode(request, Session())
        self.assertEqual(result.status, "complete")
        self.assertEqual(result.metadata["execution_backend"]["name"], "dynsyslab-native")
        self.assertTrue(result.metadata["execution_backend"]["compiled"])

    @unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
    def test_solve_ode_dispatches_dimensionless_rpe_and_returns_physical_coordinates(self):
        class Session:
            request = type("Request", (), {"history": HistoryOptions()})()
            def telemetry(self, *args, **kwargs):
                pass

        system = load_system("rpe", preset="air_water_20c")
        scale = system.dimensionless_parameters.time_scale_s
        initial = np.asarray((1.15 * system.equilibrium_radius_m, -0.02 * system.equilibrium_radius_m / scale))
        common = dict(
            system=system,
            kind="ode",
            time_span=(0.0, 0.25 * scale),
            initial_state=initial,
        )
        native = solve_ode(SolveRequest(
            **common,
            options=SolverOptions(
                method="RK4", step_size=0.0005 * scale,
                use_dimensionless_form=True, execution_backend="native",
            ),
        ), Session())
        python = solve_ode(SolveRequest(
            **common,
            options=SolverOptions(
                method="RK4", step_size=0.0005 * scale,
                use_dimensionless_form=True, execution_backend="python",
            ),
        ), Session())
        self.assertEqual(native.status, "complete")
        self.assertEqual(native.metadata["integration_form"], "dimensionless")
        self.assertEqual(native.metadata["execution_backend"]["name"], "dynsyslab-native")
        np.testing.assert_array_equal(native.time, python.time)
        np.testing.assert_allclose(native.values, python.values, rtol=5e-12, atol=5e-13)


if __name__ == "__main__":
    unittest.main()
