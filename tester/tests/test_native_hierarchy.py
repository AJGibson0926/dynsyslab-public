import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np

from dynsyslib import load_system
from hierarchy import System
from solver import solve_hierarchy
from solver.models import SolverOptions


NATIVE_AVAILABLE = importlib.util.find_spec("_dynsyslab_native") is not None


def _ensemble(root: str, label: str) -> System:
    target = System(
        load_system("duffing", preset="double_well_reference"),
        time_span=(0.0, 0.2),
        uniform_time_steps=True,
        delta_t=0.01,
        label=label,
        archive_root=root,
    )
    target.set_initial_conditions(np.asarray(((0.1, 0.0), (0.2, -0.1), (-0.3, 0.2))))
    return target


@unittest.skipUnless(NATIVE_AVAILABLE, "optional native wheel is not installed")
class NativeHierarchyTests(unittest.TestCase):
    def test_native_cluster_uses_one_kernel_call_and_preserves_archives(self):
        with tempfile.TemporaryDirectory() as directory:
            target = _ensemble(directory, "native-ensemble")
            processor = Mock(side_effect=lambda path: path)
            processor.__name__ = "record"
            from solver.native_ode import integrate_fixed_native as actual_integrate

            with patch("solver.hierarchy.integrate_fixed_native", wraps=actual_integrate) as kernel:
                result = solve_hierarchy(
                    target,
                    options=SolverOptions(
                        method="RK4",
                        step_size=0.001,
                        execution_backend="native",
                    ),
                    postprocessors=(processor,),
                )

            self.assertEqual(result.status, "complete")
            self.assertEqual(kernel.call_count, 1)
            self.assertEqual(processor.call_count, 3)
            self.assertEqual(len(result.trajectory_results), 3)
            self.assertTrue(Path(result.archive_path).is_relative_to(Path(directory)))
            self.assertTrue(all(Path(item.archive_path).is_file() for item in result.trajectory_results))
            self.assertTrue(all(Path(item.archive_path).is_relative_to(Path(directory)) for item in result.trajectory_results))
            self.assertTrue(all(item.metadata["batch_execution"]["enabled"] for item in result.trajectory_results))
            self.assertTrue(all(item.metadata["execution_backend"]["compiled"] for item in result.trajectory_results))
            self.assertEqual(target.time_series.values.shape, (21, 3, 2))

    def test_native_cluster_matches_python_hierarchy(self):
        with tempfile.TemporaryDirectory() as native_root, tempfile.TemporaryDirectory() as python_root:
            native_target = _ensemble(native_root, "native-parity")
            python_target = _ensemble(python_root, "python-parity")
            native = solve_hierarchy(
                native_target,
                options=SolverOptions(method="RK4", step_size=0.001, execution_backend="native"),
            )
            reference = solve_hierarchy(
                python_target,
                options=SolverOptions(method="RK4", step_size=0.001, execution_backend="python"),
            )
            self.assertEqual(len(native.trajectory_results), len(reference.trajectory_results))
            for actual, expected in zip(native.trajectory_results, reference.trajectory_results):
                np.testing.assert_array_equal(actual.time, expected.time)
                np.testing.assert_allclose(actual.values, expected.values, rtol=3e-12, atol=3e-12)


class HierarchyFallbackTests(unittest.TestCase):
    def test_python_preference_keeps_scalar_solver_route(self):
        with tempfile.TemporaryDirectory() as directory:
            target = _ensemble(directory, "python-scalar")
            with patch("solver.hierarchy.integrate_fixed_native") as kernel:
                result = solve_hierarchy(
                    target,
                    options=SolverOptions(method="RK4", step_size=0.01, execution_backend="python"),
                )
            kernel.assert_not_called()
            self.assertTrue(all("batch_execution" not in item.metadata for item in result.trajectory_results))
            self.assertTrue(all(item.metadata["execution_backend"]["name"] == "python-reference" for item in result.trajectory_results))


if __name__ == "__main__":
    unittest.main()
