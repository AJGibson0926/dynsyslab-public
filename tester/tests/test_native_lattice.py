import importlib.util
import unittest
from unittest.mock import patch

import numpy as np

from dynsyslib.lattice import IsingModel, PottsModel
from solver.lattice import solve_lattice
from solver.models import SolveRequest, SolverOptions


class LatticeReferenceTests(unittest.TestCase):
    def test_potts_energy_counts_both_positive_direction_bonds(self):
        system = PottsModel(size=4, states=3, coupling=1.0, field=0.0, seed=1)
        system.spins.fill(0)
        self.assertEqual(system.energy() / system.spins.size, -2.0)


@unittest.skipUnless(importlib.util.find_spec("_dynsyslab_native"), "optional native wheel is not installed")
class NativeLatticeTests(unittest.TestCase):
    def test_native_ising_is_reproducible_and_updates_model_state(self):
        first = IsingModel(size=8, seed=123)
        second = IsingModel(size=8, seed=123)
        options = SolverOptions(execution_backend="native")
        first_result = solve_lattice(SolveRequest(first, "lattice", time_span=(0.0, 20.0), options=options), None)
        second_result = solve_lattice(SolveRequest(second, "lattice", time_span=(0.0, 20.0), options=options), None)
        self.assertEqual(first_result.status, "complete")
        np.testing.assert_array_equal(first_result.values, second_result.values)
        np.testing.assert_array_equal(first.spins, second.spins)
        self.assertTrue(np.all(np.isin(first.spins, (-1, 1))))
        self.assertAlmostEqual(first_result.values[-1, 0], first.magnetization())
        self.assertAlmostEqual(first_result.values[-1, 1], first.energy() / first.spins.size)
        self.assertEqual(first_result.metadata["lattice_backend"]["random_engine"], "std::mt19937_64")
        self.assertEqual(first_result.metadata["lattice_backend"]["trajectory_parity"], "statistical-not-bitwise")

    def test_auto_does_not_change_the_established_numpy_random_stream(self):
        system = IsingModel(size=6, seed=7)
        with patch("solver.lattice.run_ising2d_native") as native:
            result = solve_lattice(SolveRequest(system, "lattice", time_span=(0.0, 2.0)), None)
        native.assert_not_called()
        self.assertEqual(result.metadata["execution_backend"]["name"], "python-reference")

    def test_explicit_native_rejects_models_without_an_adapter(self):
        request = SolveRequest(PottsModel(size=6, seed=1), "lattice", time_span=(0.0, 2.0),
                               options=SolverOptions(execution_backend="native"))
        with self.assertRaisesRegex(RuntimeError, "no native lattice adapter"):
            solve_lattice(request, None)


if __name__ == "__main__":
    unittest.main()
