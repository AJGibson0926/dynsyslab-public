from __future__ import annotations

import unittest

import numpy as np

from solver import HistoryOptions, SolveRequest, SolverOptions, solve_ode_in_memory


class _Decay:
    name="in-memory-decay"
    initial_state=np.asarray((1.,))
    def rhs(self,t,state): return -np.asarray(state)


class InMemoryODESolverTests(unittest.TestCase):
    def test_runs_canonical_kernel_without_archive_publication(self):
        result=solve_ode_in_memory(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),
            options=SolverOptions(method="DOP853",rtol=1e-10,atol=1e-12)))
        self.assertEqual(result.status,"complete")
        self.assertIsNone(result.archive_path)
        self.assertIsNone(result.history_path)
        self.assertAlmostEqual(float(result.values[-1,0]),float(np.exp(-1.)),places=9)

    def test_rejects_durable_histories_and_postprocessors(self):
        with self.assertRaisesRegex(ValueError,"non-tracking"):
            solve_ode_in_memory(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),
                history=HistoryOptions(mode="sampled",track_during_solve=True)))
        with self.assertRaisesRegex(ValueError,"postprocessors"):
            solve_ode_in_memory(SolveRequest(system=_Decay(),kind="ode",time_span=(0.,1.),
                postprocessors=(lambda path:path,)))


if __name__=="__main__":
    unittest.main()
