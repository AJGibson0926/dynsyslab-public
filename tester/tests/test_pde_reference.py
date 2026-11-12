import unittest

import numpy as np

from solver.models import HistoryOptions, SolveRequest, SolverOptions
from solver.pde import solve_pde


class _DecayField:
    name = "decay-field"
    def initial_field(self, mesh): return np.arange(1.0, 7.0).reshape(2, 3)
    def rhs(self, time, field, mesh): return -field


class _Session:
    def __init__(self): self.checkpoints=[]
    def checkpoint(self, time, values): self.checkpoints.append((time,values))


class PDEReferenceTests(unittest.TestCase):
    def test_rk4_fixed_step_preserves_field_shape_and_accuracy(self):
        request=SolveRequest(_DecayField(),"pde",time_span=(0.0,1.0),mesh=object(),
                             options=SolverOptions(method="RK4",step_size=0.01,execution_backend="python"))
        result=solve_pde(request,_Session())
        self.assertEqual(result.status,"complete")
        self.assertEqual(result.values.shape,(101,2,3))
        np.testing.assert_allclose(result.values[-1],_DecayField().initial_field(None)*np.exp(-1),rtol=2e-9,atol=2e-9)

    def test_fixed_step_history_emits_one_checkpoint_per_completed_step(self):
        session=_Session()
        request=SolveRequest(_DecayField(),"pde",time_span=(0.0,0.2),mesh=object(),
                             options=SolverOptions(method="Euler",step_size=0.05,execution_backend="python"),
                             history=HistoryOptions(track_during_solve=True))
        result=solve_pde(request,session)
        self.assertEqual(result.status,"complete")
        self.assertEqual(len(session.checkpoints),4)
        self.assertTrue(np.all(np.diff(result.time)>0))

    def test_segmented_adaptive_history_has_no_duplicate_boundaries(self):
        session=_Session()
        request=SolveRequest(_DecayField(),"pde",time_span=(0.0,0.2),mesh=object(),
                             options=SolverOptions(method="RK45",execution_backend="python"),
                             history=HistoryOptions(track_during_solve=True,checkpoint_interval=0.05))
        result=solve_pde(request,session)
        self.assertEqual(result.status,"complete")
        self.assertTrue(np.all(np.diff(result.time)>0))
        self.assertEqual(len(session.checkpoints),4)

    def test_fixed_step_requires_a_positive_step_size(self):
        request=SolveRequest(_DecayField(),"pde",time_span=(0.0,1.0),mesh=object(),
                             options=SolverOptions(method="RK4",execution_backend="python"))
        result=solve_pde(request,_Session())
        self.assertEqual(result.status,"failed")
        self.assertIn("step_size",result.message)


if __name__=="__main__":
    unittest.main()
