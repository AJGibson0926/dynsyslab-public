from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from campaigner.rpe_krom import design_rpe_krom_training_plan
from operations.rpe_krom_screening import (RPEScreeningPolicy,
    publish_screened_rpe_krom_campaign, save_rpe_krom_screening_report,
    screen_rpe_krom_training_plan)
from solver import NumericalResult
from time_series import TimeSeries


class _IdentityRPE:
    name="rpe"
    initial_state=np.asarray((1.,0.))
    dimensionless_parameters=SimpleNamespace(time_scale_s=2.)
    def to_dimensional_state(self,state): return np.asarray(state,dtype=float)
    def to_dimensionless_state(self,state): return np.asarray(state,dtype=float)
    def to_dimensional_time(self,time): return 2.*np.asarray(time,dtype=float)


def _krom():
    return SimpleNamespace(training_domain={"minimum":[.2,-30.],"maximum":[2.4,30.],
        "coordinate_system":"dimensionless"},coordinate_system="dimensionless",system=_IdentityRPE())


def _complete(request):
    initial=np.asarray(request.initial_state,dtype=float)
    values=np.vstack((initial,initial))
    return NumericalResult(TimeSeries(np.asarray((0.,1.)),values),"complete","solve completed")


class RPEKROMCampaignGateTests(unittest.TestCase):
    def _plan(self):
        with patch("campaigner.krom_training.load_krom_archive",return_value=_krom()):
            return design_rpe_krom_training_plan(source_krom="synthetic.npz",
                target_names=("rpe_multisine_equilibrium_small",),points_per_axis=(3,3),
                target_probes_per_target=9)

    def test_safe_target_probes_make_filtered_plan_publication_ready(self):
        plan=self._plan()
        with patch("operations.krom_training_screening.load_krom_archive",return_value=_krom()), \
             patch("solver.solve_ode_in_memory",side_effect=_complete):
            report=screen_rpe_krom_training_plan(plan)
        self.assertTrue(report.publication_ready)
        self.assertEqual(report.rejected_target_probe_count,0)
        self.assertEqual(report.target_summary["rpe_multisine_equilibrium_small"]["probe_count"],9)

    def test_unsafe_target_probe_blocks_publication_even_if_box_is_covered(self):
        plan=self._plan()
        def collapse_inward(request):
            initial=np.asarray(request.initial_state,dtype=float)
            final=np.asarray((.01,initial[1])) if initial[1]<0. else initial
            return NumericalResult(TimeSeries(np.asarray((0.,1.)),np.vstack((initial,final))),
                "complete","solve completed")
        with patch("operations.krom_training_screening.load_krom_archive",return_value=_krom()), \
             patch("solver.solve_ode_in_memory",side_effect=collapse_inward):
            report=screen_rpe_krom_training_plan(plan,policy=RPEScreeningPolicy(minimum_radius=.05))
        self.assertFalse(report.publication_ready)
        self.assertGreater(report.rejected_target_probe_count,0)
        with self.assertRaisesRegex(ValueError,"publication-ready"):
            publish_screened_rpe_krom_campaign(plan,report)

    def test_screening_record_is_content_addressed_and_idempotent(self):
        plan=self._plan()
        with patch("operations.krom_training_screening.load_krom_archive",return_value=_krom()), \
             patch("solver.solve_ode_in_memory",side_effect=_complete):
            report=screen_rpe_krom_training_plan(plan)
        with tempfile.TemporaryDirectory() as directory:
            first=save_rpe_krom_screening_report(plan,report,root=directory)
            second=save_rpe_krom_screening_report(plan,report,root=directory)
            self.assertEqual(first,second)
            self.assertIn(report.fingerprint,first.name)
            self.assertEqual(len(list(Path(directory).glob("*.json"))),1)


if __name__=="__main__":
    unittest.main()
