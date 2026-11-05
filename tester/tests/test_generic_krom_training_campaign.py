from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from campaigner import (Campaigner, StateEnvelopeRequirement, StateProbe,
    design_krom_training_plan, target_state_requirement)
from control import CallableTarget
from dynsyslib import ZeroDriver, load_system
from operations.krom_training_screening import (KROMTrainingScreeningPolicy,
    publish_screened_krom_training_campaign, screen_krom_training_plan)
from rom import AxisAlignedStateDomain
from solver import NumericalResult
from time_series import TimeSeries


def _krom(system, minimum, maximum, coordinates):
    return SimpleNamespace(system=system, coordinate_system=coordinates,
        training_domain={"minimum": list(minimum), "maximum": list(maximum),
            "coordinate_system": coordinates})


def _complete(request):
    initial = np.asarray(request.initial_state, dtype=float)
    return NumericalResult(TimeSeries(np.asarray((0., 1.)), np.vstack((initial, initial))),
        "complete", "solve completed")


class GenericKROMTrainingCampaignTests(unittest.TestCase):
    def test_duffing_and_lorenz_use_the_same_domain_planner(self):
        duffing = load_system("duffing", driver=ZeroDriver(), preset="double_well_reference")
        duffing_krom = _krom(duffing, (-.5, -.5), (.5, .5), "dimensionless")
        with patch("campaigner.krom_training.load_krom_archive", return_value=duffing_krom):
            plan = design_krom_training_plan(source_krom="duffing.npz",
                required_domain=AxisAlignedStateDomain((-2., -1.), (2., 1.), "dimensionless"),
                points_per_axis=(3, 3), preset="double_well_reference")
        self.assertEqual(plan.system_name, "duffing")
        self.assertEqual(plan.coverage.planned.dimension, 2)
        self.assertEqual(plan.physical_time_scale, 1.)

        lorenz = load_system("lorenz", driver=ZeroDriver(), preset="lorenz_1963")
        lorenz_krom = _krom(lorenz, (-5., -5., 0.), (5., 5., 15.), "dimensional")
        with patch("campaigner.krom_training.load_krom_archive", return_value=lorenz_krom):
            plan = design_krom_training_plan(source_krom="lorenz.npz",
                required_domain=AxisAlignedStateDomain((-20., -25., 0.), (20., 25., 50.),
                    "dimensional"), points_per_axis=(3, 3, 3), preset="lorenz_1963")
        self.assertEqual(plan.system_name, "lorenz")
        self.assertEqual(plan.coverage.planned.dimension, 3)
        self.assertNotIn("rpe", str(plan.summary()).lower())

    def test_arbitrary_target_becomes_an_envelope_anchor_and_probe_set(self):
        target = CallableTarget(lambda time: np.asarray((time, -time, time * time)))
        requirement = target_state_requirement("flight-path", target,
            sample_parameters=(0., 1., 2.), state_dimension=3,
            coordinate_system="dimensional")
        self.assertEqual(requirement.domain.minimum, (0., -2., 0.))
        self.assertEqual(requirement.domain.maximum, (2., 0., 4.))
        self.assertEqual(len(requirement.probes), 3)
        self.assertEqual(requirement.anchors, ((0., 0., 0.),))

    def test_generic_screen_and_publish_preserve_system_parameters(self):
        system = load_system("lorenz", driver=ZeroDriver(), preset="lorenz_1963")
        source = _krom(system, (-2., -2., 0.), (2., 2., 5.), "dimensional")
        requirement = StateEnvelopeRequirement("test-region",
            AxisAlignedStateDomain((-3., -3., 0.), (3., 3., 6.), "dimensional"),
            anchors=((3., 3., 6.),), probes=(StateProbe("test-region", 0., (3., 3., 6.)),))
        with patch("campaigner.krom_training.load_krom_archive", return_value=source):
            plan = design_krom_training_plan(source_krom="lorenz.npz",
                requirements=(requirement,), points_per_axis=(3, 3, 3), preset="lorenz_1963")
        with patch("operations.krom_training_screening.load_krom_archive", return_value=source), \
             patch("solver.solve_ode_in_memory", side_effect=_complete):
            report = screen_krom_training_plan(plan,
                policy=KROMTrainingScreeningPolicy(maximum_absolute_state=100.))
        self.assertTrue(report.publication_ready)
        with tempfile.TemporaryDirectory() as directory, \
             patch("operations.krom_training_screening.load_krom_archive", return_value=source), \
             patch("campaigner.designer.Campaigner", side_effect=lambda: Campaigner(directory)):
            campaign, collection = publish_screened_krom_training_campaign(plan, report)
        spec = collection.jobs[0].systems[0]
        self.assertEqual(spec["system_name"], "lorenz")
        self.assertEqual(spec["preset"], "lorenz_1963")
        self.assertEqual(spec["parameters"]["rho"], 28.)
        self.assertEqual(campaign.metadata["purpose"], "krom-training-domain-expansion")


if __name__ == "__main__":
    unittest.main()
