from __future__ import annotations

import unittest

import numpy as np

from control import rpe_control_target
from rom import AxisAlignedStateDomain, plan_state_domain_expansion


class KROMTrainingCoverageTests(unittest.TestCase):
    def test_harmonic_target_envelope_contains_dense_samples(self):
        target=rpe_control_target("rpe_multisine_offset_4_periodic")
        minimum,maximum=target.conservative_state_envelope()
        time=np.linspace(0.,20.,10001)
        values=np.vstack([target.evaluate(value,2) for value in time])
        self.assertTrue(np.all(values>=minimum-1e-12))
        self.assertTrue(np.all(values<=maximum+1e-12))

    def test_expansion_margin_applies_only_to_axes_that_need_expansion(self):
        current=AxisAlignedStateDomain((.2,-30.),(2.3,30.),"dimensionless")
        required=AxisAlignedStateDomain((.5,-20.),(7.,20.),"dimensionless")
        plan=plan_state_domain_expansion(current,required,points_per_axis=(5,5),margin_fraction=.1)
        self.assertEqual(plan.planned.minimum,current.minimum)
        self.assertEqual(plan.planned.maximum[1],current.maximum[1])
        self.assertGreater(plan.planned.maximum[0],required.maximum[0])
        self.assertTrue(plan.summary()["required_covered"])
        points=np.asarray(plan.initial_conditions)
        self.assertTrue(np.all(~current.contains(points)))

    def test_anchors_are_preserved_even_when_already_in_current_domain(self):
        current=AxisAlignedStateDomain((0.,0.),(1.,1.))
        required=AxisAlignedStateDomain((0.,0.),(2.,2.))
        plan=plan_state_domain_expansion(current,required,points_per_axis=(2,2),anchors=((.5,.5),))
        self.assertIn((.5,.5),plan.initial_conditions)

    def test_cartesian_explosion_is_gated(self):
        domain=AxisAlignedStateDomain((0.,0.,0.),(1.,1.,1.))
        with self.assertRaisesRegex(ValueError,"100,000"):
            plan_state_domain_expansion(domain,domain,points_per_axis=(100,100,100),outside_current_only=False)


if __name__=="__main__":
    unittest.main()
