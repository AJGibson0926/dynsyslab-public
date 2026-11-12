from __future__ import annotations

import unittest

import numpy as np

from scripts.generate_rpe_reference_training import (
    RPE_REFERENCE_ATOL,
    RPE_REFERENCE_DTAU,
    RPE_REFERENCE_HORIZON_TAU,
    RPE_REFERENCE_INITIAL_STATE_DIMENSIONLESS,
    RPE_REFERENCE_RTOL,
    build_request,
)


class RPEReferenceTrainingTests(unittest.TestCase):
    def test_reference_request_is_exact_dimensionless_bdf_configuration(self):
        request=build_request()
        system=request.system
        scale=system.dimensionless_parameters.time_scale_s
        self.assertEqual(system.preset.kappa,4/3)
        np.testing.assert_array_equal(RPE_REFERENCE_INITIAL_STATE_DIMENSIONLESS,
                                      np.asarray((2.3333865986017348,0.)))
        self.assertEqual(request.options.method,"BDF")
        self.assertEqual(request.options.rtol,RPE_REFERENCE_RTOL)
        self.assertEqual(request.options.atol,RPE_REFERENCE_ATOL)
        self.assertTrue(request.options.use_dimensionless_form)
        self.assertEqual(request.time_span,(0.,RPE_REFERENCE_HORIZON_TAU*scale))
        self.assertEqual(request.interpolation.time_step,RPE_REFERENCE_DTAU*scale)
        self.assertEqual(request.interpolation.method,"pchip")
        np.testing.assert_array_equal(system.to_dimensionless_state(request.initial_state),
                                      RPE_REFERENCE_INITIAL_STATE_DIMENSIONLESS)
        np.testing.assert_array_equal(system.dimensionless_initial_state,np.asarray((1.,0.)))


if __name__=="__main__":
    unittest.main()
