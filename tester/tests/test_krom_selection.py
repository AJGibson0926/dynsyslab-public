from __future__ import annotations

import unittest

import numpy as np

from postprocessing.koopman.spectrum import (DynamicalSystemReference,
    KoopmanEigenfunction,KoopmanEigenvalue,KoopmanSpectrum)
from rom import KROM, ObservableModeSelectionPolicy, select_observable_modes


class _ControlledPlane:
    name="controlled-plane"
    initial_state=np.zeros(2)
    state_dimension=2
    def input_matrix(self,state): return np.asarray(((0.,),(1.,)))


def _complex_coordinate_krom() -> KROM:
    system=_ControlledPlane(); reference=DynamicalSystemReference(None,system.name,{})
    values=[KoopmanEigenvalue.from_continuous(-.1+1j,.01,reference),
            KoopmanEigenvalue.from_continuous(-.1-1j,.01,reference)]
    positive=np.asarray((0.,1.,1j)); negative=np.conjugate(positive)
    functions=[
        KoopmanEigenfunction(values[0],"polynomial",coefficients=positive,
            state_callable=lambda x:complex(x[0],x[1]),gradient_callable=lambda x:np.asarray((1.,1j))),
        KoopmanEigenfunction(values[1],"polynomial",coefficients=negative,
            state_callable=lambda x:complex(x[0],-x[1]),gradient_callable=lambda x:np.asarray((1.,-1j))),
    ]
    metadata={"method":"selection-test","dictionary":{"encoding":"polynomial-dictionary/v1","degree":1,"include_constant":True},
        "state_preprocessing":{"center":[0.,0.],"scale":[1.,1.],"multiplier":[1.,1.]},
        "preprocessing":{"coordinate_system":"dimensional"}}
    return KROM(KoopmanSpectrum(reference,values,functions,[],metadata,dynamical_system=system),system)


class KROMSelectionTests(unittest.TestCase):
    def test_complex_mode_is_selected_as_one_group_and_reconstructs_both_states(self):
        axis=np.linspace(-2.,2.,41)
        states=np.asarray([(x,y) for x in axis for y in axis[::4]],dtype=float)
        report=select_observable_modes(_complex_coordinate_krom(),states,
            ObservableModeSelectionPolicy(r2_threshold=.999999,max_modes=3,control_authority_weight=0.))
        self.assertEqual(report.selected_retained_indices,(0,))
        self.assertEqual(report.selected_full_spectrum_indices,(0,))
        np.testing.assert_allclose(report.final_r2_by_state,np.ones(2),atol=1e-13)
        self.assertEqual(report.stop_reason,"r2_threshold")
        self.assertEqual(report.steps[0].continuous_eigenvalue,-.1+1j)

    def test_selection_policy_validation_and_report_are_serializable(self):
        with self.assertRaises(ValueError): ObservableModeSelectionPolicy(r2_threshold=0.)
        with self.assertRaises(ValueError): ObservableModeSelectionPolicy(max_modes=0)
        states=np.asarray(((-1.,-2.),(0.,0.),(1.,2.),(2.,-1.)))
        summary=select_observable_modes(_complex_coordinate_krom(),states).summary()
        self.assertEqual(summary["method"],"grouped-polynomial-state-observability/v1")
        self.assertEqual(summary["candidate_count"],1)
        self.assertIsInstance(summary["steps"][0]["continuous_eigenvalue"],dict)


if __name__=="__main__":
    unittest.main()
