from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from dynsyslib.observables import Observable
from estimation.koopman import (KLQEEstimator, design_klqe, design_klqe_from_krom,
                                load_klqe_design, load_klqe_result, run_klqe,
                                save_klqe_design, save_klqe_result)
from rom.observation import ObservableReconstruction


class _RealEigenfunction:
    is_real=True


class _ScalarKROM:
    lifted_dimension=1
    coordinate_system="dimensionless"
    eigenfunctions=(_RealEigenfunction(),)
    system=SimpleNamespace(initial_state=np.array((0.,)))

    def A(self,*,output="complex"):
        return np.array([[-1.]]) if output=="split" else np.array([[-1.+0j]])

    def psi(self,state,*,output="complex",state_coordinates="model"):
        value=float(np.asarray(state).reshape(-1)[0])
        return np.array((value,)) if output=="split" else np.array((value+0j,))


def _reconstruction(observable):
    return ObservableReconstruction(observable,np.array(((1.+0j,),)),np.zeros(1),
        {"passed":True},{"method":"unit-test"})


class KoopmanKLQETests(unittest.TestCase):
    def test_observable_driven_design_is_stable_and_satisfies_care(self):
        rom=_ScalarKROM(); observable=Observable.identity(1)
        design=design_klqe(rom,_reconstruction(observable),_reconstruction(observable),
                           np.array([[1.]]),np.array([[.25]]))
        self.assertTrue(design.detectable)
        self.assertLess(design.care_residual_norm,1e-10)
        self.assertTrue(np.all(np.real(design.observer_poles)<0.))
        self.assertEqual(design.observation_model.measurement_observable.name,"state")

    def test_exact_zoh_propagation_matches_closed_form_constant_measurement(self):
        rom=_ScalarKROM(); gain=np.array([[2.]]); covariance=np.array([[1.]])
        estimator=KLQEEstimator(rom,np.array([[1.]]),np.array([[1.]]),np.array([[1.]]),
            np.array((0.,)),lambda z:z,gain_override=gain,covariance_override=covariance)
        estimator.step_measurement(0.,np.array((2.,)))
        estimator.step_measurement(.5,np.array((2.,)))
        expected=4./3.*(1.-np.exp(-1.5))
        np.testing.assert_allclose(estimator.lifted_state,np.array((expected,)),rtol=1e-13)

    def test_offline_runner_records_estimate_measurement_innovation_and_truth(self):
        rom=_ScalarKROM(); observable=Observable.identity(1)
        design=design_klqe(rom,_reconstruction(observable),_reconstruction(observable),
                           np.array([[1.]]),np.array([[.25]]))
        time=np.linspace(0.,1.,11); truth=np.exp(-time)[:,None]
        result=run_klqe(design.estimator(np.array((0.,))),time,truth,true_states=truth)
        self.assertEqual(result.state_estimate.values.shape,(11,1))
        self.assertEqual(result.innovation.values.shape,(11,1))
        self.assertIsNotNone(result.state_error)
        with tempfile.TemporaryDirectory() as directory:
            loaded=load_klqe_result(save_klqe_result(Path(directory)/"klqe.npz",result))
        np.testing.assert_allclose(loaded.state_estimate.values,result.state_estimate.values)
        self.assertEqual(loaded.summary()["sample_count"],11)

    def test_design_artifact_verifies_fingerprint_and_reconstructed_matrices(self):
        rom=_ScalarKROM(); observable=Observable.identity(1); reconstruction=_reconstruction(observable)
        rom.observable_reconstructions={"measurement":reconstruction,"state":reconstruction}
        rom.observable_map=lambda name:rom.observable_reconstructions[name]
        design=design_klqe_from_krom(rom,"measurement","state",np.array([[1.]]),np.array([[.25]]))
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); krom_archive=root/"krom.npz"; krom_archive.write_bytes(b"frozen-krom")
            artifact=save_klqe_design(root/"design.npz",design,krom_archive=krom_archive,
                                      measurement_name="measurement",state_output_name="state")
            with patch("rom.load_krom_archive",return_value=rom):
                restored,metadata=load_klqe_design(artifact)
        np.testing.assert_allclose(restored.gain,design.gain)
        self.assertEqual(metadata["measurement_name"],"measurement")


if __name__=="__main__":
    unittest.main()
