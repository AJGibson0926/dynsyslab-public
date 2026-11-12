from __future__ import annotations

import unittest
import json
import tempfile
from pathlib import Path

import numpy as np

from postprocessing.koopman.spectrum import (DynamicalSystemReference,
    KoopmanEigenfunction,KoopmanEigenvalue,KoopmanSpectrum)
from rom import (KROM,ObservableReconstruction,ObservationMapQualityPolicy,fit_state_observation_map,
                 load_krom_archive,split_lifted_coordinates)
from dynsyslib.observables import Observable


class _Plane:
    name="observation-plane"
    initial_state=np.zeros(2)
    state_dimension=2
    def input_matrix(self,state): return np.asarray(((0.,),(1.,)))


def _complex_krom(*,observable_velocity: bool=True) -> KROM:
    system=_Plane(); reference=DynamicalSystemReference(None,system.name,{})
    values=[KoopmanEigenvalue.from_continuous(-.1+1j,.01,reference),
            KoopmanEigenvalue.from_continuous(-.1-1j,.01,reference)]
    positive=np.asarray((0.,1.,1j if observable_velocity else 0j)); negative=np.conjugate(positive)
    def positive_psi(x): return complex(x[0],x[1] if observable_velocity else 0.)
    def negative_psi(x): return np.conjugate(positive_psi(x))
    functions=[KoopmanEigenfunction(values[0],"polynomial",coefficients=positive,
            state_callable=positive_psi,gradient_callable=lambda x:np.asarray((1.,1j if observable_velocity else 0j))),
        KoopmanEigenfunction(values[1],"polynomial",coefficients=negative,
            state_callable=negative_psi,gradient_callable=lambda x:np.asarray((1.,-1j if observable_velocity else 0j)))]
    return KROM(KoopmanSpectrum(reference,values,functions,[],{"method":"observation-test"},
        dynamical_system=system),system)


class KROMObservationMapTests(unittest.TestCase):
    def setUp(self):
        axis=np.linspace(-2.,2.,31)
        self.states=np.asarray([(x,y) for x in axis for y in axis[::3]],dtype=float)

    def test_complex_split_encoding_is_explicit(self):
        model=_complex_krom()
        lifted=np.asarray(((1.+2j,),(3.-4j,)))
        np.testing.assert_array_equal(split_lifted_coordinates(lifted,model.eigenfunctions),
            np.asarray(((1.,2.),(3.,-4.))))

    def test_fit_validates_and_enables_state_reconstruction(self):
        source=_complex_krom()
        observation=fit_state_observation_map(source,self.states,
            quality_policy=ObservationMapQualityPolicy(minimum_r2=.999999,
                maximum_normalized_rmse=1e-10,maximum_condition_number=1e6))
        self.assertTrue(observation.validated)
        model=source.with_observation_map(observation)
        self.assertEqual(model.state_representation,"learned-state-observation-map")
        self.assertTrue(model.capabilities()["state_reconstruction"])
        self.assertTrue(model.capabilities()["validated_state_reconstruction"])
        state=np.asarray((.75,-1.25)); lifted=model.psi(state)
        np.testing.assert_allclose(model.reconstruct(lifted),state,atol=2e-14)
        model.require_capability("validated_state_reconstruction")
        with self.assertRaisesRegex(ValueError,"no complex Koopman mode matrix"):
            model.mode_matrix()

    def test_failed_map_requires_an_explicit_diagnostic_override(self):
        source=_complex_krom(observable_velocity=False)
        with self.assertRaisesRegex(ValueError,"failed held-out quality"):
            fit_state_observation_map(source,self.states)
        observation=fit_state_observation_map(source,self.states,require_quality=False)
        self.assertFalse(observation.validated)
        model=source.with_observation_map(observation)
        self.assertTrue(model.capabilities()["state_reconstruction"])
        self.assertFalse(model.capabilities()["validated_state_reconstruction"])
        with self.assertRaisesRegex(ValueError,"validated_state_reconstruction"):
            model.require_capability("validated_state_reconstruction")

    def test_mode_extraction_invalidates_a_decoder_that_no_longer_aligns(self):
        source=_complex_krom()
        observed=source.with_observation_map(fit_state_observation_map(source,self.states))
        extracted=observed.extract([0],include_conjugate_duals=False)
        self.assertIsNone(extracted.observation_map)
        self.assertFalse(extracted.capabilities()["state_reconstruction"])

    def test_named_observable_reconstructions_coexist_and_follow_extraction(self):
        source=_complex_krom()
        state=Observable.identity(2)
        state_map=ObservableReconstruction(state,np.asarray(((.5-.0j,),(-.0-.5j,))),
            np.zeros(2),{"passed":True},{"method":"unit-test"})
        sum_observable=Observable("x_plus_y",1,2,lambda x:np.asarray((x[0]+x[1],)),
            lambda x:np.asarray(((1.,1.),)),{"kind":"callable"})
        sum_map=ObservableReconstruction(sum_observable,np.asarray(((.5-.5j,),)),
            np.zeros(1),{"passed":True},{"method":"unit-test"})
        model=(source.with_observable_reconstruction(state_map,name="state")
                     .with_observable_reconstruction(sum_map,name="x_plus_y",make_default=False))
        self.assertEqual(sorted(model.observable_reconstructions),["state","x_plus_y"])
        self.assertIs(model.observable_map("state"),model.observable_reconstruction)
        np.testing.assert_allclose(model.C("state"),np.eye(2))
        np.testing.assert_allclose(model.C("x_plus_y"),np.asarray(((1.,1.),)))
        extracted=model.extract([0],include_conjugate_duals=False)
        self.assertEqual(sorted(extracted.observable_reconstructions),["state","x_plus_y"])
        np.testing.assert_allclose(extracted.C("x_plus_y"),model.C("x_plus_y"))

    def test_archived_observation_map_round_trips_as_an_executable_capability(self):
        from dynsyslib import load_system
        from archives.store import serializable_system
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); primary=root/"solution.npz"; archive=root/"observed-krom.npz"
            solution_metadata={"format":"dynsyslab.solution/v1","system":"rayleigh-plesset",
                "system_parameters":serializable_system(load_system("rpe",preset="air_water_20c")),
                "delta_t":.01,"integration_form":"dimensionless"}
            np.savez(primary,metadata_json=np.asarray(json.dumps(solution_metadata)))
            quality={"method":"held-out-affine-state-observation/v1","passed":True,
                "validation":{"r2_by_state":[1.,1.]}}
            observation={"encoding":"real-imaginary-split-affine/v1","coordinate_system":"model",
                "lifted_feature_dimension":2,"state_dimension":2,"validated":True,"quality":quality}
            def reconstruction_summary(observable):
                return {"encoding":"nonredundant-complex-koopman-modes/v1",
                    "observable":observable.summary(),"observable_dimension":observable.output_dimension,
                    "lifted_dimension":1,"validated":True,"quality":{"passed":True},
                    "provenance":{"method":"unit-test"}}
            radius=Observable.component(2,0,name="radius"); state=Observable.identity(2)
            named={"radius":reconstruction_summary(radius),"state":reconstruction_summary(state)}
            analysis={"method":"koopman-with-sindy","primary_source":str(primary),
                "dictionary":{"encoding":"polynomial-dictionary/v1","degree":1,"include_constant":True},
                "state_preprocessing":{"center":[0.,0.],"scale":[1.,1.],"multiplier":[1.,1.]},
                "preprocessing":{"coordinate_system":"dimensionless"},"observation_map":observation,
                "observable_reconstruction":named["radius"],
                "krom":{"delta_t":.01,"observable_reconstruction":named["radius"],
                        "observable_reconstructions":named,
                        "default_observable_reconstruction":"radius"}}
            discrete=np.exp(np.asarray((-.1+1j,-.1-1j))*.01)
            coefficients=np.asarray(((0.,1.,1j),(0.,1.,-1j)),dtype=complex)
            states=np.asarray(((-1.,-2.),(0.,0.),(1.,2.)))
            np.savez(archive,eigenvalues=discrete,coefficients=coefficients,training_states=states,
                observation_coefficients=np.eye(2),observation_intercept=np.zeros(2),
                observable_reconstruction_coefficients=np.asarray(((.5+0j,),)),
                observable_reconstruction_intercept=np.zeros(1),
                observable_reconstruction_000_coefficients=np.asarray(((.5+0j,),)),
                observable_reconstruction_000_intercept=np.zeros(1),
                observable_reconstruction_001_coefficients=np.asarray(((.5+0j,),(-.5j,))),
                observable_reconstruction_001_intercept=np.zeros(2),
                metadata_json=np.asarray(json.dumps({"format":"dynsyslab.derived/v1","analysis_metadata":analysis})))
            model=load_krom_archive(archive)
            self.assertTrue(model.capabilities()["validated_state_reconstruction"])
            probe=np.asarray((.75,-1.25))
            np.testing.assert_allclose(model.reconstruct(model.psi(probe)),probe,atol=1e-14)
            self.assertEqual(sorted(model.observable_reconstructions),["radius","state"])
            np.testing.assert_allclose(model.C("radius"),np.asarray(((1.,0.),)))
            np.testing.assert_allclose(model.C("state"),np.eye(2))


if __name__=="__main__":
    unittest.main()
