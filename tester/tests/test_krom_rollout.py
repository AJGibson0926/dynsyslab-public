import importlib
import unittest
from types import SimpleNamespace

import numpy as np

from postprocessing.koopman.spectrum import (
    DynamicalSystemReference,
    KoopmanEigenfunction,
    KoopmanEigenvalue,
    KoopmanSpectrum,
)
from rom import KROM, extract_krom
from control import ConstantTarget, KoopmanLQRController
from archives.store import serializable_system


class _LinearSystem:
    initial_state = np.asarray((1.0, 1.0))
    state_dimension = 2

    def __init__(self,input_sign=1.): self.input_sign=float(input_sign)

    def input_matrix(self, state):
        return self.input_sign*np.ones((2, 1))


def _krom() -> KROM:
    system = _LinearSystem()
    reference = DynamicalSystemReference(None, "test-linear", {})
    eigenvalues = [
        KoopmanEigenvalue.from_continuous(-1.0, 0.1, reference),
        KoopmanEigenvalue.from_continuous(-2.0, 0.1, reference),
    ]
    functions = [
        KoopmanEigenfunction(eigenvalues[0], "coordinate", state_callable=lambda state: state[0]),
        KoopmanEigenfunction(eigenvalues[1], "coordinate", state_callable=lambda state: state[1]),
    ]
    spectrum = KoopmanSpectrum(
        reference,
        eigenvalues,
        functions,
        [np.asarray((1.0, 0.0)), np.asarray((0.0, 1.0))],
        {"method": "test"},
        dynamical_system=system,
    )
    return KROM(spectrum, system)


class _ScaledLinearSystem(_LinearSystem):
    initial_state=np.asarray((2.,4.))
    dimensionless_parameters=SimpleNamespace(time_scale_s=.5)
    def to_dimensionless_state(self,state): return np.asarray(state,dtype=float)/np.asarray((2.,4.))
    def to_dimensional_state(self,state): return np.asarray(state,dtype=float)*np.asarray((2.,4.))
    def dimensionless_input_matrix(self,state): return np.asarray(((.25,),(.5,)))
    def input_matrix(self,state): raise AssertionError("dimensionless KROM must not use the physical-state input matrix")


class _UnactuatedSystem:
    initial_state=np.asarray((1.,1.))
    state_dimension=2
    def input_matrix(self,state): return None


def _scaled_krom(*,delay_modes:bool=False) -> KROM:
    system=_ScaledLinearSystem(); reference=DynamicalSystemReference(None,"scaled-linear",{})
    eigenvalues=[KoopmanEigenvalue.from_continuous(-1.,.1,reference),KoopmanEigenvalue.from_continuous(-2.,.1,reference)]
    functions=[KoopmanEigenfunction(eigenvalues[0],"coordinate",state_callable=lambda state:state[0],gradient_callable=lambda state:np.asarray((1.,0.))),
               KoopmanEigenfunction(eigenvalues[1],"coordinate",state_callable=lambda state:state[1],gradient_callable=lambda state:np.asarray((0.,1.)))]
    modes=[np.asarray((1.,0.,0.,0.)),np.asarray((0.,1.,0.,0.))] if delay_modes else [np.asarray((1.,0.)),np.asarray((0.,1.))]
    metadata={"method":"test","preprocessing":{"coordinate_system":"dimensionless"}}
    if delay_modes: metadata["delay_embedding"]={"type":"hankel"}
    return KROM(KoopmanSpectrum(reference,eigenvalues,functions,modes,metadata,dynamical_system=system),system)


def _conjugate_krom() -> KROM:
    system=_LinearSystem(); reference=DynamicalSystemReference(None,"complex-linear",{})
    continuous=(-.2+0j,-.5+2j,-.5-2j)
    eigenvalues=[KoopmanEigenvalue.from_continuous(value,.1,reference) for value in continuous]
    functions=[KoopmanEigenfunction(eigenvalues[0],"coordinate",state_callable=lambda state:state[0]),
        KoopmanEigenfunction(eigenvalues[1],"coordinate",values=np.asarray((1+1j,2+2j)),state_callable=lambda state:state[1]+1j*state[0]),
        KoopmanEigenfunction(eigenvalues[2],"coordinate",values=np.asarray((1-1j,2-2j)),state_callable=lambda state:state[1]-1j*state[0])]
    modes=[np.asarray((1.,0.)),np.asarray((0.,1.)),np.asarray((0.,1.))]
    return KROM(KoopmanSpectrum(reference,eigenvalues,functions,modes,{"method":"complex-test"},dynamical_system=system),system)


class KROMRolloutTests(unittest.TestCase):
    def test_analytic_control_system_rebinding_is_independent_and_provenance_recorded(self):
        source=_krom(); rebound=source.with_analytic_control_system(_LinearSystem(-1.),reason="correct input sign")
        self.assertIsNot(rebound,source)
        self.assertEqual(source.system.input_sign,1.)
        self.assertEqual(rebound.system.input_sign,-1.)
        binding=rebound.spectrum.metadata["control_system_binding"]
        self.assertEqual(binding["reason"],"correct input sign")
        self.assertTrue(binding["spectral_source_unchanged"])
        source.control_mode="linear"; source.B1=np.ones((2,1)); source._validate_control_representation()
        with self.assertRaisesRegex(ValueError,"analytic-control"):
            source.with_analytic_control_system(_LinearSystem(-1.),reason="unsafe fitted matrix")

    def test_extract_krom_expands_conjugate_dual_and_preserves_source(self):
        source=_conjugate_krom(); source_ids=[id(value) for value in source.spectrum.eigenvalues]
        positive=next(index for index,value in enumerate(source.eigenvalues) if value.is_complex)
        derived=source.extract([positive])
        self.assertEqual(len(derived.spectrum.eigenvalues),2)
        self.assertEqual(derived.lifted_dimension,1)
        self.assertIsNot(derived.spectrum.eigenvalues[0],source.spectrum.eigenvalues[1])
        self.assertEqual([id(value) for value in source.spectrum.eigenvalues],source_ids)
        self.assertIs(source.spectrum.eigenvalues[1].dual,source.spectrum.eigenvalues[2])
        self.assertEqual(derived.spectrum.metadata["extraction"]["selected_source_indices"],[1,2])
        derived.spectrum.eigenfunctions[0].values[0]=99.
        self.assertNotEqual(source.spectrum.eigenfunctions[1].values[0],99.)

    def test_extract_krom_accepts_exact_eigenfunction_and_can_omit_dual(self):
        source=_conjugate_krom(); selected=source.spectrum.eigenfunctions[2]
        derived=extract_krom(source,[selected],include_conjugate_duals=False)
        self.assertEqual(len(derived.spectrum.eigenvalues),1)
        self.assertEqual(derived.spectrum.eigenvalues[0].continuous_time,selected.eigenvalue.continuous_time)
        self.assertIsNone(derived.spectrum.eigenvalues[0].dual)

    def test_extract_krom_subsets_linear_control_representation(self):
        source=_krom(); source.control_mode="linear"; source.B1=np.asarray(((10.,),(20.,))); source._validate_control_representation()
        derived=source.extract([1])
        self.assertEqual(derived.control_mode,"linear")
        np.testing.assert_array_equal(derived.B1,np.asarray(((20.,),)))
        source.B1[1,0]=30.
        self.assertEqual(derived.B1[0,0],20.)

    def test_extract_krom_inherits_independent_closed_bilinear_restriction(self):
        source=_krom(); source.control_mode="bilinear"
        source.B1=np.asarray(((1.,2.),(3.,4.))); source.B2=np.asarray(((5.,),(6.,)))
        source._validate_control_representation()
        probe=np.asarray((2.,3.))
        derived=source.extract([1])
        self.assertEqual(derived.control_mode,"bilinear")
        np.testing.assert_array_equal(derived.B1,np.asarray(((4.,),)))
        np.testing.assert_array_equal(derived.B2,np.asarray(((6.,),)))
        restricted=derived.control_representations["bilinear"]
        self.assertEqual(restricted.operators.shape,(1,1,1))
        np.testing.assert_array_equal(restricted.operators[0],np.asarray(((4.,),)))
        np.testing.assert_array_equal(restricted.offsets,np.asarray(((6.,),)))
        np.testing.assert_allclose(
            derived.B_psi(probe,output="complex",representation="bilinear"),
            restricted.evaluate(derived.psi(probe,output="complex")),
        )
        source.B1[1,1]=99.; source.B2[1,0]=98.
        np.testing.assert_array_equal(restricted.operators[0],np.asarray(((4.,),)))
        np.testing.assert_array_equal(restricted.offsets,np.asarray(((6.,),)))
        nested=derived.extract([0])
        np.testing.assert_array_equal(nested.B1,np.asarray(((4.,),)))
        np.testing.assert_array_equal(nested.B2,np.asarray(((6.,),)))

    def test_extract_krom_rejects_empty_foreign_and_out_of_range_selections(self):
        source=_krom()
        with self.assertRaises(ValueError): source.extract([])
        with self.assertRaises(IndexError): source.extract([2])
        with self.assertRaises(ValueError): source.extract([object()])

    def test_unactuated_system_can_construct_and_evolve_an_autonomous_krom(self):
        system=_UnactuatedSystem(); reference=DynamicalSystemReference(None,"unactuated",{})
        eigenvalues=[KoopmanEigenvalue.from_continuous(-1.,.1,reference),KoopmanEigenvalue.from_continuous(-2.,.1,reference)]
        functions=[KoopmanEigenfunction(eigenvalues[0],"coordinate",state_callable=lambda state:state[0]),
                   KoopmanEigenfunction(eigenvalues[1],"coordinate",state_callable=lambda state:state[1])]
        model=KROM(KoopmanSpectrum(reference,eigenvalues,functions,
            [np.asarray((1.,0.)),np.asarray((0.,1.))],{"method":"test"},dynamical_system=system),system)
        self.assertTrue(model.capabilities()["lifted_evolution"])
        self.assertTrue(model.capabilities()["state_reconstruction"])
        self.assertFalse(model.capabilities()["analytic_control_coupling"])
        result=model.rollout(np.asarray((1.,1.)),np.linspace(0.,.2,3),execution_backend="python")
        self.assertEqual(result.reconstructed.values.shape,(3,2))
        with self.assertRaisesRegex(ValueError,"no compatible input matrix"):
            _=model.input_dimension

    def test_numpy_rollout_matches_exact_diagonal_flow(self):
        model = _krom()
        time = np.linspace(3.0, 4.0, 101)
        result = model.rollout(np.asarray((2.0, -3.0)), time, execution_backend="python")
        expected = np.column_stack((2.0 * np.exp(-(time - 3.0)), -3.0 * np.exp(-2.0 * (time - 3.0))))
        np.testing.assert_allclose(result.reconstructed.values, expected, rtol=2e-14, atol=2e-14)
        np.testing.assert_allclose(result.lifted.values, expected[:, ::-1], rtol=2e-14, atol=2e-14)
        self.assertEqual(result.metadata["execution_backend"]["name"], "numpy-reference")

    def test_reconstruct_accepts_multiple_lifted_rows(self):
        values = _krom().reconstruct(np.asarray(((2.0, 1.0), (4.0, 3.0))))
        np.testing.assert_array_equal(values, np.asarray(((1.0, 2.0), (3.0, 4.0))))

    def test_rollout_rejects_invalid_time(self):
        with self.assertRaises(ValueError):
            _krom().rollout(np.asarray((1.0, 1.0)), np.asarray((0.0, 0.0)))

    def test_auto_keeps_benchmark_admitted_numpy_backend(self):
        result = _krom().rollout(np.asarray((1.0, 1.0)), np.linspace(0.0, 1.0, 11))
        self.assertEqual(result.metadata["execution_backend"]["name"], "numpy-reference")

    def test_native_rollout_matches_numpy_when_capability_is_installed(self):
        try:
            module = importlib.import_module("_dynsyslab_native")
        except ImportError:
            self.skipTest("optional native wheel is not installed")
        if not module.capabilities().get("krom_rollout", False):
            self.skipTest("installed native wheel predates KROM rollout")
        model = _krom(); time = np.linspace(0.0, 2.0, 4001); initial = np.asarray((0.25, -0.5))
        expected = model.rollout(initial, time, execution_backend="python")
        actual = model.rollout(initial, time, execution_backend="native")
        np.testing.assert_allclose(actual.lifted.values, expected.lifted.values, rtol=2e-14, atol=2e-14)
        np.testing.assert_allclose(actual.reconstructed.values, expected.reconstructed.values, rtol=2e-14, atol=2e-14)
        self.assertTrue(actual.metadata["execution_backend"]["compiled"])

    def test_dimensionless_krom_converts_physical_state_time_and_input_units(self):
        model=_scaled_krom()
        np.testing.assert_allclose(model.psi(np.asarray((2.,4.)),state_coordinates="physical"),np.asarray((1.,1.)))
        np.testing.assert_allclose(model.A(time_coordinates="physical"),2.*model.A(time_coordinates="model"))
        expected=np.asarray(((1.,),(.5,)))
        np.testing.assert_allclose(model.B_psi(np.asarray((2.,4.)),output="complex",state_coordinates="physical",time_coordinates="physical"),expected)
        self.assertEqual(model.capabilities()["coordinate_system"],"dimensionless")
        self.assertTrue(model.capabilities()["state_reconstruction"])

    def test_delay_observable_modes_cannot_masquerade_as_physical_state_modes(self):
        model=_scaled_krom(delay_modes=True)
        self.assertEqual(model.state_representation,"delay-observable")
        self.assertFalse(model.capabilities()["state_reconstruction"])
        with self.assertRaisesRegex(ValueError,"delay-observable"):
            model.mode_matrix(physical_state_only=True)

    def test_koopman_lqr_inherits_dimensionless_krom_coordinates(self):
        model=_scaled_krom(); target=ConstantTarget(np.asarray((1.,1.)))
        controller=KoopmanLQRController(model,np.eye(2),np.eye(1),target=target)
        output=controller.value(0.,np.asarray((1.,1.)))
        self.assertEqual(output.shape,(1,))
        self.assertTrue(np.all(np.isfinite(output)))
        self.assertEqual(controller.actuation,"system_input")
        self.assertEqual(controller.coordinate_system,"dimensionless")
        self.assertEqual(controller.time_units,"dimensionless")
        self.assertEqual(controller.state_units,"dimensionless")
        np.testing.assert_allclose(controller.target.evaluate(0.,2),np.asarray((1.,1.)))
        np.testing.assert_allclose(controller.display_target.evaluate(0.,2),np.asarray((1.,1.)))
        metadata=serializable_system(controller)
        self.assertEqual(metadata["controller"],"koopman-lqr")
        self.assertEqual(metadata["rom"]["coordinate_system"],"dimensionless")

    def test_koopman_lqr_applies_explicit_system_input_limit(self):
        model=_scaled_krom(); target=ConstantTarget(np.asarray((10.,10.)))
        controller=KoopmanLQRController(model,np.eye(2),np.eye(1),target=target,control_limit=.01)
        output=controller.value(0.,np.asarray((2.,4.)))
        self.assertLessEqual(abs(float(output[0])),.01)
        np.testing.assert_allclose(controller.archive_metadata()["control_limit"],np.asarray((.01,)))

    def test_koopman_controller_gates_training_domain_extrapolation(self):
        model=_scaled_krom(); model.spectrum.metadata["training_domain"]={
            "encoding":"axis-aligned-training-domain/v1","minimum":[0.,0.],"maximum":[2.,2.],
            "sample_count":10,"coordinate_system":"dimensionless","source":"test"}
        inside=model.domain_report(np.asarray((2.,4.)),state_coordinates="physical")
        outside=model.domain_report(np.asarray((3.,1.)),state_coordinates="model")
        self.assertTrue(inside["inside"]); self.assertFalse(outside["inside"])
        target=ConstantTarget(np.asarray((3.,1.)))
        controller=KoopmanLQRController(model,np.eye(2),np.eye(1),target=target)
        with self.assertRaisesRegex(ValueError,"outside the KROM training domain"):
            controller.value(0.,np.asarray((2.,4.)))
        explicit=KoopmanLQRController(model,np.eye(2),np.eye(1),target=target,target_domain_policy="ignore")
        self.assertTrue(np.all(np.isfinite(explicit.value(0.,np.asarray((2.,4.))))))

    def test_control_coupling_auto_dispatches_but_explicit_analytic_is_strict(self):
        model=_krom()
        with self.assertRaisesRegex(NotImplementedError,"no analytic derivative"):
            model.B_psi(np.asarray((1.,1.)),gradient_backend="analytic")
        automatic=model.B_psi(np.asarray((1.,1.)),output="complex")
        numerical=model.B_psi(np.asarray((1.,1.)),output="complex",gradient_backend="finite_difference")
        np.testing.assert_allclose(automatic,numerical)
        self.assertEqual(numerical.shape,(2,1))


if __name__ == "__main__":
    unittest.main()
