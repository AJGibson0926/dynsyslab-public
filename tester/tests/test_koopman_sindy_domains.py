from __future__ import annotations

import unittest
from unittest.mock import patch
import json
from pathlib import Path
import tempfile

import numpy as np

from dynsyslib import load_system
from postprocessing.koopman.koopman_sindy import (
    _matlab_single_absolute_null_vector,
    _eigenvalue_transfer,
    _sparse_null_vector,
    _sparse_null_vector_from_gram,
    koopman_with_sindy,
)
from postprocessing.koopman.dictionary import PolynomialDictionary
from postprocessing.koopman.spectrum import DynamicalSystemReference, KoopmanEigenvalue
from postprocessing.koopman.storage import sindy_archive_payload


class KoopmanSINDyDomainTests(unittest.TestCase):
    def test_matlab_single_absolute_policy_is_one_unscaled_svd_then_absolute_threshold(self):
        rng=np.random.default_rng(29)
        matrix=rng.normal(size=(40,6))+1j*rng.normal(size=(40,6))
        _,_,vh=np.linalg.svd(matrix,full_matrices=False)
        expected=vh.conj().T[:,-1]
        expected[np.abs(expected)<.3]=0.
        diagnostics={}
        actual=_matlab_single_absolute_null_vector(matrix,.3,diagnostics=diagnostics)
        np.testing.assert_allclose(actual,expected)
        self.assertEqual(diagnostics["svd_subproblems"],1)

    def test_explicit_domain_points_are_retained_exactly(self):
        system=load_system("rpe")
        points=np.asarray(((.91,-.04),(.97,.02),(1.03,-.01),(1.09,.05)))
        result=koopman_with_sindy(eigenvalues=np.asarray((-.2+.3j,)),eigenvalue_domain="continuous",
            delta_t=.01,domain_points=points,system=system,dictionary=PolynomialDictionary(2,False),
            normalize=False,threshold=0.,reuse_conjugate_fits=False,
            directional_derivative_backend="dictionary-gradient",
            sparsification_policy="matlab-single-absolute")
        np.testing.assert_array_equal(result.training_states,points)
        self.assertEqual(result.metadata["input"],"domain_points")
        self.assertEqual(result.metadata["sparsification_policy"],"matlab-single-absolute")

    def test_matlab_sparsification_rejects_gram_backend(self):
        system=load_system("rpe")
        with self.assertRaisesRegex(ValueError,"requires null_space_backend='svd'"):
            koopman_with_sindy(eigenvalues=np.asarray((-.2+.3j,)),eigenvalue_domain="continuous",
                delta_t=.01,domain=(np.asarray((.9,-.1)),np.asarray((1.1,.1))),system=system,
                domain_density=3,null_space_backend="gram-eigh",
                sparsification_policy="matlab-single-absolute")

    def test_reddy_sindy_coordinates_and_control_chain_rule(self):
        system=load_system("reddy_szeri",preset="matlab_cleanup_reference")
        base=system.initial_state.copy(); points=[]
        for radius,velocity in ((.8,-.06),(.9,.03),(1.,0.),(1.1,.05),(1.2,-.02)):
            state=base.copy(); state[0]=radius; state[1]=velocity; points.append(state)
        points=np.asarray(points)
        result=koopman_with_sindy(eigenvalues=np.asarray((-.15+.4j,)),
            eigenvalue_domain="continuous",delta_t=.01,domain_points=points,system=system,
            coordinate_observable="c_k5",coordinate_system="dimensionless",
            dictionary=PolynomialDictionary(1,True),normalize=True,threshold=0.,
            reuse_conjugate_fits=False,directional_derivative_backend="dictionary-gradient",
            compute_training_values=False)
        coordinates=system.observable_coordinates("c_k5",coordinate_system="dimensionless")
        np.testing.assert_allclose(result.training_coordinates,coordinates.evaluate(points))
        self.assertEqual(result.training_states.shape,(5,12))
        self.assertEqual(result.training_coordinates.shape,(5,2))
        self.assertEqual(result.krom.coordinate_system,"dimensionless")
        function=result.krom.eigenfunctions[0]; probe=points[1]
        manual=(function.nabla_psi(coordinates.evaluate(probe))@coordinates.jacobian(probe)
                @system.dimensionless_input_matrix(probe))
        np.testing.assert_allclose(result.krom.B_psi(probe,output="complex"),manual.reshape(1,1))
        self.assertEqual(function.metadata["eigenfunction_coordinate_dimension"],2)
        np.testing.assert_allclose(result.evaluate_training_eigenfunctions()[:,0],
            np.asarray([function.psi(row) for row in result.training_coordinates]))
        residual=function.residuals(points,system=system)
        self.assertTrue(np.isfinite(residual["rms"]))
        from control import KoopmanLQRController
        controller=KoopmanLQRController(result.krom,np.eye(2),np.eye(1),
            target_domain_policy="ignore")
        control=controller.value(0.,probe)
        self.assertEqual(control.shape,(1,))
        self.assertTrue(np.all(np.isfinite(control)))
        from archives.store import serializable_system
        from rom.storage import load_krom_archive
        arrays,analysis=sindy_archive_payload(result,storage_policy="coefficient-first")
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); primary=root/"reddy-training.npz"; archive=root/"reddy-krom.npz"
            solution_metadata={"format":"dynsyslab.solution/v1","system":system.name,
                "system_parameters":serializable_system(system),"delta_t":.01,
                "integration_form":"dimensionless"}
            np.savez(primary,metadata_json=np.asarray(json.dumps(solution_metadata,default=str)))
            analysis={**analysis,"primary_source":str(primary)}
            outer={"format":"dynsyslab.derived/v1","source_archive":str(primary),
                   "analysis_metadata":analysis}
            np.savez(archive,**arrays,metadata_json=np.asarray(json.dumps(outer,default=str)))
            restored=load_krom_archive(archive)
            np.testing.assert_allclose(restored.psi(probe),result.krom.psi(probe))
            np.testing.assert_allclose(restored.B_psi(probe,output="complex"),
                                       result.krom.B_psi(probe,output="complex"))

    def test_gram_eigh_matches_exact_svd_for_a_separated_smallest_direction(self):
        rng=np.random.default_rng(7)
        left,_=np.linalg.qr(rng.normal(size=(120,5)))
        right,_=np.linalg.qr(rng.normal(size=(5,5)))
        matrix=left@np.diag((10.,7.,4.,2.,.1))@right.T
        exact=_sparse_null_vector(matrix,0.,2)
        diagnostics={}
        accelerated=_sparse_null_vector_from_gram(matrix.T@matrix,0.,2,matrix_factory=lambda: matrix,diagnostics=diagnostics)
        self.assertGreater(abs(np.vdot(exact,accelerated)),1.-1e-10)
        self.assertEqual(diagnostics.get("gram_eigh_fallbacks",0),0)
        self.assertGreater(diagnostics["gram_eigh_subproblems"],0)

    def test_gram_eigh_falls_back_for_an_unseparated_smallest_subspace(self):
        matrix=np.diag((4.,2.,1.,1.))
        diagnostics={}
        result=_sparse_null_vector_from_gram(matrix.T@matrix,0.,2,matrix_factory=lambda: matrix,diagnostics=diagnostics)
        self.assertEqual(result.shape,(4,))
        self.assertGreater(diagnostics["gram_eigh_fallbacks"],0)
        self.assertGreater(diagnostics["svd_subproblems"],0)

    def test_discrete_source_is_validated_and_declared_continuous_for_generator(self):
        reference=DynamicalSystemReference(None,"test",{})
        value=KoopmanEigenvalue.from_discrete(.98+.04j,.001,reference)
        metadata=_eigenvalue_transfer([value])
        self.assertEqual(metadata["source_original_domains"],{"discrete":1,"continuous":0})
        self.assertEqual(metadata["generator_eigenvalue_domain"],"continuous")
        self.assertEqual(metadata["display_eigenvalue_domain"],"discrete")
        self.assertLess(metadata["max_log_conversion_error"],1e-14)

    def test_inconsistent_discrete_continuous_pair_is_rejected(self):
        reference=DynamicalSystemReference(None,"test",{})
        value=KoopmanEigenvalue.from_discrete(.98+.04j,.001,reference)
        value.continuous_time+=1.
        with self.assertRaisesRegex(ValueError,r"log\(discrete_time\)/delta_t"):
            _eigenvalue_transfer([value])

    def test_conjugate_modes_reuse_one_exact_fit_and_batch_values(self):
        system=load_system("rpe")
        bounds=(np.array((.9,-.1)),np.array((1.1,.1)))
        eigenvalues=np.array((-.2+1.3j,-.2-1.3j))
        with patch("postprocessing.koopman.koopman_sindy._sparse_null_vector",wraps=_sparse_null_vector) as fitted:
            result=koopman_with_sindy(eigenvalues=eigenvalues,eigenvalue_domain="continuous",delta_t=.01,
                domain=bounds,system=system,domain_density=4,polynomial_degree=2,reuse_conjugate_fits=True)
        self.assertEqual(fitted.call_count,1)
        np.testing.assert_allclose(result.coefficients[1],np.conjugate(result.coefficients[0]))
        np.testing.assert_allclose(result.training_eigenfunction_values[:,1],
                                   np.conjugate(result.training_eigenfunction_values[:,0]))
        acceleration=result.metadata["fit_acceleration"]
        self.assertEqual(acceleration["direct_eigenvalue_fits"],1)
        self.assertEqual(acceleration["reused_conjugate_fits"],1)
        self.assertEqual(acceleration["eigenfunction_evaluation"],"batched-matrix-multiply")

    def test_conjugate_fit_reuse_can_be_disabled(self):
        system=load_system("rpe")
        bounds=(np.array((.9,-.1)),np.array((1.1,.1)))
        with patch("postprocessing.koopman.koopman_sindy._sparse_null_vector",wraps=_sparse_null_vector) as fitted:
            result=koopman_with_sindy(eigenvalues=np.array((-.2+1.3j,-.2-1.3j)),eigenvalue_domain="continuous",delta_t=.01,
                domain=bounds,system=system,domain_density=3,polynomial_degree=2,reuse_conjugate_fits=False)
        self.assertEqual(fitted.call_count,2)
        self.assertFalse(result.metadata["fit_acceleration"]["conjugate_fit_reuse_active"])

    def test_end_to_end_gram_backend_agrees_projectively_with_svd(self):
        system=load_system("rpe"); bounds=(np.array((.9,-.1)),np.array((1.1,.1)))
        settings=dict(eigenvalues=np.array((-.3+.7j,)),eigenvalue_domain="continuous",delta_t=.01,
                      domain=bounds,system=system,domain_density=5,polynomial_degree=2,
                      threshold=.01,reuse_conjugate_fits=False)
        exact=koopman_with_sindy(**settings,null_space_backend="svd")
        accelerated=koopman_with_sindy(**settings,null_space_backend="gram-eigh")
        self.assertGreater(abs(np.vdot(exact.coefficients[0],accelerated.coefficients[0])),1.-1e-7)
        backend=accelerated.metadata["fit_acceleration"]["null_space_backend"]
        self.assertEqual(backend["requested"],"gram-eigh")
        self.assertGreater(backend["gram_eigh_subproblems"],0)

    def test_coefficient_first_archive_omits_reconstructable_dense_values(self):
        system=load_system("rpe"); bounds=(np.array((.9,-.1)),np.array((1.1,.1)))
        result=koopman_with_sindy(eigenvalues=np.array((-.2+.7j,-.2-.7j)),
            eigenvalue_domain="continuous",delta_t=.01,domain=bounds,system=system,
            domain_density=3,polynomial_degree=2,reuse_conjugate_fits=True)
        compact,metadata=sindy_archive_payload(result,storage_policy="coefficient-first")
        self.assertEqual(set(compact),{"coefficients","training_states","training_coordinates","eigenvalues"})
        self.assertFalse(metadata["storage"]["full_coordinate_series_saved"])
        virtual=metadata["storage"]["virtual_arrays"]["eigenfunction_values"]
        self.assertEqual(virtual["encoding"],"polynomial-dictionary-evaluation/v1")
        self.assertEqual(virtual["states"],"training_coordinates")
        self.assertEqual(virtual["shape"],list(result.training_eigenfunction_values.shape))
        full,full_metadata=sindy_archive_payload(result,storage_policy="full")
        np.testing.assert_allclose(full["eigenfunction_values"],result.training_eigenfunction_values)
        self.assertTrue(full_metadata["storage"]["full_coordinate_series_saved"])

    def test_coefficient_only_execution_preserves_callable_krom_and_lazy_archive_values(self):
        system=load_system("rpe"); bounds=(np.array((.9,-.1)),np.array((1.1,.1)))
        settings=dict(eigenvalues=np.array((-.2+.7j,-.2-.7j)),eigenvalue_domain="continuous",
            delta_t=.01,domain=bounds,system=system,domain_density=4,polynomial_degree=2,
            reuse_conjugate_fits=True)
        materialized=koopman_with_sindy(**settings,compute_training_values=True)
        deferred=koopman_with_sindy(**settings,compute_training_values=False)
        self.assertIsNone(deferred.training_eigenfunction_values)
        self.assertTrue(all(function.values is None for function in deferred.spectrum.eigenfunctions))
        self.assertEqual(deferred.metadata["fit_acceleration"]["eigenfunction_evaluation"],
                         "deferred-polynomial-evaluation")
        np.testing.assert_allclose(deferred.coefficients,materialized.coefficients)
        np.testing.assert_allclose(deferred.evaluate_training_eigenfunctions(columns=[1,0]),
                                   materialized.training_eigenfunction_values[:,[1,0]])
        probe=np.array((1.02,.03))
        np.testing.assert_allclose(deferred.krom.psi(probe),materialized.krom.psi(probe))
        compact,metadata=sindy_archive_payload(deferred,storage_policy="coefficient-first")
        self.assertNotIn("eigenfunction_values",compact)
        self.assertEqual(metadata["storage"]["virtual_arrays"]["eigenfunction_values"]["shape"],
                         [deferred.training_states.shape[0],deferred.coefficients.shape[0]])
        full,_=sindy_archive_payload(deferred,storage_policy="full")
        np.testing.assert_allclose(full["eigenfunction_values"],materialized.training_eigenfunction_values)


if __name__=="__main__":
    unittest.main()
