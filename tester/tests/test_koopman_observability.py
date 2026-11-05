from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from rom.observability import (
    analyze_observability,
    classical_observability_matrix,
    continuous_observability_gramian,
    discrete_observability_gramian,
    load_observability_analysis,
    save_observability_analysis,
)


class KoopmanObservabilityTests(unittest.TestCase):
    def test_classical_matrix_uses_right_singular_vectors_as_column_directions(self):
        A=np.diag((.5,.2)); C=np.array(((1.,0.),))
        matrix=classical_observability_matrix(A,C,horizon=2)
        np.testing.assert_allclose(matrix,np.array(((1.,0.),(.5,0.))))
        analysis=analyze_observability(A,C,method="matrix",horizon=2)
        self.assertEqual(analysis.numerical_rank,1)
        self.assertEqual(analysis.directions.shape,(2,2))
        np.testing.assert_allclose(np.abs(analysis.directions[:,0]),np.array((1.,0.)))
        self.assertTrue(analysis.summary()["directions_are_columns"])

    def test_discrete_gramian_equals_classical_matrix_normal_product(self):
        A=np.array(((.7,.1),(0.,.4))); C=np.array(((1.,-.25),))
        matrix=classical_observability_matrix(A,C,horizon=5)
        gramian=discrete_observability_gramian(A,C,horizon=5)
        np.testing.assert_allclose(gramian,matrix.T@matrix,rtol=1e-13,atol=1e-14)

    def test_continuous_finite_horizon_gramian_has_known_diagonal_solution(self):
        A=np.diag((-1.,-2.)); C=np.eye(2); duration=1.5
        expected=np.diag(((1.-np.exp(-2.*duration))/2.,(1.-np.exp(-4.*duration))/4.))
        np.testing.assert_allclose(
            continuous_observability_gramian(A,C,horizon=duration),expected,
            rtol=1e-12,atol=1e-13)

    def test_analysis_artifact_round_trip_preserves_directions_and_metadata(self):
        analysis=analyze_observability(np.diag((-.5,-1.)),np.eye(2),
            method="continuous-gramian",horizon=2.,coordinate_groups=((0,),(1,)),
            metadata={"observable":{"name":"echo"}})
        with tempfile.TemporaryDirectory() as directory:
            path=save_observability_analysis(Path(directory)/"analysis.npz",analysis)
            loaded=load_observability_analysis(path)
        np.testing.assert_allclose(loaded.gramian,analysis.gramian)
        np.testing.assert_allclose(loaded.directions,analysis.directions)
        np.testing.assert_allclose(loaded.mode_participation,analysis.mode_participation)
        self.assertEqual(loaded.metadata["observable"]["name"],"echo")


if __name__=="__main__":
    unittest.main()
