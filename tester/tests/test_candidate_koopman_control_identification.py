from __future__ import annotations

import unittest

import numpy as np

from rom.candidates import (
    LiftedControlledDataset,
    LiftedControlledTrajectory,
    compare_control_approximations_on_grid,
    fit_aligned_lifted_residual_control_candidate,
    fit_bilinear_continuous_candidate,
    fit_bilinear_control_map_candidate,
    fit_constant_multistep_output_candidate,
    fit_observable_projection_candidate,
)


class CandidateKoopmanControlIdentificationTests(unittest.TestCase):
    def test_aligned_lifted_residual_candidate_recovers_constant_and_bilinear_forms(self):
        rng=np.random.default_rng(31)
        samples=300; modes=3
        A=np.diag((-.2+.3j,-.1-.4j,-.05))
        z=rng.normal(size=(samples,modes))+1j*rng.normal(size=(samples,modes))
        u=rng.uniform(-1.,1.,size=(samples,1))
        operator=np.asarray(((.2+.1j,-.1,.05j),(.3,-.2j,.1),(-.05,.08,.12j)))
        offset=np.asarray((.1-.02j,-.03+.04j,.05))
        dz=z@A.T+u*(z@operator.T+offset)
        result=fit_aligned_lifted_residual_control_candidate(z,dz,u,A)
        np.testing.assert_allclose(result.bilinear_operators[0],operator,atol=3e-14)
        np.testing.assert_allclose(result.bilinear_offsets[:,0],offset,atol=3e-14)
        self.assertLess(result.diagnostics["bilinear"]["relative_residual"],1e-13)

    def test_observable_projection_recovers_full_complex_C_without_promotion(self):
        rng = np.random.default_rng(4)
        lifted = rng.normal(size=(100, 3)) + 1j * rng.normal(size=(100, 3))
        expected = np.asarray(((1 + 2j, -0.5j, 0.25), (-1j, 0.75, 2 - 0.25j)))
        observed = lifted @ expected.T
        candidate = fit_observable_projection_candidate(lifted, observed)
        np.testing.assert_allclose(candidate.coefficients, expected, atol=2e-14)
        self.assertTrue(candidate.diagnostics["candidate_only"])
        self.assertEqual(candidate.diagnostics["promotion"], "none; caller must explicitly adapt this full-spectrum candidate")

    def test_multistep_output_candidate_recovers_discrete_B(self):
        rng = np.random.default_rng(7)
        A = np.diag((0.96 + 0.05j, 0.91 - 0.02j))
        C = np.asarray(((1.0, 0.4), (-0.25j, 1.0)))
        expected = np.asarray(((0.03 + 0.01j,), (-0.02 + 0.04j,)))
        trajectories = []
        for _ in range(4):
            controls = rng.uniform(-1.0, 1.0, size=(16, 1))
            lifted = [rng.normal(size=2) + 1j * rng.normal(size=2)]
            for control in controls[:, 0]:
                lifted.append(A @ lifted[-1] + expected[:, 0] * control)
            lifted = np.asarray(lifted)
            outputs = lifted @ C.T
            trajectories.append(
                LiftedControlledTrajectory(
                    np.arange(17) * 0.1, outputs, lifted, controls
                )
            )
        result = fit_constant_multistep_output_candidate(
            A, C, LiftedControlledDataset(trajectories)
        )
        np.testing.assert_allclose(result.B_discrete, expected, rtol=1e-11, atol=1e-12)
        self.assertLess(result.diagnostics["relative_residual"], 1e-11)

    def test_recursive_multistep_design_matches_explicit_convolution_for_multiple_inputs(self):
        rng = np.random.default_rng(71)
        A = np.asarray(((0.91, 0.04), (-0.02, 0.95)), dtype=complex)
        C = np.asarray(((1.0, -0.25), (0.3j, 0.8)), dtype=complex)
        expected = np.asarray(((0.03, -0.01j), (0.02j, -0.04)), dtype=complex)
        controls = rng.uniform(-1.0, 1.0, size=(9, 2))
        lifted = [np.asarray((0.2 + 0.1j, -0.3))]
        for control in controls:
            lifted.append(A @ lifted[-1] + expected @ control)
        lifted = np.asarray(lifted)
        trajectory = LiftedControlledTrajectory(
            np.arange(10) * 0.05, lifted @ C.T, lifted, controls
        )
        result = fit_constant_multistep_output_candidate(
            A, C, LiftedControlledDataset((trajectory,))
        )

        powers = [np.eye(2, dtype=complex)]
        for _ in range(9):
            powers.append(powers[-1] @ A)
        explicit = []
        for step in range(1, 10):
            block = np.zeros((2, 4), dtype=complex)
            for control_index in range(step):
                block += np.kron(
                    controls[control_index][None, :],
                    C @ powers[step - control_index - 1],
                )
            explicit.append(block)
        np.testing.assert_allclose(result.design, np.vstack(explicit), atol=2e-15)
        np.testing.assert_allclose(result.B_discrete, expected, atol=2e-13)
        self.assertEqual(
            result.diagnostics["design_assembly"],
            "recursive-control-convolution/O(horizon)",
        )

    def test_bilinear_candidate_recovers_matlab_B1_B2_form(self):
        rng = np.random.default_rng(11)
        A = np.asarray(((-0.2 + 0.4j, 0.0), (0.0, -0.1 - 0.3j)))
        operator = np.asarray(((0.4 + 0.1j, -0.2), (0.3j, -0.5 + 0.2j)))
        offset = np.asarray((0.15 - 0.1j, -0.05 + 0.08j))
        count = 120
        lifted = rng.normal(size=(count + 1, 2)) + 1j * rng.normal(size=(count + 1, 2))
        controls = rng.uniform(-1.0, 1.0, size=(count, 1))
        derivatives = np.empty((count, 2), dtype=complex)
        for index in range(count):
            derivatives[index] = A @ lifted[index] + (
                operator @ lifted[index] + offset
            ) * controls[index, 0]
        trajectory = LiftedControlledTrajectory(
            np.arange(count + 1) * 0.01,
            np.zeros((count + 1, 1)),
            lifted,
            controls,
            derivatives,
        )
        result = fit_bilinear_continuous_candidate(
            A, LiftedControlledDataset((trajectory,))
        )
        np.testing.assert_allclose(result.operators[0], operator, atol=2e-14)
        np.testing.assert_allclose(result.offsets[:, 0], offset, atol=2e-14)
        self.assertEqual(result.diagnostics["derivative_sources"], ["supplied"])

    def test_error_fields_use_matrix_norm_and_bilinear_can_be_exact(self):
        x = np.linspace(-1.0, 1.0, 5)
        y = np.linspace(-2.0, 2.0, 7)
        grid_x, grid_y = np.meshgrid(x, y)
        lifted = np.stack((grid_x + 1j * grid_y, grid_x * grid_y), axis=-1)
        operator = np.asarray(((((0.5, 0.25j), (-0.1, 0.2)),)))
        offset = np.asarray(((0.1,), (-0.05j,)))
        exact = np.einsum("imn,yxn->yxmi", operator, lifted) + offset
        fields = compare_control_approximations_on_grid(
            x,
            y,
            exact,
            lifted_values=lifted,
            constant_matrix=np.zeros((2, 1)),
            bilinear_operators=operator,
            bilinear_offsets=offset,
        )
        np.testing.assert_allclose(fields["bilinear"].absolute_error, 0.0, atol=1e-15)
        expected_constant = np.linalg.norm(exact, axis=(-2, -1))
        np.testing.assert_allclose(fields["constant"].absolute_error, expected_constant)

    def test_pointwise_control_map_projection_recovers_bilinear_truth(self):
        rng = np.random.default_rng(19)
        lifted = rng.normal(size=(80, 2)) + 1j * rng.normal(size=(80, 2))
        operators = np.asarray(
            (
                ((0.25 + 0.1j, -0.3), (0.2j, 0.45)),
                ((-0.1, 0.15j), (0.35, -0.2j)),
            )
        )
        offsets = np.asarray(((0.2, -0.1j), (-0.05 + 0.02j, 0.3)))
        exact = np.einsum("imn,sn->smi", operators, lifted) + offsets
        result = fit_bilinear_control_map_candidate(lifted, exact)
        np.testing.assert_allclose(result.operators, operators, atol=2e-14)
        np.testing.assert_allclose(result.offsets, offsets, atol=2e-14)
        self.assertLess(result.diagnostics["relative_residual"], 1e-13)


if __name__ == "__main__":
    unittest.main()
