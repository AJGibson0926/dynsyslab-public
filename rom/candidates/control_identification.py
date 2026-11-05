"""Candidate algorithms for identifying Koopman control operators.

The routines are independent of KROM mutation and retain their design matrices
so parity experiments can be audited numerically.
"""
from __future__ import annotations

from typing import Any, Sequence
import numpy as np
from scipy.linalg import expm

from .control_identification_contracts import (
    AlignedLiftedResidualControlFitCandidate, BilinearControlFitCandidate,
    ConstantControlFitCandidate, LiftedControlledDataset,
    LiftedControlledTrajectory, _readonly,
)


def fit_aligned_lifted_residual_control_candidate(
    lifted_states: np.ndarray,
    lifted_derivatives: np.ndarray,
    inputs: np.ndarray,
    A_continuous: np.ndarray,
    *,
    ridge: float=0.0,
    rcond: float=1e-15,
    metadata: dict[str,Any] | None=None,
) -> AlignedLiftedResidualControlFitCandidate:
    """Fit constant and bilinear control maps from aligned complex arrays.

    Sample-major input represents
    ``dPsi-A*Psi = B*u`` and
    ``dPsi-A*Psi = sum_j (M_j*Psi+b_j)*u_j``.
    No derivative, conjugation, or state-coordinate convention is inferred.
    """
    z=np.asarray(lifted_states,dtype=complex)
    dz=np.asarray(lifted_derivatives,dtype=complex)
    control=np.asarray(inputs,dtype=float)
    generator=np.asarray(A_continuous,dtype=complex)
    if control.ndim==1: control=control[:,None]
    if z.ndim!=2 or dz.shape!=z.shape or control.ndim!=2 or control.shape[0]!=z.shape[0]:
        raise ValueError("Psi, dPsi/dt, and inputs must be aligned sample-major arrays.")
    modes=z.shape[1]; input_dimension=control.shape[1]
    if generator.shape!=(modes,modes):
        raise ValueError("A_continuous must be square and aligned with Psi.")
    if z.shape[0]<3 or input_dimension<1 or any(np.any(~np.isfinite(item)) for item in (z,dz,control,generator)):
        raise ValueError("Aligned lifted-residual regression inputs must be finite and nonempty.")
    if ridge<0. or not np.isfinite(ridge) or rcond<=0. or not np.isfinite(rcond):
        raise ValueError("ridge and rcond settings are invalid.")
    target=dz-z@generator.T

    def solve(design: np.ndarray) -> np.ndarray:
        if ridge==0.:
            return np.linalg.pinv(design,rcond=rcond)@target
        gram=design.conj().T@design+ridge*np.eye(design.shape[1])
        return np.linalg.pinv(gram,rcond=rcond)@design.conj().T@target

    constant_design=control.astype(complex)
    constant_coefficients=solve(constant_design)
    constant=constant_coefficients.T
    constant_error=constant_design@constant_coefficients-target
    bilinear_design=np.hstack([
        part for input_index in range(input_dimension)
        for part in (control[:,input_index,None]*z,control[:,input_index,None])
    ])
    coefficients=solve(bilinear_design)
    operators=np.empty((input_dimension,modes,modes),dtype=complex)
    offsets=np.empty((modes,input_dimension),dtype=complex)
    for input_index in range(input_dimension):
        block=coefficients[input_index*(modes+1):(input_index+1)*(modes+1)]
        operators[input_index]=block[:modes].T
        offsets[:,input_index]=block[modes]
    bilinear_error=bilinear_design@coefficients-target
    common={"method":"candidate-aligned-lifted-residual-control-regression/v1","candidate_only":True,
        "sample_count":int(z.shape[0]),"lifted_dimension":int(modes),"input_dimension":int(input_dimension),
        "ridge":float(ridge),"rcond":float(rcond),"metadata":dict(metadata or {}),"promotion":"none"}
    diagnostics={**common,
        "constant":{"design_shape":list(constant_design.shape),"design_rank":int(np.linalg.matrix_rank(constant_design)),
            "relative_residual":float(np.linalg.norm(constant_error)/max(np.linalg.norm(target),np.finfo(float).tiny)),
            "complex_rmse":float(np.sqrt(np.mean(np.abs(constant_error)**2)))},
        "bilinear":{"design_shape":list(bilinear_design.shape),"design_rank":int(np.linalg.matrix_rank(bilinear_design)),
            "relative_residual":float(np.linalg.norm(bilinear_error)/max(np.linalg.norm(target),np.finfo(float).tiny)),
            "complex_rmse":float(np.sqrt(np.mean(np.abs(bilinear_error)**2)))}}
    return AlignedLiftedResidualControlFitCandidate(constant,operators,offsets,target,
        constant_design,bilinear_design,diagnostics)

def discrete_input_to_continuous(
    A_continuous: np.ndarray,
    B_discrete: np.ndarray,
    delta_t: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Solve ``Bd = integral_0^dt exp(Ac s) ds Bc`` by least squares."""
    generator = np.asarray(A_continuous, dtype=complex)
    discrete = np.asarray(B_discrete, dtype=complex)
    if generator.ndim != 2 or generator.shape[0] != generator.shape[1]:
        raise ValueError("A_continuous must be square.")
    if discrete.ndim != 2 or discrete.shape[0] != generator.shape[0]:
        raise ValueError("B_discrete rows must align with A_continuous.")
    if delta_t <= 0.0 or not np.isfinite(delta_t):
        raise ValueError("delta_t must be finite and positive.")
    dimension = generator.shape[0]
    augmented = np.zeros((2 * dimension, 2 * dimension), dtype=complex)
    augmented[:dimension, :dimension] = generator
    augmented[:dimension, dimension:] = np.eye(dimension)
    transition = expm(augmented * delta_t)
    integral = transition[:dimension, dimension:]
    continuous, _, rank, singular_values = np.linalg.lstsq(integral, discrete, rcond=1e-12)
    residual = integral @ continuous - discrete
    diagnostics = {
        "method": "exact-zoh-integral-least-squares/v1",
        "integral_rank": int(rank),
        "integral_condition_number": float(np.linalg.cond(integral)),
        "integral_singular_values": singular_values.tolist(),
        "conversion_relative_residual": float(
            np.linalg.norm(residual) / max(np.linalg.norm(discrete), np.finfo(float).tiny)
        ),
    }
    return continuous, integral, diagnostics

def fit_constant_multistep_output_candidate(
    A_discrete: np.ndarray,
    C: np.ndarray,
    dataset: LiftedControlledDataset,
    *,
    intercept: np.ndarray | None = None,
    A_continuous: np.ndarray | None = None,
    ridge: float = 0.0,
    rcond: float = 1e-12,
) -> ConstantControlFitCandidate:
    """Reproduce Korda--Mezic's multistep output regression for ``Bd``."""
    transition = np.asarray(A_discrete, dtype=complex)
    observation = np.asarray(C, dtype=complex)
    modes = dataset.lifted_dimension
    outputs = dataset.output_dimension
    inputs = dataset.input_dimension
    if transition.shape != (modes, modes):
        raise ValueError("A_discrete does not align with the lifted dataset.")
    if observation.shape != (outputs, modes):
        raise ValueError("C does not align with output and lifted dimensions.")
    if ridge < 0.0 or not np.isfinite(ridge):
        raise ValueError("ridge must be finite and nonnegative.")
    if rcond <= 0.0 or not np.isfinite(rcond):
        raise ValueError("rcond must be finite and positive.")
    bias = np.zeros(outputs, dtype=complex) if intercept is None else np.asarray(intercept, dtype=complex).reshape(-1)
    if bias.size != outputs:
        raise ValueError("Observable intercept does not align with C.")

    blocks: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    for trajectory in dataset.trajectories:
        free_state = trajectory.lifted_states[0].copy()
        control_response = np.zeros((outputs, inputs, modes), dtype=complex)
        for step in range(1, trajectory.lifted_states.shape[0]):
            control_response = np.einsum(
                "oim,mn->oin", control_response, transition, optimize=True
            ) + observation[:, None, :] * trajectory.inputs[step - 1][None, :, None]
            free_state = transition @ free_state
            blocks.append(control_response.reshape(outputs, inputs * modes))
            targets.append(trajectory.outputs[step] - bias - observation @ free_state)
    design = np.vstack(blocks)
    target = np.concatenate(targets)
    if ridge == 0.0:
        vector = np.linalg.pinv(design, rcond=rcond) @ target
        solve_method = "svd-pseudoinverse"
    else:
        gram = design.conj().T @ design + ridge * np.eye(design.shape[1])
        vector = np.linalg.pinv(gram, rcond=rcond) @ design.conj().T @ target
        solve_method = "ridge-normal-equations"
    discrete = vector.reshape((modes, inputs), order="F")
    predicted = design @ vector
    residual = predicted - target
    continuous = None
    conversion = None
    if A_continuous is not None:
        continuous, _, conversion = discrete_input_to_continuous(
            A_continuous, discrete, dataset.delta_t
        )
    diagnostics = {
        "method": "candidate-korda-mezic-multistep-output/v1",
        "candidate_only": True,
        "design_assembly": "recursive-control-convolution/O(horizon)",
        "solve_method": solve_method,
        "ridge": float(ridge),
        "rcond": float(rcond),
        "design_shape": list(design.shape),
        "design_rank": int(np.linalg.matrix_rank(design)),
        "design_condition_number": float(np.linalg.cond(design)),
        "complex_rmse": float(np.sqrt(np.mean(np.abs(residual) ** 2))),
        "relative_residual": float(
            np.linalg.norm(residual) / max(np.linalg.norm(target), np.finfo(float).tiny)
        ),
        "dataset": dataset.summary(),
        "time_domain_of_fit": "discrete",
        "conversion": conversion,
        "promotion": "none",
    }
    return ConstantControlFitCandidate(discrete, continuous, design, target, diagnostics)

def fit_bilinear_continuous_candidate(
    A_continuous: np.ndarray,
    dataset: LiftedControlledDataset,
    *,
    ridge: float = 0.0,
    rcond: float = 1e-12,
) -> BilinearControlFitCandidate:
    """Fit ``B_psi(psi)[:,j] = M[j] psi + b[:,j]`` in continuous time.

    This is the array-level form of the MATLAB experiment
    ``(DZ-A*Z)*pinv([Z.*U; U])``.  Supplied lifted derivatives are preferred;
    otherwise NumPy's second-order finite-difference gradient is used and that
    approximation is recorded.
    """
    generator = np.asarray(A_continuous, dtype=complex)
    modes = dataset.lifted_dimension
    inputs = dataset.input_dimension
    if generator.shape != (modes, modes):
        raise ValueError("A_continuous does not align with the lifted dataset.")
    if ridge < 0.0 or not np.isfinite(ridge):
        raise ValueError("ridge must be finite and nonnegative.")
    rows: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    derivative_sources: set[str] = set()
    for trajectory in dataset.trajectories:
        active = trajectory.lifted_states[:-1]
        controls = trajectory.inputs
        if trajectory.lifted_derivatives is None:
            derivative = np.gradient(
                trajectory.lifted_states, trajectory.time, axis=0, edge_order=2
            )[:-1]
            derivative_sources.add("finite-difference-second-order")
        else:
            derivative = trajectory.lifted_derivatives
            if derivative.shape[0] == trajectory.time.size:
                derivative = derivative[:-1]
            derivative_sources.add("supplied")
        features = np.hstack(
            [
                np.hstack(
                    (
                        controls[:, input_index, None] * active,
                        controls[:, input_index, None],
                    )
                )
                for input_index in range(inputs)
            ]
        )
        rows.append(features)
        targets.append(derivative - active @ generator.T)
    design = np.vstack(rows)
    target = np.vstack(targets)
    if ridge == 0.0:
        coefficients = np.linalg.pinv(design, rcond=rcond) @ target
        solve_method = "svd-pseudoinverse"
    else:
        gram = design.conj().T @ design + ridge * np.eye(design.shape[1])
        coefficients = np.linalg.pinv(gram, rcond=rcond) @ design.conj().T @ target
        solve_method = "ridge-normal-equations"
    operators = np.empty((inputs, modes, modes), dtype=complex)
    offsets = np.empty((modes, inputs), dtype=complex)
    for input_index in range(inputs):
        block = coefficients[input_index * (modes + 1) : (input_index + 1) * (modes + 1)]
        operators[input_index] = block[:modes].T
        offsets[:, input_index] = block[modes]
    residual = design @ coefficients - target
    diagnostics = {
        "method": "candidate-continuous-bilinear-lifted-residual/v1",
        "candidate_only": True,
        "source_matlab_expression": "(DZ-A*Z)*pinv([Z.*U; U])",
        "solve_method": solve_method,
        "derivative_sources": sorted(derivative_sources),
        "ridge": float(ridge),
        "rcond": float(rcond),
        "design_shape": list(design.shape),
        "design_rank": int(np.linalg.matrix_rank(design)),
        "design_condition_number": float(np.linalg.cond(design)),
        "complex_rmse": float(np.sqrt(np.mean(np.abs(residual) ** 2))),
        "relative_residual": float(
            np.linalg.norm(residual) / max(np.linalg.norm(target), np.finfo(float).tiny)
        ),
        "dataset": dataset.summary(),
        "time_domain_of_fit": "continuous",
        "promotion": "none",
    }
    return BilinearControlFitCandidate(operators, offsets, design, target, diagnostics)

def fit_bilinear_control_map_candidate(
    lifted_samples: np.ndarray,
    exact_Bpsi_samples: np.ndarray,
    *,
    ridge: float = 0.0,
    rcond: float = 1e-12,
    metadata: dict[str, Any] | None = None,
) -> BilinearControlFitCandidate:
    """Project known exact matrices onto ``B1*Psi+B2`` pointwise.

    This fit isolates representational error from derivative and autonomous
    generator error. ``exact_Bpsi_samples`` has shape
    ``(samples, modes, inputs)`` and must use the same continuous-time input
    convention as the desired approximation.
    """
    lifted = np.asarray(lifted_samples, dtype=complex)
    exact = np.asarray(exact_Bpsi_samples, dtype=complex)
    if lifted.ndim != 2 or exact.ndim != 3:
        raise ValueError("Lifted samples and exact B_psi samples have invalid rank.")
    if exact.shape[:2] != lifted.shape:
        raise ValueError("Exact B_psi samples must align with sample and lifted dimensions.")
    if np.any(~np.isfinite(lifted)) or np.any(~np.isfinite(exact)):
        raise ValueError("Control-map projection samples must be finite.")
    if ridge < 0.0 or not np.isfinite(ridge):
        raise ValueError("ridge must be finite and nonnegative.")
    if rcond <= 0.0 or not np.isfinite(rcond):
        raise ValueError("rcond must be finite and positive.")
    design = np.column_stack((lifted, np.ones(lifted.shape[0], dtype=complex)))
    modes = lifted.shape[1]
    inputs = exact.shape[2]
    operators = np.empty((inputs, modes, modes), dtype=complex)
    offsets = np.empty((modes, inputs), dtype=complex)
    predicted = np.empty_like(exact)
    if ridge == 0.0:
        inverse = np.linalg.pinv(design, rcond=rcond)
        solve_method = "svd-pseudoinverse"
        for input_index in range(inputs):
            coefficients = inverse @ exact[:, :, input_index]
            operators[input_index] = coefficients[:modes].T
            offsets[:, input_index] = coefficients[modes]
            predicted[:, :, input_index] = design @ coefficients
    else:
        gram = design.conj().T @ design + ridge * np.diag(np.r_[np.ones(modes), 0.0])
        inverse = np.linalg.pinv(gram, rcond=rcond) @ design.conj().T
        solve_method = "ridge-normal-equations"
        for input_index in range(inputs):
            coefficients = inverse @ exact[:, :, input_index]
            operators[input_index] = coefficients[:modes].T
            offsets[:, input_index] = coefficients[modes]
            predicted[:, :, input_index] = design @ coefficients
    residual = predicted - exact
    diagnostics = {
        "method": "candidate-pointwise-bilinear-control-map-projection/v1",
        "candidate_only": True,
        "purpose": "isolate representation capacity from generator/derivative error",
        "solve_method": solve_method,
        "ridge": float(ridge),
        "rcond": float(rcond),
        "design_shape": list(design.shape),
        "design_rank": int(np.linalg.matrix_rank(design)),
        "design_condition_number": float(np.linalg.cond(design)),
        "complex_rmse": float(np.sqrt(np.mean(np.abs(residual) ** 2))),
        "relative_residual": float(
            np.linalg.norm(residual) / max(np.linalg.norm(exact), np.finfo(float).tiny)
        ),
        "metadata": dict(metadata or {}),
        "time_domain_of_fit": "continuous",
        "promotion": "none",
    }
    return BilinearControlFitCandidate(operators, offsets, design, exact, diagnostics)
