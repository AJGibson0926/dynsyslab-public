"""Deterministic state sampling and virtual KROM candidate screening."""
from __future__ import annotations

from itertools import combinations_with_replacement
from typing import Iterable, Mapping, Sequence
import numpy as np


def trajectory_tube_states(
    states: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    radii: Sequence[float] = (0.025, 0.05),
) -> np.ndarray:
    """Build deterministic axis-aligned samples around a reference path.

    This is a sampling utility for local control-map identification, not a
    trajectory generator. It evaluates no dynamics and cannot silently expand
    or replace the campaign's retained training datasets.
    """
    path = np.asarray(states, dtype=float)
    lo = np.asarray(lower, dtype=float).reshape(-1)
    hi = np.asarray(upper, dtype=float).reshape(-1)
    if (
        path.ndim != 2
        or path.shape[0] < 3
        or path.shape[1] != lo.size
        or hi.shape != lo.shape
        or np.any(~np.isfinite(path))
        or np.any(~np.isfinite(lo))
        or np.any(~np.isfinite(hi))
        or np.any(hi <= lo)
    ):
        raise ValueError("A finite sample-major path and valid bounding box are required.")
    widths = tuple(float(value) for value in radii)
    if any(value <= 0.0 or not np.isfinite(value) for value in widths):
        raise ValueError("Tube radii must be finite and positive.")
    samples = [path]
    for radius in widths:
        for axis in range(path.shape[1]):
            offset = np.zeros(path.shape[1], dtype=float)
            offset[axis] = radius
            samples.extend((path + offset, path - offset))
    return np.clip(np.vstack(samples), lo, hi)

def resample_trajectory_by_arclength(
    states: np.ndarray, sample_count: int
) -> np.ndarray:
    """Resample a state path uniformly in Euclidean chord length."""
    path = np.asarray(states, dtype=float)
    count = int(sample_count)
    if path.ndim != 2 or path.shape[0] < 2 or np.any(~np.isfinite(path)):
        raise ValueError("states must be a finite sample-major path.")
    if count < 2:
        raise ValueError("sample_count must be at least two.")
    segment_lengths = np.linalg.norm(np.diff(path, axis=0), axis=1)
    keep = np.r_[True, segment_lengths > np.finfo(float).eps]
    unique_path = path[keep]
    if unique_path.shape[0] < 2:
        raise ValueError("The path must contain at least two distinct states.")
    cumulative = np.r_[
        0.0, np.cumsum(np.linalg.norm(np.diff(unique_path, axis=0), axis=1))
    ]
    targets = np.linspace(0.0, cumulative[-1], count)
    return np.column_stack(
        [np.interp(targets, cumulative, unique_path[:, axis]) for axis in range(path.shape[1])]
    )

def rank_krom_candidates(
    parent,
    selections: Iterable[Sequence[int]],
    initial_state: np.ndarray,
    target_state: np.ndarray,
    *,
    structural_closure_quality: Mapping[tuple[int, ...], float] | None = None,
) -> list[KROMCandidate]:
    """Rank virtual subsets without fitting a control representation.

    Screening uses only existing exact eigenfunctions: lifted target separation,
    exact local input authority, local linear controllability, modal frequency
    spread, and an optional analytical (not data-regressed) closure quality.
    It does not count as a campaign attempt because no fitted KROM or controller
    solve is created.
    """
    from .control_campaign import KROMCandidate
    initial = np.asarray(initial_state, dtype=float).reshape(-1)
    target = np.asarray(target_state, dtype=float).reshape(-1)
    closure = structural_closure_quality or {}
    candidates: list[KROMCandidate] = []
    seen: set[tuple[int, ...]] = set()
    for selection in selections:
        indices = tuple(sorted(int(item) for item in selection))
        if indices in seen:
            continue
        seen.add(indices)
        model = parent.extract(indices)
        generator = model.A(output="split")
        dimension = generator.shape[0]
        probes = (initial, 0.5 * (initial + target), target)
        ratios: list[float] = []
        authorities: list[float] = []
        ranks: list[int] = []
        for state in probes:
            input_matrix = model.B_psi(state, output="split", representation="exact")
            controllability = np.column_stack(
                [np.linalg.matrix_power(generator, power) @ input_matrix for power in range(dimension)]
            )
            singular = np.linalg.svd(controllability, compute_uv=False)
            ranks.append(int(np.linalg.matrix_rank(controllability, tol=1e-10)))
            ratios.append(float(singular[-1] / max(singular[0], np.finfo(float).tiny)))
            authorities.append(float(np.linalg.norm(input_matrix)))
        displacement = float(
            np.linalg.norm(
                model.psi(initial, output="split") - model.psi(target, output="split")
            )
        )
        frequencies = np.asarray(
            [abs(item.continuous_time.imag) for item in model.eigenvalues], dtype=float
        )
        frequency_spread = float(np.ptp(frequencies)) if frequencies.size > 1 else 0.0
        closure_quality = float(closure.get(indices, 1.0))
        worst_ratio = min(ratios)
        minimum_authority = min(authorities)
        # Fractional controllability weighting avoids numerically tiny singular
        # ratios dominating useful target separation while still penalizing
        # poorly actuated subsets.  Every component is retained in the ledger.
        score = float(
            displacement
            * (max(worst_ratio, np.finfo(float).eps) ** 0.25)
            * np.log1p(minimum_authority)
            * max(closure_quality, np.finfo(float).eps)
            * (1.0 + 0.05 * frequency_spread)
            / (len(indices) ** 0.15)
        )
        candidates.append(
            KROMCandidate(
                indices,
                {
                    "selection_score": score,
                    "lifted_target_displacement": displacement,
                    "worst_normalized_controllability": worst_ratio,
                    "minimum_exact_input_authority": minimum_authority,
                    "frequency_spread": frequency_spread,
                    "structural_closure_quality": closure_quality,
                    "minimum_controllability_rank": float(min(ranks)),
                    "split_dimension": float(dimension),
                },
                "virtual exact-map screen; no control approximation was fitted",
            )
        )
    return sorted(
        candidates,
        key=lambda item: (
            -item.scores["selection_score"],
            len(item.retained_indices),
            item.retained_indices,
        ),
    )

def polynomial_bilinear_closure_quality(
    parent,
    selections: Iterable[Sequence[int]],
    *,
    constancy_tolerance: float = 1e-12,
) -> dict[tuple[int, ...], float]:
    """Measure analytical closure of ``nabla(Psi) B`` in ``span(Psi, 1)``.

    This coefficient-space calculation applies to analytic polynomial KROMs
    with a constant input matrix and identity coordinate observable.  It is a
    structural screen, not a fit to controlled data.  A quality of one means
    every selected exact input-map row is exactly representable by an affine
    bilinear map of the selected lifted coordinates; zero means the worst row
    has unit or larger relative coefficient residual.
    """
    functions = list(parent.eigenfunctions)
    if not functions or any(item.coefficients is None for item in functions):
        raise ValueError("Polynomial closure screening requires coefficient-backed eigenfunctions.")
    dictionary = dict(parent.spectrum.metadata.get("dictionary") or {})
    if dictionary.get("encoding") != "polynomial-dictionary/v1":
        raise ValueError("Polynomial closure screening requires a polynomial dictionary.")
    state_dimension = int(np.asarray(parent.system.initial_state).size)
    degree = int(dictionary.get("degree"))
    include_constant = bool(dictionary.get("include_constant", True))
    exponents: list[tuple[int, ...]] = []
    if include_constant:
        exponents.append((0,) * state_dimension)
    for order in range(1, degree + 1):
        for indices in combinations_with_replacement(range(state_dimension), order):
            exponent = [0] * state_dimension
            for index in indices:
                exponent[index] += 1
            exponents.append(tuple(exponent))
    coefficients = np.vstack([np.asarray(item.coefficients, dtype=complex) for item in functions])
    if coefficients.shape[1] != len(exponents):
        raise ValueError("Eigenfunction coefficients do not align with the polynomial dictionary.")
    preprocessing = dict(parent.spectrum.metadata.get("state_preprocessing") or {})
    scale = np.asarray(preprocessing.get("scale"), dtype=float)
    multiplier = np.asarray(preprocessing.get("multiplier"), dtype=float)
    if scale.shape != (state_dimension,) or multiplier.shape != scale.shape or np.any(scale <= 0.0):
        raise ValueError("Polynomial closure screening requires archived affine state preprocessing.")
    probes = [np.zeros(state_dimension), np.ones(state_dimension), -np.ones(state_dimension)]
    matrices = [np.asarray(parent._input_matrix(item, "model"), dtype=float) for item in probes]
    if any(matrix.shape != matrices[0].shape for matrix in matrices):
        raise ValueError("Input-matrix shape changes across structural-screen probes.")
    if any(not np.allclose(matrix, matrices[0], rtol=constancy_tolerance, atol=constancy_tolerance)
           for matrix in matrices[1:]):
        raise ValueError("Polynomial coefficient closure is only valid for a constant input matrix.")
    input_matrix = matrices[0]
    lookup = {exponent: index for index, exponent in enumerate(exponents)}
    derivatives = np.zeros(
        (len(functions), len(exponents), input_matrix.shape[1]), dtype=complex
    )
    state_scale = multiplier / scale
    for mode, coefficient in enumerate(coefficients):
        for feature, exponent in enumerate(exponents):
            for state_axis in range(state_dimension):
                power = exponent[state_axis]
                if power == 0:
                    continue
                reduced = list(exponent)
                reduced[state_axis] -= 1
                destination = lookup.get(tuple(reduced))
                if destination is None:
                    continue
                derivatives[mode, destination] += (
                    coefficient[feature]
                    * power
                    * state_scale[state_axis]
                    * input_matrix[state_axis]
                )
    constant = np.zeros(len(exponents), dtype=complex)
    zero_exponent = lookup.get((0,) * state_dimension)
    if zero_exponent is not None:
        constant[zero_exponent] = 1.0
    else:
        # Even a dictionary without a constant feature permits the explicit
        # bilinear offset b.  Append an exact constant direction virtually.
        coefficients = np.column_stack((coefficients, np.zeros(coefficients.shape[0])))
        derivatives = np.concatenate(
            (derivatives, np.zeros((derivatives.shape[0], 1, derivatives.shape[2]), dtype=complex)),
            axis=1,
        )
        constant = np.zeros(coefficients.shape[1], dtype=complex)
        constant[-1] = 1.0
    quality: dict[tuple[int, ...], float] = {}
    for selection in selections:
        indices = tuple(sorted(int(item) for item in selection))
        if not indices or any(item < 0 or item >= len(functions) for item in indices):
            raise IndexError("Polynomial closure selection is outside the parent KROM.")
        basis = np.vstack((coefficients[list(indices)], constant))
        residuals: list[float] = []
        for mode in indices:
            for input_index in range(input_matrix.shape[1]):
                target = derivatives[mode, :, input_index]
                fitted = np.linalg.lstsq(basis.T, target, rcond=None)[0] @ basis
                residuals.append(
                    float(np.linalg.norm(target - fitted) / max(np.linalg.norm(target), np.finfo(float).tiny))
                )
        quality[indices] = float(max(0.0, 1.0 - max(residuals)))
    return quality
