"""Array-level contracts for candidate Koopman control identification."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence
import numpy as np


def _readonly(value: np.ndarray, *, dtype) -> np.ndarray:
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result

@dataclass(frozen=True)
class LiftedControlledTrajectory:
    """One controlled trajectory with precomputed full-complex lifting."""

    time: np.ndarray
    outputs: np.ndarray
    lifted_states: np.ndarray
    inputs: np.ndarray
    lifted_derivatives: np.ndarray | None = None
    metadata: dict[str, Any] | None = None

    def __post_init__(self):
        time = np.asarray(self.time, dtype=float).reshape(-1)
        outputs = np.asarray(self.outputs, dtype=complex)
        lifted = np.asarray(self.lifted_states, dtype=complex)
        inputs = np.asarray(self.inputs, dtype=float)
        if outputs.ndim == 1:
            outputs = outputs[:, None]
        if lifted.ndim == 1:
            lifted = lifted[:, None]
        if inputs.ndim == 1:
            inputs = inputs[:, None]
        if time.size < 3 or np.any(~np.isfinite(time)) or np.any(np.diff(time) <= 0):
            raise ValueError("Trajectory time must be finite, increasing, and contain at least three samples.")
        if outputs.ndim != 2 or outputs.shape[0] != time.size or np.any(~np.isfinite(outputs)):
            raise ValueError("Outputs must be finite and sample-major with one row per time.")
        if lifted.ndim != 2 or lifted.shape[0] != time.size or np.any(~np.isfinite(lifted)):
            raise ValueError("Lifted states must be finite and sample-major with one row per time.")
        if inputs.ndim != 2 or inputs.shape[0] != time.size - 1 or np.any(~np.isfinite(inputs)):
            raise ValueError("Inputs must have one finite row per transition interval.")
        derivatives = None
        if self.lifted_derivatives is not None:
            derivatives = np.asarray(self.lifted_derivatives, dtype=complex)
            if derivatives.ndim == 1:
                derivatives = derivatives[:, None]
            if derivatives.shape not in {lifted.shape, (time.size - 1, lifted.shape[1])}:
                raise ValueError("Lifted derivatives must have one row per state sample or transition.")
            if np.any(~np.isfinite(derivatives)):
                raise ValueError("Lifted derivatives must be finite.")
            derivatives = _readonly(derivatives, dtype=complex)
        time = _readonly(time, dtype=float)
        outputs = _readonly(outputs, dtype=complex)
        lifted = _readonly(lifted, dtype=complex)
        inputs = _readonly(inputs, dtype=float)
        object.__setattr__(self, "time", time)
        object.__setattr__(self, "outputs", outputs)
        object.__setattr__(self, "lifted_states", lifted)
        object.__setattr__(self, "inputs", inputs)
        object.__setattr__(self, "lifted_derivatives", derivatives)
        object.__setattr__(self, "metadata", dict(self.metadata or {}))

    @property
    def delta_t(self) -> float:
        steps = np.diff(self.time)
        if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
            raise ValueError("Candidate identification requires uniform sample time.")
        return float(steps[0])

    @property
    def lifted_dimension(self) -> int:
        return int(self.lifted_states.shape[1])

    @property
    def output_dimension(self) -> int:
        return int(self.outputs.shape[1])

    @property
    def input_dimension(self) -> int:
        return int(self.inputs.shape[1])

@dataclass(frozen=True)
class LiftedControlledDataset:
    trajectories: tuple[LiftedControlledTrajectory, ...]
    metadata: dict[str, Any] | None = None

    def __init__(
        self,
        trajectories: Sequence[LiftedControlledTrajectory],
        metadata: dict[str, Any] | None = None,
    ):
        items = tuple(trajectories)
        if not items:
            raise ValueError("A lifted controlled dataset needs at least one trajectory.")
        signature = (
            items[0].lifted_dimension,
            items[0].output_dimension,
            items[0].input_dimension,
        )
        if any(
            (item.lifted_dimension, item.output_dimension, item.input_dimension) != signature
            for item in items
        ):
            raise ValueError("All lifted trajectories must share lifted, output, and input dimensions.")
        if any(not np.isclose(item.delta_t, items[0].delta_t) for item in items):
            raise ValueError("All lifted trajectories must share one sample time.")
        object.__setattr__(self, "trajectories", items)
        object.__setattr__(self, "metadata", dict(metadata or {}))

    @property
    def delta_t(self) -> float:
        return self.trajectories[0].delta_t

    @property
    def lifted_dimension(self) -> int:
        return self.trajectories[0].lifted_dimension

    @property
    def output_dimension(self) -> int:
        return self.trajectories[0].output_dimension

    @property
    def input_dimension(self) -> int:
        return self.trajectories[0].input_dimension

    def summary(self) -> dict[str, Any]:
        return {
            "encoding": "candidate-prelifted-controlled-ensemble/v1",
            "candidate_only": True,
            "trajectory_count": len(self.trajectories),
            "transition_count": int(sum(item.inputs.shape[0] for item in self.trajectories)),
            "delta_t": self.delta_t,
            "lifted_dimension": self.lifted_dimension,
            "output_dimension": self.output_dimension,
            "input_dimension": self.input_dimension,
            "metadata": dict(self.metadata),
        }

@dataclass(frozen=True)
class ConstantControlFitCandidate:
    B_discrete: np.ndarray
    B_continuous: np.ndarray | None
    design: np.ndarray
    target: np.ndarray
    diagnostics: dict[str, Any]

    def __post_init__(self):
        discrete = _readonly(self.B_discrete, dtype=complex)
        continuous = None if self.B_continuous is None else _readonly(self.B_continuous, dtype=complex)
        design = _readonly(self.design, dtype=complex)
        target = _readonly(self.target, dtype=complex)
        if discrete.ndim != 2 or design.ndim != 2 or target.ndim != 1:
            raise ValueError("Candidate constant-control arrays have invalid rank.")
        if continuous is not None and continuous.shape != discrete.shape:
            raise ValueError("Continuous and discrete candidate B matrices must align.")
        object.__setattr__(self, "B_discrete", discrete)
        object.__setattr__(self, "B_continuous", continuous)
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

@dataclass(frozen=True)
class BilinearControlFitCandidate:
    operators: np.ndarray
    offsets: np.ndarray
    design: np.ndarray
    target: np.ndarray
    diagnostics: dict[str, Any]

    def __post_init__(self):
        operators = _readonly(self.operators, dtype=complex)
        offsets = _readonly(self.offsets, dtype=complex)
        design = _readonly(self.design, dtype=complex)
        target = _readonly(self.target, dtype=complex)
        if operators.ndim != 3:
            raise ValueError("Bilinear operators must have shape (inputs, modes, modes).")
        if offsets.shape != (operators.shape[1], operators.shape[0]):
            raise ValueError("Bilinear offsets must have shape (modes, inputs).")
        object.__setattr__(self, "operators", operators)
        object.__setattr__(self, "offsets", offsets)
        object.__setattr__(self, "design", design)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    def evaluate(self, lifted_state: np.ndarray) -> np.ndarray:
        psi = np.asarray(lifted_state, dtype=complex).reshape(-1)
        if psi.size != self.operators.shape[1]:
            raise ValueError("Lifted state does not align with bilinear operators.")
        return self.offsets + np.column_stack(
            [self.operators[index] @ psi for index in range(self.operators.shape[0])]
        )

@dataclass(frozen=True)
class AlignedLiftedResidualControlFitCandidate:
    """Constant and bilinear fits from explicitly aligned lifted samples.

    This is the most literal array-level form of the cleanup-suite regression:
    the caller supplies ``Psi``, ``dPsi/dt``, and ``u`` on the same samples.
    It is intentionally a candidate so compatibility conventions cannot alter
    standard KROM behavior implicitly.
    """

    constant_matrix: np.ndarray
    bilinear_operators: np.ndarray
    bilinear_offsets: np.ndarray
    lifted_control_residual: np.ndarray
    constant_design: np.ndarray
    bilinear_design: np.ndarray
    diagnostics: dict[str, Any]

    def __post_init__(self):
        constant=_readonly(self.constant_matrix,dtype=complex)
        operators=_readonly(self.bilinear_operators,dtype=complex)
        offsets=_readonly(self.bilinear_offsets,dtype=complex)
        residual=_readonly(self.lifted_control_residual,dtype=complex)
        constant_design=_readonly(self.constant_design,dtype=complex)
        bilinear_design=_readonly(self.bilinear_design,dtype=complex)
        if constant.ndim!=2 or operators.ndim!=3 or offsets.shape!=(operators.shape[1],operators.shape[0]):
            raise ValueError("Aligned lifted-residual control matrices have incompatible shapes.")
        if constant.shape!=(operators.shape[1],operators.shape[0]):
            raise ValueError("Constant and bilinear candidates must share mode and input dimensions.")
        object.__setattr__(self,"constant_matrix",constant)
        object.__setattr__(self,"bilinear_operators",operators)
        object.__setattr__(self,"bilinear_offsets",offsets)
        object.__setattr__(self,"lifted_control_residual",residual)
        object.__setattr__(self,"constant_design",constant_design)
        object.__setattr__(self,"bilinear_design",bilinear_design)
        object.__setattr__(self,"diagnostics",dict(self.diagnostics))
