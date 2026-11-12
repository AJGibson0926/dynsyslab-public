"""Opt-in research candidates that do not participate in standard workflows.

Nothing in this package is re-exported from :mod:`rom`.  Importing a candidate
must therefore be an explicit experiment-level decision.
"""

from .control_identification import (
    AlignedLiftedResidualControlFitCandidate,
    BilinearControlFitCandidate,
    ConstantControlFitCandidate,
    LiftedControlledDataset,
    LiftedControlledTrajectory,
    discrete_input_to_continuous,
    fit_aligned_lifted_residual_control_candidate,
    fit_bilinear_continuous_candidate,
    fit_bilinear_control_map_candidate,
    fit_constant_multistep_output_candidate,
)
from .error_fields import (
    ControlApproximationErrorField,
    compare_control_approximations_on_grid,
)
from .observable_projection import (
    ObservableProjectionCandidate,
    fit_observable_projection_candidate,
)

__all__ = [
    "AlignedLiftedResidualControlFitCandidate",
    "BilinearControlFitCandidate",
    "ConstantControlFitCandidate",
    "ControlApproximationErrorField",
    "LiftedControlledDataset",
    "LiftedControlledTrajectory",
    "ObservableProjectionCandidate",
    "compare_control_approximations_on_grid",
    "discrete_input_to_continuous",
    "fit_aligned_lifted_residual_control_candidate",
    "fit_bilinear_continuous_candidate",
    "fit_bilinear_control_map_candidate",
    "fit_constant_multistep_output_candidate",
    "fit_observable_projection_candidate",
]
