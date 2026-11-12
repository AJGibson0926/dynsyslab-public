"""Typed requests and results for the standalone dynsyslab solver."""
from __future__ import annotations
from dataclasses import dataclass, field, fields, is_dataclass
from typing import Any, Callable, Literal, Protocol
import numpy as np
from time_series import TimeSeries

def _identity_estimator():
    """Construct the neutral default estimator without loading plug-ins."""
    from infrastructure.estimation import IdentityEstimator
    return IdentityEstimator()

SolverKind = Literal["ode", "dde", "pde", "lattice"]

class ControlInput(Protocol):
    """Minimal additive-input protocol accepted by legacy solver requests."""
    def value(self, t: float, state: np.ndarray) -> np.ndarray: ...

@dataclass(frozen=True)
class SolverOptions:
    """Numerical method, coordinate, control-clock, and backend policy."""
    method: str = "RK45"  # RK45, DOP853, Radau, BDF, LSODA, Euler, RK4
    rtol: float = 1e-8
    atol: float = 1e-10
    max_step: float | None = None
    step_size: float | None = None  # required by Euler/RK4; used by DDE
    # None selects a system's preferred formulation: dimensionless when the
    # library model provides one (currently RPE), dimensional otherwise.
    use_dimensionless_form: bool | None = None
    # Explicit zero-order-hold clock in the active dynamics coordinates:
    # nondimensional time for a dimensionless solve, seconds otherwise.
    control_sample_time: float | None = None
    # ``automatic`` preserves the established API: a supplied sample time
    # selects zero-order hold, while its absence evaluates a stateless control
    # law at every RHS call.  Prefer an explicit policy for scientific runs.
    control_evaluation: Literal["automatic","sampled","continuous_rhs"] = "automatic"
    # low=25% (default), medium=10%, high=5% hierarchy progress/ETA updates.
    estimate_aggression: str = "low"
    # ``solve_ivp`` is the lightest default.  ``low_level_stepper`` exposes
    # accepted adaptive steps for sampled/durable histories without restarts.
    adaptive_engine: Literal["solve_ivp","low_level_stepper"] = "solve_ivp"
    # ``auto`` chooses a supported optional native kernel when present and
    # otherwise uses the reference Python/SciPy pathway.  An explicit
    # ``native`` request never silently falls back.
    execution_backend: Literal["auto","python","native"] = "auto"

@dataclass(frozen=True)
class HistoryOptions:
    """Telemetry and opt-in checkpoint policy for a solve lifecycle."""
    # telemetry writes only tiny lifecycle/progress events; no state history.
    mode: Literal["telemetry","sampled","checkpoint"] = "telemetry"
    # Default-on wall-clock heartbeat.  It writes a few scalar fields to the
    # telemetry log and never retains numerical state/history arrays.
    telemetry_enabled: bool = True
    telemetry_heartbeat: bool = True
    track_during_solve: bool = False
    checkpoint_interval: float | None = None
    sample_interval: float | None = None
    telemetry_interval_seconds: float = 60.
    # segmented retains the earlier solve_ivp-restart implementation;
    # low_level_stepper keeps the adaptive solver object alive between marks.
    checkpoint_strategy: Literal["segmented","low_level_stepper"] = "segmented"
    keep_partial_history_on_failure: bool = True
    # Opt-in cooperative pause: stop only after a durable checkpoint segment.
    pause_at_checkpoint: bool = False

@dataclass(frozen=True)
class InterpolationOptions:
    """Optional post-solve resampling for continuous-time solver outputs.

    ``method`` defaults to PCHIP, a shape-preserving cubic interpolant that is
    generally more robust than an unconstrained cubic spline.  Supplying
    ``sample_times`` on SolveRequest is a shorthand for an explicit target
    grid and always enables resampling.
    """
    enabled: bool = False
    time_step: float | None = None
    times: np.ndarray | None = None
    method: str = "pchip"  # pchip, linear, nearest

@dataclass(frozen=True)
class SolveRequest:
    """Complete, immutable specification of one numerical solve."""
    system: Any
    kind: SolverKind
    time_span: tuple[float, float] | None = None
    spatial_domain: Any = None
    mesh: Any = None
    # PDE-only optional contract.  It is intentionally typed as ``Any`` here
    # so ODE/DDE/lattice requests do not import mesh-parallel infrastructure.
    parallel_mesh_plan: Any = None
    pde_parallel_backend: str | None = None  # local or mpi
    control: ControlInput | Callable[[float, np.ndarray], np.ndarray] | None = None
    # ``control`` remains the legacy additive-RHS hook.  The controller and
    # estimator packages use these explicit plug-in points instead.
    controller: Any = None
    # Full, exact state feedback is the explicit default.  Estimators replace
    # this only when a caller supplies one.
    estimator: Any = field(default_factory=_identity_estimator)
    options: SolverOptions = field(default_factory=SolverOptions)
    history: HistoryOptions = field(default_factory=HistoryOptions)
    interpolation: InterpolationOptions = field(default_factory=InterpolationOptions)
    sample_times: np.ndarray | None = None  # legacy shorthand for interpolation.times
    initial_state: np.ndarray | None = None  # optional per-trajectory override
    postprocessors: tuple[Callable[[str], Any], ...] = ()
    on_event: Callable[[dict[str,Any]],None] | None = None
    cancel_token: Any = None
    # A settable token checked only at checkpoint boundaries when enabled.
    pause_token: Any = None
    label: str = "solve"
    # Optional registered storage root.  None retains the primary archive.
    archive_root: str | None = None
    # Set by JobManager for traceable worker telemetry; direct solver use
    # intentionally leaves it unset.
    job_id: str | None = None
    # Small caller-owned provenance (campaign attribution, sweep coordinates,
    # UI job identifiers). It is archived beside the request but never affects
    # numerical dispatch or mutable solver state.
    metadata: dict[str, Any] = field(default_factory=dict)

@dataclass
class NumericalResult:
    """Numerical outcome and aligned state, target, and control histories."""
    time_series: TimeSeries | None
    status: str
    message: str
    archive_path: str | None = None
    history_path: str | None = None
    raw_time_series: TimeSeries | None = None
    control_time_series: TimeSeries | None = None
    raw_control_time_series: TimeSeries | None = None
    control_target_time_series: TimeSeries | None = None
    raw_control_target_time_series: TimeSeries | None = None
    control_call_time: np.ndarray | None = None
    control_call_values: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def time(self): return None if self.time_series is None else self.time_series.time
    @property
    def values(self): return np.empty(0) if self.time_series is None else self.time_series.values
    @property
    def raw_time(self): return None if self.raw_time_series is None else self.raw_time_series.time
    @property
    def raw_values(self): return None if self.raw_time_series is None else self.raw_time_series.values

def serializable_system(system: Any) -> dict[str, Any]:
    """Return the provenance-safe public representation of a system object."""
    from provenance import serializable
    return serializable(system)
