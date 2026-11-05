"""Immutable named control-study policies.

These records preserve experiment contracts without coupling the control
package to solver construction or archive I/O.
"""
from __future__ import annotations

from types import MappingProxyType


RPE_CLEAN_ISOLATION_KLQR_SETTINGS=MappingProxyType({
    "name":"rpe-clean-isolation-klqr",
    "workflow_role":"established-default-reproduction",
    "accept_initial_transient":True,
    "viewer_transient_policy":"preserve-full-signal-with-automatic-detail-inset",
    "coordinate_system":"dimensionless",
    "input_units":"p_control/p0",
    "initial_state":(1.,0.),
    "target":"rpe_multisine_offset_4_mixed",
    "spectral_selection":"first positive-frequency complex pair by original full-spectrum index (indices 1/2; retained energy-order position 999)",
    "Q_diagonal":(5000.,5000.),
    "R":1e-4,
    "control_limit":None,
    "target_domain_policy":"ignore",
    "control_evaluation":"continuous_rhs",
    "solver_method":"BDF",
    "relative_tolerance":1e-9,
    "absolute_tolerance":1e-9,
    "horizon_tau":10.,
    "output_dtau":.001,
    "feedforward":False,
    "startup_ramp":False,
})

# Deferred design record only. Importing this record changes no controller,
# target, solve, or default. A caller must deliberately implement and select
# every stage under a separately named experiment before it can affect a run.
RPE_OPTIONAL_PRODUCTION_CONTROL_SEQUENCE=MappingProxyType({
    "name":"rpe-optional-production-control-sequence",
    "workflow_role":"deferred-opt-in-add-on",
    "enabled_by_default":False,
    "opt_in_required":True,
    "control_law":"u=-K(Psi(x)-Psi(x_ref))",
    "feedforward":False,
    "command_clipping":False,
    "reference_engagement":MappingProxyType({
        "kind":"quintic-C2-state-reference-engagement",
        "tested_duration_tau":2.,
        "status":"evidence-backed candidate; not an established default",
    }),
    "sequence":(
        "target-and-actuator-feasibility preflight",
        "controlled-domain KROM training and held-out validation",
        "coordinate-invariant physically derived Q/R",
        "optional C2 reference engagement selected by explicit experiment policy",
        "actuator-dynamics augmentation and control-quality KROM selection",
        "comparison against the unchanged transient-producing baseline",
    ),
    "acceptance_gate":"reduce peak and slew without degrading post-engagement tracking, phase, or stability margins",
})


__all__=["RPE_CLEAN_ISOLATION_KLQR_SETTINGS","RPE_OPTIONAL_PRODUCTION_CONTROL_SEQUENCE"]
