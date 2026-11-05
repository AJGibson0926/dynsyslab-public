"""System-specific display defaults retained by the archive-only Viewer suite."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any

@dataclass(frozen=True)
class ViewerSystemSettings:
    """Recommended visual defaults for one dynsyslib system.

    Entries are presentation-only: they never affect a solve or archived data.
    """
    system: str
    recommended_views: tuple[str,...]
    state_components: tuple[int,...]=()
    camera: dict[str,Any] | None=None
    svg_yaw: float=-.18
    svg_depth_weight: float=.13
    notes: str=""

@dataclass(frozen=True)
class InterpolationSamplingSettings:
    """Display-density policy for one state-space interpolation representation."""
    kind: str
    max_state_samples: int | None
    recommended_domain_density: int
    notes: str=""

_OBLIQUE_3D={"center":{"x":0.,"y":0.,"z":0.},
             "eye":{"x":1.55,"y":-1.65,"z":0.82},
             "up":{"x":0.,"y":0.,"z":1.}}

VIEWER_SYSTEM_SETTINGS={
    "rpe": ViewerSystemSettings("rpe", ("time", "phase", "fourier"), (0,1), notes="Use radius/velocity axes; dimensional labels follow archive metadata."),
    "rayleigh-plesset": ViewerSystemSettings("rayleigh-plesset", ("time", "phase", "fourier"), (0,1)),
    "duffing": ViewerSystemSettings("duffing", ("time", "phase", "fourier"), (0,1)),
    "van_der_pol": ViewerSystemSettings("van_der_pol", ("time", "phase"), (0,1)),
    "nonlinear_pendulum": ViewerSystemSettings("nonlinear_pendulum", ("time", "phase"), (0,1)),
    "henon_heiles": ViewerSystemSettings("henon_heiles", ("phase", "time"), (0,1)),
    "double_pendulum": ViewerSystemSettings("double_pendulum", ("time", "phase"), (0,1)),
    "spring_coupled_pendulums": ViewerSystemSettings("spring_coupled_pendulums", ("time", "phase"), (0,1)),
    "lorenz": ViewerSystemSettings("lorenz", ("three_d", "time", "koopman_spectrum"), (0,1,2), _OBLIQUE_3D, -.18, .13, "Centered oblique camera exposes both lobes."),
    "rossler": ViewerSystemSettings("rossler", ("three_d", "time", "koopman_spectrum"), (0,1,2), _OBLIQUE_3D),
    "chua": ViewerSystemSettings("chua", ("three_d", "time", "koopman_spectrum"), (0,1,2), _OBLIQUE_3D),
    "hindmarsh_rose": ViewerSystemSettings("hindmarsh_rose", ("three_d", "time"), (0,1,2), _OBLIQUE_3D),
    "hodgkin_huxley": ViewerSystemSettings("hodgkin_huxley", ("time", "phase"), (0,1)),
    "fput": ViewerSystemSettings("fput", ("time", "fourier"), (0,1)),
    "kuramoto": ViewerSystemSettings("kuramoto", ("time", "fourier"), (0,1)),
    "ising": ViewerSystemSettings("ising", ("time",), (0,)),
    "ising_3d": ViewerSystemSettings("ising_3d", ("time",), (0,)),
    "potts": ViewerSystemSettings("potts", ("time",), (0,)),
    "mackey_glass": ViewerSystemSettings("mackey_glass", ("time", "phase"), (0,1)),
}

# These limits affect display evaluation only; they never change archived data,
# fitted decompositions, or solver sample density.
VIEWER_INTERPOLATION_SAMPLING={
    "pod": InterpolationSamplingSettings("pod", 800, 120, "POD trajectory-amplitude inverse-distance interpolation."),
    "havok": InterpolationSamplingSettings("havok", 320, 120, "HAVOK PROM modes interpolate a long Hankel-coordinate trajectory."),
    "hdmd": InterpolationSamplingSettings("hdmd", None, 160, "H-DMD eigenfunction trajectory interpolation; retain the full archived time resolution by default."),
    "sindy": InterpolationSamplingSettings("sindy", None, 200, "SINDy eigenfunction expressions retain the full archived time resolution by default."),
}

def get_interpolation_sampling(kind: str | None) -> InterpolationSamplingSettings:
    """Return a conservative default; unknown interpolants use POD density."""
    return VIEWER_INTERPOLATION_SAMPLING.get((kind or "").lower(), VIEWER_INTERPOLATION_SAMPLING["pod"])

def get_viewer_settings(system: str | None) -> ViewerSystemSettings | None:
    """Return retained defaults without forcing an unrecognized system to fit one."""
    return VIEWER_SYSTEM_SETTINGS.get((system or "").lower())
