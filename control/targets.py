"""Canonical, reusable control targets for dynsyslab systems."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping
import numpy as np


@dataclass(frozen=True)
class RPEHarmonicTarget:
    """Dimensionless RPE radius target and its analytic radial velocity.

    ``angular_frequencies`` are radians per dimensionless time unit.  The
    returned state is ``(r, dr/dtau)``; no numerical differentiation is used.
    """

    name: str
    offset: float
    amplitudes: tuple[float, ...]
    angular_frequencies: tuple[float, ...]
    phases: tuple[float, ...] = ()
    description: str = ""

    def __post_init__(self) -> None:
        phases=self.phases or (0.,)*len(self.amplitudes)
        if not self.amplitudes or len(self.amplitudes)!=len(self.angular_frequencies):
            raise ValueError("RPE harmonic amplitudes and angular frequencies must be nonempty and equally sized.")
        if len(phases)!=len(self.amplitudes):
            raise ValueError("RPE harmonic phases must match the number of amplitudes.")
        values=np.asarray((self.offset,*self.amplitudes,*self.angular_frequencies,*phases),dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("RPE harmonic target parameters must be finite.")
        object.__setattr__(self,"phases",tuple(float(value) for value in phases))

    def evaluate(self, t: float, state_dimension: int) -> np.ndarray:
        if state_dimension!=2:
            raise ValueError("An RPE harmonic target is defined only for state (r, dr/dtau).")
        amplitude=np.asarray(self.amplitudes,dtype=float)
        omega=np.asarray(self.angular_frequencies,dtype=float)
        angle=omega*float(t)+np.asarray(self.phases,dtype=float)
        return np.asarray((self.offset+np.sum(amplitude*np.sin(angle)),
                           np.sum(amplitude*omega*np.cos(angle))),dtype=float)

    def conservative_state_envelope(self) -> tuple[np.ndarray,np.ndarray]:
        """Return analytic bounds containing the complete target waveform."""
        amplitude=np.abs(np.asarray(self.amplitudes,dtype=float))
        velocity_bound=float(np.sum(amplitude*np.abs(np.asarray(self.angular_frequencies,dtype=float))))
        radius_bound=float(np.sum(amplitude))
        return (np.asarray((self.offset-radius_bound,-velocity_bound),dtype=float),
                np.asarray((self.offset+radius_bound,velocity_bound),dtype=float))


_RPE_CONTROL_TARGETS={
    "rpe_multisine_equilibrium_small": RPEHarmonicTarget(
        "rpe_multisine_equilibrium_small",1.,(.06,.08),(4.,3.),(0.,np.pi),
        description="Small two-tone evaluation target that starts at RPE equilibrium."),
    "rpe_multisine_equilibrium_moderate_identification": RPEHarmonicTarget(
        "rpe_multisine_equilibrium_moderate_identification",1.,(.09,.12),(4.,3.),(0.,np.pi),
        description="Moderate two-tone B_psi-identification target that starts at RPE equilibrium."),
    "rpe_multisine_offset_3p5_periodic": RPEHarmonicTarget(
        "rpe_multisine_offset_3p5_periodic",3.5,(1.5,1/3,1.),(np.pi,2*np.pi,4*np.pi),
        description="Three-tone periodic radius target."),
    "rpe_multisine_offset_4_mixed": RPEHarmonicTarget(
        "rpe_multisine_offset_4_mixed",4.,(1.5,1.),(4.,3.),
        description="Two-tone radius target with angular frequencies 4 and 3."),
    "rpe_multisine_offset_4_periodic": RPEHarmonicTarget(
        "rpe_multisine_offset_4_periodic",4.,(1.5,1/3,1.),(np.pi,2*np.pi,4*np.pi),
        description="Three-tone periodic radius target."),
    "rpe_multisine_offset_2p85_mixed": RPEHarmonicTarget(
        "rpe_multisine_offset_2p85_mixed",2.85,(1.5,1.),(4.,3.),
        description="Two-tone radius target with angular frequencies 4 and 3."),
    "rpe_single_sine_offset_2_experimental": RPEHarmonicTarget(
        "rpe_single_sine_offset_2_experimental",2.,(1.5,),(4.59333,),
        description="Legacy single-tone experimental radius target."),
}

RPE_CONTROL_TARGETS: Mapping[str,RPEHarmonicTarget]=MappingProxyType(_RPE_CONTROL_TARGETS)


def rpe_control_target(name: str) -> RPEHarmonicTarget:
    """Return one immutable canonical RPE target by stable name."""
    try: return RPE_CONTROL_TARGETS[name]
    except KeyError as error:
        raise KeyError(f"Unknown RPE control target {name!r}; choose from {sorted(RPE_CONTROL_TARGETS)}.") from error
