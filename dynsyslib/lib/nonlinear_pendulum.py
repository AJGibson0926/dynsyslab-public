"""Driven, damped nonlinear pendulum ODE."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

NONLINEAR_PENDULUM_PRESETS={
    "earth_small_angle": SystemPreset("earth_small_angle","One-metre laboratory pendulum near the harmonic regime.",{"length_m":1.,"mass_kg":1.,"gravity_m_s2":9.80665,"damping_n_m_s":.02},(.2,0.),{"frequency_hz":.5,"amplitude":0.},"Standard terrestrial pendulum baseline."),
    "large_angle_low_damping": SystemPreset("large_angle_low_damping","Large-amplitude pendulum emphasizing nonlinear period shift.",{"length_m":1.,"mass_kg":1.,"gravity_m_s2":9.80665,"damping_n_m_s":.003},(2.4,0.),{"frequency_hz":.5,"amplitude":0.},"Useful free nonlinear-oscillation benchmark."),
    "driven_damped": SystemPreset("driven_damped","Moderately driven pendulum for resonance and separatrix studies.",{"length_m":1.,"mass_kg":1.,"gravity_m_s2":9.80665,"damping_n_m_s":.08},(.4,0.),{"frequency_hz":.5,"amplitude":.7},"Torque amplitude is in N m."),
    "holmes_chaotic": SystemPreset("holmes_chaotic","Canonical strongly driven, damped pendulum regime used for nonlinear and chaotic-dynamics benchmarks.",{"length_m":1.,"mass_kg":1.,"gravity_m_s2":9.80665,"damping_n_m_s":.5},(-3.141592653589793,0.),{"frequency_hz":(2./3.)*(9.80665**.5)/(2.*np.pi),"amplitude":1.2},"In natural-time units this is theta'' + 0.5 theta' + sin(theta) = 1.2 cos((2/3) tau); driver values are torque N m. Holmes-style periodically forced pendulum benchmark."),
}

@dataclass(frozen=True)
class NonlinearPendulumSystem(ODESystem):
    """State ``(theta, omega)``; driver is an applied pivot-axis torque."""
    driver: SinusoidalDriver
    length_m: float=1.; mass_kg: float=1.; gravity_m_s2: float=9.80665; damping_n_m_s: float=.02
    initial_state_override: tuple[float,float]|None=None; name: str="nonlinear_pendulum"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.asarray(self.initial_state_override if self.initial_state_override is not None else (.2,0.),dtype=float)
    @property
    def inertia(self): return self.mass_kg*self.length_m**2
    @property
    def control_matrix(self): return np.array(((0.,),(1/self.inertia,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,torque_n_m):
        theta,omega=np.asarray(x,dtype=float)
        torque=float(np.asarray(torque_n_m,dtype=float).reshape(-1)[0])
        return np.array((omega,-self.gravity_m_s2/self.length_m*np.sin(theta)-self.damping_n_m_s/self.inertia*omega))+self.input_matrix(x)@np.asarray((torque,))
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x):
        theta,_=np.asarray(x,dtype=float)
        return np.array(((0.,1.),(-self.gravity_m_s2/self.length_m*np.cos(theta),-self.damping_n_m_s/self.inertia)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
