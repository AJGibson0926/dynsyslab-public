"""Forced van der Pol oscillator system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

VAN_DER_POL_PRESETS = {
    "weakly_nonlinear": SystemPreset("weakly_nonlinear", "Near-harmonic van der Pol limit cycle.", {"mu": .1}, (.1, 0.), {"frequency_hz": 1., "amplitude": .1}, "van der Pol weak-nonlinearity regime."),
    "relaxation_oscillation": SystemPreset("relaxation_oscillation", "Slow-fast relaxation-oscillation regime.", {"mu": 8.}, (.1, 0.), {"frequency_hz": 1., "amplitude": .1}, "van der Pol relaxation regime."),
}

@dataclass(frozen=True)
class VanDerPolSystem(ODESystem):
    driver: SinusoidalDriver; mu: float=8.; initial_state_override: tuple[float,...]|None=None; name: str="van-der-pol"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (.1,0.))
    @property
    def control_matrix(self): return np.array(((0.,),(1.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        q,v=x; return np.array((v,self.mu*(1-q*q)*v-q))+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x): q,v=x; return np.array(((0.,1.),(-2*self.mu*q*v-1,self.mu*(1-q*q))))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
