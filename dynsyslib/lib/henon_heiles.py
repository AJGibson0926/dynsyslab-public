"""Forced Hénon--Heiles Hamiltonian oscillator system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

HENON_HEILES_PRESETS = {
    "low_energy": SystemPreset("low_energy", "Low-energy bounded Hénon--Heiles orbit (unforced).", {"lam": 1.}, (0., .1, 0., .1), {"frequency_hz": 1., "amplitude": 0.}, "Standard Hénon--Heiles Hamiltonian with unit coupling."),
    "near_escape_energy": SystemPreset("near_escape_energy", "Bounded initial condition close to the E=1/6 escape threshold.", {"lam": 1.}, (0., 0., 0., .55), {"frequency_hz": 1., "amplitude": 0.}, "Useful for mixed regular/chaotic Hénon--Heiles phase-space studies."),
}

@dataclass(frozen=True)
class HenonHeilesSystem(ODESystem):
    """Two-degree-of-freedom Hamiltonian oscillator; driver acts in x."""
    driver: SinusoidalDriver; lam: float=1.; initial_state_override: tuple[float,...]|None=None; name: str="henon-heiles"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (.0,.1,.0,.1))
    @property
    def control_matrix(self): return np.array(((0.,),(0.,),(1.,),(0.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        q,y,pq,py=x; return np.array((pq,py,-q-2*self.lam*q*y,-y-self.lam*(q*q-y*y)))+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x): q,y,_,_=x; return np.array(((0.,0.,1.,0.),(0.,0.,0.,1.),(-1.-2*self.lam*y,-2*self.lam*q,0.,0.),(-2*self.lam*q,-1.+2*self.lam*y,0.,0.)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
