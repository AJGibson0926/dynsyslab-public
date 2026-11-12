"""Dimensionless forced Chua circuit system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

CHUA_PRESETS = {
    "canonical_double_scroll": SystemPreset("canonical_double_scroll", "Dimensionless double-scroll Chua circuit.", {"alpha": 15.6, "beta": 28., "m0": -1.143, "m1": -.714}, (.7,0.,0.), {"frequency_hz": 1., "amplitude": 0.}, "Canonical dimensionless Chua double-scroll benchmark."),
    "double_scroll_seed_left": SystemPreset("double_scroll_seed_left", "Same canonical circuit with a seed near the opposite lobe.", {"alpha": 15.6, "beta": 28., "m0": -1.143, "m1": -.714}, (-.7,0.,0.), {"frequency_hz": 1., "amplitude": 0.}, "Initial-condition companion to the canonical double scroll."),
}

@dataclass(frozen=True)
class ChuaCircuit(ODESystem):
    driver: SinusoidalDriver; alpha: float=15.6; beta: float=28.; m0: float=-1.143; m1: float=-.714; initial_state_override: tuple[float,...]|None=None; name: str="chua"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (.1,0.,0.))
    def _h(self,x): return self.m1*x+.5*(self.m0-self.m1)*(abs(x+1)-abs(x-1))
    def _dh(self,x): return self.m0 if abs(x)<1 else self.m1
    @property
    def control_matrix(self): return np.array(((1.,),(0.,),(0.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        q,y,z=x; return np.array((self.alpha*(y-q-self._h(q)),q-y+z,-self.beta*y))+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x): q,_,_=x; return np.array(((-self.alpha*(1+self._dh(q)),self.alpha,0.),(1.,-1.,1.),(0.,-self.beta,0.)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
