"""Forced Hindmarsh--Rose neuron system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

HINDMARSH_ROSE_PRESETS = {
    "canonical_bursting": SystemPreset("canonical_bursting", "Slow-fast bursting parameterization used by the original library model.", {"a":1., "b":3., "c":1., "d":5., "r":.006, "s":4., "x_r":-1.6, "current":3.25}, (-.5,-1.,3.), {"frequency_hz": 1., "amplitude": 0.}, "Hindmarsh--Rose bursting-model regime."),
    "standard_spiking": SystemPreset("standard_spiking", "Frequently used tonic-spiking parameterization.", {"a":1., "b":2.7, "c":1., "d":5., "r":.01, "s":4., "x_r":-1.6, "current":4.}, (-1.,-10.,1.), {"frequency_hz": 1., "amplitude": 0.}, "Hindmarsh & Rose (1982) style parameter set."),
}

@dataclass(frozen=True)
class HindmarshRoseSystem(ODESystem):
    """Three-state slow-fast neuron model; driver is an applied-current term."""
    driver: SinusoidalDriver; a: float=1.; b: float=3.; c: float=1.; d: float=5.; r: float=.006; s: float=4.; x_r: float=-1.6; current: float=3.25; initial_state_override: tuple[float,...]|None=None; name: str="hindmarsh-rose"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (-.5,-1.,3.))
    @property
    def control_matrix(self): return np.array(((1.,),(0.,),(0.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        q,y,z=x; return np.array((y-self.a*q**3+self.b*q*q-z+self.current,self.c-self.d*q*q-y,self.r*(self.s*(q-self.x_r)-z)))+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x): q,_,_=x; return np.array(((-3*self.a*q*q+2*self.b*q,1.,-1.),(-2*self.d*q,-1.,0.),(self.r*self.s,0.,-self.r)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
