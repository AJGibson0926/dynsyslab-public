"""Frozen pre-nondimensionalization dimensional RPE implementation."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .rpe import BubblePreset, BUBBLE_PRESETS

@dataclass(frozen=True)
class DimensionalBaselineRayleighPlessetSystem(ODESystem):
    driver: SinusoidalDriver
    preset: BubblePreset=BUBBLE_PRESETS["air_water_20c"]
    equilibrium_radius_m: float=10e-6
    ambient_pressure_pa: float=101325.
    initial_state_override: tuple[float,float]|None=None
    name: str="rayleigh-plesset-dimensional-baseline"
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (self.equilibrium_radius_m,0.))
    @property
    def period(self): return self.driver.period
    @property
    def pg0(self): return self.ambient_pressure_pa-self.preset.vapor_pressure_pa+2*self.preset.sigma_n_m/self.equilibrium_radius_m
    def rhs(self,t,x):
        r,v=x;p=self.preset;pg=self.pg0*(self.equilibrium_radius_m/r)**(3*p.kappa)
        jump=pg+p.vapor_pressure_pa-(self.ambient_pressure_pa-self.driver.value(t))-2*p.sigma_n_m/r-4*p.mu_pa_s*v/r
        return np.array((v,(jump/p.rho_kg_m3-1.5*v*v)/r))
    def jacobian(self,t,x):
        r,v=x;p=self.preset;pg=self.pg0*(self.equilibrium_radius_m/r)**(3*p.kappa);jump=pg+p.vapor_pressure_pa-(self.ambient_pressure_pa-self.driver.value(t))-2*p.sigma_n_m/r-4*p.mu_pa_s*v/r;q=jump/p.rho_kg_m3-1.5*v*v
        return np.array(((0.,1.),((r*(-3*p.kappa*pg/r+2*p.sigma_n_m/r**2+4*p.mu_pa_s*v/r**2)/p.rho_kg_m3-q)/r**2,(-4*p.mu_pa_s/(p.rho_kg_m3*r)-3*v)/r)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or x[0]<=.01*self.equilibrium_radius_m
