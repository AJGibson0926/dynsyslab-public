"""Forced Kuramoto phase-oscillator system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

KURAMOTO_PRESETS = {
    "weak_coupling": SystemPreset("weak_coupling", "Weakly coupled, broadly detuned 32-oscillator population.", {"oscillators":32, "coupling":.1, "frequency_spread":.5}, None, {"frequency_hz": 1., "amplitude": 0.}, "Kuramoto mean-field weak-coupling regime."),
    "synchronizing": SystemPreset("synchronizing", "Strong-coupling 32-oscillator synchronization regime.", {"oscillators":32, "coupling":5., "frequency_spread":.5}, None, {"frequency_hz": 1., "amplitude": 0.}, "Kuramoto mean-field synchronization regime."),
}

@dataclass(frozen=True)
class KuramotoSystem(ODESystem):
    """All-to-all coupled phase oscillators with a global sinusoidal drive."""
    driver: SinusoidalDriver; oscillators: int=16; coupling: float=1.; frequency_spread: float=.5; initial_state_override: tuple[float,...]|None=None; name: str="kuramoto"
    @property
    def period(self): return self.driver.period
    @property
    def natural_frequencies(self): return np.linspace(-self.frequency_spread,self.frequency_spread,self.oscillators)
    @property
    def initial_state(self): return np.array(self.initial_state_override) if self.initial_state_override is not None else np.linspace(0.,2*np.pi,self.oscillators,endpoint=False)
    @property
    def control_matrix(self): return np.ones((self.oscillators,1))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        differences=x[None,:]-x[:,None]
        return self.natural_frequencies+self.coupling*np.sin(differences).mean(axis=1)+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x):
        differences=x[None,:]-x[:,None]; c=self.coupling*np.cos(differences)/self.oscillators
        j=c.copy(); np.fill_diagonal(j,0.); np.fill_diagonal(j,-j.sum(axis=1)); return j
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e12
