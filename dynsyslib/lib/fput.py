"""Forced Fermi--Pasta--Ulam--Tsingou lattice system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

FPUT_PRESETS = {
    "alpha_fput_32": SystemPreset("alpha_fput_32", "32-site quadratic-nonlinearity FPUT chain.", {"sites":32, "alpha":.25, "beta":0.}, None, {"frequency_hz": 1., "amplitude": 0.}, "Alpha-FPUT benchmark; fixed-end chain."),
    "beta_fput_32": SystemPreset("beta_fput_32", "32-site quartic-nonlinearity FPUT chain.", {"sites":32, "alpha":0., "beta":.25}, None, {"frequency_hz": 1., "amplitude": 0.}, "Beta-FPUT benchmark; fixed-end chain."),
}

@dataclass(frozen=True)
class FPUTSystem(ODESystem):
    """Fixed-end alpha-beta Fermi--Pasta--Ulam--Tsingou lattice."""
    driver: SinusoidalDriver; sites: int=16; alpha: float=.25; beta: float=.25; initial_state_override: tuple[float,...]|None=None; name: str="fput"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self):
        if self.initial_state_override is not None: return np.array(self.initial_state_override)
        q=np.zeros(self.sites); q[0]=.1
        return np.r_[q,np.zeros(self.sites)]
    def _vprime(self,r): return r+self.alpha*r*r+self.beta*r*r*r
    def _vsecond(self,r): return 1+2*self.alpha*r+3*self.beta*r*r
    @property
    def control_matrix(self):
        matrix=np.zeros((2*self.sites,1)); matrix[self.sites,0]=1.; return matrix
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        q,p=x[:self.sites],x[self.sites:]; padded=np.r_[0.,q,0.]; bonds=np.diff(padded)
        force=self._vprime(bonds[1:])-self._vprime(bonds[:-1])
        return np.r_[p,force]+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x):
        n=self.sites; q=x[:n]; padded=np.r_[0.,q,0.]; k=self._vsecond(np.diff(padded)); j=np.zeros((2*n,2*n)); j[:n,n:]=np.eye(n)
        for i in range(n):
            j[n+i,i]=-(k[i]+k[i+1])
            if i: j[n+i,i-1]=k[i]
            if i<n-1: j[n+i,i+1]=k[i+1]
        return j
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
