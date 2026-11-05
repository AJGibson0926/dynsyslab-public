"""Forced Lorenz system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

LORENZ_PRESETS = {
    "lorenz_1963": SystemPreset("lorenz_1963", "The canonical Lorenz-1963 chaotic-flow parameterization.", {"sigma": 10., "rho": 28., "beta": 8/3}, (1.,1.,1.), {"frequency_hz": 1., "amplitude": 0.}, "Lorenz, Deterministic Nonperiodic Flow (1963)."),
    "arbabi_mezic_hdmd_2017": SystemPreset("arbabi_mezic_hdmd_2017", "Unforced Lorenz setup used for the Arbabi--Mezić Hankel-DMD illustration: observe z1 at 0.01-second intervals with a 10,000-by-500 Hankel matrix.", {"sigma": 10., "rho": 28., "beta": 8/3}, None, {"frequency_hz": 1., "amplitude": 0.}, "Arbabi & Mezić, Ergodic theory, Dynamic Mode Decomposition and Computation of Spectral Properties of the Koopman Operator (2017), Sec. 4. The paper specifies a random trajectory, not an initial state.", {"hdmd": {"settings": {"matrix_size": 501, "delay_stride": 1, "sample_stride": 1}, "observable": {"state_index": 0, "name": "z1"}, "sampling_interval": 0.01, "delay_spacing": 0.01, "delay_window": 5.0, "hankel_matrix_paper_shape": [10000, 500], "hankel_matrix_dynsyslab_shape": [501, 10000], "source_sample_count": 10500, "orientation_note": "The paper's tall Hankel matrix uses time samples as rows and f through U^500 f as columns. dynsyslab uses the transposed delay-by-snapshot orientation, so matrix_size=501 preserves the paper's 500 time-delay intervals."}}),
    "near_classic_transition": SystemPreset("near_classic_transition", "Near the classic Hopf-transition region for the nonzero equilibria.", {"sigma": 10., "rho": 24.06, "beta": 8/3}, (1.,1.,1.), {"frequency_hz": 1., "amplitude": 0.}, "Useful continuation/bifurcation starting point."),
}

@dataclass(frozen=True)
class LorenzSystem(ODESystem):
    driver: SinusoidalDriver; sigma: float=10.; rho: float=28.; beta: float=8/3; initial_state_override: tuple[float,...]|None=None; name: str="lorenz"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (1.,1.,1.))
    @property
    def control_matrix(self): return np.array(((1.,),(0.,),(0.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        q,y,z=x; return np.array((self.sigma*(y-q),q*(self.rho-z)-y,q*y-self.beta*z))+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x): q,y,z=x; return np.array(((-self.sigma,self.sigma,0.),(self.rho-z,-1.,-q),(y,q,-self.beta)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
