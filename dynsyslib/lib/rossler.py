"""Forced Rössler system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

ROSSLER_PRESETS = {
    "classic_chaos": SystemPreset("classic_chaos", "Rössler's widely used spiral-attractor parameterization.", {"a": .2, "b": .2, "c": 5.7}, (.1,.1,.1), {"frequency_hz": 1., "amplitude": 0.}, "Rössler (1976): a=b=0.2, c=5.7."),
    "common_alternative": SystemPreset("common_alternative", "Commonly used alternative Rössler parameterization.", {"a": .1, "b": .1, "c": 14.}, (.1,.1,.1), {"frequency_hz": 1., "amplitude": 0.}, "Frequently used Rössler alternative."),
    "brunton_shavok_2021": SystemPreset("brunton_shavok_2021", "Rossler trajectory used for Brunton-coauthored structured-HAVOK comparisons.", {"a": .1, "b": .1, "c": 14.}, (1.,1.,1.), {"frequency_hz": 1., "amplitude": 0.}, "Hirsh, Ichinaga, Brunton, Kutz & Brunton, Structured time-delay models for dynamical systems with connections to Frenet-Serret frame (2021), doi:10.1098/rspa.2021.0097.", {"havok": {"settings": {"delays": 101, "rank": 15, "observable": 0}, "observable": {"state_index": 0, "name": "x(t)"}, "sampling": {"delta_t": .001, "sample_count": 70000, "time_span": [0., 69.999], "integration_method": "adaptive DOP853 sampled onto the stated uniform grid"}, "delay_embedding": {"delay_count": 101, "delay_spacing": .001, "window_duration": .1, "rank": 15, "setting_origin": "dynsyslab reproducible implementation setting; the paper reports sampling but not these specific HAVOK truncation choices"}, "paper_reported": {"parameters": {"a": .1, "b": .1, "c": 14.}, "initial_state": [1.,1.,1.], "sample_count": 70000, "delta_t": .001, "observable": "x(t)"}}}),
    "funnel_chaos_2021": SystemPreset("funnel_chaos_2021", "Funnel-regime chaotic Rossler parameterization.", {"a": .2925, "b": .1, "c": 8.5}, (1.,1.,1.), {"frequency_hz": 1., "amplitude": 0.}, "Kraemer et al., New Journal of Physics 23, 033017 (2021): funnel regime a=0.2925, b=0.1, c=8.5. The paper samples initial conditions uniformly on [0,2]; dynsyslab fixes (1,1,1) for reproducibility."),
    "b045_chaos_2013": SystemPreset("b045_chaos_2013", "Chaotic RÃ¶ssler regime obtained by increasing b from the standard spiral setup.", {"a": .2, "b": .45, "c": 5.7}, (.1,.1,.1), {"frequency_hz": 1., "amplitude": 0.}, "Non-Linear Dynamical Classification of Short Time Series of the RÃ¶ssler System in High Noise Regimes (2013): a=0.2, c=5.7, b=0.45 described as chaotic; b=1.0 is a simple attractor."),
}

@dataclass(frozen=True)
class RosslerSystem(ODESystem):
    driver: SinusoidalDriver; a: float=.2; b: float=.2; c: float=5.7; initial_state_override: tuple[float,...]|None=None; name: str="rossler"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (.1,.1,.1))
    @property
    def control_matrix(self): return np.array(((1.,),(0.,),(0.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        q,y,z=x; return np.array((-y-z,q+self.a*y,self.b+z*(q-self.c)))+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x): q,_,z=x; return np.array(((0.,-1.,-1.),(1.,self.a,0.),(z,0.,q-self.c)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
