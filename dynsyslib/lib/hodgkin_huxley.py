"""Classic forced Hodgkin--Huxley neuron system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

HODGKIN_HUXLEY_PRESETS = {
    "squid_axon_rest": SystemPreset("squid_axon_rest", "Classic squid-giant-axon conductances at zero applied bias.", {"capacitance":1., "g_na":120., "g_k":36., "g_l":.3, "e_na":50., "e_k":-77., "e_l":-54.387, "bias_current":0.}, (-65.,.0529,.596,.3177), {"frequency_hz": 100., "amplitude": 0.}, "Hodgkin--Huxley 1952 parameterization."),
    "tonic_firing_10": SystemPreset("tonic_firing_10", "Classic squid-axon parameters with a 10 microamp/cm^2 DC bias.", {"capacitance":1., "g_na":120., "g_k":36., "g_l":.3, "e_na":50., "e_k":-77., "e_l":-54.387, "bias_current":10.}, (-65.,.0529,.596,.3177), {"frequency_hz": 100., "amplitude": 0.}, "Classic tonic-firing benchmark."),
}

@dataclass(frozen=True)
class HodgkinHuxleySystem(ODESystem):
    """Classic squid-axon model (V in mV, time in seconds)."""
    driver: SinusoidalDriver; capacitance: float=1.; g_na: float=120.; g_k: float=36.; g_l: float=.3
    e_na: float=50.; e_k: float=-77.; e_l: float=-54.387; bias_current: float=0.; initial_state_override: tuple[float,...]|None=None; name: str="hodgkin-huxley"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (-65.,.0529,.596,.3177))
    @staticmethod
    def _ratio(x,scale,coefficient): return coefficient*scale if abs(x)<1e-8 else coefficient*x/(-np.expm1(-x/scale))
    @staticmethod
    def _rate_derivative(rate,voltage):
        step=1e-4
        return (rate(voltage+step)-rate(voltage-step))/(2*step)
    def _rates(self,v):
        am=lambda z:self._ratio(z+40.,10.,.1); an=lambda z:self._ratio(z+55.,10.,.01)
        bm=lambda z:4*np.exp(-(z+65.)/18.); bn=lambda z:.125*np.exp(-(z+65.)/80.)
        ah=lambda z:.07*np.exp(-(z+65.)/20.); bh=lambda z:1/(1+np.exp(-(z+35.)/10.))
        return am(v),bm(v),ah(v),bh(v),an(v),bn(v),(am,bm,ah,bh,an,bn)
    @property
    def control_matrix(self): return np.array(((1000/self.capacitance,),(0.,),(0.,),(0.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,u):
        v,m,h,n=x; am,bm,ah,bh,an,bn,_=self._rates(v)
        ionic=self.g_na*m**3*h*(v-self.e_na)+self.g_k*n**4*(v-self.e_k)+self.g_l*(v-self.e_l)
        base=1000*np.array(((self.bias_current-ionic)/self.capacitance,am*(1-m)-bm*m,ah*(1-h)-bh*h,an*(1-n)-bn*n))
        return base+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x):
        v,m,h,n=x; am,bm,ah,bh,an,bn,rate_fns=self._rates(v); dam,dbm,dah,dbh,dan,dbn=[self._rate_derivative(f,v) for f in rate_fns]
        j=np.zeros((4,4)); scale=1000/self.capacitance
        j[0]=(-scale*(self.g_na*m**3*h+self.g_k*n**4+self.g_l),-scale*3*self.g_na*m*m*h*(v-self.e_na),-scale*self.g_na*m**3*(v-self.e_na),-scale*4*self.g_k*n**3*(v-self.e_k))
        j[1]=(1000*(dam*(1-m)-dbm*m),1000*(-am-bm),0.,0.)
        j[2]=(1000*(dah*(1-h)-dbh*h),0.,1000*(-ah-bh),0.)
        j[3]=(1000*(dan*(1-n)-dbn*n),0.,0.,1000*(-an-bn))
        return j
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or abs(x[0])>1e5 or np.any(x[1:]<-1) or np.any(x[1:]>2)
