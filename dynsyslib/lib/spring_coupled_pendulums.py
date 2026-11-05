"""Two nonlinear pendula coupled by an extensible spring between their bobs."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

SPRING_COUPLED_PENDULUM_PRESETS={
    "symmetric_weak_coupling": SystemPreset("symmetric_weak_coupling","Identical pendula with weak bob-to-bob spring coupling.",{"mass1_kg":1.,"mass2_kg":1.,"length1_m":1.,"length2_m":1.,"pivot_separation_m":1.,"spring_constant_n_m":.5,"spring_rest_length_m":1.,"damping1_n_m_s":.01,"damping2_n_m_s":.01},(.45,0.,-.2,0.),{"frequency_hz":.5,"amplitude":0.},"Useful for in-phase and anti-phase normal-mode studies."),
    "strong_nonlinear_coupling": SystemPreset("strong_nonlinear_coupling","Large angles and stiff spring, retaining exact bob geometry.",{"mass1_kg":1.,"mass2_kg":1.,"length1_m":1.,"length2_m":1.,"pivot_separation_m":1.,"spring_constant_n_m":20.,"spring_rest_length_m":1.,"damping1_n_m_s":.005,"damping2_n_m_s":.005},(1.1,0.,-.9,0.),{"frequency_hz":.7,"amplitude":.2},"Not the small-angle coupled-pendulum approximation."),
    "detuned_energy_exchange": SystemPreset("detuned_energy_exchange","Unequal pendula for beats and nonlinear energy-transfer experiments.",{"mass1_kg":.8,"mass2_kg":1.2,"length1_m":.85,"length2_m":1.15,"pivot_separation_m":1.,"spring_constant_n_m":3.,"spring_rest_length_m":1.,"damping1_n_m_s":.01,"damping2_n_m_s":.015},(.7,0.,.05,0.),{"frequency_hz":.55,"amplitude":0.},"Geometric coupling produces state-dependent effective stiffness."),
}

@dataclass(frozen=True)
class SpringCoupledPendulumsSystem(ODESystem):
    """State ``(theta1, omega1, theta2, omega2)``; driver torque acts on bob 1."""
    driver: SinusoidalDriver
    mass1_kg: float=1.; mass2_kg: float=1.; length1_m: float=1.; length2_m: float=1.; pivot_separation_m: float=1.
    spring_constant_n_m: float=1.; spring_rest_length_m: float=1.; gravity_m_s2: float=9.80665; damping1_n_m_s: float=.01; damping2_n_m_s: float=.01
    initial_state_override: tuple[float,float,float,float]|None=None; name: str="spring_coupled_pendulums"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.asarray(self.initial_state_override if self.initial_state_override is not None else (.45,0.,-.2,0.),dtype=float)
    @property
    def control_matrix(self): return np.array(((0.,),(1/(self.mass1_kg*self.length1_m**2),),(0.,),(0.,)))
    control_affine: bool=True
    def _spring_gradients(self,a,b):
        l1,l2=self.length1_m,self.length2_m
        difference=np.array((self.pivot_separation_m+l2*np.sin(b)-l1*np.sin(a),l1*np.cos(a)-l2*np.cos(b)))
        distance=float(np.linalg.norm(difference))
        if distance<1e-12: raise FloatingPointError("Spring bobs coincide; spring direction is undefined.")
        derivative1=np.array((-l1*np.cos(a),-l1*np.sin(a)))
        derivative2=np.array((l2*np.cos(b),l2*np.sin(b)))
        scale=self.spring_constant_n_m*(distance-self.spring_rest_length_m)/distance
        return scale*np.dot(difference,derivative1),scale*np.dot(difference,derivative2)
    def rhs_with_input(self,t,x,torque_n_m):
        a,w1,b,w2=np.asarray(x,dtype=float); grad1,grad2=self._spring_gradients(a,b)
        torque=float(np.asarray(torque_n_m,dtype=float).reshape(-1)[0])
        dd1=-self.gravity_m_s2/self.length1_m*np.sin(a)-(self.damping1_n_m_s*w1+grad1)/(self.mass1_kg*self.length1_m**2)
        dd2=-self.gravity_m_s2/self.length2_m*np.sin(b)-(self.damping2_n_m_s*w2+grad2)/(self.mass2_kg*self.length2_m**2)
        return np.array((w1,dd1,w2,dd2))+self.input_matrix(x)@np.asarray((torque,))
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x):
        state=np.asarray(x,dtype=float); step=1e-6*np.maximum(1.,np.abs(state)); matrix=np.empty((4,4))
        for index,h in enumerate(step):
            plus=state.copy(); minus=state.copy(); plus[index]+=h; minus[index]-=h
            matrix[:,index]=(self.rhs(t,plus)-self.rhs(t,minus))/(2*h)
        return matrix
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
