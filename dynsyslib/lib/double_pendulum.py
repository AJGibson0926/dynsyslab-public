"""Driven, damped planar double pendulum ODE."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from .presets import SystemPreset

DOUBLE_PENDULUM_PRESETS={
    "equal_arm_demo": SystemPreset("equal_arm_demo","Equal one-metre, one-kilogram arms with light viscous joints.",{"mass1_kg":1.,"mass2_kg":1.,"length1_m":1.,"length2_m":1.,"damping1_n_m_s":.01,"damping2_n_m_s":.01},(.5,0.,.7,0.),{"frequency_hz":.5,"amplitude":0.},"Standard educational double-pendulum baseline."),
    "classic_chaotic_release": SystemPreset("classic_chaotic_release","Large release angles in the familiar sensitive double-pendulum regime.",{"mass1_kg":1.,"mass2_kg":1.,"length1_m":1.,"length2_m":1.,"damping1_n_m_s":.001,"damping2_n_m_s":.001},(2.1,0.,1.2,0.),{"frequency_hz":.5,"amplitude":0.},"Finite damping is retained for numerical regularity."),
    "unequal_lab_rig": SystemPreset("unequal_lab_rig","Unequal masses and lengths, representative of a benchtop rig.",{"mass1_kg":.8,"mass2_kg":1.2,"length1_m":.75,"length2_m":1.1,"damping1_n_m_s":.015,"damping2_n_m_s":.02},(.9,0.,-.6,.2),{"frequency_hz":.8,"amplitude":.15},"Driver torque is applied at the first joint."),
}

@dataclass(frozen=True)
class DoublePendulumSystem(ODESystem):
    """State ``(theta1, omega1, theta2, omega2)`` with torque at joint 1."""
    driver: SinusoidalDriver
    mass1_kg: float=1.; mass2_kg: float=1.; length1_m: float=1.; length2_m: float=1.; gravity_m_s2: float=9.80665
    damping1_n_m_s: float=.01; damping2_n_m_s: float=.01
    initial_state_override: tuple[float,float,float,float]|None=None; name: str="double_pendulum"
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.asarray(self.initial_state_override if self.initial_state_override is not None else (.5,0.,.7,0.),dtype=float)
    @property
    def control_matrix(self): return np.array(((0.,),(1/(self.mass1_kg*self.length1_m**2),),(0.,),(0.,)))
    control_affine: bool=True
    def rhs_with_input(self,t,x,torque_n_m):
        a,w1,b,w2=np.asarray(x,dtype=float); m1,m2,l1,l2,g=self.mass1_kg,self.mass2_kg,self.length1_m,self.length2_m,self.gravity_m_s2
        delta=a-b; denominator=2*m1+m2-m2*np.cos(2*delta)
        # Standard point-mass double pendulum equations, with independent
        # viscous joint torques and an external torque at the upper joint.
        dd1=(-g*(2*m1+m2)*np.sin(a)-m2*g*np.sin(a-2*b)-2*np.sin(delta)*m2*(w2*w2*l2+w1*w1*l1*np.cos(delta)))/(l1*denominator)
        dd2=(2*np.sin(delta)*(w1*w1*l1*(m1+m2)+g*(m1+m2)*np.cos(a)+w2*w2*l2*m2*np.cos(delta)))/(l2*denominator)
        torque=float(np.asarray(torque_n_m,dtype=float).reshape(-1)[0])
        dd1+=-self.damping1_n_m_s*w1/(m1*l1*l1)
        dd2-=self.damping2_n_m_s*w2/(m2*l2*l2)
        return np.array((w1,dd1,w2,dd2))+self.input_matrix(x)@np.asarray((torque,))
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x):
        """Central-difference Jacobian; equations retain their full nonlinear form."""
        state=np.asarray(x,dtype=float); step=1e-6*np.maximum(1.,np.abs(state)); matrix=np.empty((4,4))
        for index,h in enumerate(step):
            plus=state.copy(); minus=state.copy(); plus[index]+=h; minus[index]-=h
            matrix[:,index]=(self.rhs(t,plus)-self.rhs(t,minus))/(2*h)
        return matrix
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
