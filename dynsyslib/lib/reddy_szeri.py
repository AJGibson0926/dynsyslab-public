"""Twelve-state Reddy--Szeri bubble model retained from MATLAB Bubble.m."""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from ..observables import Observable
from ..analytic_expressions import state_variables
from .presets import SystemPreset
from .rpe import BUBBLE_PRESETS, BubblePreset


REDDY_SZERI_STATE_NAMES=("radius","radial_velocity","translation","translation_velocity",
    "shape_2","shape_2_velocity","shape_3","shape_3_velocity",
    "shape_4","shape_4_velocity","shape_5","shape_5_velocity")


@dataclass(frozen=True)
class ReddySzeriDimensionlessParameters:
    time_scale_s: float
    natural_frequency_hz: float
    liquid_mach: float
    gas_mach: float
    euler: float
    vapor_pressure_ratio: float
    capillarity: float
    viscosity: float
    gas_pressure_ratio: float
    pressure_amplitude: float
    pressure_derivative_amplitude: float
    angular_frequency: float
    wave_number: float
    legacy_acoustic_pressure_amplitude: float
    convention: str="matlab-bubble-reddy-szeri/v1"


REDDY_SZERI_PRESETS={
    "matlab_cleanup_reference":SystemPreset(
        "matlab_cleanup_reference",
        "Reddy--Szeri translation and n=2..5 shape-mode setup in generate_training.m.",
        {"equilibrium_radius_m":2e-6,"translation_enabled":True,
         "acoustic_pressure_amplitude":0.,"drive_frequency_hz":850.75e3,
         "liquid_sound_speed_m_s":0.,"gas_sound_speed_m_s":0.},
        (.524736182055959,0.,0.,.1,.1,0.,.1,0.,.1,0.,.1,0.),
        {"frequency_hz":850.75e3,"amplitude":0.},
        "MATLAB cleanup suite: Bubble.m and generate_training.m.",
        {"training":{"trajectory_count":1,"time_span":[0.,10.],"delta_tau":.001,
         "coordinate_system":"dimensionless","initial_condition_source":"explicit"}}),
    "matlab_cleanup_control":SystemPreset(
        "matlab_cleanup_control",
        "Lower-amplitude initial shape/translation state used by control_dynamics.m.",
        {"equilibrium_radius_m":2e-6,"translation_enabled":True,
         "acoustic_pressure_amplitude":0.,"drive_frequency_hz":850.75e3,
         "liquid_sound_speed_m_s":0.,"gas_sound_speed_m_s":0.},
        (1.,0.,0.,.01,.01,0.,.01,0.,.01,0.,.01,0.),
        {"frequency_hz":850.75e3,"amplitude":50e3},
        "MATLAB cleanup suite: control_dynamics.m; the 50 kPa forcing call in generate_training.m is commented out."),
}


@dataclass(frozen=True)
class ReddySzeriSystem(ODESystem):
    driver: SinusoidalDriver
    preset: BubblePreset=BUBBLE_PRESETS["air_water_20c"]
    equilibrium_radius_m: float=2e-6
    ambient_pressure_pa: float=100e3
    liquid_sound_speed_m_s: float=0.
    gas_sound_speed_m_s: float=0.
    translation_enabled: bool=True
    acoustic_pressure_amplitude: float=0.
    drive_frequency_hz: float=850.75e3
    initial_state_override: tuple[float,...]|None=None
    name: str="reddy-szeri"
    prefer_dimensionless_form: bool=True
    control_affine: bool=True

    @property
    def period(self): return 1/self.drive_frequency_hz
    @property
    def initial_state(self):
        default=(.524736182055959,0.,0.,.1,.1,0.,.1,0.,.1,0.,.1,0.)
        state=np.asarray(self.initial_state_override if self.initial_state_override is not None else default,dtype=float)
        if state.shape!=(12,): raise ValueError("Reddy--Szeri initial state must contain 12 components.")
        return state
    @property
    def dimensionless_initial_state(self): return self.initial_state.copy()
    @property
    def pg0(self):
        p=self.preset; return self.ambient_pressure_pa-p.vapor_pressure_pa+2*p.sigma_n_m/self.equilibrium_radius_m
    @property
    def dimensionless_parameters(self):
        p=self.preset; r0=self.equilibrium_radius_m
        stiffness=3*p.kappa*self.pg0-2*p.sigma_n_m/r0
        if stiffness<=0: raise ValueError("Reddy--Szeri natural-frequency stiffness must be positive.")
        frequency=np.sqrt(stiffness)/(2*np.pi*r0*np.sqrt(p.rho_kg_m3)); T=1/frequency
        ml=0. if self.liquid_sound_speed_m_s==0 else r0/(self.liquid_sound_speed_m_s*T)
        mg=0. if self.gas_sound_speed_m_s==0 else r0/(3*self.gas_sound_speed_m_s*T)
        eu=self.ambient_pressure_pa*T*T/(p.rho_kg_m3*r0*r0)
        omega=2*np.pi*self.drive_frequency_hz*T
        pressure=-eu*ml*self.acoustic_pressure_amplitude
        pressure_dt=-omega*eu*ml*self.acoustic_pressure_amplitude
        wavelength=0. if self.liquid_sound_speed_m_s==0 else self.liquid_sound_speed_m_s/self.drive_frequency_hz
        wave_number=0. if wavelength==0 else 2*np.pi*r0/wavelength
        return ReddySzeriDimensionlessParameters(T,frequency,ml,mg,eu,
            p.vapor_pressure_pa/self.ambient_pressure_pa,
            2*p.sigma_n_m/(self.ambient_pressure_pa*r0),4*p.mu_pa_s/(self.ambient_pressure_pa*T),
            self.pg0/self.ambient_pressure_pa,pressure,pressure_dt,omega,wave_number,
            self.acoustic_pressure_amplitude)

    @staticmethod
    def _effective_kappa(radius,kappa): return np.where(np.asarray(radius)>1.,1.,kappa)
    def _volume_terms(self,state):
        x=np.asarray(state,dtype=float); R=x[0]; V=x[1]; p=self.preset; g=self.dimensionless_parameters
        k=float(self._effective_kappa(R,p.kappa))
        t1=-1.5*V*V
        t2=g.euler*g.gas_pressure_ratio*(1-3*k*g.liquid_mach*V)*R**(-3*k)
        t3=3*k*g.euler*g.gas_pressure_ratio*g.gas_mach*V*((1.5*(k+1)*g.liquid_mach*V)-1)*R**(-1.5*(k+1))
        laplace=g.euler*g.capillarity*(g.liquid_mach*V-1)/R
        viscous=g.euler*g.viscosity*(V/R)*(1-g.liquid_mach*V)
        denominator=R+g.euler*g.liquid_mach*g.viscosity+3*k*g.euler*g.gas_pressure_ratio*g.gas_mach*g.liquid_mach*R**(-(3*k+1)/2)
        balance=t1+t2+t3-g.euler-viscous+laplace
        return t1,t2,t3,laplace,viscous,denominator,balance,balance/denominator

    def dimensionless_input_matrix(self,state):
        x=np.asarray(state,dtype=float); R=x[0]; p=self.preset; g=self.dimensionless_parameters
        k=float(self._effective_kappa(R,p.kappa))
        b=R+g.euler*g.liquid_mach*g.viscosity+3*k*g.euler*g.gas_pressure_ratio*g.gas_mach*g.liquid_mach*R**(-(3*k+1)/2)
        d=R*R+g.euler*g.liquid_mach*g.viscosity*R+3*k*g.euler*g.gas_pressure_ratio*g.gas_mach*g.liquid_mach*R**(-(3*k-1)/2)
        return -np.asarray((0.,g.euler/b,0.,0.,0.,g.euler*x[4]/d,0.,2*g.euler*x[6]/d,
                           0.,3*g.euler*x[8]/d,0.,4*g.euler*x[10]/d)).reshape(12,1)
    def input_matrix(self,state): return self.dimensionless_input_matrix(state)

    def dimensionless_rhs_with_input(self,tau,state,control_input):
        x=np.asarray(state,dtype=float); R,V=x[:2]; g=self.dimensionless_parameters
        u=float(np.asarray(control_input,dtype=float).reshape(-1)[0]); dx=np.zeros(12)
        dx[0]=V; dx[1]=self._volume_terms(x)[-1]+self.dimensionless_input_matrix(x)[1,0]*u
        if self.translation_enabled:
            position,velocity=x[2],x[3]
            relative=velocity+g.pressure_amplitude*np.sin(g.angular_frequency*tau)
            acceleration=-(3*V/R+4.5*g.euler*g.viscosity/R**2)*relative+3*g.pressure_derivative_amplitude*np.cos(g.angular_frequency*tau)
            relative_prime=-(g.pressure_derivative_amplitude*np.cos(g.angular_frequency*tau)+acceleration)
            dx[2]=velocity; dx[3]=acceleration
        else:
            relative=0.; relative_prime=0.
        c=3*V/R
        def k(n): return (.5*g.euler*g.capillarity*(n-1)*(n+1)*(n+2)/R**3
                          -(n-1)*dx[1]/R-4.5*n*(n+1)**2/((2*n-1)*(2*n+3))*(relative/R)**2)
        def d(n): return 9*n*(n-2)*(n-1)*(n+1)/(4*(2*n-3)*(2*n-1))*(relative/R)**2
        def f(n): return 3*(n+1)*(2*(n-1)*relative*V+n*R*relative_prime)/(2*(2*n-1)*R**2)
        def h(n): return 3*(n+1)**2*(-2*relative*V+R*relative_prime)/(2*(2*n+3)*R**2)
        def j(n): return 9*(n+1)**2*(n+2)*(n+3)/(4*(2*n+3)*(2*n+5))*(relative/R)**2
        def q(n): return 1.5*(n+1)*relative/R
        dx[4]=x[5]
        dx[5]=-k(2)*x[4]-c*x[5]+h(2)*x[6]+j(2)*x[8]-2.25*relative**2/R
        dx[6]=x[7]
        dx[7]=f(3)*x[4]+q(3)*x[5]-k(3)*x[6]-c*x[7]+h(3)*x[8]+j(3)*x[10]
        dx[8]=x[9]
        dx[9]=d(4)*x[4]+f(4)*x[6]+q(4)*x[7]-k(4)*x[8]-c*x[9]+h(4)*x[10]
        dx[10]=x[11]
        dx[11]=d(5)*x[6]+f(5)*x[8]+q(5)*x[9]-k(5)*x[10]-c*x[11]
        return dx
    def dimensionless_rhs(self,tau,state):
        time=float(self.to_dimensional_time(tau)); return self.dimensionless_rhs_with_input(tau,state,self.driver.value(time)/self.ambient_pressure_pa)
    # Bubble.m stores and advances this model natively in nondimensional form.
    def rhs_with_input(self,t,state,u): return self.dimensionless_rhs_with_input(t,state,u)
    def rhs(self,t,state): return self.dimensionless_rhs(t,state)
    def to_dimensionless_time(self,time): return np.asarray(time)/self.dimensionless_parameters.time_scale_s
    def to_dimensional_time(self,time): return np.asarray(time)*self.dimensionless_parameters.time_scale_s
    def to_dimensionless_state(self,state): return np.asarray(state,dtype=float)
    def to_dimensional_state(self,state): return np.asarray(state,dtype=float)
    def to_dimensionless_input(self,value): return np.asarray(value,dtype=float)/self.ambient_pressure_pa
    def to_dimensional_input(self,value): return np.asarray(value,dtype=float)*self.ambient_pressure_pa
    def control_input_units(self,coordinate_system): return "p_control/p0" if coordinate_system=="dimensionless" else "native Reddy-Szeri input"
    def control_input_scale_to_dimensional(self,coordinate_system): return self.ambient_pressure_pa if coordinate_system=="dimensionless" else 1.

    def jacobian(self,t,state):
        x=np.asarray(state,dtype=float); result=np.empty((12,12))
        for index in range(12):
            step=1e-6*max(1.,abs(x[index])); plus=x.copy(); minus=x.copy(); plus[index]+=step; minus[index]-=step
            result[:,index]=(self.rhs(t,plus)-self.rhs(t,minus))/(2*step)
        return result
    dimensionless_jacobian=jacobian
    def invalid(self,t,state):
        x=np.asarray(state,dtype=float); return not np.all(np.isfinite(x)) or x[0]<=.01 or np.linalg.norm(x)>1e10

    def state_component_names(self,coordinate_system=None): return REDDY_SZERI_STATE_NAMES
    def observable_aliases(self,coordinate_system=None):
        return {"r":"radius","R":"radius","rdot":"radial_velocity","Rdot":"radial_velocity",
                "x":"translation","xdot":"translation_velocity","Rdotdot":"radial_acceleration"}
    def state_component_labels(self,coordinate_system=None):
        return ("r = R/R0","dr/dtau","x/R0","dx/dtau","a2","da2/dtau",
                "a3","da3/dtau","a4","da4/dtau","a5","da5/dtau")
    def observable_coordinate_presets(self,coordinate_system=None):
        selected=coordinate_system or "dimensionless"
        if selected not in {"dimensionless","native"}:
            raise ValueError("Reddy--Szeri Koopman coordinates are natively dimensionless.")
        return {
            "radial_state":("radius","radial_velocity"),
            "r_rdot":("radius","radial_velocity"),
            "c_k2":("c","k2"),
            "c_k5":("c","k5"),
        }
    def _legacy_observable_terms(self,state):
        x=np.asarray(state,dtype=float); R,V=x[:2]; p=self.preset; g=self.dimensionless_parameters
        kappa=p.kappa
        t3=3*kappa*g.euler*g.gas_pressure_ratio*g.gas_mach*V*((1.5*(kappa+1)*g.liquid_mach*V)-1)*R**(-1.5*(kappa+1))
        # Bubble.m's observable definitions use the base kappa and the
        # spherical unforced acceleration, even where the EOM switches kappa.
        t1=-1.5*V*V
        t2=g.euler*g.gas_pressure_ratio*(1-3*kappa*g.liquid_mach*V)*R**(-3*kappa)
        lp=g.euler*g.capillarity*(g.liquid_mach*V-1)/R
        vc=g.euler*g.viscosity*(V/R)*(1-g.liquid_mach*V)
        D=R+g.euler*g.liquid_mach*g.viscosity+3*kappa*g.euler*g.gas_pressure_ratio*g.gas_mach*g.liquid_mach*R**(-(3*kappa+1)/2)
        rdd=(t1+t2+t3-g.euler-vc+lp)/D
        c=3*V/R; c_dt=3*(rdd/R-(V/R)**2); c_r=-3*V/R**2; c_v=3/R
        def kval(n): return .5*g.euler*g.capillarity*(n-1)*(n+1)*(n+2)/R**3-(n-1)*rdd/R
        return t3/R,rdd,c,c_dt,c_r,c_v,kval

    def _legacy_observable_expressions(self):
        """Exact expression-tree counterparts of Bubble.m's Reddy observables."""
        variables=state_variables(12); R,V=variables[:2]; p=self.preset; g=self.dimensionless_parameters
        kappa=float(p.kappa)
        t1=-1.5*V**2
        t2=g.euler*g.gas_pressure_ratio*(1-3*kappa*g.liquid_mach*V)*R**(-3*kappa)
        t3=3*kappa*g.euler*g.gas_pressure_ratio*g.gas_mach*V*((1.5*(kappa+1)*g.liquid_mach*V)-1)*R**(-1.5*(kappa+1))
        laplace=g.euler*g.capillarity*(g.liquid_mach*V-1)/R
        viscous=g.euler*g.viscosity*(V/R)*(1-g.liquid_mach*V)
        denominator=R+g.euler*g.liquid_mach*g.viscosity+3*kappa*g.euler*g.gas_pressure_ratio*g.gas_mach*g.liquid_mach*R**(-(3*kappa+1)/2)
        rdd=(t1+t2+t3-g.euler-viscous+laplace)/denominator
        c=3*V/R
        def kval(n):
            return .5*g.euler*g.capillarity*(n-1)*(n+1)*(n+2)*R**-3-(n-1)*rdd/R
        return {"compressibility_loss":t3/R,"radial_acceleration":rdd,"c":c,
                "k2":kval(2),"k3":kval(3),"k4":kval(4),"k5":kval(5),"g1":R*V}

    def model_observables(self,coordinate_system=None):
        selected=coordinate_system or "dimensionless"
        if selected not in {"dimensionless","native"}:
            raise ValueError("The retained Reddy--Szeri model and its observables are natively dimensionless.")
        return dict(self._dimensionless_model_observables)

    @cached_property
    def _dimensionless_model_observables(self):
        validation=tuple(self.initial_state); expressions=self._legacy_observable_expressions()
        axis_names=REDDY_SZERI_STATE_NAMES
        labels={"compressibility_loss":"Compressibility loss","radial_acceleration":"Radial acceleration",
                "c":"c","k2":"k2","k3":"k3","k4":"k4","k5":"k5","g1":"g1 = r rdot"}
        def scalar(key):
            return Observable.from_expressions(key,(expressions[key],),state_dimension=12,
                descriptor={"model_kind":"reddy-szeri-observable","label":labels[key],
                    "units":"dimensionless","state_coordinate_system":"dimensionless",
                    "matlab_model":"Reddy"},validation_state=validation,
                derivative_axis_names=axis_names)
        items={key:scalar(key) for key in labels}
        R,V=state_variables(12)[:2]; rdd=expressions["radial_acceleration"]
        for key in ("c","k2","k5"):
            expression=expressions[key]
            time_derivative=expression.derivative(0)*V+expression.derivative(1)*rdd
            items[f"{key}_time_derivative"]=Observable.from_expressions(
                f"{key}_time_derivative",(time_derivative,),state_dimension=12,
                descriptor={"model_kind":"reddy-szeri-observable-derivative","label":f"d{key}/dtau",
                    "units":"dimensionless","state_coordinate_system":"dimensionless",
                    "source_observable":key},validation_state=validation,derivative_axis_names=axis_names)
            partials=items[key].derivative_observables[0]
            items[f"{key}_radius_derivative"]=partials[0]
            items[f"{key}_velocity_derivative"]=partials[1]
        return items


__all__=["REDDY_SZERI_PRESETS","REDDY_SZERI_STATE_NAMES","ReddySzeriDimensionlessParameters","ReddySzeriSystem"]
