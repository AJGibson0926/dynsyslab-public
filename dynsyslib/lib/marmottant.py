"""Two-state Marmottant coated-bubble model from the retained MATLAB suite.

The continuous state remains ``(R, Rdot)``. Shell rupture/reformation is the
same discrete hysteresis latch used by the MATLAB implementation; it is not a
third state and can be reset explicitly between integrations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from ..observables import Observable
from .presets import SystemPreset
from .rpe import BUBBLE_PRESETS, BubblePreset


@dataclass(frozen=True)
class MarmottantDimensionlessParameters:
    time_scale_s: float
    natural_frequency_hz: float
    liquid_mach: float
    gas_mach: float
    euler: float
    vapor_pressure_ratio: float
    capillarity: float
    viscosity: float
    gas_pressure_ratio: float
    shell_viscosity: float
    buckling_radius: float
    rupture_radius: float
    break_radius: float
    jacobian_stiffness: float
    jacobian_damping: float
    convention: str="matlab-bubble-marmottant/v1"


MARMOTTANT_PRESETS={
    "matlab_cleanup_linear":SystemPreset(
        "matlab_cleanup_linear",
        "Linear-shell Marmottant setup used by generate_training.m.",
        {"equilibrium_radius_m":.975e-6,"ambient_pressure_pa":101e3,
         "shell_viscosity_kg_s":15e-9,"buckling_radius_m":.975e-6,
         "critical_surface_tension_n_m":1.,"shell_elasticity_n_m":1.,
         "liquid_sound_speed_m_s":1480.,"gas_sound_speed_m_s":0.,
         "shell_law":"linear"},
        None,{"frequency_hz":1.,"amplitude":0.},
        "MATLAB cleanup suite: Bubble.m and generate_training.m.",
        {"training":{"seed":2141442,"trajectory_count":100,
         "initial_condition_ranges":{"radius":[0.,1.],"radial_velocity":[-10.,10.]},
         "alternate_radius_range":[0.,1.5],"time_span":[0.,2.],"delta_tau":.001,
         "coordinate_system":"dimensionless"}}),
    "matlab_cleanup_long":SystemPreset(
        "matlab_cleanup_long",
        "Long-horizon Marmottant experiment retained in generate_training.m.",
        {"equilibrium_radius_m":.82e-6,"shell_viscosity_kg_s":7.2e-9,
         "critical_surface_tension_n_m":.13,"shell_elasticity_n_m":1.,
         "shell_law":"linear"},(1.,0.),{"frequency_hz":2e6,"amplitude":200e3},
        "MATLAB cleanup suite: generate_training.m (forcing call is commented out).",
        {"training":{"trajectory_count":1,"time_span":[0.,600.],"delta_tau":.001,
         "coordinate_system":"dimensionless","forcing_status":"configured only in commented code"}}),
}


@dataclass
class MarmottantSystem(ODESystem):
    driver: SinusoidalDriver
    preset: BubblePreset=BUBBLE_PRESETS["marmottant_matlab_reference"]
    equilibrium_radius_m: float=.975e-6
    ambient_pressure_pa: float=101e3
    shell_viscosity_kg_s: float=15e-9
    buckling_radius_m: float=.975e-6
    critical_surface_tension_n_m: float=1.
    shell_elasticity_n_m: float=1.
    liquid_sound_speed_m_s: float=1480.
    gas_sound_speed_m_s: float=0.
    shell_law: Literal["linear","sigmoid","gompertz"]="linear"
    initial_state_override: tuple[float,float]|None=None
    name: str="marmottant"
    prefer_dimensionless_form: bool=True
    control_affine: bool=True
    _ruptured: bool=field(default=False,init=False,repr=False,compare=False)

    def __post_init__(self):
        positive=(self.equilibrium_radius_m,self.ambient_pressure_pa,self.shell_viscosity_kg_s,
                  self.buckling_radius_m,self.critical_surface_tension_n_m,self.shell_elasticity_n_m)
        if any(not np.isfinite(value) or value<=0 for value in positive):
            raise ValueError("Marmottant radii, pressures, shell viscosity, and elasticity must be positive.")
        if self.shell_law not in {"linear","sigmoid","gompertz"}:
            raise ValueError("shell_law must be 'linear', 'sigmoid', or 'gompertz'.")

    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self):
        state=self.initial_state_override if self.initial_state_override is not None else (1.,0.)
        return np.asarray(state,dtype=float)
    @property
    def dimensionless_initial_state(self): return np.asarray(self.initial_state,dtype=float)

    def reset_hysteresis(self, *, ruptured: bool=False) -> None:
        """Reset the discrete shell latch without changing the two-state model."""
        self._ruptured=bool(ruptured)

    @property
    def hysteresis_ruptured(self) -> bool: return bool(self._ruptured)

    @property
    def pg_reference_pa(self) -> float:
        p=self.preset
        return self.ambient_pressure_pa-p.vapor_pressure_pa+2*p.sigma_n_m/self.equilibrium_radius_m

    @property
    def dimensionless_parameters(self) -> MarmottantDimensionlessParameters:
        p=self.preset; r0=self.equilibrium_radius_m
        stiffness=3*p.kappa*self.pg_reference_pa-2*p.sigma_n_m/r0
        if stiffness<=0: raise ValueError("Marmottant natural-frequency stiffness must be positive.")
        frequency=.5*np.sqrt(stiffness)/(2*np.pi*r0*np.sqrt(p.rho_kg_m3)); T=1/frequency
        sigma0=self.shell_elasticity_n_m*((r0/self.buckling_radius_m)**2-1.)
        pg0=self.ambient_pressure_pa+2*sigma0/r0
        ml=0. if self.liquid_sound_speed_m_s==0 else r0/(self.liquid_sound_speed_m_s*T)
        mg=0. if self.gas_sound_speed_m_s==0 else r0/(3*self.gas_sound_speed_m_s*T)
        eu=self.ambient_pressure_pa*T*T/(p.rho_kg_m3*r0*r0)
        shell=4*self.shell_viscosity_kg_s/(self.ambient_pressure_pa*r0*T)
        rrupt=self.buckling_radius_m*np.sqrt(1+p.sigma_n_m/self.shell_elasticity_n_m)
        rbreak=self.buckling_radius_m*np.sqrt(1+self.critical_surface_tension_n_m/self.shell_elasticity_n_m)
        j1=(3*p.kappa*pg0+4*self.shell_elasticity_n_m/r0)*T*T/(p.rho_kg_m3*r0*r0)
        j2=(3*p.kappa*pg0*r0*(0. if self.liquid_sound_speed_m_s==0 else 1/self.liquid_sound_speed_m_s)
            +4*p.mu_pa_s+4*self.shell_viscosity_kg_s/r0)*T/(p.rho_kg_m3*r0*r0)
        return MarmottantDimensionlessParameters(T,frequency,ml,mg,eu,
            p.vapor_pressure_pa/self.ambient_pressure_pa,
            2*p.sigma_n_m/(self.ambient_pressure_pa*r0),4*p.mu_pa_s/(self.ambient_pressure_pa*T),
            pg0/self.ambient_pressure_pa,shell,self.buckling_radius_m/r0,
            rrupt/r0,rbreak/r0,j1,j2)

    def surface_tension(self, radius: float, *, update_hysteresis: bool=True) -> float:
        """Return shell tension for nondimensional radius, preserving MATLAB branches."""
        r=float(radius); g=self.dimensionless_parameters
        if update_hysteresis:
            if r<g.rupture_radius: self._ruptured=False
            elif r>=g.break_radius: self._ruptured=True
        if r<g.buckling_radius:
            return 0.
        if g.buckling_radius<=r<=g.rupture_radius:
            if self.shell_law=="linear":
                return self.shell_elasticity_n_m*((r/g.buckling_radius)**2-1.)
            if self.shell_law=="sigmoid":
                scaled=(r-g.buckling_radius)/(g.rupture_radius-g.buckling_radius)
                c=3.; f=lambda value,mid:value**c/mid**(c-1)
                shaped=f(scaled,.5) if scaled<=.5 else 1-f(1-scaled,.5)
                return self.preset.sigma_n_m*shaped
            sigma0=self.preset.sigma_n_m/2
            c=(2*self.shell_elasticity_n_m*np.e/self.preset.sigma_n_m)*np.sqrt(1+self.preset.sigma_n_m/(2*self.shell_elasticity_n_m))
            b=-np.log(sigma0/self.preset.sigma_n_m)/np.exp(c*(1-1/g.buckling_radius))
            return self.preset.sigma_n_m*np.exp(-b*np.exp(c*(1-r/g.buckling_radius)))
        if g.rupture_radius<r<g.break_radius:
            if self.shell_law=="linear":
                return self.preset.sigma_n_m if self._ruptured else self.shell_elasticity_n_m*((r/g.buckling_radius)**2-1.)
        if r>g.break_radius:
            return self.preset.sigma_n_m
        # This equality branch is retained from Bubble.m rather than silently
        # changing the legacy equation at the switching surface.
        return 0.

    def _dimensionless_terms(self,state,*,update_hysteresis: bool):
        r,v=np.asarray(state,dtype=float); p=self.preset; g=self.dimensionless_parameters
        sigma_r=self.surface_tension(r,update_hysteresis=update_hysteresis)
        inertial=-1.5*v*v
        gas=(1-3*p.kappa*g.liquid_mach*v)*g.euler*g.gas_pressure_ratio*r**(-3*p.kappa)
        laplace=g.euler*(2*sigma_r/(self.ambient_pressure_pa*self.equilibrium_radius_m))/r
        viscous=g.euler*g.viscosity*v/r
        shell=g.euler*g.shell_viscosity*v/(r*r)
        momentum=inertial+gas-(g.euler+laplace+viscous+shell)
        return sigma_r,shell,inertial,gas,laplace,viscous,momentum

    def dimensionless_input_matrix(self,state):
        radius=float(np.asarray(state)[0]); return np.asarray(((0.,),(-self.dimensionless_parameters.euler/radius,)))
    def control_matrix(self,state):
        radius=float(np.asarray(state)[0]); return np.asarray(((0.,),(-1/(self.preset.rho_kg_m3*radius),)))
    def dimensionless_rhs_with_input(self,tau,state,control_input):
        r,v=np.asarray(state,dtype=float); momentum=self._dimensionless_terms(state,update_hysteresis=True)[-1]
        u=float(np.asarray(control_input,dtype=float).reshape(-1)[0])
        return np.asarray((v,momentum/r))+self.dimensionless_input_matrix(state)@np.asarray((u,))
    def dimensionless_rhs(self,tau,state):
        time=float(self.to_dimensional_time(tau)); u=self.driver.value(time)/self.ambient_pressure_pa
        return self.dimensionless_rhs_with_input(tau,state,u)
    def rhs_with_input(self,t,state,pressure_pa):
        R,V=np.asarray(state,dtype=float); z=self.to_dimensionless_state((R,V)); tau=float(self.to_dimensionless_time(t))
        dz=self.dimensionless_rhs_with_input(tau,z,float(np.asarray(pressure_pa).reshape(-1)[0])/self.ambient_pressure_pa)
        T=self.dimensionless_parameters.time_scale_s; r0=self.equilibrium_radius_m
        return np.asarray((dz[0]*r0/T,dz[1]*r0/(T*T)))
    def rhs(self,t,state): return self.rhs_with_input(t,state,self.driver.value(t))

    def to_dimensionless_time(self,time): return np.asarray(time)/self.dimensionless_parameters.time_scale_s
    def to_dimensional_time(self,time): return np.asarray(time)*self.dimensionless_parameters.time_scale_s
    def to_dimensionless_state(self,state):
        values=np.asarray(state,dtype=float); T=self.dimensionless_parameters.time_scale_s; r0=self.equilibrium_radius_m
        return np.asarray((values[0]/r0,values[1]*T/r0))
    def to_dimensional_state(self,state):
        values=np.asarray(state,dtype=float); T=self.dimensionless_parameters.time_scale_s; r0=self.equilibrium_radius_m
        return np.asarray((values[0]*r0,values[1]*r0/T))
    def to_dimensionless_input(self,value): return np.asarray(value,dtype=float)/self.ambient_pressure_pa
    def to_dimensional_input(self,value): return np.asarray(value,dtype=float)*self.ambient_pressure_pa
    def control_input_units(self,coordinate_system): return "p_control/p0" if coordinate_system=="dimensionless" else "Pa"
    def control_input_scale_to_dimensional(self,coordinate_system): return self.ambient_pressure_pa if coordinate_system=="dimensionless" else 1.

    def _finite_jacobian(self,function,state):
        x=np.asarray(state,dtype=float); result=np.empty((2,2)); held=self._ruptured
        try:
            for index in range(2):
                step=1e-6*max(1.,abs(x[index])); plus=x.copy(); minus=x.copy(); plus[index]+=step; minus[index]-=step
                self._ruptured=held; fp=function(plus); self._ruptured=held; fm=function(minus)
                result[:,index]=(fp-fm)/(2*step)
        finally: self._ruptured=held
        return result
    def jacobian(self,t,state): return self._finite_jacobian(lambda x:self.rhs(t,x),state)
    def dimensionless_jacobian(self,tau,state): return self._finite_jacobian(lambda x:self.dimensionless_rhs(tau,x),state)
    def invalid(self,t,state):
        x=np.asarray(state,dtype=float)
        return not np.all(np.isfinite(x)) or x[0]<=0.

    def state_component_names(self,coordinate_system=None): return ("radius","radial_velocity")
    def observable_aliases(self,coordinate_system=None): return {"r":"radius","R":"radius","rdot":"radial_velocity","Rdot":"radial_velocity"}
    def state_component_labels(self,coordinate_system=None):
        return ("R (m)","dR/dt (m/s)") if coordinate_system=="dimensional" else ("r = R/R0","dr/dtau")
    def model_observables(self,coordinate_system=None):
        selected=coordinate_system or "dimensionless"; validation=(self.equilibrium_radius_m,0.) if selected=="dimensional" else (1.,0.)
        def zstate(x): return self.to_dimensionless_state(x) if selected=="dimensional" else np.asarray(x,dtype=float)
        def terms(x): return self._dimensionless_terms(zstate(x),update_hysteresis=False)
        def item(key,label,index,units="dimensionless"):
            return Observable(key,1,2,lambda x,i=index:np.asarray((terms(x)[i],)),None,
                {"kind":"marmottant-observable","label":label,"units":units,
                 "state_coordinate_system":selected,"matlab_model":"Marmottant"},validation)
        return {
            "surface_tension":item("surface_tension","Surface tension",0,"N/m"),
            "shell_viscosity":item("shell_viscosity","Shell viscosity",1),
            "inertial_kinetic_energy":item("inertial_kinetic_energy","Inertial kinetic energy",2),
            "internal_gas_pressure":item("internal_gas_pressure","Internal gas pressure",3),
            "laplace_pressure":item("laplace_pressure","Laplace pressure",4),
            "viscous_damping":item("viscous_damping","Viscous damping",5),
            "momentum_balance":item("momentum_balance","Momentum balance",6),
        }


__all__=["MARMOTTANT_PRESETS","MarmottantDimensionlessParameters","MarmottantSystem"]
