"""Dimensional and nondimensional Rayleigh--Plesset equations."""
from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Literal
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from ..observables import Observable

@dataclass(frozen=True)
class BubblePreset:
    key: str; description: str; temperature_c: float
    rho_kg_m3: float; mu_pa_s: float; sigma_n_m: float; vapor_pressure_pa: float
    gas: str="air"; kappa: float=1.4; caution: str=""
    def construction_kwargs(self) -> dict[str,object]: return {"preset":self}

BUBBLE_PRESETS={
 "air_water_0c":BubblePreset("air_water_0c","Air in ultrapure water",0,999.8,1.792e-3,.0756,611),
 "air_water_20c":BubblePreset("air_water_20c","Air in ultrapure water",20,998.,1.e-3,.0725,2330,kappa=4/3),
 "air_water_37c":BubblePreset("air_water_37c","Air in ultrapure water",37,993.3,.691e-3,.0696,6280),
 "air_saline_37c":BubblePreset("air_saline_37c","Air in 0.9% saline",37,1004.6,.75e-3,.0730,6280),
 "air_blood_37c":BubblePreset("air_blood_37c","Air in blood/plasma, Newtonian approximation",37,1025.,3.5e-3,.055,6280,caution="Blood is non-Newtonian; use only as a baseline."),
 "air_glycerol_20c":BubblePreset("air_glycerol_20c","Air in anhydrous glycerol",20,1261.,1.49,.0634,.01,caution="Very viscous: incompressible RPE may be inadequate."),
 "air_glycerol_water_50pct_20c":BubblePreset("air_glycerol_water_50pct_20c","Air in 50 wt% glycerol/water",20,1126.,.006,.067,1800,caution="Mixture values are nominal."),
 "air_silicone_oil_100cst_25c":BubblePreset("air_silicone_oil_100cst_25c","Air in 100 cSt PDMS oil",25,965.,.0965,.020,1.,caution="Air-oil interfacial tension is formulation-dependent."),
 "co2_water_20c":BubblePreset("co2_water_20c","CO2 in water",20,998.2,1.002e-3,.0728,2339,gas="carbon dioxide",kappa=1.30,caution="Gas dissolution/thermal effects are omitted."),
 "co2_honey_25c":BubblePreset("co2_honey_25c","CO2 in honey",25,1420.,10.,.055,3200,gas="carbon dioxide",kappa=1.30,caution="Order-of-magnitude Newtonian surrogate; honey is non-Newtonian and RPE is not validated here."),
 "marmottant_matlab_reference":BubblePreset("marmottant_matlab_reference","Marmottant cleanup-suite liquid/gas preset",20,1000.,1.e-3,.0725,0.,kappa=1.095,caution="Retained exactly from the MATLAB Bubble.m Marmottant preset."),
}

@dataclass(frozen=True)
class DimensionalRPEParameters:
    rho_kg_m3: float; mu_pa_s: float; sigma_n_m: float; vapor_pressure_pa: float; kappa: float
    equilibrium_radius_m: float; ambient_pressure_pa: float; initial_gas_pressure_pa: float

@dataclass(frozen=True)
class DimensionlessRPEParameters:
    """Groups for ``r=R/R0`` and ``tau=f0*t`` using the natural period ``1/f0``."""
    time_scale_s: float; natural_frequency_hz: float; pressure_scale_pa: float
    gas_pressure: float; vapor_pressure: float; ambient_pressure: float
    capillarity: float; viscosity: float
    control_pressure_scale_pa: float
    convention: str="rpe-natural-period/v1"
    pressure_input_convention: str="external-acoustic-pressure/v2"
    control_input_convention: str="p_control/p0"

RPEPressureInputConvention = Literal[
    "external-acoustic-pressure/v2",
    "legacy-positive-acceleration/v1",
]

@dataclass(frozen=True)
class RayleighPlessetSystem(ODESystem):
    """RPE with equivalent dimensional and dimensionless RHS/Jacobian calls."""
    driver: SinusoidalDriver
    preset: BubblePreset=BUBBLE_PRESETS["air_water_20c"]
    equilibrium_radius_m: float=10e-6
    ambient_pressure_pa: float=100000.
    initial_state_override: tuple[float,float]|None=None
    pressure_input_convention: RPEPressureInputConvention="external-acoustic-pressure/v2"
    name: str="rayleigh-plesset"
    prefer_dimensionless_form: bool=True
    def __post_init__(self):
        if self.pressure_input_convention not in {"external-acoustic-pressure/v2","legacy-positive-acceleration/v1"}:
            raise ValueError("Unknown RPE pressure_input_convention.")
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (self.equilibrium_radius_m,0.))
    # Positive external acoustic pressure compresses the bubble. The legacy
    # convention is retained only for exact reconstruction of older archives.
    control_affine: bool=True
    @property
    def pressure_input_sign(self) -> float:
        return -1. if self.pressure_input_convention=="external-acoustic-pressure/v2" else 1.
    def control_matrix(self,x):
        radius=float(np.asarray(x,dtype=float)[0])
        return np.array(((0.,),(self.pressure_input_sign/(self.preset.rho_kg_m3*radius),)))
    @property
    def pg0(self):
        p=self.preset; return self.ambient_pressure_pa-p.vapor_pressure_pa+2*p.sigma_n_m/self.equilibrium_radius_m
    @property
    def dimensional_parameters(self):
        p=self.preset; return DimensionalRPEParameters(p.rho_kg_m3,p.mu_pa_s,p.sigma_n_m,p.vapor_pressure_pa,p.kappa,self.equilibrium_radius_m,self.ambient_pressure_pa,self.pg0)
    @property
    def dimensionless_parameters(self):
        d=self.dimensional_parameters
        stiffness=3*d.kappa*d.initial_gas_pressure_pa-2*d.sigma_n_m/d.equilibrium_radius_m
        if stiffness<=0: raise ValueError("RPE natural-frequency stiffness must be positive.")
        frequency=np.sqrt(stiffness)/(2*np.pi*d.equilibrium_radius_m*np.sqrt(d.rho_kg_m3)); T=1/frequency
        pressure_scale=d.rho_kg_m3*(d.equilibrium_radius_m/T)**2
        return DimensionlessRPEParameters(T,frequency,pressure_scale,
            d.initial_gas_pressure_pa/pressure_scale,d.vapor_pressure_pa/pressure_scale,d.ambient_pressure_pa/pressure_scale,
            2*d.sigma_n_m/(d.equilibrium_radius_m*pressure_scale),4*d.mu_pa_s/(T*pressure_scale),d.ambient_pressure_pa,
            pressure_input_convention=self.pressure_input_convention)
    def with_dimensional_parameters(self,**changes):
        """Return a new system and therefore freshly derived dimensionless groups."""
        system_fields={"equilibrium_radius_m","ambient_pressure_pa","initial_state_override","pressure_input_convention"}; preset_fields={"rho_kg_m3","mu_pa_s","sigma_n_m","vapor_pressure_pa","kappa","gas","temperature_c","description","caution","key"}
        unknown=set(changes)-system_fields-preset_fields
        if unknown: raise ValueError(f"Unknown dimensional RPE parameters: {sorted(unknown)}")
        preset_changes={k:v for k,v in changes.items() if k in preset_fields}; system_changes={k:v for k,v in changes.items() if k in system_fields}
        return replace(self,preset=replace(self.preset,**preset_changes),**system_changes)
    def to_dimensionless_time(self,t): return np.asarray(t)/self.dimensionless_parameters.time_scale_s
    def to_dimensional_time(self,tau): return np.asarray(tau)*self.dimensionless_parameters.time_scale_s
    def to_dimensionless_state(self,x):
        x=np.asarray(x,dtype=float); T=self.dimensionless_parameters.time_scale_s; return np.array((x[0]/self.equilibrium_radius_m,x[1]*T/self.equilibrium_radius_m))
    def to_dimensional_state(self,x):
        x=np.asarray(x,dtype=float); T=self.dimensionless_parameters.time_scale_s; return np.array((x[0]*self.equilibrium_radius_m,x[1]*self.equilibrium_radius_m/T))
    def to_dimensionless_input(self,pressure_pa):
        """Convert the RPE control pressure to ``u*=p_control/p0``."""
        return np.asarray(pressure_pa,dtype=float)/self.ambient_pressure_pa
    def to_dimensional_input(self,control_input):
        """Convert ``u*=p_control/p0`` back to pressure in pascals."""
        return np.asarray(control_input,dtype=float)*self.ambient_pressure_pa
    def control_input_units(self,coordinate_system: str) -> str:
        if coordinate_system=="dimensionless": return "p_control/p0"
        if coordinate_system=="dimensional": return "Pa"
        raise ValueError("coordinate_system must be 'dimensionless' or 'dimensional'.")
    def control_input_scale_to_dimensional(self,coordinate_system: str) -> float:
        if coordinate_system=="dimensionless": return float(self.ambient_pressure_pa)
        if coordinate_system=="dimensional": return 1.
        raise ValueError("coordinate_system must be 'dimensionless' or 'dimensional'.")
    @property
    def dimensionless_initial_state(self): return self.to_dimensionless_state(self.initial_state)
    def state_component_names(self,coordinate_system=None): return ("radius","radial_velocity")
    def observable_aliases(self,coordinate_system=None): return {"r":"radius","R":"radius","rdot":"radial_velocity","Rdot":"radial_velocity"}
    def state_component_labels(self,coordinate_system=None):
        if coordinate_system=="dimensional": return ("R (m)","dR/dt (m/s)")
        return ("r = R/R0","dr/dtau")
    def model_observables(self,coordinate_system=None):
        """MATLAB Bubble/RPE physical-balance observables.

        The legacy catalog is defined in nondimensional bubble coordinates.
        For a dimensional trajectory the input state is converted first, so
        the reported quantities retain exactly one canonical interpretation.
        """
        selected=coordinate_system or "dimensionless"; g=self.dimensionless_parameters; p=self.preset
        validation=(self.equilibrium_radius_m,0.) if selected=="dimensional" else (1.,0.)
        def zstate(x):
            values=np.asarray(x,dtype=float)
            return self.to_dimensionless_state(values) if selected=="dimensional" else values
        def terms(x):
            r,v=zstate(x); inertial=1.5*v*v
            gas=g.gas_pressure*r**(-3*p.kappa)
            laplace=g.capillarity/r; viscous=g.viscosity*v/r
            momentum=gas-inertial+g.vapor_pressure-(g.ambient_pressure+laplace+viscous)
            echo=self.equilibrium_radius_m*g.pressure_scale_pa*r*(2*v*v+momentum)
            return inertial,gas,laplace,viscous,momentum,echo
        def item(key,label,index,units="dimensionless"):
            return Observable(key,1,2,lambda x,i=index:np.asarray((terms(x)[i],)),None,
                {"kind":"rpe-observable","label":label,"units":units,
                 "state_coordinate_system":selected,"matlab_model":"Rayleigh-Plesset"},validation)
        return {
            "inertial_kinetic_energy":item("inertial_kinetic_energy","Inertial kinetic energy",0),
            "internal_gas_pressure":item("internal_gas_pressure","Internal gas pressure",1),
            "laplace_pressure":item("laplace_pressure","Laplace pressure",2),
            "viscous_damping":item("viscous_damping","Viscous damping",3),
            "momentum_balance":item("momentum_balance","Momentum balance",4),
            "far_field_acoustic_echo":item("far_field_acoustic_echo","Far-field acoustic echo",5,"Pa m"),
        }
    def rhs_with_input(self,t,x,drive_pa):
        r,v=x;p=self.preset;pg=self.pg0*(self.equilibrium_radius_m/r)**(3*p.kappa)
        drive=float(np.asarray(drive_pa,dtype=float).reshape(-1)[0])
        jump=pg+p.vapor_pressure_pa-self.ambient_pressure_pa-2*p.sigma_n_m/r-4*p.mu_pa_s*v/r
        return np.array((v,(jump/p.rho_kg_m3-1.5*v*v)/r))+self.input_matrix(x)@np.asarray((drive,))
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    def jacobian(self,t,x):
        r,v=x;p=self.preset;pg=self.pg0*(self.equilibrium_radius_m/r)**(3*p.kappa);jump=pg+p.vapor_pressure_pa-self.ambient_pressure_pa+self.pressure_input_sign*self.driver.value(t)-2*p.sigma_n_m/r-4*p.mu_pa_s*v/r;q=jump/p.rho_kg_m3-1.5*v*v
        return np.array(((0.,1.),((r*(-3*p.kappa*pg/r+2*p.sigma_n_m/r**2+4*p.mu_pa_s*v/r**2)/p.rho_kg_m3-q)/r**2,(-4*p.mu_pa_s/(p.rho_kg_m3*r)-3*v)/r)))
    def dimensionless_control_matrix(self,x):
        radius=float(np.asarray(x,dtype=float)[0])
        return np.array(((0.,),(self.pressure_input_sign*self.dimensionless_parameters.ambient_pressure/radius,)))
    def dimensionless_input_matrix(self,x):
        """Return ``d(dz/dtau)/d(u*)`` for ``u*=p_control/p0``."""
        return self.dimensionless_control_matrix(x)
    def dimensionless_rhs_with_input(self,tau,x,control_input):
        r,u=x;p=self.preset;g=self.dimensionless_parameters
        drive=float(np.asarray(control_input,dtype=float).reshape(-1)[0])
        numerator=g.gas_pressure*r**(-3*p.kappa)+g.vapor_pressure-g.ambient_pressure-g.capillarity/r-g.viscosity*u/r-1.5*u*u
        return np.array((u,numerator/r))+self.dimensionless_input_matrix(x)@np.asarray((drive,))
    def dimensionless_rhs(self,tau,x):
        pressure_pa=self.driver.value(float(self.to_dimensional_time(tau)))
        return self.dimensionless_rhs_with_input(tau,x,self.to_dimensionless_input(pressure_pa))
    def dimensionless_jacobian(self,tau,x):
        r,u=x;p=self.preset;g=self.dimensionless_parameters;A=g.gas_pressure*r**(-3*p.kappa)+g.vapor_pressure-g.ambient_pressure+self.pressure_input_sign*self.driver.value(float(self.to_dimensional_time(tau)))/g.pressure_scale_pa-g.capillarity/r-g.viscosity*u/r
        dA=-3*p.kappa*g.gas_pressure*r**(-3*p.kappa-1)+g.capillarity/r**2+g.viscosity*u/r**2
        return np.array(((0.,1.),((r*dA-A+1.5*u*u)/r**2,(-g.viscosity/r-3*u)/r)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or x[0]<=.01*self.equilibrium_radius_m
