"""Forced Duffing oscillator system definition."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ..drivers import SinusoidalDriver
from ..interfaces import ODESystem
from ..observables import Observable
from .presets import SystemPreset

DUFFING_PRESETS = {
    "double_well_reference": SystemPreset("double_well_reference", "Canonical weakly damped double-well Duffing oscillator.", {"delta": .2, "alpha": 1., "beta": 1.}, (.1, 0.), {"frequency_hz": 1.2, "amplitude": .3}, "Canonical nondimensional Duffing benchmark."),
    "hardening_single_well": SystemPreset("hardening_single_well", "Single-well, hardening-restoring-force variant in this sign convention.", {"delta": .2, "alpha": -1., "beta": 1.}, (.2, 0.), {"frequency_hz": 1., "amplitude": .2}, "Standard Duffing sign-convention variant."),
    "gibson_undamped_klqr_2022": SystemPreset("gibson_undamped_klqr_2022", "Undamped double-well Duffing trajectory used in the Gibson--Calvisi--Yee complex-Koopman-LQR example.", {"delta": 0., "alpha": 1., "beta": 1.}, (0., -2.8), {"frequency_hz": 1., "amplitude": 0.}, "Gibson, Calvisi & Yee, Koopman Linear Quadratic Regulator Using Complex Eigenfunctions for Nonlinear Dynamical Systems, SIADS 21 (2022), Sec. 4.1.2, doi:10.1137/21M1456078.", {"control": {"initial_state": [0., -2.8], "activation_time": 3.85, "target_state": [1., 0.], "paper_parameters": {"alpha": 1., "beta": -1., "delta": 0.}, "dynsyslab_mapping": "paper alpha -> dynsyslab beta; paper beta -> -dynsyslab alpha"}, "hdmd": {"observable": {"state_index": 0, "name": "x"}, "training_initial_conditions": "Three random trajectories were used for the damping comparison; their coordinates are not reported."}}),
    "gibson_undamped_vdp_tracking_2022": SystemPreset("gibson_undamped_vdp_tracking_2022", "Undamped double-well Duffing system used by Gibson--Calvisi--Yee for tracking a separately simulated Van der Pol target trajectory.", {"delta": 0., "alpha": 1., "beta": 1.}, None, {"frequency_hz": 1., "amplitude": 0.}, "Gibson, Calvisi & Yee, Koopman Linear Quadratic Regulator Using Complex Eigenfunctions for Nonlinear Dynamical Systems, SIADS 21 (2022), Sec. 4.1.4, doi:10.1137/21M1456078. The paper restates the physical parameters but does not state a new Duffing initial condition.", {"control": {"target_source": "separately simulated Van der Pol trajectory", "eigenfunction_representation": "one complex-conjugate pair", "paper_parameters": {"alpha": 1., "beta": -1., "delta": 0.}, "initial_condition_reported": False, "dynsyslab_mapping": "paper alpha -> dynsyslab beta; paper beta -> -dynsyslab alpha"}}),
    "gibson_damped_stabilization_2022": SystemPreset("gibson_damped_stabilization_2022", "Damped double-well Duffing case used by Gibson--Calvisi--Yee to stabilize a non-equilibrium target in the opposite potential well.", {"delta": .5, "alpha": 1., "beta": 4.}, (.5, 0.), {"frequency_hz": 1., "amplitude": 0.}, "Gibson, Calvisi & Yee, Koopman Linear Quadratic Regulator Using Complex Eigenfunctions for Nonlinear Dynamical Systems, SIADS 21 (2022), Sec. 4.1.5, doi:10.1137/21M1456078.", {"control": {"initial_state": [.5, 0.], "active_time_window": [1., 10.], "target_state": [-.25, 0.], "paper_parameters": {"alpha": 4., "beta": -1., "delta": .5}, "dynsyslab_mapping": "paper alpha -> dynsyslab beta; paper beta -> -dynsyslab alpha", "source_note": "The retained MATLAB cleanup experiment activates control on [1,10]; this supersedes the earlier [1,5] transcription."}, "koopman": {"eigenvalue_method": "Korda et al. non-recurrent-set method", "eigenfunction_method": "SINDy with monomial basis", "selected_complex_eigenvalues": 8, "training_ensemble": {"trajectory_count": 30, "forcing": "unforced", "initial_conditions": "randomly generated around the unit circle"}}}),
    "korda_mezic_optimal_2020": SystemPreset("korda_mezic_optimal_2020", "Canonical DynSysLab form exactly equivalent to the Duffing dynamics in Korda--Mezic's March 2020 optimal-eigenfunction code.", {"delta": .5, "alpha": 1., "beta": 4., "input_gain": .5}, None, {"frequency_hz": 1., "amplitude": 0.}, "Korda & Mezic, Optimal construction of Koopman eigenfunctions for prediction and control, IEEE TAC; companion Duffing.m dated March 2020.", {"free_training": {"seed": 2141442, "trajectory_count": 100, "initial_conditions": "normalized two-dimensional randn draws", "delta_t": .01, "sample_count": 800, "final_sample_time": 7.99}, "controlled_training": {"top_level_seed": 2141442, "exact_control_stream": "depends on intervening MATLAB RNG use; export Ubig for parity", "trajectory_count": 100, "input_distribution": "independent uniform [-1,1] in source convention", "dynsyslab_equivalent_input_distribution": "independent uniform [-0.5,0.5] when input_gain=1", "delta_t": .01, "transition_count": 200, "hold": "one value per fixed RK4 interval"}, "equation_reconciliation": {"canonical": "qdot=v; vdot=-0.5*v+q-4*q^3+0.5*u", "autonomous_equivalent": True, "input_gain_difference_from_standard_dynsyslab": .5}}),
}

@dataclass(frozen=True)
class DuffingDimensionlessParameters:
    """Identity scaling for the canonically nondimensional Duffing model."""
    time_scale_s: float=1.
    state_scale: tuple[float,float]=(1.,1.)
    input_scale: float=1.

@dataclass(frozen=True)
class DuffingSystem(ODESystem):
    driver: SinusoidalDriver; delta: float=.2; alpha: float=1.; beta: float=1.; input_gain: float=1.; initial_state_override: tuple[float,...]|None=None; name: str="duffing"
    def __post_init__(self):
        if not np.isfinite(self.input_gain): raise ValueError("Duffing input_gain must be finite.")
    @property
    def period(self): return self.driver.period
    @property
    def initial_state(self): return np.array(self.initial_state_override if self.initial_state_override is not None else (0.,0.))
    @property
    def control_matrix(self): return np.array(((0.,),(self.input_gain,)))
    control_affine: bool=True
    # Preserve the established native/default solver route.  Callers may set
    # use_dimensionless_form=True when they need explicit coordinate metadata;
    # the equations and transforms are identity-equivalent either way.
    prefer_dimensionless_form: bool=False
    def rhs_with_input(self,t,x,u):
        q,v=x; return np.array((v,-self.delta*v+self.alpha*q-self.beta*q**3))+self.input_matrix(x)@np.asarray(u,dtype=float).reshape(-1)
    def rhs(self,t,x): return self.rhs_with_input(t,x,self.driver.value(t))
    @property
    def dimensionless_parameters(self): return DuffingDimensionlessParameters()
    @property
    def dimensionless_initial_state(self): return np.asarray(self.initial_state,dtype=float)
    def to_dimensionless_state(self,state): return np.asarray(state,dtype=float)
    def to_dimensional_state(self,state): return np.asarray(state,dtype=float)
    def to_dimensionless_time(self,time): return np.asarray(time,dtype=float)
    def to_dimensional_time(self,time): return np.asarray(time,dtype=float)
    def dimensionless_input_matrix(self,state): return self.input_matrix(state)
    def dimensionless_rhs_with_input(self,t,x,u): return self.rhs_with_input(t,x,u)
    def dimensionless_rhs(self,t,x): return self.rhs(t,x)
    def control_input_units(self,coordinate_system): return "dimensionless forcing" if coordinate_system=="dimensionless" else "system input"
    def control_input_scale_to_dimensional(self,coordinate_system): return 1.
    def state_component_names(self,coordinate_system=None): return ("position","velocity")
    def observable_aliases(self,coordinate_system=None): return {"x":"position","xdot":"velocity"}
    def state_component_labels(self,coordinate_system=None): return ("x","dx/dt")
    def model_observables(self,coordinate_system=None):
        selected=coordinate_system or "dimensionless"
        def scalar(key,label,function,jacobian=None,*,value_type="real"):
            return Observable(key,1,2,lambda x:np.asarray((function(np.asarray(x)),)),
                None if jacobian is None else lambda x:np.asarray((jacobian(np.asarray(x)),)),
                {"kind":"duffing-observable","label":label,"units":"dimensionless",
                 "state_coordinate_system":selected,"matlab_model":"Duffing"},
                (0.,0.),value_type)
        kinetic=lambda x:.5*x[1]**2
        potential=lambda x:-.5*self.alpha*x[0]**2+.25*self.beta*x[0]**4
        random_observation=lambda x:np.exp(1j*x[0])/np.cosh(x[1])
        return {
            "kinetic_energy":scalar("kinetic_energy","Kinetic energy",kinetic,lambda x:(0.,x[1])),
            "potential_energy":scalar("potential_energy","Potential energy",potential,
                lambda x:(-self.alpha*x[0]+self.beta*x[0]**3,0.)),
            "total_energy":scalar("total_energy","Total energy",lambda x:kinetic(x)+potential(x),
                lambda x:(-self.alpha*x[0]+self.beta*x[0]**3,x[1])),
            "random_observation":scalar("random_observation","Random observation",
                random_observation,lambda x:(1j*random_observation(x),-random_observation(x)*np.tanh(x[1])),
                value_type="complex"),
        }
    def jacobian(self,t,x): q,_=x; return np.array(((0.,1.),(self.alpha-3*self.beta*q*q,-self.delta)))
    def invalid(self,t,x): return not np.all(np.isfinite(x)) or np.linalg.norm(x)>1e8
