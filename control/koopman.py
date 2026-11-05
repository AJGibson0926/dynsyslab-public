"""ROM-facing Koopman controller contracts; no dependency on postprocessing internals."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import warnings
import numpy as np
from .core import LQRController, MPCController, Target
from rom import KROM

def _rom_boundary(rom: KROM, coordinate_system: str) -> tuple[str,str]:
    if coordinate_system not in {"dimensionless","dimensional"}:
        raise ValueError("coordinate_system must be 'dimensionless' or 'dimensional'.")
    if coordinate_system=="dimensionless" and rom.coordinate_system!="dimensionless":
        raise ValueError("A dimensional KROM cannot control dimensionless dynamics without an explicit model-coordinate transformation.")
    return ("model","model") if coordinate_system=="dimensionless" else ("physical","physical")

@dataclass
class KoopmanLiftedTarget:
    """Lift a model-state target into a KROM's nonredundant real coordinates."""
    rom: KROM
    state_target: Target
    coordinate_system: str="dimensionless"
    def bind_coordinate_system(self,coordinate_system: str) -> None:
        _rom_boundary(self.rom,coordinate_system); self.coordinate_system=coordinate_system
    def evaluate(self,t: float,state_dimension: int) -> np.ndarray:
        state_coordinates,time_coordinates=_rom_boundary(self.rom,self.coordinate_system)
        model_time=float(t) if time_coordinates=="model" else float(np.asarray(self.rom._model_time(np.asarray(t),"physical")))
        physical_dimension=int(np.asarray(self.rom.system.initial_state).size)
        model_state=np.asarray(self.state_target.evaluate(model_time,physical_dimension),dtype=float)
        lifted=np.asarray(self.rom.psi(model_state,output="split",state_coordinates="model"),dtype=float)
        if lifted.size!=state_dimension: raise ValueError("Lifted target dimension does not match the KROM controller state.")
        return lifted

@dataclass
class KoopmanDisplayTarget:
    """Expose a model-state target in physical coordinates for archive/Viewer use."""
    rom: KROM
    state_target: Target
    coordinate_system: str="dimensionless"
    def bind_coordinate_system(self,coordinate_system: str) -> None:
        _rom_boundary(self.rom,coordinate_system); self.coordinate_system=coordinate_system
    def evaluate(self,t: float,state_dimension: int) -> np.ndarray:
        _,time_coordinates=_rom_boundary(self.rom,self.coordinate_system)
        model_time=float(t) if time_coordinates=="model" else float(np.asarray(self.rom._model_time(np.asarray(t),"physical")))
        model_state=np.asarray(self.state_target.evaluate(model_time,state_dimension),dtype=float)
        if self.coordinate_system=="dimensional" and self.rom.coordinate_system=="dimensionless":
            return np.asarray(self.rom.system.to_dimensional_state(model_state),dtype=float)
        return model_state

class KoopmanLQRController(LQRController):
    """State-aware KLQR law over an executable KROM and control map."""
    def __init__(self,rom: KROM,Q,R,target: Target|None=None,active_time_window: tuple[float,float]|None=None,*,target_space: str="state",control_limit=None,target_domain_policy: str="error",control_representation: str="default"):
        self.rom=rom
        if target_space not in {"state","lifted"}: raise ValueError("target_space must be 'state' or 'lifted'.")
        self.state_target=target if target_space=="state" else None
        initial_coordinates="dimensionless" if rom.coordinate_system=="dimensionless" else "dimensional"
        lifted_target=KoopmanLiftedTarget(rom,target,initial_coordinates) if target is not None and target_space=="state" else target
        self.display_target=KoopmanDisplayTarget(rom,target,initial_coordinates) if target is not None and target_space=="state" else None
        self.control_limit=None if control_limit is None else np.asarray(control_limit,dtype=float).reshape(-1)
        if self.control_limit is not None and (not self.control_limit.size or np.any(~np.isfinite(self.control_limit)) or np.any(self.control_limit<=0)):
            raise ValueError("control_limit must contain positive finite bounds in the controller's active coordinates.")
        if target_domain_policy not in {"error","warn","ignore"}:
            raise ValueError("target_domain_policy must be 'error', 'warn', or 'ignore'.")
        self.target_domain_policy=target_domain_policy; self._domain_warnings=set()
        selected=rom.default_control_representation if control_representation=="default" else control_representation
        if selected not in rom.available_control_representations():
            raise ValueError(f"KLQR control representation {selected!r} is unavailable; available={rom.available_control_representations()}.")
        self.control_representation=selected
        super().__init__(lambda x:rom.A(output="split",time_coordinates="model"),lambda x:rom.B_psi(x,output="split",state_coordinates="model",time_coordinates="model",representation=self.control_representation),Q,R,lifted_target,active_time_window=active_time_window)
        self.actuation="system_input"
        self.target_space=target_space
        self.bind_coordinate_system(initial_coordinates,system=rom.system)
    def bind_coordinate_system(self,coordinate_system: str,*,system=None) -> None:
        if system is not None and system is not self.rom.system:
            raise ValueError("A Koopman controller must be bound to the same analytic system as its KROM.")
        _rom_boundary(self.rom,coordinate_system)
        self.coordinate_system=coordinate_system
        self.time_units=coordinate_system; self.state_units=coordinate_system
        if isinstance(self.target,KoopmanLiftedTarget): self.target.bind_coordinate_system(coordinate_system)
        if self.display_target is not None: self.display_target.bind_coordinate_system(coordinate_system)
        if getattr(self.rom.system,"name",None)=="rayleigh-plesset":
            self.control_units="p_control/p0" if coordinate_system=="dimensionless" else "Pa"
        else: self.control_units="dimensionless system input" if coordinate_system=="dimensionless" else "system input"
    def archive_metadata(self) -> dict[str,Any]:
        target=self.state_target if self.state_target is not None else self.target
        target_name=getattr(target,"name",None)
        return {"controller":"koopman-lqr","Q":np.asarray(self.Q,dtype=float),"R":np.asarray(self.R,dtype=float),
            "target":target_name or (None if target is None else type(target).__name__),"target_space":self.target_space,
            "active_time_window":self.active_time_window,"actuation":self.actuation,"coordinate_system":self.coordinate_system,"time_units":self.time_units,
            "state_units":self.state_units,"control_limit":self.control_limit,"control_units":self.control_units,
            "control_representation":self.control_representation,
            "target_domain_policy":self.target_domain_policy,"training_domain":self.rom.training_domain,
            "rom":{"method":self.rom.method,"source_archive":self.rom.spectrum.metadata.get("loaded_from_archive"),
            "coordinate_system":self.rom.coordinate_system,"retained_eigenvalues":[complex(value.continuous_time) for value in self.rom.eigenvalues],
            "control_system_binding":self.rom.spectrum.metadata.get("control_system_binding"),
            "pressure_input_convention":getattr(self.rom.system,"pressure_input_convention",None)}}
    def value(self,t,state):
        active=np.asarray(state,dtype=float)
        self._enforce_domain(active,t)
        state_coordinates,time_coordinates=_rom_boundary(self.rom,self.coordinate_system)
        lifted=np.asarray(self.rom.psi(active,output="split",state_coordinates=state_coordinates),dtype=float)
        reference=np.zeros_like(lifted) if self.target is None else self.target.evaluate(t,lifted.size)
        A=self.rom.A(output="split",time_coordinates=time_coordinates); B=self.rom.B_psi(active,output="split",state_coordinates=state_coordinates,time_coordinates=time_coordinates,representation=self.control_representation)
        P=__import__("scipy.linalg",fromlist=["solve_continuous_are"]).solve_continuous_are(A,B,self.Q,self.R)
        gain=np.linalg.solve(self.R,B.T@P); self.last_gain=gain
        u=-gain@(lifted-reference)
        if self.control_limit is not None:
            limits=self.control_limit
            if limits.size==1: limits=np.full(u.size,limits.item())
            if limits.size!=u.size: raise ValueError("control_limit must be scalar or match the system input dimension.")
            u=np.clip(u,-limits,limits)
        return np.asarray(u,dtype=float)

    def _enforce_domain(self,state: np.ndarray,t: float) -> None:
        if self.target_domain_policy=="ignore" or self.rom.training_domain is None: return
        state_coordinates,time_coordinates=_rom_boundary(self.rom,self.coordinate_system)
        reports=[("state",self.rom.domain_report(state,state_coordinates=state_coordinates))]
        if self.state_target is not None:
            model_time=float(t) if time_coordinates=="model" else float(np.asarray(self.rom._model_time(np.asarray(t),"physical")))
            dimension=int(np.asarray(self.rom.system.initial_state).size)
            target_state=np.asarray(self.state_target.evaluate(model_time,dimension),dtype=float)
            reports.append(("target",self.rom.domain_report(target_state,state_coordinates="model")))
        violations=[(name,report) for name,report in reports if report.get("available") and not report.get("inside")]
        if not violations: return
        detail="; ".join(f"{name} axes={report['violating_axes']} excursion={report['maximum_relative_excursion']:.6g}" for name,report in violations)
        message=f"KLQR evaluation is outside the KROM training domain: {detail}."
        if self.target_domain_policy=="error": raise ValueError(message+" Set target_domain_policy='warn' or 'ignore' only for an explicit extrapolation study.")
        key=tuple(name for name,_ in violations)
        if key not in self._domain_warnings:
            warnings.warn(message,RuntimeWarning,stacklevel=2); self._domain_warnings.add(key)

class KoopmanMPCController(MPCController):
    """Finite-horizon Koopman MPC using a selected KROM input map."""
    def __init__(self,rom: KROM,Q,R,target: Target|None=None,horizon:int=20,sample_time:float=1e-3,control_limit=None,active_time_window: tuple[float,float]|None=None,*,target_space: str="state",target_domain_policy: str="error",control_representation: str="default"):
        self.rom=rom
        if target_space not in {"state","lifted"}: raise ValueError("target_space must be 'state' or 'lifted'.")
        self.state_target=target if target_space=="state" else None
        initial_coordinates="dimensionless" if rom.coordinate_system=="dimensionless" else "dimensional"
        lifted_target=KoopmanLiftedTarget(rom,target,initial_coordinates) if target is not None and target_space=="state" else target
        self.display_target=KoopmanDisplayTarget(rom,target,initial_coordinates) if target is not None and target_space=="state" else None
        if target_domain_policy not in {"error","warn","ignore"}: raise ValueError("target_domain_policy must be 'error', 'warn', or 'ignore'.")
        self.target_domain_policy=target_domain_policy; self._domain_warnings=set()
        selected=rom.default_control_representation if control_representation=="default" else control_representation
        if selected not in rom.available_control_representations():
            raise ValueError(f"KMPC control representation {selected!r} is unavailable; available={rom.available_control_representations()}.")
        self.control_representation=selected
        super().__init__(lambda x:rom.A(output="split",time_coordinates="model"),lambda x:rom.B_psi(x,output="split",state_coordinates="model",time_coordinates="model",representation=self.control_representation),Q,R,lifted_target,active_time_window=active_time_window,horizon=horizon,sample_time=sample_time,control_limit=control_limit)
        self.actuation="system_input"
        self.target_space=target_space
        KoopmanLQRController.bind_coordinate_system(self,initial_coordinates,system=rom.system)
    def bind_coordinate_system(self,coordinate_system: str,*,system=None) -> None:
        KoopmanLQRController.bind_coordinate_system(self,coordinate_system,system=system)
    def archive_metadata(self) -> dict[str,Any]:
        metadata=KoopmanLQRController.archive_metadata(self); metadata["controller"]="koopman-mpc"
        metadata.update({"horizon":self.horizon,"sample_time":self.sample_time,"control_limit":self.control_limit})
        return metadata
    def value(self,t,state):
        active=np.asarray(state,dtype=float); KoopmanLQRController._enforce_domain(self,active,t)
        state_coordinates,time_coordinates=_rom_boundary(self.rom,self.coordinate_system)
        lifted=np.asarray(self.rom.psi(active,output="split",state_coordinates=state_coordinates),dtype=float)
        reference=np.zeros_like(lifted) if self.target is None else self.target.evaluate(t,lifted.size)
        A=self.rom.A(output="split",time_coordinates=time_coordinates); B=self.rom.B_psi(active,output="split",state_coordinates=state_coordinates,time_coordinates=time_coordinates,representation=self.control_representation)
        Ad=np.eye(A.shape[0])+self.sample_time*A; Bd=self.sample_time*B; P=np.asarray(self.Q,dtype=float).copy()
        for _ in range(self.horizon):
            gain=np.linalg.solve(self.R+Bd.T@P@Bd,Bd.T@P@Ad); P=np.asarray(self.Q,dtype=float)+Ad.T@P@(Ad-Bd@gain)
        self.last_gain=gain; u=-gain@(lifted-reference)
        if self.control_limit is not None:
            limits=np.asarray(self.control_limit,dtype=float).reshape(-1)
            if limits.size==1: limits=np.full(u.size,limits.item())
            u=np.clip(u,-limits,limits)
        return np.asarray(u,dtype=float)
