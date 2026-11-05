"""ODE pathway: scipy adaptive methods plus transparent Euler and RK4."""
from __future__ import annotations
from time import monotonic
import numpy as np
from scipy.integrate import solve_ivp, RK23, RK45, DOP853, Radau, BDF, LSODA
from typing import Any
from provenance import serializable
from .models import NumericalResult, SolveRequest
from .backends import select_ode_backend
from .native_ode import integrate_fixed_native, native_fixed_eligibility
from .sundials import adaptive_native_eligibility, integrate_adaptive_sundials
from time_series import TimeSeries

class _ControlClock:
    """Explicit sample-and-hold controller clock, independent of RHS density."""
    def __init__(self,step: float,origin: float=0.):
        if step<=0: raise ValueError("options.control_sample_time must be positive.")
        self.step=float(step); self.origin=float(origin); self.next_tick=0; self.held=None
    def value(self,t: float,compute):
        tolerance=16.*np.finfo(float).eps*max(1.,abs(float(t)),abs(self.origin))
        scheduled=self.origin+self.next_tick*self.step
        if self.held is None or float(t)+tolerance>=scheduled:
            self.held=compute()
            elapsed=max(0.,float(t)-self.origin)
            self.next_tick=int(np.floor((elapsed+tolerance)/self.step))+1
        return self.held

def _control_evaluation_policy(request: SolveRequest) -> str:
    """Resolve and validate the controller evaluation policy once per solve."""
    policy=request.options.control_evaluation
    if policy not in {"automatic","sampled","continuous_rhs"}:
        raise ValueError("options.control_evaluation must be 'automatic', 'sampled', or 'continuous_rhs'.")
    if policy=="automatic":
        policy="sampled" if request.options.control_sample_time is not None else "continuous_rhs"
    if policy=="sampled" and request.options.control_sample_time is None:
        raise ValueError("sampled control evaluation requires options.control_sample_time.")
    if policy=="continuous_rhs" and request.options.control_sample_time is not None:
        raise ValueError("continuous_rhs control evaluation cannot also declare control_sample_time.")
    if policy=="continuous_rhs":
        controller_safe=bool(getattr(request.controller,"continuous_rhs_safe",not getattr(request.controller,"requires_sample_clock",False)))
        estimator_safe=not bool(getattr(request.estimator,"requires_sample_clock",False))
        if not controller_safe or not estimator_safe:
            unsafe=[]
            if not controller_safe: unsafe.append(type(request.controller).__name__)
            if not estimator_safe: unsafe.append(type(request.estimator).__name__)
            raise ValueError("continuous_rhs requires stateless RHS-safe control and estimation: "+", ".join(unsafe))
    return policy

def _bind_controller_coordinates(request: SolveRequest, *, dimensionless: bool) -> str | None:
    """Bind control to the exact coordinates used by the active dynamics.

    Controllers without an explicit coordinate contract inherit the solver
    coordinates. Controllers that advertise a binding hook must perform any
    necessary boundary conversion themselves; silently feeding physical-state
    controllers dimensionless states (or conversely) is forbidden.
    """
    controller=request.controller
    if controller is None: return None
    coordinates="dimensionless" if dimensionless else "dimensional"
    binder=getattr(controller,"bind_coordinate_system",None)
    if callable(binder):
        binder(coordinates,system=request.system)
        return coordinates
    declarations={getattr(controller,name,None) for name in ("coordinate_system","time_units","state_units")}
    declarations.discard(None); declarations.discard("inherit")
    normalized={"dimensional" if value=="physical" else value for value in declarations}
    if normalized and normalized!={coordinates}:
        raise ValueError(f"Controller coordinates {sorted(normalized)} do not match the active {coordinates} dynamics. Provide bind_coordinate_system() for an explicit conversion boundary.")
    return coordinates

def _control_input_metadata(request: SolveRequest,coordinate_system: str | None) -> tuple[str|None,float|None]:
    if request.controller is None or coordinate_system is None: return None,None
    if getattr(request.controller,"actuation","additive_state")=="system_input":
        units=getattr(request.system,"control_input_units",None)
        scale=getattr(request.system,"control_input_scale_to_dimensional",None)
        if callable(units):
            return str(units(coordinate_system)),float(scale(coordinate_system)) if callable(scale) else None
    fallback=getattr(request.controller,"control_units",None)
    return (str(fallback) if fallback is not None else ("dimensionless system input" if coordinate_system=="dimensionless" else "system input")),None

def _controller_signal(request: SolveRequest,t: float,y: np.ndarray,records: list[tuple[float,np.ndarray]]|None=None,clock: _ControlClock|None=None,*,dimensionless: bool=False):
    """Evaluate the configured estimator then controller at one solver state."""
    controller=request.controller
    if controller is None: return None
    controller_time=float(t)
    controller_state=np.asarray(y,dtype=float)
    window=getattr(controller,"active_time_window",None)
    if window is not None and not (float(window[0])<=controller_time<=float(window[1])): return None
    def compute():
        estimate=request.estimator.estimate(controller_time,controller_state)
        computed=np.asarray(controller.value(controller_time,np.asarray(estimate,dtype=float)),dtype=float)
        if records is not None: records.append((float(t),computed.reshape(-1).copy()))
        return computed
    signal=clock.value(t,compute) if clock is not None else compute()
    return signal

def _apply_control(base: np.ndarray,request: SolveRequest,t: float,y: np.ndarray,*,dimensionless: bool=False,records: list[tuple[float,np.ndarray]]|None=None,clock: _ControlClock|None=None):
    """Apply the single explicit controller, or retain the legacy RHS hook."""
    if request.controller is not None:
        signal=_controller_signal(request,t,y,records,clock,dimensionless=dimensionless)
        if signal is None: return base
        if getattr(request.controller,"actuation","additive_state")=="system_input":
            name="dimensionless_rhs_with_input" if dimensionless else "rhs_with_input"
            method=getattr(request.system,name,None)
            if method is not None: return np.asarray(method(t,y,signal),dtype=float)
            # Generic control-affine fallback for a plug-in that exposes B or
            # B(x), but has not implemented a custom input channel.
            if not dimensionless and hasattr(request.system,"input_matrix"):
                matrix=request.system.input_matrix(y)
                if matrix is not None: return base+np.asarray(matrix,dtype=float)@signal.reshape(-1)
            raise TypeError(f"{type(request.system).__name__} does not implement {name} or a usable control matrix; use additive_state actuation or add a system input channel.")
        extra=signal
    elif request.control is not None:
        control=request.control
        extra=control.value(t,y) if hasattr(control,"value") else control(t,y)
    else:
        return base
    extra=np.asarray(extra,dtype=float)
    if extra.shape!=base.shape: raise ValueError("Additive control must return one value per state variable.")
    return base+extra

def _controlled_rhs(request: SolveRequest,records: list[tuple[float,np.ndarray]]|None=None,clock: _ControlClock|None=None):
    system=request.system
    def rhs(t,y):
        if request.cancel_token is not None and request.cancel_token.is_set(): raise RuntimeError("solve cancelled")
        base=np.asarray(system.rhs(t,y),dtype=float)
        return _apply_control(base,request,t,y,records=records,clock=clock)
    return rhs

def _fixed_step(rhs,t0,tf,y0,h,method,session):
    times=[t0]; values=[np.asarray(y0,dtype=float)]; t=t0; y=np.asarray(y0,dtype=float)
    while t < tf:
        dt=min(h,tf-t)
        if method=="Euler": y=y+dt*rhs(t,y)
        else:
            k1=rhs(t,y); k2=rhs(t+dt/2,y+dt*k1/2); k3=rhs(t+dt/2,y+dt*k2/2); k4=rhs(t+dt,y+dt*k3)
            y=y+dt*(k1+2*k2+2*k3+k4)/6
        t+=dt; times.append(t); values.append(y.copy())
        if session.request.history.track_during_solve: session.checkpoint(np.asarray((t,)),np.asarray((y,)))
    return np.asarray(times),np.asarray(values)

def _control_history(time: np.ndarray, state_values: np.ndarray, records: list[tuple[float,np.ndarray]], request: SolveRequest, *, dimensionless: bool, policy: str) -> TimeSeries | None:
    """Build the accepted-solution control series.

    Sampled control is reconstructed from the held-command call log. Continuous
    control is evaluated on accepted solver states; the separate call log
    retains every internal RHS evaluation, including rejected trial states.
    """
    if request.controller is None: return None
    output_time=np.asarray(time,dtype=float)
    width=(1 if getattr(request.controller,"actuation","additive_state")=="system_input" else int(np.asarray(request.system.initial_state).size))
    control_values=np.zeros((output_time.size,width),dtype=float)
    if policy=="continuous_rhs":
        sampled=[]
        for item_time,item_state in zip(output_time,np.asarray(state_values,dtype=float)):
            signal=_controller_signal(request,float(item_time),item_state,dimensionless=dimensionless)
            sampled.append(np.zeros(width,dtype=float) if signal is None else np.asarray(signal,dtype=float).reshape(-1))
        if sampled: control_values=np.vstack(sampled)
    elif records:
        record_time=np.asarray([item[0] for item in records],dtype=float); record_values=np.vstack([item[1] for item in records])
        width=record_values.shape[1]; control_values=np.zeros((output_time.size,width),dtype=float)
        unique,indices=np.unique(record_time,return_index=True); samples=record_values[indices]
        for index in range(width): control_values[:,index]=np.interp(output_time,unique,samples[:,index],left=0.,right=float(samples[-1,index]))
        window=getattr(request.controller,"active_time_window",None)
        if window is not None:
            inactive=(output_time<float(window[0])) | (output_time>float(window[1]))
            control_values[inactive]=0.
    return TimeSeries(output_time,control_values)

def _control_target_history(time: np.ndarray, request: SolveRequest, state_dimension: int) -> TimeSeries | None:
    """Sample the configured reference on the solution grid for durable display.

    This intentionally archives values, not executable target/controller code,
    so Viewer can remain read-only and deterministic.
    """
    controller=request.controller
    target=None if controller is None else (getattr(controller,"display_target",None) or getattr(controller,"target",None))
    if target is None: return None
    integration_time=np.asarray(time,dtype=float)
    values=np.vstack([np.asarray(target.evaluate(float(item),state_dimension),dtype=float).reshape(-1) for item in integration_time])
    if values.shape!=(integration_time.size,state_dimension):
        raise ValueError("Controller target history must return one value per controlled state.")
    if not np.all(np.isfinite(values)):
        raise ValueError("Controller target history contains non-finite values.")
    return TimeSeries(integration_time,values)

def _low_level_adaptive(rhs, t0, tf, y0, request: SolveRequest, session: Any, max_step: float | None):
    """Adaptive stepping without reinitialization, for optional histories.

    It is intentionally opt-in: ordinary ``solve_ivp`` remains faster/simpler
    when only final output is needed.  Samples/checkpoints are emitted only at
    coarse requested intervals, never after every RHS evaluation.
    """
    classes={"RK23":RK23,"RK45":RK45,"DOP853":DOP853,"Radau":Radau,"BDF":BDF,"LSODA":LSODA}
    method=request.options.method
    if method not in classes: raise ValueError(f"low_level_stepper does not support {method!r}.")
    solver=classes[method](rhs,t0,np.asarray(y0,dtype=float),tf,rtol=request.options.rtol,atol=request.options.atol,max_step=max_step or np.inf)
    time=[float(t0)]; values=[np.asarray(y0,dtype=float).copy()]
    history=request.history; sample_every=history.sample_interval; checkpoint_every=history.checkpoint_interval
    next_sample=None if sample_every is None else t0+sample_every; next_checkpoint=None if checkpoint_every is None else t0+checkpoint_every
    checkpoint_path=None
    while solver.status=="running":
        solver.step()
        if solver.status=="failed": raise RuntimeError("low-level adaptive solver failed")
        time.append(float(solver.t)); values.append(np.asarray(solver.y,dtype=float).copy())
        if history.mode=="sampled" and next_sample is not None and solver.t>=next_sample:
            checkpoint_path=session.checkpoint(np.asarray((solver.t,)),np.asarray((solver.y,)),note="sampled-history")
            next_sample+=sample_every
        if history.mode=="checkpoint" and next_checkpoint is not None and solver.t>=next_checkpoint:
            checkpoint_path=session.checkpoint(np.asarray((solver.t,)),np.asarray((solver.y,)),note="durable-checkpoint")
            next_checkpoint+=checkpoint_every
            if history.pause_at_checkpoint and request.pause_token is not None and request.pause_token.is_set():
                return np.asarray(time),np.asarray(values),True,checkpoint_path
    return np.asarray(time),np.asarray(values),False,checkpoint_path

class _TelemetryHeartbeat:
    """Wall-clock-throttled scalar progress events for long ODE solves."""
    def __init__(self, request: SolveRequest, session: Any, t0: float, tf: float):
        self.request,self.session,self.t0,self.tf=request,session,float(t0),float(tf)
        self.started=self.last=monotonic()
        self.enabled=bool(request.history.telemetry_enabled and request.history.telemetry_heartbeat)
        self.interval=max(0.1,float(request.history.telemetry_interval_seconds))
    def report(self, t: float, *, accepted_steps: int | None=None) -> None:
        if not self.enabled: return
        now=monotonic()
        if now-self.last<self.interval: return
        span=self.tf-self.t0
        fraction=1. if span==0 else min(1.,max(0.,(float(t)-self.t0)/span))
        details={"elapsed_wall_seconds":round(now-self.started,6),"simulation_time":float(t),"progress_fraction":fraction}
        if accepted_steps is not None: details["accepted_steps"]=int(accepted_steps)
        self.session.telemetry("heartbeat",**details)
        if self.request.on_event is not None:
            self.request.on_event({"event":"heartbeat",**details})
        self.last=now

def _resolved_initial_state(request: SolveRequest) -> tuple[np.ndarray,str]:
    """Resolve free and controlled initial-state policies at the solver gate."""
    explicit=None if request.initial_state is None else np.asarray(request.initial_state,dtype=float)
    if request.controller is not None and getattr(request.system,"name",None)=="rayleigh-plesset":
        radius=float(request.system.equilibrium_radius_m)
        equilibrium=np.asarray((radius,0.),dtype=float)
        if explicit is not None: return explicit,"controlled-rpe-request-override"
        return equilibrium,"controlled-rpe-equilibrium-default"
    if explicit is not None: return explicit,"request-override"
    return np.asarray(request.system.initial_state,dtype=float),"system-default"

def solve_ode(request: SolveRequest, session: Any) -> NumericalResult:
    if request.time_span is None: raise ValueError("ODE solve requires time_span.")
    t0,tf=request.time_span
    if request.control is not None and request.controller is not None:
        raise ValueError("Provide either legacy control or controller, not both.")
    physical_initial,initial_state_policy=_resolved_initial_state(request)
    requested_form=request.options.use_dimensionless_form
    has_dimensionless=(hasattr(request.system,"dimensionless_rhs") and hasattr(request.system,"dimensionless_initial_state"))
    dimensionless=(has_dimensionless and bool(getattr(request.system,"prefer_dimensionless_form",True))) if requested_form is None else requested_form
    controller_coordinates=_bind_controller_coordinates(request,dimensionless=dimensionless)
    control_policy=None if request.controller is None else _control_evaluation_policy(request)
    window=None if request.controller is None else getattr(request.controller,"active_time_window",None)
    if window is not None and (len(window)!=2 or float(window[1])<float(window[0])):
        raise ValueError("controller.active_time_window must be a (start, end) pair with start <= end.")
    records: list[tuple[float,np.ndarray]]=[]
    clock=None
    if request.on_event is not None: request.on_event({"event":"started","kind":"ode","time_span":request.time_span})
    if dimensionless:
        if request.control is not None: raise ValueError("Legacy external control is dimensional; use a controller or solve dimensionally.")
        required=("dimensionless_rhs","dimensionless_initial_state","to_dimensionless_time","to_dimensional_time","to_dimensional_state")
        if not all(hasattr(request.system,name) for name in required): raise TypeError(f"{type(request.system).__name__} does not provide a dimensionless ODE form.")
        t0,tf=request.system.to_dimensionless_time(t0),request.system.to_dimensionless_time(tf)
        time_scale=request.system.dimensionless_parameters.time_scale_s
        sample_time=request.options.control_sample_time
        clock=_ControlClock(sample_time,t0) if control_policy=="sampled" else None
        def rhs(t,y):
            base=np.asarray(request.system.dimensionless_rhs(t,y),dtype=float)
            if request.cancel_token is not None and request.cancel_token.is_set(): raise RuntimeError("solve cancelled")
            return _apply_control(base,request,t,y,dimensionless=True,records=records,clock=clock)
        y0=np.asarray(request.system.to_dimensionless_state(physical_initial),dtype=float)
    else:
        time_scale=1.
        clock=_ControlClock(request.options.control_sample_time,t0) if control_policy=="sampled" else None
        rhs=_controlled_rhs(request,records,clock); y0=physical_initial
    method=request.options.method
    if method in {"Euler","RK4"}:
        native_eligible,native_reason=native_fixed_eligibility(request,dimensionless=dimensionless)
    else:
        native_eligible,native_reason=adaptive_native_eligibility(request)
    backend=select_ode_backend(request.options.execution_backend,method,native_eligible=native_eligible,ineligibility_reason=native_reason)
    max_step=None if request.options.max_step is None else request.options.max_step/time_scale
    step_size=None if request.options.step_size is None else request.options.step_size/time_scale
    heartbeat=_TelemetryHeartbeat(request,session,t0,tf)
    raw_rhs=rhs
    def rhs(t,y):
        heartbeat.report(float(t))
        return raw_rhs(t,y)
    paused=False; checkpoint_path=None
    try:
        if method in {"Euler","RK4"}:
            if not step_size: raise ValueError(f"{method} requires options.step_size.")
            if backend.compiled:
                time,values=integrate_fixed_native(
                    request.system,t0,tf,y0,step_size,method,
                    coordinate_system="dimensionless" if dimensionless else "dimensional",
                )
            else:
                time,values=_fixed_step(rhs,t0,tf,y0,step_size,method,session)
        elif backend.name=="sundials-cvode":
            time,values,adaptive_metadata=integrate_adaptive_sundials(rhs,t0,tf,y0,method=method,rtol=request.options.rtol,atol=request.options.atol,max_step=max_step)
        elif request.history.mode in {"sampled","checkpoint"} and (request.history.checkpoint_strategy=="low_level_stepper" or request.options.adaptive_engine=="low_level_stepper"):
            time,values,paused,checkpoint_path=_low_level_adaptive(rhs,t0,tf,y0,request,session,max_step)
        elif request.history.track_during_solve:
            interval=request.history.checkpoint_interval or (tf-t0)/20
            time_parts=[np.asarray((t0,))]; value_parts=[y0[None,:]]; current_t,current_y=t0,y0
            while current_t < tf:
                end=min(tf,current_t+interval)
                sol=solve_ivp(rhs,(current_t,end),current_y,method=method,rtol=request.options.rtol,atol=request.options.atol,
                              max_step=max_step or np.inf)
                if not sol.success: raise RuntimeError(sol.message)
                chunk_t=sol.t; chunk_y=sol.y.T; checkpoint_path=session.checkpoint(chunk_t,chunk_y)
                # ``solve_ivp`` repeats the segment's initial state; retain
                # it only once so checkpoint-paused series stay strictly
                # increasing and directly appendable on resume.
                time_parts.append(chunk_t[1:]); value_parts.append(chunk_y[1:]); current_t,current_y=end,sol.y[:,-1]
                if request.history.pause_at_checkpoint and request.pause_token is not None and request.pause_token.is_set():
                    paused=True; break
            time=np.concatenate(time_parts); values=np.vstack(value_parts)
        else:
            sol=solve_ivp(rhs,(t0,tf),y0,method=method,rtol=request.options.rtol,atol=request.options.atol,max_step=max_step or np.inf)
            if not sol.success: raise RuntimeError(sol.message)
            time,values=sol.t,sol.y.T
        control_series=_control_history(time,values,records,request,dimensionless=dimensionless,policy=control_policy or "continuous_rhs")
        target_series=_control_target_history(time,request,int(y0.size))
        control_call_time=np.asarray([item[0] for item in records],dtype=float) if records else None
        control_call_values=np.vstack([item[1] for item in records]) if records else None
        if dimensionless:
            time=np.asarray(request.system.to_dimensional_time(time),dtype=float)
            values=np.asarray([request.system.to_dimensional_state(row) for row in values],dtype=float)
            if control_series is not None: control_series=TimeSeries(np.asarray(request.system.to_dimensional_time(control_series.time),dtype=float),control_series.values)
            if control_call_time is not None: control_call_time=np.asarray(request.system.to_dimensional_time(control_call_time),dtype=float)
            if target_series is not None:
                target_time=np.asarray(request.system.to_dimensional_time(target_series.time),dtype=float)
                target_values=np.asarray([request.system.to_dimensional_state(row) for row in target_series.values],dtype=float)
                target_series=TimeSeries(target_time,target_values)
        status="paused" if paused else "complete"; message="paused at durable checkpoint" if paused else "solve completed"
        control_units,control_scale=_control_input_metadata(request,controller_coordinates)
        result=NumericalResult(TimeSeries(time,values),status,message,control_time_series=control_series,control_target_time_series=target_series,control_call_time=control_call_time,control_call_values=control_call_values,metadata={"state_dimension":int(y0.size),"initial_state_policy":initial_state_policy,"integration_form":"dimensionless" if dimensionless else "dimensional","nondimensionalization":serializable(request.system.dimensionless_parameters) if dimensionless else None,"controller_coordinate_system":controller_coordinates,"control_input_coordinate_system":controller_coordinates,"control_input_units":control_units,"control_input_scale_to_dimensional":control_scale,"control_sample_time":request.options.control_sample_time,"control_evaluation":control_policy,"controller_present":request.controller is not None,"control_target_present":target_series is not None,"controller_call_count":len(records),"control_history_sample_count":None if control_series is None else control_series.sample_count,"checkpoint_pause":paused,"checkpoint_path":None if checkpoint_path is None else str(checkpoint_path),"adaptive_backend":locals().get("adaptive_metadata"),"execution_backend":backend.provenance(requested=request.options.execution_backend)})
        if request.on_event is not None: request.on_event({"event":"paused" if paused else "completed","samples":result.time_series.sample_count,"checkpoint_path":None if checkpoint_path is None else str(checkpoint_path)})
        return result
    except Exception as exc:
        if request.on_event is not None: request.on_event({"event":"failed","message":str(exc)})
        return NumericalResult(None,"failed",str(exc),metadata={"state_dimension":int(y0.size),"execution_backend":backend.provenance(requested=request.options.execution_backend)})
