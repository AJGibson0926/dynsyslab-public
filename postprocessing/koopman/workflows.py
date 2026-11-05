"""Explicit standard operating procedures for Koopman training workflows."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from dynsyslib import Observable, ZeroDriver, load_system
from provenance import serializable
from time_series import TimeSeries

from .dmd import hankel_dmd_ensemble
from .koopman_sindy import koopman_with_sindy
from postprocessing.core import Trajectory
from rom import ControlledDataset, ControlledTrajectory


DUFFING_REFERENCE_IC_PATH=(Path(__file__).resolve().parents[2]/"data"/"duffing_korda_mezic_initial_conditions.csv")
DUFFING_REFERENCE_IC_NUMERIC_SHA256="decabdb25e40936dec8ab5281c4ea299d474ef0483b9b4fcfe66c0f974713cc7"


@dataclass(frozen=True)
class DuffingKoopmanSOP:
    """Reproducible free-data H-DMD-to-SINDy Duffing pathway."""

    random_seed: int=2141442
    trajectory_count: int=30
    free_final_time: float=8.0
    free_delta_t: float=.01
    free_solver: str="BDF"
    free_rtol: float=1e-10
    free_atol: float=1e-12
    hdmd_scalar_state_index: int=0
    hdmd_matrix_layout: str="window-snapshots"
    hdmd_mode_sort: str="terminal-amplitude-energy"
    hdmd_rank: int=20
    # The historical global default is 2,000, but an 8-second, 0.01-sampled
    # Duffing member has only 801 samples.  This per-experiment geometry is
    # explicit rather than silently reducing the global default.
    hdmd_matrix_size: int=400
    hdmd_sampling_window: int=400
    sindy_coordinate_observable: str="full-state-identity"
    sindy_polynomial_degree: int=4
    sindy_threshold: float=.05
    sindy_normalize: bool=True
    sindy_null_space_backend: str="gram-eigh"

    def metadata(self) -> dict[str,Any]:
        return {"encoding":"duffing-hdmd-sindy-sop/v1",**asdict(self),
            "observable_policy":{"hdmd":"scalar first state x","sindy":"full dynamical state (x, xdot)"},
            "hdmd_geometry_note":"2,000-delay global default is infeasible for an 801-sample member; this SOP declares 400 rather than silently clipping it."}


@dataclass(frozen=True)
class DuffingControlledTrainingSOP:
    """Korda--Mezic controlled-data geometry in canonical coordinates.

    The source command enters its unusual equation with a factor of ``0.5``.
    DynSysLab keeps the canonical equation and records that conversion instead
    of adopting the source's rescaled state-equation notation.
    """

    source_seed: int=2141442
    trajectory_count: int=100
    transition_count: int=200
    delta_t: float=.01
    source_input_low: float=-1.
    source_input_high: float=1.
    source_to_dynsyslab_input_scale: float=.5
    state_rejection_magnitude: float=3.

    def metadata(self) -> dict[str,Any]:
        return {"encoding":"duffing-korda-mezic-controlled-data/v1",**asdict(self),
            "integration":"fixed classical RK4; input held over each interval",
            "canonical_equation":"xdot=v; vdot=-0.5*v+x-4*x^3+u_dynsyslab",
            "source_input_conversion":"u_dynsyslab=0.5*u_source",
            "rng_warning":"Exact source controls must be supplied. MATLAB's Ubig depends on the RNG state after earlier script sections, including an optional rng(31) reset."}


@dataclass(frozen=True)
class RPEControlledTrainingSOP:
    """Designed persistent-excitation data for RPE control-map regression.

    Commands are exact interval-held values in ``p_control/p0``.  The input is
    deliberately zero mean, contains several hold times and a deterministic
    multisine component, and is independent of the stratified initial states.
    This avoids the unavailable adaptive-call command tape in the legacy RPE
    spreadsheet while exciting both the target tube and nearby free dynamics.
    """

    random_seed: int=742391
    trajectory_count: int=24
    transition_count: int=6000
    delta_tau: float=.001
    initial_radius_bounds: tuple[float,float]=(.78,1.22)
    initial_velocity_bounds: tuple[float,float]=(-.8,.8)
    input_amplitudes: tuple[float,...]=(.45,.60,.75,.90)
    hold_steps: tuple[int,...]=(10,20,40,80)
    multisine_fraction: float=.28
    multisine_angular_frequencies: tuple[float,...]=(1.5,3.,4.8,7.)
    radius_rejection_bounds: tuple[float,float]=(.25,2.)
    velocity_rejection_magnitude: float=20.

    def __post_init__(self) -> None:
        if self.trajectory_count<2 or self.transition_count<2 or self.delta_tau<=0.:
            raise ValueError("RPE controlled training needs at least two trajectories and transitions and positive delta_tau.")
        if not self.input_amplitudes or any(value<=0. for value in self.input_amplitudes):
            raise ValueError("RPE controlled-training input amplitudes must be positive.")
        if not self.hold_steps or any(int(value)!=value or value<1 for value in self.hold_steps):
            raise ValueError("RPE controlled-training hold steps must be positive integers.")
        if not 0.<=self.multisine_fraction<=1.:
            raise ValueError("multisine_fraction must lie in [0, 1].")

    def metadata(self) -> dict[str,Any]:
        return {"encoding":"rpe-designed-persistent-excitation/v1",**asdict(self),
            "integration":"fixed classical RK4; exact dimensionless pressure held over every output interval",
            "state_coordinates":"dimensionless (r, dr/dtau)","input_coordinates":"p_control/p0",
            "excitation":"balanced multilevel random telegraph plus random-phase multisine",
            "identification_purpose":"constant and bilinear B_psi regression"}


@dataclass(frozen=True)
class RPESymmetricPulseSOP:
    """Paired local experiments for data-driven pointwise ``B_psi`` samples.

    Every anchor launches otherwise identical ``+amplitude`` and
    ``-amplitude`` trajectories.  Their first-step lifted difference cancels
    autonomous evolution to leading order and isolates the input direction
    without consulting the analytic ``B_psi`` implementation.
    """

    random_seed: int=742391
    target_anchor_count: int=400
    rectangle_anchor_count: int=200
    delta_tau: float=.001
    transition_count: int=2
    pulse_magnitude: float=1.
    target_name: str="rpe_multisine_equilibrium_moderate_identification"
    target_offset: float=1.
    target_amplitudes: tuple[float,float]=(.09,.12)
    target_angular_frequencies: tuple[float,float]=(4.,3.)
    target_phases: tuple[float,float]=(0.,float(np.pi))
    target_tau_span: tuple[float,float]=(0.,20.)
    tube_radius_standard_deviation: float=.025
    tube_velocity_standard_deviation: float=.08
    rectangle_radius_bounds: tuple[float,float]=(.76,1.24)
    rectangle_velocity_bounds: tuple[float,float]=(-.9,.9)
    anchor_radius_clip: tuple[float,float]=(.72,1.28)

    def __post_init__(self) -> None:
        if self.target_anchor_count<1 or self.rectangle_anchor_count<1:
            raise ValueError("Symmetric-pulse training needs target-tube and rectangle anchors.")
        if self.delta_tau<=0. or self.transition_count<2 or self.pulse_magnitude<=0.:
            raise ValueError("Symmetric-pulse time, transition count, and magnitude must be positive.")

    def metadata(self) -> dict[str,Any]:
        return {"encoding":"rpe-symmetric-pulse-identification/v1",**asdict(self),
            "integration":"fixed classical RK4; pulse held over each interval",
            "state_coordinates":"dimensionless (r, dr/dtau)","input_coordinates":"p_control/p0",
            "pairing":"same anchor and opposite commands; pair_id and polarity retained",
            "estimator":"(Psi(x_plus(dt))-Psi(x_minus(dt)))/((u_plus-u_minus)*dt)",
            "analytic_Bpsi_used_to_generate_data":False}


def duffing_reference_initial_conditions(*,count: int|None=None) -> np.ndarray:
    """Load the exact retained MATLAB/Korda unit-circle IC table as rows."""
    values=np.loadtxt(DUFFING_REFERENCE_IC_PATH,delimiter=",").T
    digest=sha256(np.ascontiguousarray(values.T,dtype="<f8").tobytes()).hexdigest()
    if values.shape!=(100,2) or digest!=DUFFING_REFERENCE_IC_NUMERIC_SHA256:
        raise RuntimeError("Retained Duffing initial-condition table failed its numeric checksum.")
    if not np.allclose(np.linalg.norm(values,axis=1),1.,rtol=0.,atol=1e-13):
        raise RuntimeError("Retained Duffing initial conditions are no longer on the unit circle.")
    if count is None: return values
    if int(count)!=count or not 1<=count<=values.shape[0]:
        raise ValueError("count must select between 1 and 100 retained Duffing ICs.")
    return values[:int(count)].copy()


def generate_duffing_free_training(sop: DuffingKoopmanSOP|None=None):
    """Integrate the canonical Duffing equations on the exact retained ICs."""
    policy=sop or DuffingKoopmanSOP(); system=load_system("duffing",driver=ZeroDriver(),
        delta=.5,alpha=1.,beta=4.,input_gain=1.)
    initial=duffing_reference_initial_conditions(count=policy.trajectory_count)
    time=np.arange(0.,policy.free_final_time+.5*policy.free_delta_t,policy.free_delta_t)
    members=[]; started=perf_counter(); evaluations=0
    for state in initial:
        solution=solve_ivp(system.rhs,(time[0],time[-1]),state,method=policy.free_solver,
            t_eval=time,rtol=policy.free_rtol,atol=policy.free_atol)
        if not solution.success: raise RuntimeError(f"Duffing free-training integration failed: {solution.message}")
        members.append(solution.y.T); evaluations+=int(solution.nfev)
    elapsed=perf_counter()-started
    series=TimeSeries(time,np.stack(members,axis=1),pooled_members=len(members))
    metadata={"uniform_time_steps":True,"system":system.name,"system_parameters":serializable(system),
        "integration_form":"dimensionless","postprocessing_preprocessing":{"coordinate_system":"dimensionless"},
        "training":{"kind":"free","initial_condition_source":"retained exact Korda/cleanup table",
        "initial_condition_numeric_sha256":DUFFING_REFERENCE_IC_NUMERIC_SHA256,
        "trajectory_count":len(members),"solver":policy.free_solver,"rtol":policy.free_rtol,"atol":policy.free_atol,
        "elapsed_seconds":elapsed,"rhs_evaluations":evaluations},"sop":policy.metadata()}
    return series,system,metadata


def sample_duffing_controlled_inputs(
    sop: DuffingControlledTrainingSOP|None=None,
    *,
    seed: int|None=None,
) -> np.ndarray:
    """Sample a fresh source-convention control cube with explicit provenance.

    This is useful for new DynSysLab experiments, but it is intentionally not
    advertised as reproducing MATLAB's random stream.  Exact source parity uses
    the ``Ubig`` array exported from the source run and passes it to
    :func:`generate_duffing_controlled_training`.
    """
    policy=sop or DuffingControlledTrainingSOP()
    generator=np.random.default_rng(policy.source_seed if seed is None else int(seed))
    return generator.uniform(policy.source_input_low,policy.source_input_high,
        size=(policy.trajectory_count,policy.transition_count,1))


def generate_duffing_controlled_training(
    source_controls: np.ndarray,
    *,
    initial_conditions: np.ndarray|None=None,
    sop: DuffingControlledTrainingSOP|None=None,
) -> tuple[ControlledDataset,Any]:
    """Reproduce the controlled fixed-RK4 experiment from supplied controls.

    ``source_controls`` is normally shaped ``(trajectory, transition, input)``.
    For direct use of MATLAB's ``Ubig`` it may instead be the two-dimensional
    source layout ``(transition, trajectory)``.  Stored trajectory inputs are
    canonical DynSysLab commands, so constant/bilinear regressions are aligned
    with the KROM system's ``input_gain=1`` contract.
    """
    policy=sop or DuffingControlledTrainingSOP()
    raw=np.asarray(source_controls,dtype=float)
    if raw.ndim==2 and raw.shape==(policy.transition_count,policy.trajectory_count):
        raw=raw.T[:,:,None]
    elif raw.ndim==2 and policy.trajectory_count==1 and raw.shape==(policy.transition_count,1):
        raw=raw.T[:,:,None]
    if raw.shape!=(policy.trajectory_count,policy.transition_count,1):
        raise ValueError("source_controls must have shape (trajectory_count, transition_count, 1) or MATLAB Ubig shape (transition_count, trajectory_count).")
    if np.any(~np.isfinite(raw)):
        raise ValueError("source_controls must be finite.")
    initial=(duffing_reference_initial_conditions(count=policy.trajectory_count)
        if initial_conditions is None else np.asarray(initial_conditions,dtype=float))
    if initial.shape!=(policy.trajectory_count,2) or np.any(~np.isfinite(initial)):
        raise ValueError("initial_conditions must have shape (trajectory_count, 2) and be finite.")

    system=load_system("duffing",driver=ZeroDriver(),delta=.5,alpha=1.,beta=4.,input_gain=1.)
    effective=policy.source_to_dynsyslab_input_scale*raw
    time=policy.delta_t*np.arange(policy.transition_count+1,dtype=float)
    trajectories=[]; accepted=[]

    def rhs(state: np.ndarray, command: float) -> np.ndarray:
        x,v=state
        return np.asarray((v,-.5*v+x-4.*x**3+command),dtype=float)

    for member in range(policy.trajectory_count):
        states=np.empty((policy.transition_count+1,2),dtype=float); states[0]=initial[member]
        for step in range(policy.transition_count):
            command=float(effective[member,step,0]); state=states[step]; dt=policy.delta_t
            k1=rhs(state,command)
            k2=rhs(state+.5*dt*k1,command)
            k3=rhs(state+.5*dt*k2,command)
            k4=rhs(state+dt*k3,command)
            states[step+1]=state+(dt/6.)*(k1+2.*k2+2.*k3+k4)
        valid=bool(np.all(np.isfinite(states)) and np.all(np.abs(states)<=policy.state_rejection_magnitude))
        accepted.append(valid)
        if valid:
            trajectories.append(ControlledTrajectory(time,states,effective[member],{
                "source_member_index":member,"input_coordinates":"canonical DynSysLab",
                "source_to_dynsyslab_input_scale":policy.source_to_dynsyslab_input_scale}))
    if not trajectories:
        raise RuntimeError("Every controlled Duffing trajectory failed the retained state-quality gate.")
    raw_digest=sha256(np.ascontiguousarray(raw,dtype="<f8").tobytes()).hexdigest()
    dataset=ControlledDataset(trajectories,{**policy.metadata(),
        "source_control_numeric_sha256":raw_digest,"candidate_trajectory_count":policy.trajectory_count,
        "accepted_trajectory_count":len(trajectories),"accepted_member_mask":accepted,
        "stored_input_coordinates":"canonical DynSysLab; multiply source commands by 0.5"})
    return dataset,system


def _rk4_control_step(system,state: np.ndarray,command: float,time: float,step: float) -> np.ndarray:
    """Advance one dimensionless, interval-held control step."""
    rhs=lambda value: system.dimensionless_rhs_with_input(time,value,command)
    k1=rhs(state)
    k2=rhs(state+.5*step*k1)
    k3=rhs(state+.5*step*k2)
    k4=rhs(state+step*k3)
    return state+(step/6.)*(k1+2.*k2+2.*k3+k4)


def rpe_controlled_initial_conditions(
    sop: RPEControlledTrainingSOP|None=None,
    *,
    seed: int|None=None,
) -> np.ndarray:
    """Return stratified, independently permuted RPE identification states."""
    policy=sop or RPEControlledTrainingSOP()
    generator=np.random.default_rng(policy.random_seed if seed is None else int(seed))
    count=policy.trajectory_count
    radius_order=generator.permutation(count)
    velocity_order=generator.permutation(count)
    radius_fraction=(radius_order+generator.uniform(.1,.9,count))/count
    velocity_fraction=(velocity_order+generator.uniform(.1,.9,count))/count
    r0,r1=policy.initial_radius_bounds; v0,v1=policy.initial_velocity_bounds
    return np.column_stack((r0+(r1-r0)*radius_fraction,v0+(v1-v0)*velocity_fraction))


def sample_rpe_controlled_inputs(
    sop: RPEControlledTrainingSOP|None=None,
    *,
    seed: int|None=None,
) -> np.ndarray:
    """Generate exact, bounded and persistently exciting RPE command tapes."""
    policy=sop or RPEControlledTrainingSOP()
    generator=np.random.default_rng(policy.random_seed if seed is None else int(seed))
    controls=np.empty((policy.trajectory_count,policy.transition_count,1),dtype=float)
    sample_time=(np.arange(policy.transition_count,dtype=float)+.5)*policy.delta_tau
    base_levels=np.asarray((-1.,-.5,.5,1.),dtype=float)
    frequencies=np.asarray(policy.multisine_angular_frequencies,dtype=float)
    for member in range(policy.trajectory_count):
        amplitude=float(policy.input_amplitudes[member%len(policy.input_amplitudes)])
        hold=int(policy.hold_steps[(member//len(policy.input_amplitudes))%len(policy.hold_steps)])
        block_count=int(np.ceil(policy.transition_count/hold))
        levels=np.resize(base_levels,block_count)
        generator.shuffle(levels)
        telegraph=np.repeat(levels,hold)[:policy.transition_count]
        phases=generator.uniform(0.,2.*np.pi,frequencies.size)
        multisine=np.sum(np.sin(sample_time[:,None]*frequencies[None,:]+phases[None,:]),axis=1)
        multisine/=max(float(np.max(np.abs(multisine))),np.finfo(float).tiny)
        fraction=float(policy.multisine_fraction)
        command=amplitude*((1.-fraction)*telegraph+fraction*multisine)
        command-=np.mean(command)
        peak=float(np.max(np.abs(command)))
        if peak>amplitude: command*=amplitude/peak
        controls[member,:,0]=command
    return controls


def generate_rpe_controlled_training(
    controls: np.ndarray|None=None,
    *,
    initial_conditions: np.ndarray|None=None,
    sop: RPEControlledTrainingSOP|None=None,
) -> tuple[ControlledDataset,Any]:
    """Integrate exact-command RPE persistent-excitation trajectories."""
    policy=sop or RPEControlledTrainingSOP()
    commands=sample_rpe_controlled_inputs(policy) if controls is None else np.asarray(controls,dtype=float)
    if commands.ndim==2: commands=commands[:,:,None]
    expected=(policy.trajectory_count,policy.transition_count,1)
    if commands.shape!=expected or np.any(~np.isfinite(commands)):
        raise ValueError(f"RPE controls must be finite with shape {expected}.")
    initial=(rpe_controlled_initial_conditions(policy) if initial_conditions is None
        else np.asarray(initial_conditions,dtype=float))
    if initial.shape!=(policy.trajectory_count,2) or np.any(~np.isfinite(initial)) or np.any(initial[:,0]<=0.):
        raise ValueError("RPE identification initial conditions must be finite, positive-radius rows.")
    system=load_system("rpe",driver=ZeroDriver(),preset="air_water_20c")
    time=policy.delta_tau*np.arange(policy.transition_count+1,dtype=float)
    trajectories=[]; accepted=[]; rejection_reasons=[]
    for member in range(policy.trajectory_count):
        states=np.empty((policy.transition_count+1,2),dtype=float); states[0]=initial[member]
        reason=None
        for step in range(policy.transition_count):
            states[step+1]=_rk4_control_step(system,states[step],float(commands[member,step,0]),
                float(time[step]),policy.delta_tau)
            radius,velocity=states[step+1]
            if not np.all(np.isfinite(states[step+1])):
                reason="nonfinite-state"; break
            if not policy.radius_rejection_bounds[0]<=radius<=policy.radius_rejection_bounds[1]:
                reason="radius-quality-gate"; break
            if abs(velocity)>policy.velocity_rejection_magnitude:
                reason="velocity-quality-gate"; break
        valid=reason is None; accepted.append(valid); rejection_reasons.append(reason)
        if valid:
            amplitude=float(np.max(np.abs(commands[member,:,0])))
            trajectories.append(ControlledTrajectory(time,states,commands[member],{
                "source_member_index":member,"input_coordinates":"p_control/p0",
                "command_tape_exact":True,"command_peak":amplitude,
                "hold_steps":int(policy.hold_steps[(member//len(policy.input_amplitudes))%len(policy.hold_steps)])}))
    if not trajectories: raise RuntimeError("Every designed RPE controlled trajectory failed its state-quality gate.")
    digest=sha256(np.ascontiguousarray(commands,dtype="<f8").tobytes()).hexdigest()
    metadata={**policy.metadata(),"command_numeric_sha256":digest,"command_tape_exact":True,
        "candidate_trajectory_count":policy.trajectory_count,"accepted_trajectory_count":len(trajectories),
        "accepted_member_mask":accepted,"rejection_reasons":rejection_reasons,
        "initial_conditions":initial.tolist(),"physical_parameters":serializable(system.dimensional_parameters),
        "dimensionless_groups":serializable(system.dimensionless_parameters)}
    return ControlledDataset(trajectories,metadata),system


def rpe_symmetric_pulse_anchors(
    sop: RPESymmetricPulseSOP|None=None,
    *,
    seed: int|None=None,
) -> tuple[np.ndarray,np.ndarray]:
    """Sample target-tube and rectangular anchors with retained cohort labels."""
    policy=sop or RPESymmetricPulseSOP()
    generator=np.random.default_rng(policy.random_seed if seed is None else int(seed))
    tau=np.linspace(policy.target_tau_span[0],policy.target_tau_span[1],policy.target_anchor_count,endpoint=False)
    amplitude=np.asarray(policy.target_amplitudes); omega=np.asarray(policy.target_angular_frequencies)
    phase=np.asarray(policy.target_phases)
    angle=tau[:,None]*omega[None,:]+phase[None,:]
    target=np.column_stack((policy.target_offset+np.sum(amplitude[None,:]*np.sin(angle),axis=1),
        np.sum((amplitude*omega)[None,:]*np.cos(angle),axis=1)))
    target[:,0]+=generator.normal(0.,policy.tube_radius_standard_deviation,policy.target_anchor_count)
    target[:,1]+=generator.normal(0.,policy.tube_velocity_standard_deviation,policy.target_anchor_count)
    rectangle=np.column_stack((generator.uniform(*policy.rectangle_radius_bounds,policy.rectangle_anchor_count),
        generator.uniform(*policy.rectangle_velocity_bounds,policy.rectangle_anchor_count)))
    anchors=np.vstack((target,rectangle)); anchors[:,0]=np.clip(anchors[:,0],*policy.anchor_radius_clip)
    cohort=np.asarray(["target-tube"]*policy.target_anchor_count+["support-rectangle"]*policy.rectangle_anchor_count)
    return anchors,cohort


def generate_rpe_symmetric_pulse_training(
    *,
    anchors: np.ndarray|None=None,
    sop: RPESymmetricPulseSOP|None=None,
) -> tuple[ControlledDataset,Any]:
    """Generate paired opposite-command rollouts around prescribed anchors."""
    policy=sop or RPESymmetricPulseSOP()
    if anchors is None: initial,cohort=rpe_symmetric_pulse_anchors(policy)
    else:
        initial=np.asarray(anchors,dtype=float); cohort=np.asarray(["provided"]*initial.shape[0])
    if initial.ndim!=2 or initial.shape[1]!=2 or initial.shape[0]<2 or np.any(~np.isfinite(initial)) or np.any(initial[:,0]<=0.):
        raise ValueError("Symmetric-pulse anchors must be finite positive-radius state rows.")
    system=load_system("rpe",driver=ZeroDriver(),preset="air_water_20c")
    time=policy.delta_tau*np.arange(policy.transition_count+1,dtype=float)
    trajectories=[]
    for pair_id,(anchor,label) in enumerate(zip(initial,cohort)):
        for polarity in (-1,1):
            command=polarity*policy.pulse_magnitude
            inputs=np.full((policy.transition_count,1),command,dtype=float)
            states=np.empty((policy.transition_count+1,2),dtype=float); states[0]=anchor
            for step in range(policy.transition_count):
                states[step+1]=_rk4_control_step(system,states[step],command,float(time[step]),policy.delta_tau)
            if np.any(~np.isfinite(states)) or np.any(states[:,0]<=0.):
                raise RuntimeError(f"Symmetric-pulse pair {pair_id} produced an invalid RPE state.")
            trajectories.append(ControlledTrajectory(time,states,inputs,{"pair_id":pair_id,
                "polarity":polarity,"cohort":str(label),"command_tape_exact":True,
                "input_coordinates":"p_control/p0"}))
    metadata={**policy.metadata(),"anchor_count":int(initial.shape[0]),"trajectory_count":len(trajectories),
        "cohort_counts":{"target-tube":int(np.sum(cohort=="target-tube")),
            "support-rectangle":int(np.sum(cohort=="support-rectangle")),
            "provided":int(np.sum(cohort=="provided"))},"command_tape_exact":True,
        "physical_parameters":serializable(system.dimensional_parameters),
        "dimensionless_groups":serializable(system.dimensionless_parameters)}
    return ControlledDataset(trajectories,metadata),system


def duffing_hdmd_sindy(sop: DuffingKoopmanSOP|None=None):
    """Run the retained scalar-H-DMD/full-state-SINDy operating procedure."""
    from .duffing_pipeline import duffing_hdmd_sindy as run_pipeline
    return run_pipeline(sop)


__all__=["DUFFING_REFERENCE_IC_NUMERIC_SHA256","DUFFING_REFERENCE_IC_PATH",
    "DuffingControlledTrainingSOP","DuffingKoopmanSOP","RPEControlledTrainingSOP","RPESymmetricPulseSOP","duffing_hdmd_sindy",
    "duffing_reference_initial_conditions","generate_duffing_controlled_training",
    "generate_duffing_free_training","sample_duffing_controlled_inputs"]
__all__ += ["generate_rpe_controlled_training","generate_rpe_symmetric_pulse_training",
    "rpe_controlled_initial_conditions","rpe_symmetric_pulse_anchors","sample_rpe_controlled_inputs"]
