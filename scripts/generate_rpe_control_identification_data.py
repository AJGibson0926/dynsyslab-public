"""Generate and audit RPE controlled data designed for B_psi identification."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np


ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from postprocessing.core import integration_coordinates, load_solution
from postprocessing.koopman import (FreeTrajectoryDataset, RPEControlledTrainingSOP, RPESymmetricPulseSOP,
    generate_rpe_controlled_training, generate_rpe_symmetric_pulse_training)
from rom import (ControlledDataset, SymmetricPulseControlMapSamples,
    ConstantControlMatrix, estimate_symmetric_pulse_control_map,
    fit_bilinear_control_map, fit_bilinear_control_map_samples,
    fit_bilinear_control_matrix, fit_constant_control_map_samples,
    fit_conjugate_bilinear_control_map, fit_constant_control_matrix,
    load_krom_archive, save_controlled_dataset, save_krom_error_evaluations)


KROM_ARCHIVE=ROOT/"archives/derived/koopman-krom-extract/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-krom-extract__20260827T183646.410439Z__0dea0923.npz"
DEFAULT_OUTPUT=ROOT/"archives/experiments/rpe_control_identification_designed"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _subset_samples(samples: SymmetricPulseControlMapSamples,index: np.ndarray,role: str) -> SymmetricPulseControlMapSamples:
    return SymmetricPulseControlMapSamples(samples.states[index],samples.matrices[index],
        {**samples.diagnostics,"subset_role":role,"subset_sample_count":int(np.sum(index))})


def _target_states(count: int=4001) -> np.ndarray:
    tau=np.linspace(0.,10.,count); amplitude=np.asarray((.09,.12)); omega=np.asarray((4.,3.)); phase=np.asarray((0.,np.pi))
    angle=tau[:,None]*omega[None,:]+phase[None,:]
    return np.column_stack((1.+np.sum(amplitude[None,:]*np.sin(angle),axis=1),
        np.sum((amplitude*omega)[None,:]*np.cos(angle),axis=1)))


def _representation_error(krom,representation,states: np.ndarray) -> dict[str,float]:
    exact=np.asarray([krom.B_psi(state,output="complex",representation="exact") for state in states])
    predicted=np.asarray([representation.evaluate(krom.psi(state,output="complex",state_coordinates="model")) for state in states])
    difference=np.linalg.norm((predicted-exact).reshape(states.shape[0],-1),axis=1)
    magnitude=np.linalg.norm(exact.reshape(states.shape[0],-1),axis=1)
    relative=difference/np.maximum(magnitude,np.finfo(float).tiny)
    return {"rmse":float(np.sqrt(np.mean(difference*difference))),
        "relative_frobenius":float(np.linalg.norm(predicted-exact)/np.linalg.norm(exact)),
        "median_pointwise_relative":float(np.median(relative)),
        "p95_pointwise_relative":float(np.percentile(relative,95)),
        "maximum_pointwise_relative":float(np.max(relative))}


def _dataset_diagnostics(dataset: ControlledDataset,krom) -> dict[str,object]:
    states=np.vstack([item.states for item in dataset.trajectories])
    inputs=np.vstack([item.inputs for item in dataset.trajectories])
    # A one-in-ten audit sample is enough to diagnose persistent-excitation
    # rank and conditioning without turning a data-generation command into a
    # second full regression.
    sampled_states=np.vstack([item.states[:-1:10] for item in dataset.trajectories])
    sampled_inputs=np.vstack([item.inputs[::10] for item in dataset.trajectories])
    lifted=np.vstack([krom.psi(state,output="complex",state_coordinates="model") for state in sampled_states])
    bilinear=np.column_stack((sampled_inputs[:,0,None]*lifted,sampled_inputs[:,0]))
    return {"state_minimum":np.min(states,axis=0).tolist(),"state_maximum":np.max(states,axis=0).tolist(),
        "state_p01":np.percentile(states,1,axis=0).tolist(),"state_p99":np.percentile(states,99,axis=0).tolist(),
        "input_minimum":float(np.min(inputs)),"input_maximum":float(np.max(inputs)),
        "input_mean":float(np.mean(inputs)),"input_rms":float(np.sqrt(np.mean(inputs*inputs))),
        "constant_design_rank":int(np.linalg.matrix_rank(sampled_inputs)),
        "bilinear_design_rank":int(np.linalg.matrix_rank(bilinear)),
        "bilinear_design_shape":list(bilinear.shape),"bilinear_design_condition_number":float(np.linalg.cond(bilinear)),
        "audit_stride":10}


def _sample_estimation_error(krom,samples: SymmetricPulseControlMapSamples) -> dict[str,float]:
    exact=np.asarray([krom.B_psi(state,output="complex",representation="exact") for state in samples.states])
    residual=samples.matrices-exact
    return {"rmse":float(np.sqrt(np.mean(np.abs(residual)**2))),
        "relative_frobenius":float(np.linalg.norm(residual)/np.linalg.norm(exact))}


def generate(*,output: Path,quick: bool=False,skip_dynamics_regression: bool=False) -> Path:
    """Generate designed RPE excitation data for control-map identification."""
    output=output.expanduser().resolve(); output.mkdir(parents=True,exist_ok=True)
    persistent_policy=RPEControlledTrainingSOP()
    pulse_policy=RPESymmetricPulseSOP()
    if quick:
        persistent_policy=replace(persistent_policy,trajectory_count=8,transition_count=1000)
        pulse_policy=replace(pulse_policy,target_anchor_count=80,rectangle_anchor_count=40)
    krom=load_krom_archive(KROM_ARCHIVE)

    started=perf_counter(); persistent,_=generate_rpe_controlled_training(sop=persistent_policy)
    persistent_seconds=perf_counter()-started
    started=perf_counter(); pulses,_=generate_rpe_symmetric_pulse_training(sop=pulse_policy)
    pulse_seconds=perf_counter()-started
    persistent_archive,persistent_manifest=save_controlled_dataset(persistent,output/"rpe-persistent-excitation.npz")
    pulse_archive,pulse_manifest=save_controlled_dataset(pulses,output/"rpe-symmetric-pulses.npz")

    pulse_samples=estimate_symmetric_pulse_control_map(krom,pulses)
    validation=np.arange(pulse_samples.states.shape[0])%5==0
    train_samples=_subset_samples(pulse_samples,~validation,"training")
    validation_samples=_subset_samples(pulse_samples,validation,"held-out")
    constant=fit_constant_control_map_samples(krom,train_samples)
    bilinear=fit_bilinear_control_map_samples(krom,train_samples,ridge=1e-12)
    conjugate=fit_bilinear_control_map_samples(krom,train_samples,ridge=1e-12,conjugate_complete=True)
    target=_target_states()
    approximations={
        "pulse_constant":{"diagnostics":constant.diagnostics,
            "held_out_error":_representation_error(krom,constant.representation,validation_samples.states),
            "target_error":_representation_error(krom,constant.representation,target)},
        "pulse_bilinear":{"diagnostics":bilinear.diagnostics,
            "held_out_error":_representation_error(krom,bilinear.representation,validation_samples.states),
            "target_error":_representation_error(krom,bilinear.representation,target)},
        "pulse_conjugate_bilinear":{"diagnostics":conjugate.diagnostics,
            "held_out_error":_representation_error(krom,conjugate.representation,validation_samples.states),
            "target_error":_representation_error(krom,conjugate.representation,target)},
    }
    # These exact-map projections are audit-only representation floors.  They
    # are fitted after the controlled-data candidates and are never substituted
    # into those candidates or their archives.
    oracle_exact=np.asarray([krom.B_psi(state,output="complex",representation="exact")
        for state in train_samples.states])
    oracle_constant=ConstantControlMatrix(np.mean(oracle_exact,axis=0))
    oracle_bilinear=fit_bilinear_control_map(krom,train_samples.states,ridge=1e-12).representation
    oracle_conjugate=fit_conjugate_bilinear_control_map(krom,train_samples.states,ridge=1e-12).representation
    oracle_floor={name:{"held_out_error":_representation_error(krom,representation,validation_samples.states),
            "target_error":_representation_error(krom,representation,target)}
        for name,representation in (("constant",oracle_constant),("bilinear",oracle_bilinear),
            ("conjugate_bilinear",oracle_conjugate))}
    arrays={"pulse_sample_states":pulse_samples.states,"pulse_Bpsi_estimates":pulse_samples.matrices,
        "constant_matrix":constant.representation.matrix,"bilinear_operators":bilinear.representation.operators,
        "bilinear_offsets":bilinear.representation.offsets,
        "conjugate_bilinear_operators":conjugate.representation.operators,
        "conjugate_bilinear_conjugate_operators":conjugate.representation.conjugate_operators,
        "conjugate_bilinear_offsets":conjugate.representation.offsets}

    dynamics={}
    if not skip_dynamics_regression:
        train_members=[item for index,item in enumerate(persistent.trajectories) if index%4!=0]
        validation_members=[item for index,item in enumerate(persistent.trajectories) if index%4==0]
        train_dataset=ControlledDataset(train_members,{**persistent.metadata,"subset_role":"training"})
        started=perf_counter(); dynamics_constant=fit_constant_control_matrix(krom,train_dataset,
            method="continuous-residual",ridge=1e-10); constant_seconds=perf_counter()-started
        started=perf_counter(); dynamics_bilinear=fit_bilinear_control_matrix(krom,train_dataset,ridge=1e-10)
        bilinear_seconds=perf_counter()-started
        validation_states=np.vstack([item.states[::10] for item in validation_members])
        dynamics={"constant":{"seconds":constant_seconds,"diagnostics":dynamics_constant.diagnostics,
            "held_out_error":_representation_error(krom,dynamics_constant.representation,validation_states),
            "target_error":_representation_error(krom,dynamics_constant.representation,target)},
            "bilinear":{"seconds":bilinear_seconds,"diagnostics":dynamics_bilinear.diagnostics,
            "held_out_error":_representation_error(krom,dynamics_bilinear.representation,validation_states),
            "target_error":_representation_error(krom,dynamics_bilinear.representation,target)}}
        arrays["dynamics_constant_matrix"]=dynamics_constant.representation.matrix
        arrays["dynamics_bilinear_operators"]=dynamics_bilinear.representation.operators
        arrays["dynamics_bilinear_offsets"]=dynamics_bilinear.representation.offsets

    # Persist autonomous and controlled-map errors as first-class records. The
    # candidate maps coexist here; none becomes the KROM's default controller.
    free_source=integration_coordinates(load_solution(krom.system_reference.source_archive))
    free_dataset=FreeTrajectoryDataset((free_source.time_series,),label="eigenfunction free training source",
        role="training",coordinate_system="model",
        provenance={"source_archive":str(Path(krom.system_reference.source_archive).resolve())})
    krom.evaluate_eigenfunction_errors(free_dataset,horizons=(.001,.01,.1,.5,1.,5.,10.))
    bounds=((pulse_policy.rectangle_radius_bounds[0],pulse_policy.rectangle_velocity_bounds[0]),
        (pulse_policy.rectangle_radius_bounds[1],pulse_policy.rectangle_velocity_bounds[1]))
    constant.evaluate(krom,pulses,bounds,label="pulse-data constant",grid_shape=(61,61),
        metadata={"fit_subset":"four-of-five pulse anchors"})
    bilinear.evaluate(krom,pulses,bounds,label="pulse-data complex bilinear",grid_shape=(61,61),
        metadata={"fit_subset":"four-of-five pulse anchors"})
    conjugate.evaluate(krom,pulses,bounds,label="pulse-data conjugate-complete bilinear",grid_shape=(61,61),
        metadata={"fit_subset":"four-of-five pulse anchors"})
    if not skip_dynamics_regression:
        dynamics_constant.evaluate(krom,train_dataset,bounds,label="persistent-data constant",
            grid_shape=(61,61))
        dynamics_bilinear.evaluate(krom,train_dataset,bounds,label="persistent-data bilinear",
            grid_shape=(61,61))
    error_bundle=save_krom_error_evaluations(output/"rpe-minimal-krom-error-evaluations.npz",krom,
        metadata={"purpose":"compare controlled-data B_psi fits while retaining autonomous eigenfunction error"})

    candidate_path=output/"rpe-control-map-candidates.npz"; np.savez_compressed(candidate_path,**arrays)
    payload={"schema":"dynsyslab.rpe-designed-control-identification/v1","quick":quick,
        "krom_archive":str(KROM_ARCHIVE.resolve()),"krom_sha256":_sha256(KROM_ARCHIVE),
        "persistent_policy":asdict(persistent_policy),"pulse_policy":asdict(pulse_policy),
        "persistent_dataset":str(persistent_archive),"persistent_manifest":str(persistent_manifest),
        "persistent_sha256":_sha256(persistent_archive),"persistent_generation_seconds":persistent_seconds,
        "persistent_diagnostics":_dataset_diagnostics(persistent,krom),
        "pulse_dataset":str(pulse_archive),"pulse_manifest":str(pulse_manifest),
        "pulse_sha256":_sha256(pulse_archive),"pulse_generation_seconds":pulse_seconds,
        "pulse_estimator":pulse_samples.diagnostics,"pulse_estimator_exact_validation":_sample_estimation_error(krom,validation_samples),
        "target_contract":{"name":"rpe_multisine_equilibrium_moderate_identification",
            "radius":"1 + 0.09 sin(4 tau) - 0.12 sin(3 tau)",
            "velocity":"0.36 cos(4 tau) - 0.36 cos(3 tau)","horizon_tau":10.},
        "approximations":approximations,"dynamics_residual_approximations":dynamics,
        "oracle_representation_floor":oracle_floor,
        "candidate_arrays":str(candidate_path),"candidate_arrays_sha256":_sha256(candidate_path),
        "error_evaluation_bundle":str(error_bundle),"error_evaluation_bundle_sha256":_sha256(error_bundle),
        "interpretation":{"training":"No exact B_psi samples enter either controlled dataset or pulse-derived fits.",
            "evaluation":"The analytic exact B_psi is used only after fitting, to quantify held-out and target-path error and an oracle representation-capacity floor.",
            "decision_gate":"Use the pulse-derived data as the primary map-identification source; retain long rollouts for dynamics-regression comparison and future noise studies."}}
    summary=output/"summary.json"; summary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return summary


def main() -> int:
    """Plan or explicitly generate an RPE control-identification dataset."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument("--quick",action="store_true",help="Generate a small smoke-test dataset.")
    parser.add_argument("--skip-dynamics-regression",action="store_true")
    args=parser.parse_args(); summary=generate(output=args.output,quick=args.quick,
        skip_dynamics_regression=args.skip_dynamics_regression)
    print(summary)
    return 0


if __name__=="__main__": raise SystemExit(main())
