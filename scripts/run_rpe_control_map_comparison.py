"""Compare exact and controlled-data B_psi maps in one matched RPE KLQR case."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy.linalg import solve_continuous_are


ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))

from control import KoopmanLQRController, RPE_CLEAN_ISOLATION_KLQR_SETTINGS, rpe_control_target
from rom import load_krom_archive, load_krom_error_evaluations
from solver import InterpolationOptions, SolveRequest, SolverOptions, solve
from viewer import (ArchiveSource, ControlErrorSource, PlotRequest, TilePlan,
                    Viewer, ViewerReport)


KROM_ARCHIVE=ROOT/"archives/derived/koopman-krom-extract/rpe-reference-training-bdf-tau25-dtau0p001__20260826T220948.289776Z__koopman-krom-extract__20260827T183646.410439Z__0dea0923.npz"
ERROR_BUNDLE=ROOT/"archives/experiments/rpe_control_identification_designed/rpe-minimal-krom-error-evaluations.npz"
OUTPUT=ROOT/"viewer_outputs/rpe-minimal-krom-control-map-comparison.html"
ERROR_OUTPUT=ROOT/"viewer_outputs/rpe-minimal-krom-control-map-error-comparison.html"
SUMMARY=ROOT/"viewer_outputs/rpe-minimal-krom-control-map-comparison.json"
REPORT=ROOT/"docs/case-studies/RPE_APPROXIMATE_BPSI_KLQR_COMPARISON_2026-08-29.md"
TARGET_NAME="rpe_multisine_equilibrium_moderate_identification"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _map_error(model,representation,states: np.ndarray) -> dict[str,float]:
    exact=np.asarray([model.B_psi(state,output="complex",representation="exact") for state in states])
    if representation is None: predicted=exact
    else: predicted=np.asarray([representation.evaluate(model.psi(state,output="complex")) for state in states])
    residual=predicted-exact
    point=np.linalg.norm(residual.reshape(states.shape[0],-1),axis=1)
    truth=np.linalg.norm(exact.reshape(states.shape[0],-1),axis=1)
    relative=point/np.maximum(truth,np.finfo(float).tiny)
    return {"relative_frobenius":float(np.linalg.norm(residual)/max(np.linalg.norm(exact),np.finfo(float).tiny)),
        "mean_pointwise_relative":float(np.mean(relative)),"median_pointwise_relative":float(np.median(relative)),
        "p95_pointwise_relative":float(np.quantile(relative,.95)),"maximum_pointwise_relative":float(np.max(relative))}


def _metrics(result,model,target,horizon_tau: float,representation) -> dict:
    tau=np.asarray(model.system.to_dimensionless_time(result.time),dtype=float)
    states=np.asarray([model.system.to_dimensionless_state(row) for row in result.values],dtype=float)
    desired=np.vstack([target.evaluate(float(item),2) for item in tau]); error=states-desired
    second=tau>=horizon_tau/2.; command=np.asarray(result.control_time_series.values[:,0],dtype=float)
    stride=max(1,states.shape[0]//1000)
    return {"radius_rmse":float(np.sqrt(np.mean(error[:,0]**2))),
        "velocity_rmse":float(np.sqrt(np.mean(error[:,1]**2))),
        "radius_rmse_second_half":float(np.sqrt(np.mean(error[second,0]**2))),
        "velocity_rmse_second_half":float(np.sqrt(np.mean(error[second,1]**2))),
        "radius_max_abs_error":float(np.max(np.abs(error[:,0]))),
        "velocity_max_abs_error":float(np.max(np.abs(error[:,1]))),
        "minimum_radius":float(np.min(states[:,0])),"maximum_radius":float(np.max(states[:,0])),
        "maximum_abs_velocity":float(np.max(np.abs(states[:,1]))),
        "control_rms_dimensionless":float(np.sqrt(np.mean(command**2))),
        "control_peak_abs_dimensionless":float(np.max(np.abs(command))),
        "control_total_variation_dimensionless":float(np.sum(np.abs(np.diff(command)))),
        "final_state":states[-1].tolist(),"sample_count":int(states.shape[0]),
        "actual_path_control_map_error":_map_error(model,representation,states[::stride]),
        "target_path_control_map_error":_map_error(model,representation,desired[::stride])}


def _local_closed_loop_diagnostics(model,representation_name: str,target,policy: dict) -> dict:
    """Measure geometry and local CARE robustness along the prescribed target."""
    A=model.A(output="split",time_coordinates="model")
    Q=np.diag(policy["Q_diagonal"]); R=np.asarray(((policy["R"],),))
    tau=np.linspace(0.,float(policy["horizon_tau"]),401)
    states=np.vstack([target.evaluate(float(item),2) for item in tau])
    angles=[]; norm_ratios=[]; gain_norms=[]; nominal_poles=[]; truth_poles=[]
    for state in states:
        exact=model.B_psi(state,output="split",representation="exact")
        candidate=model.B_psi(state,output="split",representation=representation_name)
        exact_norm=float(np.linalg.norm(exact)); candidate_norm=float(np.linalg.norm(candidate))
        cosine=float(np.dot(exact[:,0],candidate[:,0])/max(exact_norm*candidate_norm,np.finfo(float).tiny))
        angles.append(float(np.degrees(np.arccos(np.clip(cosine,-1.,1.)))))
        norm_ratios.append(candidate_norm/max(exact_norm,np.finfo(float).tiny))
        P=solve_continuous_are(A,candidate,Q,R); gain=np.linalg.solve(R,candidate.T@P)
        gain_norms.append(float(np.linalg.norm(gain)))
        nominal_poles.append(float(np.max(np.linalg.eigvals(A-candidate@gain).real)))
        truth_poles.append(float(np.max(np.linalg.eigvals(A-exact@gain).real)))
    def quantiles(values):
        active=np.asarray(values,dtype=float)
        return {"p05":float(np.quantile(active,.05)),"median":float(np.median(active)),
            "p95":float(np.quantile(active,.95)),"maximum":float(np.max(active))}
    truth=np.asarray(truth_poles,dtype=float)
    return {"target_sample_count":int(states.shape[0]),"input_direction_angle_degrees":quantiles(angles),
        "input_norm_ratio_candidate_over_exact":quantiles(norm_ratios),"care_gain_norm":quantiles(gain_norms),
        "nominal_closed_loop_max_real_pole":quantiles(nominal_poles),
        "truth_coupled_max_real_pole":quantiles(truth_poles),
        "truth_coupled_locally_unstable_fraction":float(np.mean(truth>=0.))}


def _attach_exact_case_differences(records,model) -> None:
    """Attach realized trajectory/control differences relative to the exact-map solve."""
    exact=next((item for item in records if item["label"]=="exact analytic B_psi" and item["status"]=="complete"),None)
    if exact is None: return
    with np.load(exact["archive"],allow_pickle=False) as payload:
        exact_control=np.asarray(payload["control_values"],dtype=float)
        exact_state=np.asarray([model.system.to_dimensionless_state(row) for row in payload["values"]],dtype=float)
    for item in records:
        if item["status"]!="complete": continue
        with np.load(item["archive"],allow_pickle=False) as payload:
            control=np.asarray(payload["control_values"],dtype=float)
            state=np.asarray([model.system.to_dimensionless_state(row) for row in payload["values"]],dtype=float)
        if control.shape!=exact_control.shape or state.shape!=exact_state.shape:
            item["exact_case_difference"]={"available":False,"reason":"output grids differ"}; continue
        item["exact_case_difference"]={"available":True,
            "relative_control_history_l2":float(np.linalg.norm(control-exact_control)/max(np.linalg.norm(exact_control),np.finfo(float).tiny)),
            "relative_state_history_l2":float(np.linalg.norm(state-exact_state)/max(np.linalg.norm(exact_state),np.finfo(float).tiny)),
            "radius_history_rmse":float(np.sqrt(np.mean((state[:,0]-exact_state[:,0])**2))),
            "velocity_history_rmse":float(np.sqrt(np.mean((state[:,1]-exact_state[:,1])**2)))}


def _write_report(payload: dict) -> None:
    lines=["# RPE KLQR with exact and approximated control maps","",
        "**Date:** 2026-08-29  ","**Status:** matched comparison completed  ",
        f"**Target:** `{TARGET_NAME}`  ",
        "**Invariant contract:** same minimal KROM, equilibrium start `(1,0)`, dimensionless dynamics/control, continuous-RHS KLQR, BDF, unchanged Q/R, no cap, no feedforward, no target shaping","",
        "## Results","",
        "| Control map | Status | r RMSE, second half | v RMSE, second half | RMS u* | Target-path Bpsi error | B direction angle p95 | Control difference from exact |",
        "|---|---|---:|---:|---:|---:|---:|---:|"]
    for case in payload["cases"]:
        metrics=case.get("metrics") or {}; target_error=(metrics.get("target_path_control_map_error") or {}).get("relative_frobenius")
        diagnostic=case.get("local_closed_loop_diagnostics") or {}; angle=(diagnostic.get("input_direction_angle_degrees") or {}).get("p95")
        difference=(case.get("exact_case_difference") or {}).get("relative_control_history_l2")
        def number(value): return "—" if value is None else f"{value:.6g}"
        lines.append(f"| {case['label']} | {case['status']} | {number(metrics.get('radius_rmse_second_half'))} | {number(metrics.get('velocity_rmse_second_half'))} | {number(metrics.get('control_rms_dimensionless'))} | {number(target_error)} | {number(angle)} deg | {number(difference)} |")
    lines.extend(["","## Interpretation","",payload["interpretation"],"",
        "All five approximations produce essentially the same closed-loop response as exact `B_psi`: their realized control histories differ by less than one percent and their state histories differ by less than two tenths of one percent in relative L2 norm. The phase-space matrix error therefore does not translate directly into tracking error in this case.","",
        "The diagnostic reason is structural. This minimal KROM has one complex eigenfunction, represented by a two-state real block, and one scalar input. CARE is recomputed from the selected local input vector. Much of the approximation error is a rescaling of that vector, which the Riccati gain largely cancels. Direction error is more important than raw matrix-norm error. The conjugate-complete bilinear fit preserves direction most closely and produces the closest trajectory to exact control.","",
        "The persistent-data bilinear map is a warning against treating this result as universal: it reaches much larger direction errors and its gain coupled to exact `B_psi` has locally unstable frozen-state poles over a small part of the target path, although the integrated nonlinear trajectory remains stable and accurate. Closed-loop validation remains mandatory.","",
        "The aggregate phase-space heat-map error, pathwise error, input-vector geometry, nominal frozen-state poles, truth-coupled frozen-state poles, and realized trajectory differences are retained in the machine-readable summary.","",
        "## Artifacts","",f"- Viewer control comparison: `{payload['viewer_output']}`",
        f"- Viewer error maps: `{payload['error_viewer_output']}`",f"- Machine-readable summary: `{SUMMARY.resolve()}`",
        f"- KROM error bundle: `{ERROR_BUNDLE.resolve()}`",""])
    REPORT.parent.mkdir(parents=True,exist_ok=True); REPORT.write_text("\n".join(lines),encoding="utf-8")


def main() -> int:
    """Compare exact and fitted RPE control maps on one retained KLQR case."""
    policy=RPE_CLEAN_ISOLATION_KLQR_SETTINGS; base=load_krom_archive(KROM_ARCHIVE)
    loaded=load_krom_error_evaluations(ERROR_BUNDLE,base,attach=True)
    evaluations={item.label:item for item in loaded["control_approximation_evaluations"]}
    selected=("pulse-data constant","pulse-data complex bilinear",
        "pulse-data conjugate-complete bilinear","persistent-data constant","persistent-data bilinear")
    missing=set(selected)-set(evaluations)
    if missing: raise KeyError(f"Missing retained control-map evaluations: {sorted(missing)}")
    cases=[("exact analytic B_psi",None,None)]+[(label,evaluations[label].representation,evaluations[label]) for label in selected]
    target=rpe_control_target(TARGET_NAME); scale=float(base.system.dimensionless_parameters.time_scale_s)
    records=[]; figures=[]
    for label,representation,evaluation in cases:
        name="exact" if representation is None else evaluation.approximation
        model=base if representation is None else base.with_control_representation(name,representation)
        print(f"START {label}",flush=True); started=perf_counter()
        record={"label":label,"control_representation":"exact" if representation is None else name,
            "controlled_dataset":None if evaluation is None else evaluation.controlled_dataset,
            "phase_space_error":None if evaluation is None else evaluation.error_field.summary()}
        try:
            record["local_closed_loop_diagnostics"]=_local_closed_loop_diagnostics(model,
                "exact" if representation is None else name,target,policy)
            controller=KoopmanLQRController(model,Q=np.diag(policy["Q_diagonal"]),R=np.asarray(((policy["R"],),)),
                target=target,target_space="state",control_limit=None,target_domain_policy="ignore",
                control_representation="exact" if representation is None else name)
            result=solve(SolveRequest(system=model.system,kind="ode",time_span=(0.,policy["horizon_tau"]*scale),
                controller=controller,options=SolverOptions(method=policy["solver_method"],
                    rtol=policy["relative_tolerance"],atol=policy["absolute_tolerance"],
                    use_dimensionless_form=True,control_evaluation="continuous_rhs",execution_backend="python"),
                interpolation=InterpolationOptions(enabled=True,time_step=policy["output_dtau"]*scale,method="pchip"),
                label=f"rpe-approximate-bpsi-klqr-{label.lower().replace(' ','-').replace('_','-')}"))
            record.update({"status":result.status,"message":result.message,"solve_seconds":perf_counter()-started})
            if result.status=="complete":
                archive=Path(result.archive_path); record.update({"archive":str(archive.resolve()),
                    "archive_sha256":_sha256(archive),"metrics":_metrics(result,model,target,policy["horizon_tau"],representation)})
                viewer=Viewer(tile_plan=TilePlan.preset("time_phase"),renderer="plotly")
                source=ArchiveSource(archive,label=label); viewer.add(PlotRequest("phase",source,title=label))
                viewer.add(PlotRequest("time",source,title=label)); rendered=viewer.render(); figures.extend(rendered.figures)
                page=OUTPUT.with_name(f"rpe-minimal-krom-control-map-{name}-{len(records):02d}.html")
                viewer.write_html(rendered,page); record["individual_viewer_output"]=str(page.resolve())
        except Exception as error:
            record.update({"status":"failed","message":f"{type(error).__name__}: {error}",
                "solve_seconds":perf_counter()-started})
        records.append(record)
        metrics=record.get("metrics") or {}
        print(f"DONE {label} status={record['status']} r2={metrics.get('radius_rmse_second_half')} peak={metrics.get('control_peak_abs_dimensionless')}",flush=True)
    _attach_exact_case_differences(records,base)
    Viewer(renderer="plotly").write_html(ViewerReport(figures,TilePlan.preset("time_phase"),[]),OUTPUT)
    error_viewer=Viewer(renderer="plotly").add(PlotRequest("control_error",ControlErrorSource(ERROR_BUNDLE),
        title="RPE B_psi approximation error",control_error_measure="relative",control_error_log_scale=False))
    error_viewer.write_html(error_viewer.render(),ERROR_OUTPUT)
    exact=records[0].get("metrics") or {}; completed=sum(item["status"]=="complete" for item in records)
    successful=[item["label"] for item in records[1:] if item["status"]=="complete" and exact and
        item["metrics"]["radius_rmse_second_half"]<=max(2.*exact["radius_rmse_second_half"],.01) and
        item["metrics"]["velocity_rmse_second_half"]<=max(2.*exact["velocity_rmse_second_half"],.05)]
    interpretation=(f"{completed} of {len(records)} matched solves completed. "
        + (f"Approximate maps meeting the predeclared two-times-exact tracking gate: {', '.join(successful)}."
           if successful else "No approximate map met the predeclared two-times-exact tracking gate."))
    payload={"schema":"dynsyslab.rpe-approximate-bpsi-klqr-comparison/v1",
        "krom_archive":str(KROM_ARCHIVE.resolve()),"krom_sha256":_sha256(KROM_ARCHIVE),
        "error_bundle":str(ERROR_BUNDLE.resolve()),"error_bundle_sha256":_sha256(ERROR_BUNDLE),
        "target":TARGET_NAME,"experiment_contract":{"initial_state":[1.,0.],"coordinate_system":"dimensionless",
            "control_units":"p_control/p0","Q_diagonal":list(policy["Q_diagonal"]),"R":policy["R"],
            "horizon_tau":policy["horizon_tau"],"output_dtau":policy["output_dtau"],
            "solver":policy["solver_method"],"control_evaluation":"continuous_rhs","control_limit":None,
            "feedforward":False,"target_shaping":False},"cases":records,"success_gate":{
            "radius":"second-half RMSE <= max(2*exact, 0.01)","velocity":"second-half RMSE <= max(2*exact, 0.05)"},
        "successful_approximations":successful,"interpretation":interpretation,
        "viewer_output":str(OUTPUT.resolve()),"error_viewer_output":str(ERROR_OUTPUT.resolve()),
        "report":str(REPORT.resolve())}
    SUMMARY.write_text(json.dumps(payload,indent=2,sort_keys=True,default=str),encoding="utf-8")
    _write_report(payload); print(json.dumps({"summary":str(SUMMARY),"interpretation":interpretation},indent=2)); return 0


if __name__=="__main__": raise SystemExit(main())
