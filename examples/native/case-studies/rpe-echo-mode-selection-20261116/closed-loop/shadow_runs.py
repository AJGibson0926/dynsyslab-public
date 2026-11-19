"""Shadow KLQE -> KLQR runs: the plant receives the known correct drive; the echo-only KLQE -> KLQR
runs alongside it (KLQE predicts with the applied drive), and its force is archived, never applied.

Writes shadow-*.json, runs them with dynsyslab-cli and prints estimate / force errors.
"""
import json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from echo_loop import ROOT, CLI, SYSTEM, TARGETS  # noqa: E402

KROM = ROOT + "/archives/derived/"
REFERENCE_KROM = KROM + "koopman-input-map/experiment-ee-ic0p37-split-krom__input-map-constant__20261115T014750Z.npz"
ORBIT_KROM = {"mild": "20261117T002607Z", "strong": "20261117T002614Z", "hard": "20261117T002622Z"}
POLY5_KROM = (KROM + "koopman-mode-selection/rpe-viscous-echo-long-validation-five-ics-trajectory-0__~7f37764ac80c~Z"
              "__mode-selection-krom__20261116T233044Z.npz")
ESTIMATORS = {  # echo-only KLQE -> KLQR set-ups of the study (no ramp: the ideal drive has none)
    "ee-orbit": lambda t: (KROM + f"koopman-mode-selection/experiment-ee-ic0p37-split-krom__mode-selection-krom__{ORBIT_KROM[t]}.npz",
                           5000, 0.01, 1e4),
    "poly5": lambda t: (POLY5_KROM, 1, 3, 1e5),
    # SINDy psi(r, r-dot), IC 0.37, constant map: the separated-run controller (full-state R 0.01 mild, R 1 strong/hard)
    "rr": lambda t: (KROM + "koopman-input-map/experiment-rr-ic0p37-recentred-krom__input-map-constant__20261115T013524Z.npz",
                     5000, 0.01 if t == "mild" else 1, 1e4),
}


# psi(r, r-dot) on the true state cannot jump to the strong target from rest: ramp it in (as in separated-*).
RAMP = {("rr", "strong"): 20}


def target_spec(target, estimator):
    ramp = RAMP.get((estimator, target))
    return dict(TARGETS[target], ramp={"duration": ramp, "from": [1, 0]}) if ramp else TARGETS[target]


def ideal_drive(target, krom, Q, R, spec=None):
    """The known correct drive: the same KLQR evaluated on the TRUE lifted state (full-state feedback).

    KLQR has no feedforward (u = -K(z - psi_ref)), so the exact inverse drive is not what it would
    command even with a perfect state; the same controller on the true state is. With the same
    controller in the shadow, the force error is the estimation error alone, -K(z_hat - psi(x)).
    Where that full-state KLQR fails (poly5 K-ROM beyond r ~ 1.25: strong, hard) the 2026-11-15
    full-state drive is used and the force error also contains the controller difference.
    """
    if krom == POLY5_KROM and target != "mild":
        krom, Q, R = REFERENCE_KROM, 5000, 0.01
    return {"kind": "klqr", "settings": {"krom": krom, "Q": Q, "R": R, "control_limit": None,
                                         "control_representation": "constant"}, "target": spec or TARGETS[target]}, [1, 0]


def job(target, estimator):
    krom, Q, R, W = ESTIMATORS[estimator](target)
    spec = target_spec(target, estimator)
    drive, x0 = ideal_drive(target, krom, Q, R, spec)
    return {"schema": "dynsyslab.compute-job/v1", "kind": "solve", "label": f"shadow-{target}-{estimator}",
            "description": "Shadow echo-only KLQE -> KLQR: the plant receives the known correct drive "
                           "(controller.settings.applied_controller); the KLQE measures far_field_acoustic_echo alone and "
                           "predicts with the applied drive; the KLQR force is archived as commanded_input, never applied.",
            "system": SYSTEM, "time": {"start": 0, "stop": 100, "output_step": 0.01, "coordinate_system": "dimensionless"},
            "initial_conditions": {"mode": "explicit", "count": 1, "clusters": 1, "seed": 0, "values": [x0], "ranges": []},
            "controller": {"kind": "klqr", "settings": {"krom": krom, "Q": Q, "R": R, "control_representation": "constant",
                                                        "feedback": "estimate", "control_limit": None,
                                                        "applied_controller": drive},
                           "target": spec},
            "estimator": {"kind": "klqe", "krom": krom, "measurement": "far_field_acoustic_echo", "output": "state",
                          "process_covariance": W, "measurement_covariance": 1e-4, "initial_estimate": "zero",
                          "control_representation": "constant", "sample_time": 0.001},
            "solver": {"method": "BDF", "rtol": 1e-09, "atol": 1e-09, "max_step": 0.001,
                       "control_evaluation": "sampled", "control_sample_time": 0.001, "on_failure": "truncate"},
            "archive": {"enabled": True}}


def score(path, window=(40, 100)):
    z = np.load(path)
    t, x, te = z["time"], z["values"], z["estimator_time"]
    xh = z["state_estimate"]; u = z["applied_input"].ravel(); uk = z["commanded_input"].ravel()
    s = (te >= window[0]) & (te <= window[1])
    xt = np.column_stack([np.interp(te, t, x[:, k]) for k in range(2)])
    err = xh - xt
    rms = lambda a: float(np.sqrt(np.mean(a ** 2)))
    return {"r_err_rms": rms(err[s, 0]), "r_err_mean": float(err[s, 0].mean()), "r_rms": rms(xt[s, 0] - xt[s, 0].mean()),
            "rdot_err_rms": rms(err[s, 1]), "rdot_err_mean": float(err[s, 1].mean()), "rdot_rms": rms(xt[s, 1]),
            "force_err_rms": rms(uk[s] - u[s]), "force_rms": rms(u[s]), "force_rel": rms(uk[s] - u[s]) / rms(u[s]),
            "force_corr": float(np.corrcoef(uk[s], u[s])[0, 1])}


if __name__ == "__main__":
    chosen = sys.argv[1:] or list(ESTIMATORS)  # e.g. shadow_runs.py rr
    path_scores = os.path.join(HERE, "shadow-scores.json")
    results = json.load(open(path_scores)) if os.path.exists(path_scores) else {}
    for target in ("mild", "strong", "hard"):
        for estimator in chosen:
            spec = job(target, estimator)
            path = os.path.join(HERE, spec["label"] + ".json")
            json.dump(spec, open(path, "w"), indent=2)
            out = subprocess.run([CLI, "run", path, "--archive-root", ROOT + "/archives"], capture_output=True, text=True, cwd=ROOT)
            if out.returncode:
                print(spec["label"], "FAILED", out.stderr[-800:]); continue
            archive = json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]
            results[spec["label"]] = dict(archive=archive, **score(archive))
            print(spec["label"], json.dumps(results[spec["label"]]), flush=True)
    json.dump(results, open(path_scores, "w"), indent=2)
