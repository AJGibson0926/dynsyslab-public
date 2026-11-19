"""Echo-only KLQE -> KLQR closed-loop runs and scores (driver for the mode-selection case study)."""
import glob, json, subprocess, sys, os
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), *[".."] * 5)).replace(os.sep, "/")
CLI = ROOT + "/build/linear-systems/bin/dynsyslab-cli.exe"
JOBS = ROOT + "/examples/native/case-studies/rpe-echo-mode-selection-20261116/closed-loop"
SYSTEM = {"name": "rpe", "preset": "air_water_20c", "parameters": {"ambient_pressure_pa": 100000, "equilibrium_radius_m": 1e-05,
          "kappa": 1.3333333333333333, "mu_pa_s": 0.001, "pressure_input_sign": -1, "rho_kg_m3": 998, "sigma_n_m": 0.0725,
          "vapor_pressure_pa": 2330}}
TARGETS = {
    "mild": {"amplitudes": [0.05], "angular_frequencies": [1.8927362514], "kind": "waveform", "offset": [1, 0], "phases_rad": [0], "space": "state"},
    "strong": {"amplitudes": [0.4], "angular_frequencies": [1.8927362514], "kind": "waveform", "offset": [1.5, 0], "phases_rad": [0], "space": "state"},
    "hard": {"amplitudes": [0.16, 0.35, 0.07, 0.05], "angular_frequencies": [0.94636813, 2.46055713, 3.97474613, 5.48893513],
             "kind": "waveform", "offset": [1.12, 0], "phases_rad": [1.6, 5.2, 0.8, 1.3], "space": "state"},
}
REFERENCE = {"mild": "klqe-klqr-ee-ic37-mild-full-R0-01-ramp0__20261115T014805Z.npz",
             "strong": "klqe-klqr-ee-ic37-strong-full-R0-01-ramp0__20261115T014807Z.npz",
             "hard": "klqe-klqr-ee-ic37-hard-full-R0-01-ramp0__20261115T014820Z.npz"}


def job(name, krom, target, Q, R, W, V, stop=100, full_state=False, window=None, cap=None, ramp=None, ctrl_extra=None, est_krom=None, feedback="estimate"):
    j = {"schema": "dynsyslab.compute-job/v1", "kind": "solve", "label": name,
            "description": "Echo-only KLQE -> KLQR: the estimator measures far_field_acoustic_echo alone (y = e0 + C_e psi), "
                           "the controller acts on the estimate with the constant B_psi of a mode_selection K-ROM.",
            "system": SYSTEM, "time": {"start": 0, "stop": stop, "output_step": 0.01, "coordinate_system": "dimensionless"},
            "initial_conditions": {"mode": "explicit", "count": 1, "clusters": 1, "seed": 0, "values": [[1, 0]], "ranges": []},
            "controller": {"kind": "klqr", "settings": {"krom": krom, "Q": Q, "R": R, "control_representation": "constant",
                                                        "feedback": "estimate", "control_limit": None},
                           "target": TARGETS[target]},
            "estimator": {"kind": "klqe", "krom": krom, "measurement": "far_field_acoustic_echo", "output": "state",
                          "process_covariance": W, "measurement_covariance": V, "initial_estimate": "zero",
                          "control_representation": "constant", "sample_time": 0.001},
            "solver": {"method": "BDF", "rtol": 1e-09, "atol": 1e-09, "max_step": 0.001,
                       "control_evaluation": "sampled", "control_sample_time": 0.001,
                       "on_failure": "truncate"},  # keep the computed part of a run that blows up
            "archive": {"enabled": True}}
    if full_state:  # diagnostic only: KLQR on the true lifted state, no estimator
        del j["estimator"]; del j["controller"]["settings"]["feedback"]
    if window: j["controller"]["settings"]["active_time_window"] = window
    if cap is not None: j["controller"]["settings"]["control_limit"] = cap
    if ramp: j["controller"]["target"] = dict(j["controller"]["target"], ramp={"duration": ramp, "from": [1, 0]})
    if ctrl_extra: j["controller"]["settings"].update(ctrl_extra)  # e.g. integral_action, pid_correction
    if est_krom and not full_state: j["estimator"]["krom"] = est_krom  # separate estimation eigenfunctions
    if not full_state and feedback == "estimated_state":
        # KLQE (own KROM, echo only) estimates x_hat; the controller evaluates its own psi at x_hat.
        j["controller"]["settings"]["feedback"] = "estimated_state"
    elif not full_state and feedback == "passive":  # controller on the true state, KLQE observes only
        del j["controller"]["settings"]["feedback"]
    return j


def score(path, target):
    z = np.load(path); t = z["time"]; x = z["values"]; u = z["control_values"].ravel()
    ref = np.load(ROOT + "/archives/solutions/rpe/" + REFERENCE[target])
    rt, ru = ref["time"], ref["control_values"].ravel()
    tgt = ref["control_target_values"]
    if t[-1] < 40: return {"terminated_at": float(t[-1]), "finite": bool(np.all(np.isfinite(x)))}
    s = (t >= 40) & (t <= min(t[-1], rt[-1]))
    ur = np.interp(t[s], rt, ru); rr = np.interp(t[s], rt, tgt[:, 0])
    drive = np.sqrt(np.mean((u[s] - ur) ** 2) / np.mean((ur - ur.mean()) ** 2))
    radius = np.sqrt(np.mean((x[s, 0] - rr) ** 2))
    def lowpass(v, dt, cutoff=2.0):  # cycles per tau
        V = np.fft.rfft(v); f = np.fft.rfftfreq(len(v), dt); V[f > cutoff] = 0; return np.fft.irfft(V, len(v))
    dt = t[1] - t[0]; ul, url = lowpass(u[s], dt), lowpass(ur, dt)
    drive_lp = np.sqrt(np.mean((ul - url) ** 2) / np.mean((url - url.mean()) ** 2))
    est = None
    if "state_estimate" in z.files:
        te = z["estimator_time"]; xe = z["state_estimate"]; se = te >= 40
        est = np.sqrt(np.mean((xe[se, 0] - np.interp(te[se], t, x[:, 0])) ** 2))
    return {"drive_error": float(drive), "drive_error_lowpass": float(drive_lp), "radius_rmse": float(radius), "radius_estimate_rmse": None if est is None else float(est),
            "finite": bool(np.all(np.isfinite(x)))}


def run(name, krom, target, Q, R, W, V, stop=100, full_state=False, window=None, cap=None, ramp=None, ctrl_extra=None, est_krom=None, feedback="estimate"):
    os.makedirs(JOBS, exist_ok=True)
    path = f"{JOBS}/{name}.json"
    json.dump(job(name, krom, target, Q, R, W, V, stop, full_state, window, cap, ramp, ctrl_extra, est_krom, feedback), open(path, "w"), indent=2)
    out = subprocess.run([CLI, "run", path, "--archive-root", ROOT + "/archives"], capture_output=True, text=True, cwd=ROOT)
    try:
        archive = json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]
    except Exception:
        return {"error": (out.stdout + out.stderr)[-300:]}
    result = score(archive, target); result["archive"] = os.path.basename(archive)
    return result


if __name__ == "__main__":
    spec = json.loads(sys.argv[1])
    print(json.dumps(run(**spec)))
