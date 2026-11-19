"""Echo observability study: run the jobs here, summarize, and draw Viewer CLI figures.

empirical-*.json  empirical_observability of (r, r-dot) through the echo along the plant
                  trajectory of each shadow run (plant on the ideal drive);
krom-*.json       linear_systems on each KLQE K-ROM: observability gramian SVD, readout
                  signature of (r, r-dot), Kalman-predicted readout uncertainty.

Usage: python run_observability.py   (writes observability.json and viewer_outputs/echo-observability-20261118/)
"""
import glob, json, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
BIN = os.path.join(ROOT, "build", "linear-systems", "bin")
OUT = os.path.join(ROOT, "viewer_outputs", "echo-observability-20261118")
RED, BLUE, GOLD = "#D8242A", "#004DB8", "#DE9B00"
COLORS = {"mild": BLUE, "strong": RED, "hard": GOLD}


def run(job):
    out = subprocess.run([os.path.join(BIN, "dynsyslab-cli.exe"), "run", job, "--archive-root", os.path.join(ROOT, "archives")],
                         capture_output=True, text=True, cwd=ROOT, check=True)
    return json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]


def viewer(*args):
    subprocess.run([os.path.join(BIN, "DynSysLabViewerCli.exe"), *args], capture_output=True, text=True, cwd=ROOT, check=True)


def series(name, curves, axis, title, scale="log", height=700):
    """One time-series figure from {key: (array, label, color, style)} on the sample-time clock."""
    args = ["time-series", "--source", os.path.join(OUT, name + ".npz"), "--time-key", "time", "--time-label", "tau", "--width", "1400",
            "--height", str(height), "--scale", scale, "--axis", axis, "--title", title]
    for key, (_, label, color, style) in curves.items():
        args += ["--curve", f"{key}={label}", "--color", color, "--line-style", style]
    viewer(*args, "--output", os.path.join(OUT, name + ".png"))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    summary = {"empirical": {}, "krom": {}}
    emp = {}
    for target in ("mild", "strong", "hard"):
        z = np.load(run(os.path.join(HERE, f"empirical-{target}.json")))
        emp[target] = z
        rel = z["relative_estimation_bound"]
        summary["empirical"][target] = {
            "median_scaled_eigenvalues": np.median(z["scaled_eigenvalues"], axis=0).tolist(),
            "median_scaled_condition": float(np.median(z["scaled_condition_number"])),
            "median_relative_bound_r_rdot": np.median(rel, axis=0).tolist(),
            "echo_jacobian_scaled_median": np.median(np.abs(z["output_jacobian"][:, 0, :] * z["state_scale"]), axis=0).tolist()}
    t = emp["mild"]["sample_times"]

    def save(name, arrays):
        np.savez(os.path.join(OUT, name + ".npz"), time=t, **arrays)

    curves = {}
    for target, z in emp.items():
        curves[f"l1_{target}"] = (z["scaled_eigenvalues"][:, 0], f"{target}: largest", COLORS[target], "solid")
        curves[f"l2_{target}"] = (z["scaled_eigenvalues"][:, 1], f"{target}: smallest", COLORS[target], "dash")
    save("empirical_eigenvalues", {k: v[0] for k, v in curves.items()})
    series("empirical_eigenvalues", curves, "eigenvalue (orbit-scaled)",
           "Echo observability of (r, r-dot) along each orbit: empirical gramian eigenvalues over one period")
    curves = {}
    for target, z in emp.items():
        curves[f"r_{target}"] = (z["relative_estimation_bound"][:, 0], f"{target}: r", COLORS[target], "solid")
        curves[f"v_{target}"] = (z["relative_estimation_bound"][:, 1], f"{target}: r-dot", COLORS[target], "dash")
    save("estimation_bound", {k: v[0] for k, v in curves.items()})
    series("estimation_bound", curves, "Cramer-Rao std / orbit RMS",
           "Best possible estimate from one period of echo (noise intensity 1e-4), relative to the orbit size")
    curves = {}
    for target, z in emp.items():
        j = np.abs(z["output_jacobian"][:, 0, :] * z["state_scale"])
        curves[f"jr_{target}"] = (j[:, 0], f"{target}: |de/dr| x r-scale", COLORS[target], "solid")
        curves[f"jv_{target}"] = (j[:, 1], f"{target}: |de/dr-dot| x r-dot-scale", COLORS[target], "dash")
    save("echo_jacobian", {k: v[0] for k, v in curves.items()})
    series("echo_jacobian", curves, "instantaneous sensitivity",
           "What one echo sample sees: |de/dr| vs |de/dr-dot| (orbit-scaled)")
    curves = {}
    for target, z in emp.items():
        curves[f"c_{target}"] = (z["scaled_condition_number"], target, COLORS[target], "solid")
    save("condition_number", {k: v[0] for k, v in curves.items()})
    series("condition_number", curves, "condition number", "Empirical gramian condition number (orbit-scaled)")

    # K-ROM observability gramian singular values, normalized by the largest, on an index axis.
    sv = {}
    for job in sorted(glob.glob(os.path.join(HERE, "krom-*.json"))):
        name = os.path.basename(job)[5:-5]
        z = np.load(run(job)); meta = json.loads(str(z["metadata_json"]))
        s = z["observability_singular_values"]
        sv[name] = s
        summary["krom"][name] = {"singular_values": s.tolist(), **meta["readout_observability"],
                                 "readout_signature": z["readout_direction_signature"].tolist()}
    index = np.arange(1, 5, dtype=float)
    styles = {"ee-orbit-mild": (BLUE, "solid"), "ee-orbit-strong": (RED, "solid"), "ee-orbit-hard": (GOLD, "solid"),
              "poly5": (BLUE, "dash"), "rr": (RED, "dash")}
    np.savez(os.path.join(OUT, "krom_singular_values.npz"), index=index,
             **{k.replace("-", "_"): np.maximum(v / v[0], 1e-17) for k, v in sv.items()})
    args = ["time-series", "--source", os.path.join(OUT, "krom_singular_values.npz"), "--time-key", "index",
            "--time-label", "singular value index k", "--scale", "log", "--width", "1400", "--height", "700",
            "--axis", "sigma_k / sigma_1", "--title", "K-ROM observability gramian through the echo (split realization; "
            "k = 3, 4 are the duplicate conjugate copies)"]
    for k, (color, style) in styles.items():
        args += ["--curve", f"{k.replace('-', '_')}={k}", "--color", color, "--line-style", style]
    viewer(*args, "--output", os.path.join(OUT, "krom_singular_values.png"))
    json.dump(summary, open(os.path.join(HERE, "observability.json"), "w"), indent=2)
    print(json.dumps(summary, indent=1)[:3000])
