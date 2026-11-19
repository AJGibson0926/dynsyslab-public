"""Stem plots of the observability-direction magnitudes per target (Viewer CLI `--line-style stem`).

For each target (mild, strong, hard):
1. K-ROM observability gramian singular values sigma_k / sigma_1 through the echo (that target's
   psi(e, e') K-ROM, smooth H-DMD, psi(r, r-dot)).
2. Singular values of the classical observability matrix obsv(A, C) = [C; CA; CA^2; CA^3], / sigma_1.
3. Physical (empirical) observability gramian eigenvalues of (r, r-dot) at each sample point of the orbit.
Run run_observability.py first (observability.json). Writes obsv-*.json jobs and runs them.
"""
import glob, json, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
BIN = os.path.join(ROOT, "build", "linear-systems", "bin")
OUT = os.path.join(ROOT, "viewer_outputs", "echo-observability-20261118")
RED, BLUE, GOLD = "#D8242A", "#004DB8", "#DE9B00"
# One figure set per target: that target's psi(e, e') K-ROM (orbit readouts) beside the
# target-independent smooth H-DMD and psi(r, r-dot) K-ROMs, one colour per model.
TARGETS = ("mild", "strong", "hard")
MODELS = {"ee": ("psi(e, e'), {t} orbit readouts", RED), "poly5": ("smooth H-DMD", BLUE), "rr": ("psi(r, r-dot)", GOLD)}
FLOOR = 1e-17  # exact zeros (duplicate conjugate copies) drawn at this floor on the log axis


def run(job):
    out = subprocess.run([os.path.join(BIN, "dynsyslab-cli.exe"), "run", job, "--archive-root", os.path.join(ROOT, "archives")],
                         capture_output=True, text=True, cwd=ROOT, check=True)
    return json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]


def stem(name, series, xlabel, axis, title, offsets=None):
    """series: {key: (positions, values, label, color)} placed on one shared increasing clock."""
    clock = sorted({p for pos, *_ in series.values() for p in pos})
    arrays = {}
    for key, (pos, vals, *_ ) in series.items():
        a = np.full(len(clock), np.nan)
        for p, v in zip(pos, vals): a[clock.index(p)] = v
        arrays[key] = a
    src = os.path.join(OUT, name + ".npz")
    np.savez(src, k=np.array(clock), **arrays)
    args = [os.path.join(BIN, "DynSysLabViewerCli.exe"), "time-series", "--source", src, "--time-key", "k",
            "--time-label", xlabel, "--scale", "log", "--axis", axis, "--title", title, "--width", "1400", "--height", "750"]
    for key, (_, _, label, color) in series.items():
        args += ["--curve", f"{key}={label}", "--color", color, "--line-style", "stem", "--line-width", "2.5"]
    subprocess.run(args + ["--output", src[:-4] + ".png"], check=True, capture_output=True, text=True, cwd=ROOT)


def obsv_singular_values(name, summary):
    krom = json.load(open(os.path.join(HERE, f"krom-{name}.json")))["sources"][0]
    job = os.path.join(HERE, f"obsv-{name}.json")
    json.dump({"schema": "dynsyslab.compute-job/v1", "kind": "postprocess", "label": f"obsv(A,C) / {name}",
               "description": "Classical observability matrix [C; CA; CA^2; CA^3] of the K-ROM with the echo output and its SVD.",
               "sources": [krom], "methods": [{"name": "linear_analysis", "enabled": True, "settings": {
                   "measurement": "far_field_acoustic_echo", "analysis_method": "matrix", "blocks": 4}}],
               "archive": {"enabled": True}}, open(job, "w"), indent=2)
    s = np.linalg.svd(np.load(run(job))["observability_matrix"], compute_uv=False)  # plotting summary of the archived matrix
    summary["krom"][name]["obsv_singular_values"] = s.tolist()
    return s


def latest_empirical(target):
    files = sorted(glob.glob(os.path.join(ROOT, "archives", "derived", "empirical-observability",
                                          f"shadow-{target}-ee-orbit__*__empirical-observability__*.npz")))
    return np.load(files[-1])


if __name__ == "__main__":
    summary = json.load(open(os.path.join(HERE, "observability.json")))
    obsv = {n: obsv_singular_values(n, summary) for n in ("poly5", "rr")}
    for t in TARGETS:
        kroms = {"ee": f"ee-orbit-{t}", "poly5": "poly5", "rr": "rr"}
        obsv["ee"] = obsv_singular_values(kroms["ee"], summary)
        offset = {m: -0.12 + 0.12 * i for i, m in enumerate(MODELS)}
        label = lambda m, s1: MODELS[m][0].format(t=t) + f" (sigma_1 = {s1:.3g})"
        series = {}
        for m, name in kroms.items():
            s = np.array(summary["krom"][name]["singular_values"])
            series[m] = ([k + 1 + offset[m] for k in range(len(s))], np.maximum(s / s[0], FLOOR), label(m, s[0]), MODELS[m][1])
        stem(f"stem_{t}_krom_gramian", series, "observability direction k (split realization; k = 3, 4: duplicate conjugate copies)",
             "sigma_k / sigma_1", f"{t}: K-ROM observability gramian through the echo, singular value of each direction")
        series = {}
        for m in kroms:
            s = obsv[m]
            series[m] = ([k + 1 + offset[m] for k in range(len(s))], np.maximum(s / s[0], FLOOR), label(m, s[0]), MODELS[m][1])
        stem(f"stem_{t}_krom_obsv", series, "observability direction k (right singular vectors of obsv(A, C))",
             "sigma_k / sigma_1", f"{t}: classical observability matrix obsv(A, C) = [C; CA; CA^2; CA^3], singular values")
        z = latest_empirical(t)
        times, e = z["sample_times"], z["scaled_eigenvalues"]
        series = {"l1": (list(times - 0.15), e[:, 0], "direction 1 (strongest)", RED),
                  "l2": (list(times + 0.15), e[:, 1], "direction 2 (weakest, mostly r-dot)", BLUE)}
        stem(f"stem_{t}_empirical_gramian", series, "tau (sample point on the orbit; one-period window)",
             "eigenvalue (orbit-scaled)", f"{t} orbit: physical echo observability of (r, r-dot), empirical gramian eigenvalues")
    json.dump(summary, open(os.path.join(HERE, "observability.json"), "w"), indent=2)
    print("ok")
