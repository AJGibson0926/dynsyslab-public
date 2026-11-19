"""Per-eigenfunction observability of K-ROMs from the echo (linear_analysis, analysis_method
"eigenfunction-observability"), drawn as stem plots with the Viewer CLI.

For each K-ROM, over a window of `--horizon` tau (default 1.0), noise intensity V (default 1e-4)
and the KLQE process covariance W (default 1e4):
  eigobs_<name>_visibility.png  per eigenfunction: signal sqrt(G_nn) and distinguishable 1/sqrt((G^-1)_nn)
  eigobs_<name>_precision.png   per eigenfunction: Cramer-Rao std sqrt(V (G^-1)_nn), its conjugate-pair
                                worst case, and the KLQE's steady-state std

Usage: python eigenfunction_observability.py [--horizon T] [--krom name=path ...]
(default K-ROMs: the psi(r, r-dot) variants of the study).
"""
import argparse, json, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
BIN = os.path.join(ROOT, "build", "linear-systems", "bin")
OUT = os.path.join(ROOT, "viewer_outputs", "echo-observability-20261118")
D = ROOT.replace("\\", "/") + "/archives/derived/koopman-krom-extract/"
DEFAULT = {"rr-recentred": D + "experiment-rr-ic0p37-recentred-krom.npz",
           "rr-anchored": D + "experiment-rr-ic0p37-recentred-anchored-krom.npz",
           "rr-local": D + "experiment-rr-ic0p37-recentred-local-krom.npz",
           "rr-local-pair": D + "experiment-rr-ic0p37-recentred-local-pair-krom.npz"}
RED, BLUE, GOLD = "#D8242A", "#004DB8", "#DE9B00"


def run(job):
    out = subprocess.run([os.path.join(BIN, "dynsyslab-cli.exe"), "run", job, "--archive-root", os.path.join(ROOT, "archives")],
                         capture_output=True, text=True, cwd=ROOT, check=True)
    return json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]


def stems(src, curves, axis, title, out, n):
    args = [os.path.join(BIN, "DynSysLabViewerCli.exe"), "time-series", "--source", src, "--time-key", "x",
            "--time-label", "eigenfunction index", "--integer-x-ticks", "--scale", "log", "--x-range", f"0.5:{n + 0.5}",
            "--axis", axis, "--title", title, "--width", "1100", "--height", "650"]
    for key, label, color in curves: args += ["--curve", f"{key}={label}", "--color", color, "--line-style", "stem"]
    subprocess.run(args + ["--output", out], check=True, capture_output=True, text=True, cwd=ROOT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--horizon", type=float, default=1.0)
    ap.add_argument("--noise", type=float, default=1e-4)
    ap.add_argument("--process", type=float, default=1e4)
    ap.add_argument("--krom", action="append", default=[], help="name=path to a K-ROM archive")
    a = ap.parse_args()
    kroms = dict(k.split("=", 1) for k in a.krom) or DEFAULT
    summary = {}
    for name, path in kroms.items():
        job = os.path.join(HERE, f"eigobs-{name}.json")
        json.dump({"schema": "dynsyslab.compute-job/v1", "kind": "postprocess", "label": f"eigenfunction observability / {name}",
                   "description": "How well each eigenfunction can be determined from the echo over the window (finite-horizon modal gramian).",
                   "sources": [path], "methods": [{"name": "linear_analysis", "enabled": True, "settings": {
                       "measurement": "far_field_acoustic_echo", "analysis_method": "eigenfunction-observability",
                       "horizon": a.horizon, "noise_intensity": a.noise, "process_covariance": a.process}}],
                   "archive": {"enabled": True}}, open(job, "w"), indent=2)
        z = np.load(run(job))
        n = int(z["eigenfunction_index"].size)
        lam = z["continuous_eigenvalues"][:, 0] + 1j * z["continuous_eigenvalues"][:, 1]
        row = {k: z[k].tolist() for k in ("signal", "distinguishable", "isolation", "bound_std", "pair_bound_std", "kalman_std") if k in z}
        row["eigenvalues"] = [[l.real, l.imag] for l in lam]
        summary[name] = row
        x = z["eigenfunction_index"]

        def offset_npz(tag, curves):
            """Curves side by side around each integer index on one shared clock (NaN elsewhere)."""
            shifts = np.linspace(-0.15, 0.15, len(curves)) if len(curves) > 1 else [0.0]
            clock = np.sort(np.concatenate([x + d for d in shifts]))
            arrays = {}
            for (key, values), d in zip(curves, shifts):
                a = np.full(clock.size, np.nan)
                for xi, v in zip(x, values): a[np.argmin(np.abs(clock - (xi + d)))] = v
                arrays[key] = a
            path = os.path.join(OUT, f"eigobs_{name}_{tag}.npz")
            np.savez(path, x=clock, **arrays)
            return path

        lam_text = ", ".join(f"{l.real:.3g}{l.imag:+.3g}i" for l in lam[:4]) + (" ..." if n > 4 else "")
        src = offset_npz("visibility", [("signal", z["signal"]), ("distinguishable", z["distinguishable"])])
        stems(src, [("signal", "signal sqrt(G_nn)", BLUE), ("distinguishable", "distinguishable 1/sqrt((G^-1)_nn)", RED)],
              "echo RMS over the window per unit psi_n", f"{name}: how strongly each eigenfunction shows in the echo (T = {a.horizon:g}; lambda = {lam_text})",
              src[:-4] + ".png", n)
        src = offset_npz("precision", [("bound", z["bound_std"]), ("pair", z["pair_bound_std"]), ("kalman", z["kalman_std"])])
        stems(src, [("bound", "Cramer-Rao std of psi_n", RED), ("pair", "conjugate-pair worst direction", GOLD),
                    ("kalman", f"KLQE steady-state std (W = {a.process:g})", BLUE)],
              "standard deviation of psi_n-hat", f"{name}: how precisely each eigenfunction can be estimated (V = {a.noise:g})",
              src[:-4] + ".png", n)
        print(name, {k: np.round(v, 4).tolist() if isinstance(v, list) else v for k, v in row.items() if k != "eigenvalues"})
    json.dump(summary, open(os.path.join(HERE, "eigenfunction-observability.json"), "w"), indent=2)
