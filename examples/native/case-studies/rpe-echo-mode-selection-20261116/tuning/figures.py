"""Viewer CLI figures and summary table for the KLQE tuning sweep (run run_tuning.py first)."""
import json, math, os, subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
VIEWER = os.path.join(ROOT, "build", "linear-systems", "bin", "DynSysLabViewerCli.exe")
OUT = os.path.join(ROOT, "viewer_outputs", "klqe-tuning-20261119")
RED, BLUE, GOLD = "#D8242A", "#004DB8", "#DE9B00"
KROMS = {"ee-orbit": ("psi(e, e')", RED), "poly5": ("smooth H-DMD", BLUE), "rr": ("psi(r, r-dot)", GOLD)}
SCALARS = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6]


def rows():
    return [json.loads(line) for line in open(os.path.join(HERE, "tuning-results.jsonl"))]


def value(r, key):
    v = r.get(key)
    return v if isinstance(v, (int, float)) and math.isfinite(v) else float("nan")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    data = rows()
    get = {(r["target"], r["krom"], r["variant"]): r for r in data}
    x = np.log10(SCALARS)
    for target in ("mild", "strong", "hard"):
        for metric, axis in (("rdot_rel", "r-dot estimate RMS error / r-dot RMS"), ("r_rel", "r estimate RMS error / r RMS"),
                             ("force_rel", "KLQR force error / ideal drive RMS")):
            arrays, args = {}, []
            for krom, (label, color) in KROMS.items():
                key = f"{krom.replace('-', '_')}"
                arrays[key] = np.array([value(get.get((target, krom, f"W{w:g}"), {}), metric) for w in SCALARS])
                args += ["--curve", f"{key}={label}", "--color", color]
            src = os.path.join(OUT, f"{target}_{metric}.npz")
            np.savez(src, log10W=x, **arrays)
            subprocess.run([VIEWER, "time-series", "--source", src, "--time-key", "log10W", "--time-label", "log10 W (scalar process covariance)",
                            "--scale", "log", "--axis", axis, "--title", f"{target}: shadow KLQE, scalar W sweep (V = 1e-4)",
                            "--width", "1400", "--height", "700", *args, "--output", src[:-4] + ".png"],
                           check=True, capture_output=True, text=True, cwd=ROOT)
    # summary: best scalar W, residual-shaped and input-direction variants
    table = {}
    for target in ("mild", "strong", "hard"):
        for krom in KROMS:
            cand = [r for r in data if r["target"] == target and r["krom"] == krom and "error" not in r]
            def best(prefix):
                c = [r for r in cand if r["variant"].startswith(prefix) and math.isfinite(value(r, "rdot_rel"))]
                return min(c, key=lambda r: value(r, "rdot_rel")) if c else None
            entry = {}
            for prefix in ("W", "residual", "input"):
                b = best(prefix)
                if b: entry[prefix] = {k: round(value(b, k), 4) for k in ("r_rel", "rdot_rel", "force_rel", "rho_0.1", "innovation_rel")} | {"variant": b["variant"]}
            base = get.get((target, krom, "W10000" if krom != "poly5" else "W100000"))
            if base and "error" not in base:
                entry["study_setting"] = {k: round(value(base, k), 4) for k in ("r_rel", "rdot_rel", "force_rel", "rho_0.1")} | {"variant": base["variant"]}
            table[f"{target}/{krom}"] = entry
    json.dump(table, open(os.path.join(HERE, "tuning-summary.json"), "w"), indent=1)
    for k, v in table.items():
        print(k, json.dumps(v))
