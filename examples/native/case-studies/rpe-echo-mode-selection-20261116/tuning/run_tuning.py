"""KLQE tuning sweep on the shadow runs (plant on the ideal drive; echo-only KLQE -> KLQR in shadow).

For each target (mild, strong, hard) and K-ROM (psi(e, e') orbit readouts, smooth H-DMD, psi(r, r-dot)):
  scalar W = w I                         w in 1e-2 ... 1e6
  residual-shaped W (K-ROM one-step prediction residual on the 2026-11-15 full-state run of the
  same target, tau 20-100), times 1, 10, 100
  W = I + q2 B B^T (input-direction)     q2 in 1e2, 1e4, 1e6
Scores on tau 40-100: r-hat and r-dot-hat RMS error relative to the orbit RMS, shadow-force error
relative to the ideal drive, and innovation whiteness (normalized autocorrelation at 0.01, 0.1, 0.5 tau).

Usage: python run_tuning.py [--resume]   (appends results to tuning-results.jsonl)
"""
import json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LOOP = os.path.join(HERE, "..", "closed-loop")
sys.path.insert(0, LOOP)
from echo_loop import ROOT, CLI, REFERENCE  # noqa: E402
import shadow_runs as S  # noqa: E402

RESULTS = os.path.join(HERE, "tuning-results.jsonl")
JOBS = os.path.join(HERE, "jobs")
RECORDS = {t: ROOT + "/archives/solutions/rpe/" + REFERENCE[t] for t in REFERENCE}
SCALARS = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6]


def variants(target):
    for w in SCALARS:
        yield f"W{w:g}", {"process_covariance": w}
    for scale in (1, 10, 100):
        yield f"residual-x{scale}", {"process_covariance": {"from_records": [RECORDS[target]], "window": [20, 100],
                                                            "scale": scale, "floor": 1e-6}}
    for q2 in (1e2, 1e4, 1e6):
        yield f"input-q{q2:g}", {"process_covariance": 1.0, "process_covariance_input_weight": q2}


def autocorrelation(x, lag):
    x = x - x.mean()
    d = float(np.dot(x, x))
    return float(np.dot(x[:-lag], x[lag:]) / d) if d > 0 and lag < len(x) else float("nan")


def score(path):
    base = S.score(path)
    z = np.load(path)
    te = z["estimator_time"]; inn = z["innovation"].ravel(); y = z["measurement_values"].ravel()
    s = (te >= 40) & (te <= 100)
    dt = te[1] - te[0]
    base.update({"r_rel": base["r_err_rms"] / base["r_rms"], "rdot_rel": base["rdot_err_rms"] / base["rdot_rms"],
                 "innovation_rel": float(np.sqrt(np.mean(inn[s] ** 2)) / np.sqrt(np.mean((y[s] - y[s].mean()) ** 2))),
                 **{f"rho_{lag}": autocorrelation(inn[s], max(1, int(round(lag / dt)))) for lag in (0.01, 0.1, 0.5)}})
    return base


def done():
    if not os.path.exists(RESULTS): return set()
    return {(r["target"], r["krom"], r["variant"]) for r in map(json.loads, open(RESULTS))}


if __name__ == "__main__":
    os.makedirs(JOBS, exist_ok=True)
    skip = done() if "--resume" in sys.argv else set()
    if "--resume" not in sys.argv and os.path.exists(RESULTS): os.remove(RESULTS)
    for target in ("mild", "strong", "hard"):
        for krom in ("ee-orbit", "poly5", "rr"):
            for name, override in variants(target):
                if (target, krom, name) in skip: continue
                spec = S.job(target, krom)
                spec["label"] = f"tuning-{target}-{krom}-{name}"
                spec["estimator"].update(override)
                path = os.path.join(JOBS, spec["label"] + ".json")
                json.dump(spec, open(path, "w"), indent=1)
                out = subprocess.run([CLI, "run", path, "--archive-root", ROOT + "/archives"], capture_output=True, text=True, cwd=ROOT)
                row = {"target": target, "krom": krom, "variant": name}
                if out.returncode:
                    row["error"] = out.stderr.strip()[-300:]
                else:
                    archive = json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]
                    try: row.update(archive=archive, **score(archive))
                    except Exception as e: row["error"] = f"score: {e}"
                with open(RESULTS, "a") as f: f.write(json.dumps(row) + "\n")
                print(target, krom, name, {k: round(v, 4) for k, v in row.items() if k in ("r_rel", "rdot_rel", "force_rel", "rho_0.1")} or row.get("error"), flush=True)
