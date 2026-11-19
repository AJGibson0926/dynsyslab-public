"""Viewer CLI figures for the shadow KLQE -> KLQR runs (run shadow_runs.py first).

Per run: r and r-dot with target and KLQE estimate, the estimate errors, the KLQR force against the
applied (ideal) drive, and the force error.
"""
import json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from echo_loop import ROOT  # noqa: E402

VIEWER = ROOT + "/build/linear-systems/bin/DynSysLabViewerCli.exe"
OUT = ROOT + "/viewer_outputs/shadow-klqe-klqr-20261118"
RED, BLUE, GOLD = "#D8242A", "#004DB8", "#DE9B00"
NAMES = {"ee-orbit": "psi(e, e') K-ROM, orbit readouts", "poly5": "smooth H-DMD K-ROM", "rr": "psi(r, r-dot) K-ROM"}


def viewer(*args):
    out = subprocess.run([VIEWER, *args], capture_output=True, text=True, cwd=ROOT)
    if out.returncode: raise RuntimeError(out.stderr)


def figures(label, archive, start, stop):
    tag = f"{label}_t{start}-{stop}"
    target, estimator = label.split("-", 2)[1], label.split("-", 2)[2]
    title = f"{target}, {NAMES[estimator]}: "
    z = np.load(archive)
    te, t, x = z["estimator_time"], z["time"], z["values"]
    truth = np.column_stack([np.interp(te, t, x[:, k]) for k in range(2)])
    err = z["state_estimate"] - truth
    ideal, klqr = z["applied_input"].ravel(), z["commanded_input"].ravel()
    data = f"{OUT}/{label}_errors.npz"
    tg = z["control_target_values"]
    target_values = np.column_stack([np.interp(te, t, tg[:, k]) for k in range(2)])
    np.savez(data, time=te, x=truth, xhat=z["state_estimate"], target=target_values, err=err,
             ideal=ideal, klqr=klqr, force_err=klqr - ideal)
    common = ["--source", data, "--start", str(start), "--stop", str(stop), "--time-label", "tau", "--width", "1400"]
    viewer("time-series", *common, "--columns", "1", "--curve", "x=actual", "--color", RED, "--line-width", "2.7",
           "--curve", "target=target", "--color", BLUE, "--line-style", "dash",
           "--curve", "xhat=echo-only KLQE estimate", "--color", GOLD, "--line-style", "time-dash",
           "--axis", "radius r", "--axis", "radial velocity r-dot", "--title", title + "states (plant on the ideal drive)",
           "--height", "1000", "--output", f"{OUT}/{tag}_states.png")
    viewer("time-series", *common, "--columns", "1", "--curve", "err=KLQE estimate minus actual state", "--color", BLUE,
           "--axis", "r estimate error", "--axis", "r-dot estimate error", "--title", title,
           "--height", "1000", "--output", f"{OUT}/{tag}_estimate_error.png")
    viewer("time-series", *common, "--curve", "ideal=ideal drive (true-state KLQR, applied)", "--color", RED,
           "--line-width", "2.7", "--curve", "klqr=KLQR force on the echo-only estimate (not applied)", "--color", BLUE,
           "--axis", "drive p/p0", "--title", title + "KLQR force vs ideal drive",
           "--height", "700", "--output", f"{OUT}/{tag}_force.png")
    viewer("time-series", *common, "--curve", "force_err=KLQR force minus ideal drive", "--color", GOLD,
           "--axis", "force error p/p0", "--title", title + "force error",
           "--height", "600", "--output", f"{OUT}/{tag}_force_error.png")


def startup(label, archive, stop=20):
    """Start-up set: states in time, phase portrait, control signals (applied ideal drive vs KLQR force)."""
    out = f"{OUT}/startup"
    target, estimator = label.split("-", 2)[1], label.split("-", 2)[2]
    title = f"{target}, {NAMES[estimator]}, start-up: "
    z = np.load(archive)
    te, t, x, tg = z["estimator_time"], z["time"], z["values"], z["control_target_values"]
    on = lambda a: np.column_stack([np.interp(te, t, a[:, k]) for k in range(2)])
    xs, xt, xh = on(x), on(tg), z["state_estimate"]
    data = f"{out}/{label}_startup.npz"
    np.savez(data, time=te, x=xs, target=xt, xhat=xh, r=xs[:, 0], rdot=xs[:, 1], r_t=xt[:, 0], rdot_t=xt[:, 1],
             r_hat=xh[:, 0], rdot_hat=xh[:, 1], ideal=z["applied_input"].ravel(), klqr=z["commanded_input"].ravel())
    common = ["--source", data, "--start", "0", "--stop", str(stop), "--width", "1400"]
    viewer("time-series", *common, "--columns", "1", "--time-label", "tau",
           "--curve", "x=actual", "--color", RED, "--line-width", "2.7",
           "--curve", "target=target", "--color", BLUE, "--line-style", "dash",
           "--curve", "xhat=echo-only KLQE estimate", "--color", GOLD, "--line-style", "time-dash",
           "--axis", "radius r", "--axis", "radial velocity r-dot", "--title", title + "states",
           "--height", "1000", "--output", f"{out}/{label}_1_states.png")
    viewer("time-series", *common, "--curve", "r,rdot=actual", "--color", RED, "--line-width", "2.7",
           "--curve", "r_t,rdot_t=target", "--color", BLUE, "--line-style", "dash",
           "--curve", "r_hat,rdot_hat=echo-only KLQE estimate", "--color", GOLD, "--line-style", "dash",
           "--x-label", "radius r", "--value-label", "radial velocity r-dot", "--axis", "r-dot versus r",
           "--title", title + "phase portrait",
           "--height", "1000", "--output", f"{out}/{label}_2_phase.png")
    viewer("time-series", *common, "--time-label", "tau",
           "--curve", "ideal=ideal drive (true-state KLQR, applied)", "--color", RED, "--line-width", "2.7",
           "--curve", "klqr=KLQR force on the echo-only estimate (not applied)", "--color", BLUE,
           "--axis", "drive p/p0", "--title", title + "control signals",
           "--height", "700", "--output", f"{out}/{label}_3_control.png")


if __name__ == "__main__":
    os.makedirs(OUT + "/startup", exist_ok=True)
    scores = json.load(open(os.path.join(HERE, "shadow-scores.json")))
    only = [a for a in sys.argv[1:] if not a.startswith("--")]
    for label, s in scores.items():
        if only and label.split("-", 2)[2] not in only: continue
        if "--startup" in sys.argv:
            startup(label, s["archive"])
        else:
            for window in ((0, 100), (85, 100)):
                figures(label, s["archive"], *window)
        print(label, "done", flush=True)
