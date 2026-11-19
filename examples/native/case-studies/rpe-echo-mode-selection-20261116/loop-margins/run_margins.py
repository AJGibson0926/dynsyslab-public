"""Run the loop-margin jobs (linear_systems on each echo-only K-ROM) and print the robustness table.

Usage: python run_margins.py            (runs every *.json here, writes margins.json)
"""
import glob, json, os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
CLI = os.path.join(ROOT, "build", "linear-systems", "bin", "dynsyslab-cli.exe")
FIELDS = ["minimum_return_difference", "disk_gain_lower", "disk_gain_upper", "disk_phase_margin_deg",
          "gain_margin_upper", "phase_margin_deg", "sensitivity_peak", "closed_loop_stable"]


def run(job):
    out = subprocess.run([CLI, "run", job, "--archive-root", os.path.join(ROOT, "archives")],
                         capture_output=True, text=True, cwd=ROOT, check=True)
    return json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]


if __name__ == "__main__":
    import numpy as np
    table = {}
    for job in sorted(glob.glob(os.path.join(HERE, "*.json"))):
        if job.endswith("margins.json"): continue
        archive = run(job)
        z = np.load(archive); meta = json.loads(str(z["metadata_json"]))
        zeros = z["transfer_function_zeros__y0_u0"]
        entry = {"archive": archive, "zeros": [[float(a), float(b)] for a, b in zeros],
                 "transfer_function": meta["transfer_functions"][0],
                 "loops": {n: {k: meta["loop_robustness"][n].get(k) for k in FIELDS + ["ltr_weight"]}
                           for n in meta["loop_names"]}}
        name = os.path.splitext(os.path.basename(job))[0]
        table[name] = entry
        print(f"\n== {name}  zeros {np.round(zeros[:, 0] + 1j * zeros[:, 1], 3)}")
        for n, r in entry["loops"].items():
            f = lambda k: "inf" if r[k] is None else (f"{r[k]:.3g}" if isinstance(r[k], float) else str(r[k]))
            print(f"  {n:14s} alpha {f('minimum_return_difference'):>6} disk GM [{f('disk_gain_lower')}, {f('disk_gain_upper')}]"
                  f" disk PM {f('disk_phase_margin_deg')}  GM+ {f('gain_margin_upper')} PM {f('phase_margin_deg')}"
                  f" Ms {f('sensitivity_peak')} stable {f('closed_loop_stable')}" + (f"  q2={r['ltr_weight']:g}" if r.get("ltr_weight") else ""))
    json.dump(table, open(os.path.join(HERE, "margins.json"), "w"), indent=2)
