"""Observability matrix obsv(M, C) of each K-ROM in discrete modal coordinates (one complex
coordinate per eigenfunction, M = diag(exp(lambda dt))) and its SVD, drawn with the Viewer CLI:

  modal_<name>_singular_values.png   sigma_k against the integer index k (stems)
  modal_<name>_directions.png        one tile per observability direction n = column V(:, n): signed real
                                     parts over the eigenfunction index, titled with sigma_n
  modal_<name>_rows.png              one tile per row V(n, :): eigenfunction n's weight in every direction
                                     (the tiles playground.m draws), titled with sigma_n as there

K-ROMs: the psi(e, e') K-ROM the echo-only KLQE used on each target (one conjugate pair; its whole
lineage has only that pair) and the 15-mode H-DMD parent of the smooth H-DMD K-ROM.
"""
import json, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 5))
BIN = os.path.join(ROOT, "build", "linear-systems", "bin")
OUT = os.path.join(ROOT, "viewer_outputs", "echo-observability-20261118")
D = ROOT.replace("\\", "/") + "/archives/derived/"
KROMS = {
    "ee-mild": ("psi(e, e') K-ROM, mild", D + "koopman-mode-selection/experiment-ee-ic0p37-split-krom__mode-selection-krom__20261117T002607Z.npz"),
    "ee-strong": ("psi(e, e') K-ROM, strong", D + "koopman-mode-selection/experiment-ee-ic0p37-split-krom__mode-selection-krom__20261117T002614Z.npz"),
    "ee-hard": ("psi(e, e') K-ROM, hard", D + "koopman-mode-selection/experiment-ee-ic0p37-split-krom__mode-selection-krom__20261117T002622Z.npz"),
    "hdmd-parent": ("H-DMD parent K-ROM (15 modes)", D + "koopman-krom-extract/rpe-viscous-echo-long-validation-five-ics-trajectory-0__~32729aacd2d4~T194836Z__krom-extract__20261116T194904Z.npz"),
}
BLUE = "#004DB8"


def run(job):
    out = subprocess.run([os.path.join(BIN, "dynsyslab-cli.exe"), "run", job, "--archive-root", os.path.join(ROOT, "archives")],
                         capture_output=True, text=True, cwd=ROOT, check=True)
    return json.loads(out.stdout[out.stdout.index('{"archives"'):])["archives"][0]


def viewer(*args):
    subprocess.run([os.path.join(BIN, "DynSysLabViewerCli.exe"), "time-series", *args], check=True,
                   capture_output=True, text=True, cwd=ROOT)


if __name__ == "__main__":
    summary = {}
    for name, (title, krom) in KROMS.items():
        job = os.path.join(HERE, f"modal-{name}.json")
        json.dump({"schema": "dynsyslab.compute-job/v1", "kind": "postprocess", "label": f"modal obsv / {name}",
                   "description": "Observability matrix [C; C M; ...; C M^(n-1)] of the K-ROM in discrete modal coordinates "
                                  "(M = diag(exp(lambda dt)), one complex coordinate per eigenfunction) with the echo output, and its SVD.",
                   "sources": [krom], "methods": [{"name": "linear_analysis", "enabled": True, "settings": {
                       "measurement": "far_field_acoustic_echo", "analysis_method": "modal-matrix",
                       **({"sample_time": 0.001} if name == "hdmd-parent" else {})}}],  # H-DMD records delta_t only in metadata
                   "archive": {"enabled": True}}, open(job, "w"), indent=2)
        z = np.load(run(job))
        s = z["modal_singular_values"]; mag = z["modal_direction_magnitudes"]  # modes x k
        n, k = mag.shape
        meta = json.loads(str(z["metadata_json"]))
        summary[name] = {"singular_values": s.tolist(), "modes": int(n), "sample_time": meta["sample_time"]}
        src = os.path.join(OUT, f"modal_{name}.npz")
        np.savez(src, index=np.arange(1.0, len(s) + 1), sv=s, mode=np.arange(1.0, n + 1), directions=mag)
        viewer("--source", src, "--time-key", "index", "--time-label", "singular value index k", "--integer-x-ticks",
               "--curve", "sv=singular value", "--color", BLUE, "--line-style", "stem",
               "--scale", "log", "--x-range", f"0.5:{len(s) + 0.5}", "--axis", "singular value sigma_k",
               "--title", f"{title}: singular values of obsv(M, C), echo output",
               "--width", "1100", "--height", "700", "--output", src[:-4] + "_singular_values.png")
        # Singular vectors are defined up to a phase each: rotate column k so its largest entry is
        # real and positive (a reproducible convention), then plot signed real parts.
        d = z["modal_directions"]; V = d[..., 0] + 1j * d[..., 1]  # modes x k, column k = direction k
        for c in range(V.shape[1]):
            j = np.argmax(np.abs(V[:, c])); V[:, c] *= np.conj(V[j, c]) / abs(V[j, c])
        shown = min(k, 9)
        cols = 1 if shown <= 2 else 3

        def tiles(clock, label, data, axis_names, out, xlabel):
            np.savez(src[:-4] + f"_{out}.npz", **{label: clock}, values=data)
            args = ["--source", src[:-4] + f"_{out}.npz", "--time-key", label, "--time-label", xlabel, "--integer-x-ticks",
                    "--curve", "values=Re", "--color", BLUE, "--line-style", "stem",
                    "--x-range", f"0.5:{len(clock) + 0.5}", "--value-label", "real part"]
            for a in axis_names: args += ["--axis", a]
            if cols > 1: args += ["--columns", str(cols)]
            viewer(*args, "--width", "1500" if cols > 1 else "1100", "--height", "1200" if cols > 1 else "900",
                   "--output", src[:-4] + f"_{out}.png")

        prefix = "" if cols > 1 else f"{title}, "
        # columns V(:, n): the n-th observability direction, over the eigenfunctions
        tiles(np.arange(1.0, n + 1), "mode", V[:, :shown].real,
              [f"{prefix}direction {c + 1} = V(:, {c + 1}): sigma = {s[c]:.3g}" for c in range(shown)],
              "directions", "eigenfunction index")
        # rows V(n, :): eigenfunction n's weight in each direction (the playground.m tiles)
        tiles(np.arange(1.0, k + 1), "index", V[:shown, :].T.real,
              [f"{prefix}row {r + 1} = V({r + 1}, :), eigenfunction {r + 1} (title sigma = {s[r]:.3g})" for r in range(shown)],
              "rows", "singular value index k")
        print(name, n, np.array2string(s, precision=3))
    json.dump(summary, open(os.path.join(HERE, "modal-obsv.json"), "w"), indent=2)
