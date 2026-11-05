"""Present native experiment results using the existing Viewer PlotCanvas.

No fitting, projection, control, differentiation or simulation is performed
here. All numerical results come from the compiled echo-input experiment.
"""
from __future__ import annotations

import html
import json
from pathlib import Path

import numpy as np

from plot_rpe_echo_jet_pipeline_evidence import figure, render, trace

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "archives/experiments/echo-input-regression-20261006-v1"
REPORT = ROOT / "docs/reports/ECHO_HISTORY_BPSI_EXPERIMENT_2026-10-06.html"


def load(name: str, root: Path = DATA) -> dict:
    """Read an immutable native JSON result."""
    return json.loads((root / name).read_text(encoding="utf-8"))


def percent(value: float) -> str:
    """Format a native relative-error statistic."""
    return f"{100 * value:,.2f}%"


def row_table(headers: list[str], rows: list[list[str]]) -> str:
    """Produce a plain accessible result table."""
    return "<table><thead><tr>" + "".join(f"<th>{html.escape(h)}</th>" for h in headers) + (
        "</tr></thead><tbody>" + "".join(
            "<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in row) + "</tr>" for row in rows
        ) + "</tbody></table>"
    )


def plot_case(case: dict, folder: Path, target: str, tag: str) -> str:
    """Arrange archived signals with the established native Viewer schema."""
    measured = np.load(folder / f"{case['id']}.measurements.npz", allow_pickle=False)
    truth = np.load(folder / f"{case['id']}.truth.npz", allow_pickle=False)
    reference = np.load(DATA / f"target-{target}.measurements.npz", allow_pickle=False)
    reference_truth = np.load(DATA / f"target-{target}.truth.npz", allow_pickle=False)
    time = measured["time"]
    late = np.flatnonzero(time >= 54)[::4]
    full = np.arange(0, len(time), 20)
    state, desired = truth["states"], reference_truth["states"]

    def paired(title: str, label: str, y: np.ndarray, ref: np.ndarray, indices=late) -> dict:
        return figure("time_series", title, r"\tau", label, [
            trace("Target", time[indices], ref[indices], color="#606872", line_style="dash", width=1.6),
            trace("KLQR", time[indices], y[indices], color="#007f8b", width=1.8),
        ])

    figures = [
        paired("Applied input vs known sine driver", "u = p / p0", measured["input"], reference["input"]),
        paired("Measured intrinsic echo", "e", measured["echo"], reference["echo"]),
        paired("Radius (post-hoc simulator truth)", "r", state[:, 0], desired[:, 0]),
        paired("Radial velocity (post-hoc simulator truth)", r"\dot{r}", state[:, 1], desired[:, 1]),
        figure("phase_portrait", "Physical phase portrait (post-hoc only)", "r", r"\dot{r}", [
            trace("Target", desired[late, 0], desired[late, 1], color="#606872", line_style="dash"),
            trace("KLQR", state[late, 0], state[late, 1], color="#007f8b"),
        ]),
        paired("Radius over the complete run", "r", state[:, 0], desired[:, 0], full),
    ]
    name = f"echo-input-{tag}-{target}"
    render(name=name, label=f"Echo-history KLQR | {tag} | {target}",
           page_label="Exploratory best late echo error; not a success claim",
           source=folder / f"{case['id']}.measurements.npz", figures=figures, columns=2, height=1350)
    return f"../../viewer_outputs/{name}.png"


def main() -> None:
    """Build the report from finished native computation outputs."""
    broad, local, sweep = (load(p) for p in ("results.json", "refinement.json", "eigenfunction-sweep.json"))
    short_root = DATA / "short-history-0p05"
    short = load("eigenfunction-sweep.json", short_root)
    groups = [("Original", broad, DATA), ("Frequency-local", local, DATA),
              ("Long-window refit", sweep, DATA), ("Short-window refit", short, short_root)]
    controls = [(label, c, folder) for label, result, folder in groups for c in result["controls"]]
    winners = []
    panels = []
    for tag, choices in [("long-window", controls[:-len(short['controls'])]),
                         ("short-window", [v for v in controls if v[0] == "Short-window refit"])]:
        for target in ("mild", "nonlinear"):
            candidates = [v for v in choices if target in v[1]["id"] and "late" in v[1]]
            label, case, folder = min(candidates, key=lambda v: v[1]["late"]["echo_nrmse"])
            metrics = case["late"]
            winners.append([tag, target, case["id"], percent(metrics["echo_nrmse"]),
                            percent(metrics["driver_nrmse"]), percent(metrics["saturation_fraction"])])
            png = plot_case(case, folder, target, tag)
            panels.append(f'<section><h2>{html.escape(tag)} / {target}</h2><p class="mono">'
                          f'{html.escape(case["id"])}</p><p>Late echo error: {percent(metrics["echo_nrmse"])}; '
                          f'driver error: {percent(metrics["driver_nrmse"])}. Targets are grey dashed. '
                          f'The known driver is displayed only for evaluation.</p><a href="{png}">'
                          f'<img src="{png}" alt="Native Viewer control and response plots for {target}"></a></section>')
    candidate_rows = [["Original, 4.3215 τ", "2", "5", "50", "−0.060245 + 6.282949 i",
                       percent(sweep["original_free_validation"]["one_tau"]),
                       percent(sweep["original_free_validation"]["four_tau"])]]
    for label, data in [("Refit, 4.3215 τ", sweep), ("Refit, 0.05 τ", short)]:
        for choice in data["selected"]:
            candidate_rows.append([label, str(int(choice["rank"])), str(int(choice["degree"])),
                                   str(choice["fit_start"]), f'{choice["lambda_re"]:.6f} + {choice["lambda_im"]:.6f} i',
                                   percent(choice["heldout_free"]["one_tau"]), percent(choice["heldout_free"]["four_tau"])])
    all_rows = []
    for label, case, _ in controls:
        if "failure" in case:
            all_rows.append([label, case["id"], "FAILED", case["failure"], "—"])
        else:
            m = case["late"]
            all_rows.append([label, case["id"], percent(m["echo_nrmse"]), percent(m["driver_nrmse"]),
                             percent(m["saturation_fraction"])])
    map_rows = []
    for label, data in [("Long-window refit", sweep), ("Short-window refit", short)]:
        for m in data["maps"]:
            map_rows.append([label, m["id"], percent(m["fit_relative_residual"]),
                             percent(m["target_relative_residual"]), str(m["rank"]), f'{m["condition"]:.3g}'])
    failed = sum("failure" in case for _, case, _ in controls)
    storage_failed = sum("write NPZ" in case.get("failure", "") for _, case, _ in controls)
    report = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Echo-history Bψ regression and KLQR experiments</title><style>
body{{font:17px/1.6 system-ui,sans-serif;color:#253344;background:#f3f6f8;margin:0}}
main{{max-width:1280px;margin:auto;padding:40px 28px}}h1{{font-size:36px;line-height:1.2}}h2{{margin-top:0}}
section{{background:white;border:1px solid #d9e2e9;border-radius:10px;margin:24px 0;padding:26px;overflow:auto}}
table{{border-collapse:collapse;font-size:14px;width:100%}}td,th{{padding:10px;border-bottom:1px solid #dce4ea;text-align:left}}
th{{background:#e8f2f5}}.mono,code{{font:13px/1.5 ui-monospace,monospace;overflow-wrap:anywhere}}
.note{{border-left:5px solid #b77515;padding:16px 22px;background:#fff8e9}}img{{width:100%;height:auto}}
a{{color:#066778}}.flow{{font-size:20px;padding:18px;background:#e8f2f5}}summary{{cursor:pointer;font-weight:600}}
</style><main><p>DynSysLab · echo-only · 6 October 2026</p>
<h1>Echo-history Bψ regression and KLQR control</h1>
<p>Compiled C++ experiments, with the existing solver, eigenfunction fitting, numerical regression, LQR and Viewer components.</p>
<p class="note"><strong>Outcome: partial success for the mild target, no success for the nonlinear target.</strong> The best short-history, two-latent trial achieved 8.88% late-time echo error and 21.61% driver error. This is approximate tracking, not accurate recovery of the known driver, and does not establish a generally reliable control map.</p>
<section><h2>What was tested</h2><p>76 forced calibration recordings (60 broad-band combinations and 16 frequency-local runs),
9 held-out forced recordings, 8 new unforced training recordings and 6 unforced validation recordings, and 2 sine-driven targets.
The existing free-trained ψ was compared with free-only polynomial refits in two and four latent coordinates.
There were {len(controls)} KLQR attempts; {failed} did not produce a complete saved result, including {storage_failed} NPZ write failures when the disk filled. No failed attempt is omitted below. Storage failures are incomplete experiments, not evidence of controller instability.</p>
<p class="note">A low one-step prediction error is not sufficient evidence of control. Judge echo tracking, recovery of the known
driver, long-time behavior and saturation together. Best runs below are selected retrospectively from a parameter sweep;
they are exploratory comparisons, not untouched confirmation tests.</p>
{row_table(["History", "Target", "Selected trial", "Echo NRMSE", "Driver NRMSE", "Saturated"], winners)}</section>
<section><h2>The information boundary</h2><p class="flow">Echo history → ξ → ψ(ξ), B̂ψ(ξ) → existing KLQR → applied pressure</p>
<p>The simulator alone knows the RPE equations and physical state. It writes separate measurement-only archives containing
time, intrinsic echo and applied input. The regression reader rejects state-bearing archives. The history transform,
eigenfunction fitting, Bψ fitting and controller receive no r, ṙ, state Jacobians or equations of motion.
The target is supplied as an echo recording transformed through the same basis and polynomial. Its recorded input is withheld
from the control calculation and used only for the evaluation metrics and overlays.</p>
<p>Physical states in the plots are post-hoc simulator truth, not controller inputs. Driver scattering, direct acoustic
feedthrough, dipoles and higher multipoles are intentionally excluded, as agreed. No sensor noise, actuator dynamics or laboratory
uncertainty is included.</p></section>
<section><h2>Bubble, data and targets</h2><p>The source KROM fixes the air–water RPE parameters: equilibrium radius 10 μm,
ambient pressure 100 kPa, viscosity 0.001 Pa·s, density 998 kg/m³, surface tension 0.0725 N/m,
vapour pressure 2330 Pa and κ = 4/3. No parameter was retuned.</p>
<p>The broad calibration grid uses amplitudes 0.003, 0.01, 0.03, 0.08, 0.2 and 0.4 in p/p₀, frequencies
0.65, 0.85, 1.0, 1.15 and 1.35 cycles/τ, and phases 0 and π/2. Each drive lasts 16 cycles, with a
6τ unforced prefix and 6τ free tail. Nine held-out runs use new amplitudes, frequencies, phase 0.37 and initial radius 1.005.
The local calibration uses 32 cycles at each target frequency, different amplitudes and phases 0.4 and 1.1;
its fits exclude the first 6τ of the driven segment. It is frequency-local calibration, not frequency-blind generalization.</p>
<p>The mild target is 1.2 kPa at 0.93 cycles/τ (about 309 kHz), and the nonlinear target is 30 kPa at 0.83 cycles/τ
(about 276 kHz). Both start at equilibrium and run to 60τ, with driving from 6τ onward. The nonlinear reference has
r ≈ 0.573–1.478 and ṙ ≈ −3.900–3.849, a substantial nonlinear oscillation. The mild reference has r ≈ 0.967–1.033.</p>
<p>Sensor sampling is Δτ = 0.0005. Applied sine inputs are sampled and zero-order held for 0.01τ, about 100 actuator
samples per natural period. This is a numerical sampled-sine experiment, not an exact continuously varying actuator.
The final 30τ supplies the long-time error metrics. Inputs are clipped to ±0.6 p/p₀.</p></section>
<section><h2>How ψ was fitted</h2><p>The original basis has 8,644 samples spanning 4.3215τ; the original polynomial is degree five
in ξ₁–ξ₂ and was fitted over τ = 50–80. For the first refits its U and singular values were frozen.
New free training IC radii are 0.65, 0.8, 0.9, 1.05, 1.2, 1.4, 1.7 and 2.0, all with zero initial velocity.
Six separate free validation radii are 0.7, 0.85, 1.1, 1.3, 1.55 and 1.85.</p>
<p>The existing pooled polynomial EDMD routine fits snapshot pairs 0.01τ apart, with relative ridge 10⁻¹⁰.
Candidates use starts 5, 15, 35 or 50τ and degrees 1, 3 or 5 for two latents; four latents use degrees 1 or 3.
Candidate modes require −0.5 ≤ Re λ ≤ −0.001 and 4 ≤ Im λ ≤ 7. Each is stored in the normal polynomial/SINDy
eigenfunction contract; these coefficients come from EDMD, not sparse thresholded generator-PDE fitting.</p>
<p>Selection minimizes the sum of the mean held-out 1τ and 4τ semigroup relative-L₂ errors, always evaluated over τ = 5–80,
including early motion even when fitting starts later. This is validation-based selection, not an independent test of the selected ψ.
The additional 0.05τ basis uses 101 samples and is fitted only to the eight unforced records by the existing pooled Hankel routine.
Neither its basis nor any eigenfunction is fitted to forced data.</p>
{row_table(["Model", "Latents", "Degree", "Fit start", "λ", "Free 1τ error", "Free 4τ error"], candidate_rows)}</section>
<section><h2>How Bψ was regressed and used</h2>
<p>At each actuator tick the recorded echo history gives zₖ = ψ(ξₖ). With μ = exp(λh), the regression target is
qₖ = zₖ₊₁ − μzₖ. Least squares fits qₖ ≈ B̂d(ξₖ)uₖ. No numerical echo derivative or physical-state derivative is required.
Constant maps and the existing conjugate-complete bilinear family are compared:</p>
<p class="flow">B̂d(ξ) = c₀ + c₁ Re ψ(ξ) + c₂ Im ψ(ξ).</p>
<p>This is a restricted function of the latent coordinates, through ψ; it is not an unrestricted polynomial in every ξ component.
The design columns are normalized before the existing least-squares solve, with ridge 10⁻⁸. Whole runs are kept separate;
no pair crosses a recording boundary. Small-amplitude and pooled fits are distinct comparisons.</p>
<p>For the existing continuous-time KLQR solver, the fitted discrete map is converted as B̂ψ = B̂d / G, where
G = (exp(λh) − 1)/λ. For a varying B this is a frozen-map, zero-order-hold approximation over h, not an exact instantaneous Jacobian.
The maps are saved using native ConstantControlMap or BilinearControlMap array conventions.</p>
<p>The controller evaluates the current polynomial and fitted map at each tick, calls the existing
solveContinuousLQR routine, and applies u = −K[Re(z − z*), Im(z − z*)]ᵀ. There is no feedforward term and no use of the known target
driver. Q = I, and R = α·RMS(z*)²/0.1²; α is swept and recorded per trial. Stronger feedback is not assumed to be better.</p>
{row_table(["Study", "Map", "Fit residual / baseline", "Target residual / baseline", "Rank", "Scaled condition"], map_rows)}
<p>The baseline is q itself (zero fitted input response), not the full next-step z signal. A ratio near 100% means that
the map explains little of the forcing-related residual; a ratio above 100% is worse than that baseline.</p></section>
{''.join(panels)}
<section><h2>Interpretation and limits</h2><p>The original nonlinear target is outside the old two-latent training bounding box
for 93.64% of available history windows. The mild target is entirely inside that box and still exhibits substantial errors.
Thus polynomial extrapolation matters, but cannot alone explain every failure. Being inside a bounding box also does not prove
that a forced history is covered by the manifold of unforced histories.</p>
<p>A history projection can predict unforced evolution better after a ψ refit and still fail to provide a closed forced model
ż = λz + Bψ(ξ)u. Past actuation can change the history while producing similar retained coordinates. The comparisons are consistent
with that issue, but do not independently prove it: restricted map expressiveness, the unforced eigenfunction approximation,
control gain/sampling and the continuous-time conversion are additional limitations. These results must not be interpreted as
a proof that echo-only control is impossible.</p>
<p>The next scientifically discriminating test is matched-history input perturbations and conditional residual analysis: determine
whether different histories with similar ξ require different B values. If so, enrich the observation coordinates or represent
input memory before demanding more from the same regression. If not, a richer Bψ(ξ) regression may suffice.
No production controller replacement is justified by these experiments.</p></section>
<section><details><summary>All {len(controls)} KLQR attempts, including failures</summary>
{row_table(["Study", "Trial", "Late echo error", "Late driver error / failure", "Saturated"], all_rows)}</details></section>
<section><h2>Reproducibility and local assets</h2><p>All scientific computation is C++. Python only assembles these native result
tables and Viewer payloads. No existing completed archive was overwritten; large datasets and plots remain local and Git-ignored.</p>
<p>Sources: <code>native/tools/echo_input/</code>, <code>native/src/HistoryInputRegression.cpp</code> and
<code>native/include/dynsyslab/core/HistoryInputRegression.hpp</code>. The existing causal projector gained append/projectLatest
methods so acquisition can run faster than control. Production LQR.cpp and the existing physical-state regression pipeline were not modified.</p>
<p>Native executable: <code>build/native-suite-fast/bin/dynsyslab_echo_input_experiment.exe</code>.
Supply the frozen KROM path and a fresh output directory for the base run, then use <code>--refine</code>,
<code>--eigenfunctions</code>, and <code>--short-history</code> in that order. Completed results are immutable.</p>
<p><a href="../../archives/experiments/echo-input-regression-20261006-v1/results.json">Base results</a> ·
<a href="../../archives/experiments/echo-input-regression-20261006-v1/refinement.json">Longer-burst results</a> ·
<a href="../../archives/experiments/echo-input-regression-20261006-v1/eigenfunction-sweep.json">Free-only ψ sweep and control</a> ·
<a href="../../archives/experiments/echo-input-regression-20261006-v1/short-history-0p05/eigenfunction-sweep.json">Short-history comparison</a>.</p>
<p>Verification: 99/99 native tests passed, including history projection, canonical map conversion, synthetic map recovery and existing
KLQR regression tests. This establishes software checks, not scientific success or laboratory readiness.</p></section></main></html>'''
    REPORT.write_text(report, encoding="utf-8")
    print(REPORT)


if __name__ == "__main__":
    main()
