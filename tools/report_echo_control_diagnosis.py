"""Presentation only: arrange native diagnostic arrays with the existing Viewer."""
from __future__ import annotations
import html
import json
from pathlib import Path
import numpy as np
from plot_rpe_echo_jet_pipeline_evidence import figure, render, trace
from report_echo_input_experiment import row_table, percent

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "archives/experiments/echo-control-diagnosis-20261007-v2"
FIRST = ROOT / "archives/experiments/echo-control-diagnosis-20261007-v1"
REPORT = ROOT / "docs/reports/ECHO_HISTORY_CONTROL_DIAGNOSIS_2026-10-07.html"


def read(path: Path) -> dict:
    """Read the saved JSON evidence used to assemble this local report."""
    return json.loads(path.read_text(encoding="utf-8-sig"))


def panels(kind: str, rank: int = 2) -> str:
    """Render saved diagnostics for either retained latent dimension; no new solves."""
    if kind not in ("causal", "inverse") or rank not in (2, 4):
        raise ValueError("diagnostic panels require causal/inverse and two/four latents")
    figures = []
    for target in ("mild", "nonlinear"):
        for window in ("long", "short"):
            tag = f"{window}-xi{rank}-{target}-local"
            a = np.load(DATA / f"{tag}-{'paired' if kind == 'causal' else 'inverse'}.npz")
            if kind == "causal":
                x = a["time"]
                traces = [trace("Existing fitted map", x, a["fitted_B_magnitude"], color="#007f8b"),
                          trace("Matched-history sensitivity", x, a["causal_B_magnitude"], color="#b66b16", marker=True)]
                f = figure("time", f"{window} history / {target} / {rank} latents: input sensitivity", r"\tau", "|B_d|", traces)
                f["value_scale"] = "log"
            else:
                keep = (a["time"] >= 54) & (a["time"] <= 58)
                x = a["time"][keep]
                f = figure("time", f"{window} history / {target} / {rank} latents: inverse input", r"\tau", "u = p / p0", [
                    trace("Known driver (score only)", x, a["known_input"][keep], color="#626a74", line_style="dash"),
                    trace("From target echo + fitted map", x, a["inverse_input"][keep], color="#007f8b")])
                f["value_scale"] = "linear"
            figures.append(f)
    name = f"echo-control-diagnosis-{kind}-20261007" + ("-xi4" if rank == 4 else "")
    render(name=name, label=f"Echo-only control diagnosis / {kind} / {rank} latents",
           page_label=f"{rank} latents; frozen models; diagnostic only; no controller changes",
           source=DATA / "diagnosis.json", figures=figures, columns=2, height=1100)
    return f'<a href="../../viewer_outputs/{name}.png"><img src="../../viewer_outputs/{name}.png" alt="Native Viewer {kind} diagnostic plots"></a>'


def main() -> None:
    """Assemble the saved local evidence and native Viewer plots into the dated HTML report."""
    result = read(DATA / "diagnosis.json")
    resumed = read(FIRST / "diagnosis.json")["resumed_controls"]
    capacity = read(DATA / "paired-capacity/diagnosis.json")["paired_capacity"]
    causal_plot, inverse_plot = panels("causal"), panels("inverse")
    paired_rows, inverse_rows, free_rows, capacity_rows, release_rows = [], [], [], [], []
    for r in result["results"]:
        if not r["id"].endswith("local"):
            continue
        m, v = r["matched_history"], r["inverse"]
        paired_rows.append([r["id"], f'{m["zero_input_drift_over_actual_input_effect"]:,.2f}×',
                            f'{m["fitted_B_relative_error_to_causal_B"]:,.2f}',
                            percent(m["local_input_affine_residual"])])
        inverse_rows.append([r["id"], percent(v["late_inverse_driver_relative_error"]),
                             f'{v["late_inverse_driver_cosine"]:.3f}',
                             percent(v["late_perpendicular_residual_fraction"]),
                             percent(v["late_known_driver_model_residual"])])
    for r in result["free_validation"]:
        for interval in r["intervals"]:
            free_rows.append([r["case"], f'{interval["from"]:g}–{interval["to"]:g}',
                              percent(interval["pair_relative_error"]),
                              percent(interval["one_tau_relative_error"]),
                              percent(interval["four_tau_relative_error"])])
    for r in capacity:
        m = r["metrics"]
        capacity_rows.append([r["case"], r["target"], percent(m["paired_bilinear_training_relative_error"]),
                              percent(m["paired_bilinear_heldout_phase_relative_error"]),
                              f'{m["paired_bilinear_condition"]:.2f}', percent(m["paired_inverse_driver_relative_error"])])
    for r in result["release_checks"]:
        values = r["intervals"]
        release_rows.append([r["case"], r["target"], percent(values[3]["pair_relative_error"]),
                             percent(values[4]["pair_relative_error"]), percent(values[5]["pair_relative_error"])])
    resumed_rows = [[str(r["input_penalty_multiplier"]), percent(r["late"]["echo_nrmse"]),
                    percent(r["late"]["driver_nrmse"]), percent(r["late"]["saturation_fraction"])] for r in resumed]
    report = f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>Echo-history control: diagnosis and limits</title>
<style>body{{font:17px/1.65 system-ui,sans-serif;color:#253344;background:#f3f6f8;margin:0}}
main{{max-width:1280px;margin:auto;padding:40px 28px}}h1{{font-size:38px;line-height:1.2}}h2{{margin-top:0}}
section{{background:white;border:1px solid #d9e2e9;border-radius:10px;margin:24px 0;padding:28px;overflow:auto}}
table{{border-collapse:collapse;font-size:14px;width:100%}}td,th{{padding:10px;border-bottom:1px solid #dce4ea;text-align:left}}
th{{background:#e8f2f5}}code,pre{{font:14px/1.65 ui-monospace,monospace;overflow-wrap:anywhere;white-space:pre-wrap}}
.note{{border-left:5px solid #b77515;padding:16px 22px;background:#fff8e9}}.equation{{background:#eef5f7;padding:16px;overflow:auto;font:19px/1.7 Georgia,serif}}
img{{width:100%;height:auto}}a{{color:#066778}}summary{{cursor:pointer;font-weight:600}}li{{margin:10px 0}}
@media print{{body{{background:white}}section{{break-inside:avoid;border:0}}}}</style></head><body><main>
<p>DynSysLab · echo-only · 7 October 2026 · diagnosis-first investigation</p>
<h1>Why echo-history KLQR is not yet reliable</h1>
<p class="note"><strong>Recommendation: stop broad gain sweeps and validate the lifted dynamics before asking KLQR to control them.</strong>
The main obstacle is not simply a bad numerical solve or an insufficient number of sine bursts. The present fit confounds
autonomous/history drift with current-input response. A second, independent limitation is that the controller has no nominal
forcing term for persistent periodic tracking. The current frozen models cannot support a claim of reliable nonlinear control.</p>
<p>The primary pipeline was left unchanged. One isolated C++ audit executable reused existing projection, polynomial,
integration, least-squares, archive and LQR components. A small presentation script assembled these results using native Viewer.
No new SVD, eigenfunction training, controller architecture or primary-pipeline edits were made. An isolated diagnostic refit
used the existing bilinear regression routine; it was not installed or used for closed-loop control.</p>
<nav><a href="#scope">Scope</a> · <a href="#causal">Causal test</a> · <a href="#capacity">Map capacity</a> ·
<a href="#history">History problem</a> · <a href="#free">Eigenfunctions</a> · <a href="#tracking">Tracking law</a> ·
<a href="#recommendations">Recommendations</a> · <a href="#reproduce">Evidence</a></nav>

<section id="scope"><h2>1. What was kept fixed, and what was tested</h2>
<p>The system remains the viscous air-in-water RPE bubble: equilibrium radius 10 μm, ambient pressure 100 kPa,
viscosity 0.001 Pa·s, density 998 kg/m³, surface tension 0.0725 N/m, polytropic exponent 4/3 and vapor pressure 2330 Pa.
The synthetic sensor is the intrinsic, zero-input-convention far-field echo. Driver scattering and direct acoustic feedthrough
remain intentionally excluded. These are noiseless simulation tests, not laboratory validation.</p>
<p>Sensor spacing is 0.0005 τ; inputs are held for h = 0.01 τ. The old long history contains 8644 samples spanning 4.3215 τ;
the short history contains 101 samples spanning 0.05 τ. Each stored U has four columns. The tested polynomials use either the
first two or all four latent coordinates, but each model produces <strong>one complex eigenfunction</strong>—two real feedback coordinates.</p>
{row_table(['Model','Polynomial degree','Free training interval','Eigenvalue λ'],[
['Long / ξ₁–ξ₂','1','35–80 τ','−0.061506 + 6.261833 i'],['Long / ξ₁–ξ₄','3','5–80 τ','−0.034230 + 6.184099 i'],
['Short / ξ₁–ξ₂','1','50–80 τ','−0.074893 + 6.317775 i'],['Short / ξ₁–ξ₄','1','35–80 τ','−0.079402 + 6.339246 i']])}
<p>U and ψ remain trained from <strong>unforced trajectories only</strong>. Existing forced recordings train B, not ψ.
The mild target has amplitude 0.012 and frequency 0.93 cycles/τ; the nonlinear target has amplitude 0.3 and frequency 0.83 cycles/τ.
These correspond to pressure amplitudes 1.2 kPa and 30 kPa. Both start forcing at τ = 6 and end at τ = 60.</p>
<p>New work comprised: projection parity; free-motion accuracy by time interval; 24 matched-history branch points,
each tested with eight one-hold input choices; two input-release recordings; inverse-driver diagnostics;
a tiny paired-sensitivity bilinear fit; and completion of four previously storage-interrupted KLQR runs.
The 24 probe points lie along the two target histories and therefore are <strong>diagnostic interventions, not independent validation data</strong>.</p></section>

<section><h2>2. Implementation checks: what is not the main problem</h2>
<ul><li>Offline FFT-based projection and online ring-buffer projection agree to relative L² error below 2.9 × 10⁻¹⁵ on the sampled windows.
Both use oldest-to-newest history ordering and the same U / singular-value normalization. Maximum absolute discrepancy is below 4.6 × 10⁻¹⁷.</li>
<li>The integrator stores the newly selected input at each tick; that is the input driving the following interval. The regression uses that same interval-start input.
The target indexing in the feedback loop is aligned to the completed history endpoint, not its midpoint.</li>
<li>Replaying the recorded input in the probe acquisition adapter reproduced the stored echo exactly in this run (relative discrepancy zero).
Halving/doubling the probe amplitude changed the central input sensitivity by less than 1.2 × 10⁻⁶ in relative L² across these cases.</li>
<li>The existing algebraic regression test verifies exact recovery for synthetic data generated by the assumed bilinear model.
That test passes, but it does not establish that this bubble/history representation actually obeys that model.</li>
<li>The sampled experiment calls the existing continuous-time Riccati solver every hold and applies u = −K error.
It is an isolated sampled adapter, not a claim that the production continuous-RHS controller has become a discrete-time LQR.</li></ul>
<p>There is no evidence here of a projector-order, target-time or simple input-recording error causing the scientific failure.
This does not certify every implementation detail, nor eliminate laboratory timing/noise issues.</p></section>

<section id="causal"><h2>3. The decisive experiment: change only the next input</h2>
<p>Let z = ψ(ξ), μ = exp(λh). The existing fit assumes:</p>
<div class="equation">qₖ = zₖ₊₁ − μzₖ ≈ B̂<sub>d</sub>(zₖ)uₖ,<br>B̂<sub>d</sub>(z) = c₀ + c₁ Re(z) + c₂ Im(z).</div>
<p>At each probe point, the complete echo prehistory was held identical. The acquisition adapter then continued the plant for one hold with
u = 0, ±δ, ±2δ, ±δ/2, and the recorded input. Here δ = 0.01 for the mild case and 0.1 for the nonlinear case.
Only echo histories, echo continuations, time and input were passed into the diagnostic calculations. Simulator states stayed inside acquisition.</p>
<div class="equation">aₖ = zₖ₊₁(0) − μzₖ,<br>
B<sub>causal</sub>(Hₖ) ≈ [zₖ₊₁(+δ) − zₖ₊₁(−δ)] / (2δ),<br>
zₖ₊₁(u) − μzₖ = aₖ + B<sub>causal</sub>(Hₖ)u + higher-order input terms.</div>
<p><strong>aₖ does not vanish on these forced histories.</strong> No choice of a finite B can repair that discrepancy at u = 0.
The old regression tries to put this residual into B because current input is correlated with the past driven motion.</p>
<p>For the simplest constant-B illustration, if q = a + Bu, least squares through the origin gives
B̂ = B + E[u a] / E[u²]. A periodic input and its driven history can make E[u a] nonzero.
More samples of the same correlated experiment need not remove this bias. The actual bilinear regression has the same omitted-drift issue with a larger feature matrix.</p>
{row_table(['Model / target / local map','Zero-input residual / actual input effect','‖B̂ − Bcausal‖ / ‖Bcausal‖','Input-affine response residual'],paired_rows)}
<p>Ratios are norms over 12 probe phases, not pointwise maxima. The B-error column is a dimensionless relative error, not a multiplicative scale estimate;
it includes direction as well as magnitude error. The final column compares the measured finite input response with B<sub>causal</sub>u,
using the recorded input only as an evaluation condition.</p>
{causal_plot}
<p>For the long, two-latent mild model, zero-input drift is about 11,444 times the actual one-hold input effect.
The stored local B differs from measured sensitivity by relative error 8,609. This is much more than a modest fit error.
The short, two-latent mild model reduces those figures to about 22.3 and 19.6, but the confounding remains.</p></section>

<section id="capacity"><h2>4. Can the existing B-map family fit a cleaner target?</h2>
<p>To separate regression-target contamination from limited map capacity, the same native bilinear fitter was given the paired sensitivities.
Eight phase points were used for fitting and four interleaved phases were withheld. The model family remained c₀ + c₁ Re(ψ) + c₂ Im(ψ).
No new controller was run with these diagnostic maps.</p>
{row_table(['History / coordinates','Target','Fit error','Withheld-phase error','Normalized design condition','Paired inverse driver error'],capacity_rows)}
<p><strong>For the mild target, the existing map family is locally adequate:</strong> withheld-phase sensitivity errors are 0.28–2.61%.
That is promising evidence for improving the data/fit target without inventing a new map representation.
It is not evidence that the autonomous drift has disappeared, or that a map fitted on this one periodic history generalizes.</p>
<p><strong>For the nonlinear target, it is not adequate in this test:</strong> withheld-phase errors are 32.4%, 51.3%, 107.0% and 131.2%
for long-2, long-4, short-2 and short-4 respectively. All fits have full rank three; normalized design condition numbers are 1.1–6.6.
Simple numerical collinearity is therefore not the leading explanation. The tiny sample and single-orbit coverage limit the conclusion:
this is a warning about current features/representation, not a global approximation impossibility theorem.</p>
<p>The paired inverse column uses both the locally measured zero-input continuation and causal sensitivity—information unavailable from a single target echo record.
It is an optimistic diagnostic, not a deployable controller. Mild input recovery is about 0.14%; nonlinear recovery is still about 50–55%.
One nonlinear phase has very weak one-step sensitivity: in short-2, at τ = 30.9, |B<sub>causal</sub>| ≈ 0.001736 and the local inverse gives
−0.6315 instead of −0.2601. Small forward errors can become large inverse errors near weak sensitivity. This does not prove loss of multi-step controllability.</p></section>

<section id="history"><h2>5. Why a free-trained history coordinate is not automatically a state coordinate under forcing</h2>
<p>The causal coordinate is a linear filter of a whole measured history:</p>
<div class="equation">Hₖ = [eₖ₋d₊₁, …, eₖ]ᵀ,<br>ξⱼ,ₖ = Σ<sub>ℓ=0…d−1</sub> U<sub>ℓj</sub> eₖ₋d₊₁₊ℓ / σⱼ,<br>zₖ = ψ(ξ₁,ₖ, …, ξᵣ,ₖ).</div>
<p>Along an autonomous free trajectory, its history is tied to the autonomous flow. Training ψ there constrains its evolution on that collection of histories.
Under arbitrary forcing, the same current physical state can have different past echoes because different past inputs brought it there.
The stored history, and therefore ξ, need not be a function of current (r, ṙ) alone. A polynomial that is approximately an eigenfunction on free histories
need not retain the same autonomous evolution on forced histories.</p>
<p>This is the important correction to the earlier idea that <em>only the coordinate transform changes and everything after it is identical</em>.
The polynomial data classes and evaluators can remain identical. The dynamical closure assumption cannot simply be carried across.
A causal input derivative on the full history state can exist without the identity ż = λz + B(z)u being valid.</p>
<p>Over one control hold, only 20 newly measured samples enter the buffer. The long window retains most of its 8644 samples.
That does not itself quantify delay or authority—the weights and observable matter—but it explains why instantaneous sensitivity and accumulated forcing effects
are different quantities. The matched-history test measures that distinction directly.</p>
<h3>Input-release check</h3><p>Two additional recordings replayed the target to τ = 30, set input to zero, and continued for 12 τ.
The table shows ‖zₖ₊₁ − μzₖ‖ / ‖zₖ‖ in intervals measured after release.</p>
{row_table(['Model','Target history','1–4.3215 τ after release','4.3215–8 τ','8–12 τ'],release_rows)}
<p>The long-four nonlinear residual falls from 5.54% to 0.338% after the old forced samples have left the long buffer.
This is consistent with a history-manifold mismatch. Amplitude also decays during release, so this experiment does not uniquely isolate memory from nonlinear-amplitude effects.
The short-window nonlinear models remain poor after their much shorter history has washed out, demonstrating that memory is not the only problem.</p></section>

<section id="free"><h2>6. What eigenfunctions do we actually have?</h2>
<p>These are polynomial approximate eigenfunctions selected by sampled EDMD from unforced data and stored in the existing polynomial KROM representation.
They were not newly fitted by a continuous generator SINDy derivative equation in this investigation.
The representation class is shared; the identification objective differs.</p>
<p>Six free trajectories excluded from coefficient fitting were re-evaluated by interval. They were previously used for candidate selection,
so they are validation records, not an untouched final test set. Errors below pool squared norms across trajectories;
the older report averaged per-trajectory errors, so its aggregate percentages are not directly interchangeable.</p>
{row_table(['Model','Starting-time interval','One-hold error','1 τ prediction error','4 τ prediction error'],free_rows)}
<p>Long-2 is quite good late (1.91% at one τ) but much worse early (53.4%). Short-2, which produced the best mild closed loop,
has 147.1% early one-τ prediction error and 3.91% late error. Long-4 improves early one-τ accuracy to 18.1%,
but short-4 worsens it to 173.1%. <strong>More coordinates do not automatically produce a better eigenfunction.</strong></p>
<p>The dimension comparisons also change polynomial degree and selected fitting interval, so they do not isolate dimension alone.
Moreover, four latent inputs still become one complex ψ. The B map sees only its real and imaginary parts, not four independent feedback coordinates.
These results do not establish a globally injective observation or a globally accurate free eigenfunction over the nonlinear operating region.</p></section>

<section id="tracking"><h2>7. A separate limitation: the tracking law</h2>
<p>Both the sampled experiment and the existing native KLQR value calculation use:</p>
<div class="equation">uₖ = −Kₖ(zₖ − z*ₖ).</div>
<p>At exact tracking, the error is zero and the commanded input is zero. A viscous bubble’s sustained sine-driven response requires persistent nonzero forcing.
Therefore this law generally needs a persistent tracking error to generate that forcing. Better B estimation alone cannot remove this structural limitation.</p>
<p>Even with a perfect lifted model, in split-real coordinates the error obeys the following equation.
Here M is the two-by-two real representation of μ, b<sub>d</sub> is the two-by-one input vector, and η = z − z*:</p>
<div class="equation">ηₖ₊₁ = (M − b<sub>d</sub>Kₖ)ηₖ − (z*ₖ₊₁ − Mz*ₖ).</div>
<p>The final term is nonzero for a forced reference. Standard nominal-trajectory LQR includes a nominal input in addition to error feedback;
see <a href="https://underactuated.mit.edu/lqr.html">MIT’s LQR and trajectory-stabilization derivation</a>.
No nominal-input term was added here. Exact latent tracking need not imply exact echo tracking if ψ is non-injective, which makes the control interpretation harder, not easier.</p>
<h3>Can the existing fitted map recover the target input offline?</h3>
<p>With the true driver withheld from the calculation, the minimum-residual real input for the assumed one-complex-mode model is:</p>
<div class="equation">ûₖ = Re[conj(B̂<sub>d</sub>(z*ₖ)) · (z*ₖ₊₁ − μz*ₖ)] / |B̂<sub>d</sub>(z*ₖ)|².<br>
q⊥,ₖ = (z*ₖ₊₁ − μz*ₖ) − B̂<sub>d</sub>(z*ₖ)ûₖ.</div>
<p>The known driver is used only for scoring. This offline inversion uses the next target sample and is not a causal state estimator or deployed control law.
q⊥ measures the residual no real scalar input can explain under this particular fitted model.</p>
{row_table(['Model / target / local map','Inverse driver error','Driver cosine','Perpendicular residual','Residual using known driver'],inverse_rows)}
{inverse_plot}
<p>The short-2 mild controller’s prior 8.88% echo error does <em>not</em> mean its B map identifies the driver: this inverse has 128.2% driver error.
For the nonlinear target, all these local-map inverses have driver error above 100%. Long-window models leave 78–86% of the nonlinear residual perpendicular to B.
Short-window maps may have smaller perpendicular fractions yet still have completely wrong gain and input recovery.</p>
<p>Thus simply adding inverse-model feedforward would not currently solve the nonlinear problem. The inverse model fails its own prerequisite test.</p></section>

<section><h2>8. Observability, controllability and numerical conditioning are different questions</h2>
<p>Nothing here proves that echo plus its full history is intrinsically insufficient for control. We have established failure of these particular truncated coordinates,
polynomials, input maps and tracking law. We have not measured a global observability rank, found an exact indistinguishable-history counterexample,
or proved that no alternate free-trained representation works.</p>
<p>A controllability test of the fitted lifted pair alone can be misleading. For one complex λ = a + iω and a real input vector b in its two-real-coordinate realization,</p>
<div class="equation">A = [[a, −ω], [ω, a]], &nbsp; det[b, Ab] = ω‖b‖².</div>
<p>For ω ≠ 0 and b ≠ 0, that model has full two-state controllability rank. This says nothing about whether its B is the plant’s causal response.
Likewise, full regression-design rank does not eliminate bias from an incorrect drift model. A weak one-step response can make input inversion fragile
without proving the physical system is uncontrollable over a longer horizon.</p>
<p>Delay-coordinate methods and data-driven controlled Koopman models are established possibilities, but not guarantees for this chosen representation.
<a href="https://www.nature.com/articles/s41467-017-00030-8">The original HAVOK paper</a> distinguishes leading delay coordinates from unresolved forcing coordinates;
<a href="https://arxiv.org/abs/1611.03537">Korda and Mezić’s controlled Koopman work</a> and
<a href="https://arxiv.org/abs/1804.05291">Arbabi, Korda and Mezić’s sparse-measurement control study</a> provide relevant context.
These references are not evidence that our specific free-trained ψ and memoryless B must work, and no alternative method from them was implemented.</p></section>

<section><h2>9. Completion of the interrupted control runs</h2>
<p>The four short-history, four-latent, nonlinear-local-map attempts interrupted by disk exhaustion were repeated in a new output directory using their frozen models and unchanged controller.</p>
{row_table(['R multiplier','Late echo error','Late driver error','Late saturation'],resumed_rows)}
<p>None achieved nonlinear tracking. Including these results, 102 of the original 104 parameter cases now have completed saved outcomes;
two original small-penalty cases still have Riccati numerical errors. Those two were not repaired because primary changes are out of scope and other evidence already rejects the scientific model.
The earlier best mild result remains 8.88% late echo error and 21.61% driver error. The earlier best nonlinear echo error remains about 98.3%.</p>
<details><summary>Prior best mild control and response plots</summary><img src="../../viewer_outputs/echo-input-short-window-mild.png" alt="Prior best mild control: driver, echo, radius, velocity and phase plots"></details>
<details><summary>Prior nonlinear comparison plots</summary><img src="../../viewer_outputs/echo-input-long-window-nonlinear.png" alt="Prior nonlinear control comparison"></details></section>

<section id="recommendations"><h2>10. Recommended next steps, within the existing framework</h2>
<ol><li><strong>Retire the present maps as validated causal Bψ estimates.</strong> Preserve them as experiment artifacts. They are residual regressions,
often measuring forced-history correlation rather than the next-input effect. Do not use fit residual or closed-loop appearance alone as validation.</li>
<li><strong>Use paired perturbations to identify input response.</strong> Repeated, matched echo histories followed by zero and ±input continuations allow
the existing fitter to target sensitivity instead of mixed drift. Begin with repeated unforced prefixes, then assess histories reached under control.
In a laboratory these are repeated preparations, not exact state cloning; reproducibility and echo-history matching must be quantified.</li>
<li><strong>Make zero-input closure a gate before further KLQR tests.</strong> Keep U and ψ training unforced. Evaluate ψ after switching input off from driven histories.
If zₖ₊₁(0) ≈ μzₖ fails by an amount comparable to or larger than the actual input response, a B-only regression cannot repair that ψ model.
The tiny mild sensitivity-fit result is encouraging, but the drift gate still fails.</li>
<li><strong>Continue the requested two-to-four progression with controlled comparisons.</strong> Use the same free-training coverage, comparable dictionaries,
separate amplitude/time validation bands and fixed test cases. Require gains in autonomous closure and sensitivity predictability, not just Hankel energy.
If needed, the existing multi-mode KROM structures can retain more than one eigenfunction; that is a future within-framework option, not something tested here.</li>
<li><strong>Treat nonlinear map capacity as a separate gate.</strong> The tiny paired test already leaves 32–131% held-out sensitivity error.
More independent local calibration trajectories are justified, but merely repeating longer sines cannot prove a memoryless bilinear map is adequate.
Do not assume general B(ξ) has been tested: the current map is affine in one complex ψ.</li>
<li><strong>Only after model validation, revisit persistent-target tracking.</strong> A nominal-input or proper tracking formulation would require an explicitly approved
controller change. First demand reliable offline target-input reconstruction and low perpendicular residual without supplying the true target driver to the calculation.
Then retain KLQR as the feedback part; do not mistake gain tuning for solving the nominal-forcing issue.</li></ol>
<p class="note"><strong>Scope conclusion:</strong> with these frozen free-trained ψ, fixed μ and B(z)u-only lifted dynamics,
accurate general nonlinear control is outside what this particular fitted model presently supports. The observed nonzero zero-input drift cannot be corrected by B alone.
This is not a rejection of the whole DynSysLab framework or echo-only control. If closure cannot be obtained with validated free-trained coordinates/eigenfunctions
and the existing model structures, stop and seek approval before introducing history-aware drift, input-history augmentation, observers or a different controller.</p>
<p>No such expansion was implemented. Successful nonlinear control remains unresolved; the diagnosis is the completed outcome of this goal.</p></section>

<section id="reproduce"><h2>11. Evidence, reproducibility and limits</h2>
<ul><li><a href="../../tools/experiments/echo_control_diagnosis.cpp">Isolated C++ diagnostic source</a> (no primary source edits).
It links existing native DLLs and the unchanged <a href="../../native/tools/echo_input/Experiment.cpp">experiment adapter</a>.</li>
<li><a href="../../archives/experiments/echo-control-diagnosis-20261007-v2/diagnosis.json">Full native diagnostic metrics</a>;
<a href="../../archives/experiments/echo-control-diagnosis-20261007-v2/paired-capacity/diagnosis.json">paired-fit capacity metrics</a>;
<a href="../../archives/experiments/echo-control-diagnosis-20261007-v1/diagnosis.json">completed interrupted controls</a>.</li>
<li><a href="../../archives/experiments/echo-control-diagnosis-20261007-v2/mild-probes.npz">Mild measurement-only probes</a> and
<a href="../../archives/experiments/echo-control-diagnosis-20261007-v2/nonlinear-probes.npz">nonlinear measurement-only probes</a>.
Their metadata records replay error and the information boundary. Other arrays and native Viewer payloads are retained locally alongside these outputs.</li>
<li><a href="ECHO_HISTORY_BPSI_EXPERIMENT_2026-10-06.html">Prior experiment report</a> remains unchanged;
<a href="../../native/src/HistoryInputRegression.cpp">regression implementation</a>,
<a href="../../native/src/HistoryKoopman.cpp">projection implementation</a>, and
<a href="../../native/src/LQR.cpp">native LQR/controller implementation</a> are the audited sources.</li></ul>
<p>Standalone build: compile this audit source and unchanged Experiment.cpp with MSVC C++17, /O2 /MD, native/include and native/tools/echo_input includes;
link dynsyslab_foundation, dynsyslab_koopman, dynsyslab_eigenfunctions, dynsyslab_history_koopman, dynsyslab_control and dynsyslab_control_maps import libraries
from build/native-suite-fast/lib. Put build/native-suite-fast/bin on PATH when executing.</p>
<pre>diagnosis.exe ROOT NEW_OUTPUT_DIRECTORY
diagnosis.exe ROOT NEW_OUTPUT_DIRECTORY --skip-controls
diagnosis.exe ROOT NEW_OUTPUT_DIRECTORY --paired-only</pre>
<p>The paired-only mode reads the recorded v2 probes. Outputs are immutable and existing diagnosis.json is rejected.
Version 1 completed the interrupted controls; version 2 added interval-resolved release analysis; its paired-capacity subdirectory contains the final small sensitivity regression test.</p>
<p>Limitations: two target amplitudes/frequencies, one period sampled at 12 phases per target, no measurement noise, no direct scattering,
no actuator dynamics, no independent closed-loop confirmation, no proof of injectivity, and no formal nonlinear stability certification.
Finite-difference sensitivity is about u = 0 at fixed history, not an exact instantaneous state Jacobian or a global input map.
The sampled-to-continuous conversion used by prior controls freezes B over a hold; its error was not separately isolated here.</p>
<p>The final verification record is <a href="../../archives/experiments/echo-control-diagnosis-20261007-v2/verification.json">verification.json</a>:
primary-source hashes, targeted tests and rendered-report checks. No large scientific/media assets were committed or published.</p></section>
</main></body></html>'''
    REPORT.write_text(report, encoding="utf-8")
    print(REPORT)


if __name__ == "__main__":
    main()
