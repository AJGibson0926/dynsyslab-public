"""Presentation of native subset-specific B fits; no scientific fitting in Python."""
from __future__ import annotations
import html
import json
from pathlib import Path
import numpy as np
from report_echo_multimode import read, pct, table, image, get, SUBSETS, COLORS
from plot_rpe_echo_jet_pipeline_evidence import figure, render, trace

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'archives/experiments/echo-multimode-20261009-v2'
DATA = ROOT / 'archives/experiments/echo-multimode-joint-20261011-v1'
INPUTS = ROOT / 'archives/experiments/echo-input-regression-20261006-v1'
REPORT = ROOT / 'docs/reports/ECHO_MULTIMODE_JOINT_BPSI_2026-10-11.html'


def response(rows, window, target, penalty):
    """Render refitted-controller responses at one common target, window, and penalty."""
    figures = []
    with np.load(INPUTS / f'target-{target}.measurements.npz') as source:
        ref = {k: source[k] for k in ('time', 'input', 'echo')}
    with np.load(INPUTS / f'target-{target}.truth.npz') as source:
        ref['r'], ref['v'] = source['states'][:, 0], source['states'][:, 1]
    for key, label in [('input', 'u = p / p0'), ('echo', 'Echo'), ('r', 'r / R0'), ('v', 'dr / dτ')]:
        ix = np.flatnonzero((ref['time'] >= 54) & (ref['time'] <= 58))
        if key == 'input':
            ix = ix[::20]  # Display only: interval-start samples of held input.
        traces = [trace('Known driver / target', ref['time'][ix], ref[key][ix],
                        color='#58616b', line_style='dash', width=2.5)]
        for subset, color in zip(('1', '12', '123'), COLORS):
            row = get(rows, window, target, subset, penalty)
            if 'failure' in row:
                continue
            with np.load(DATA / (row['id'] + '.measurements.npz')) as source:
                t = source['time']
                if key in ('input', 'echo'):
                    y = source[key]
                else:
                    with np.load(DATA / (row['id'] + '.truth.npz')) as truth:
                        y = truth['states'][:, 0 if key == 'r' else 1]
            assert np.array_equal(t, ref['time'])
            traces.append(trace(f'ψ{subset}; refitted B; ρ={penalty:g}', t[ix], y[ix], color=color))
        f = figure('time', f'{window.capitalize()} / {target}: {label}', r'\tau', label, traces)
        f['value_scale'] = 'linear'
        figures.append(f)
    name = f'echo-multimode-joint-{window}-{target}-20261011'
    render(name=name, label=f'Subset-specific B / {window} / {target}',
           page_label=f'Same penalty ρ={penalty:g}; two latent coordinates; one, two, three eigenfunctions',
           source=DATA / 'audit.json', figures=figures, columns=2, height=1100)
    return image(name, f'{window} {target}: refitted control, echo, radius and velocity')


def sweep(rows, target):
    """Plot all saved subset and penalty comparisons for the chosen target."""
    figures = []
    for key, label in [('echo_nrmse', 'Echo relative RMS error'), ('driver_nrmse', 'Driver relative RMS error')]:
        for window in ('long', 'short'):
            traces = []
            for penalty, color in zip((100, 1, .01), COLORS):
                points = [(i + 1, get(rows, window, target, subset, penalty)) for i, subset in enumerate(SUBSETS)]
                points = [(i, r['late'][key]) for i, r in points if 'late' in r]
                traces.append(trace(f'ρ={penalty:g}', np.array([p[0] for p in points]),
                                    np.array([p[1] for p in points]), color=color, marker=True, marker_size=6))
            f = figure('time', f'{window.capitalize()} / {target}', 'Subset ID (key in report)', label, traces)
            f['value_scale'] = 'log'
            figures.append(f)
    name = f'echo-multimode-joint-{target}-sweep-20261011'
    render(name=name, label=f'Joint B refit / {target} / all subsets',
           page_label='Subset IDs: 1, 2, 3, 12, 13, 23, 123; numerical failures omitted',
           source=DATA / 'audit.json', figures=figures, columns=2, height=1100)
    return image(name, f'{target}: all seven refitted subsets and three penalties')


def main():
    """Assemble the saved local evidence and native Viewer plots into the dated HTML report."""
    audit, old_audit = read(DATA / 'audit.json'), read(OLD / 'audit.json')
    rows = {r['id']: r for r in audit['results']}
    prior = {r['id']: r for r in old_audit['results']}
    assert len(rows) == 100 and audit['completed'] == 85 and audit['failed'] == 15
    same = []
    for window, penalty in [('short', 1), ('long', 100)]:
        for subset in ('1', '12', '123'):
            a, b = get(prior, window, 'mild', subset, penalty), get(rows, window, 'mild', subset, penalty)
            same.append([window, subset, penalty, pct(a['late']['echo_nrmse']), pct(b['late']['echo_nrmse']),
                         pct(a['late']['driver_nrmse']), pct(b['late']['driver_nrmse'])])
    best = []
    for window in ('long', 'short'):
        for target in ('mild', 'nonlinear'):
            for n in (1, 2, 3):
                eligible = [r for r in rows.values() if r['id'].startswith(f'{window}-{target}-') and
                            r['modes'] == n and r['R_multiplier'] in (100, 1, .01) and 'late' in r]
                r = min(eligible, key=lambda x: x['late']['echo_nrmse'])
                best.append([window, target, r['id'].split('-m')[-1].split('-r')[0], r['R_multiplier'],
                             pct(r['late']['echo_nrmse']), pct(r['late']['driver_nrmse']),
                             pct(r['truth_late']['radius_oscillation_nrmse']), pct(r['truth_late']['velocity_nrmse'])])
    fits = []
    for window in ('long', 'short'):
        for f in read(DATA / f'{window}-joint-maps.json'):
            for j, row in enumerate(f['train']['per_mode']):
                fits.append([window, f['target'], f['subset'], row['mode'], pct(row['relative_residual']),
                             pct(f['heldout_waveforms']['per_mode'][j]['relative_residual']),
                             pct(f['target_score_only']['per_mode'][j]['relative_residual']),
                             f"{f['rank']}/{f['columns']}", f"{f['condition']:,.2g}"])
    weight = []
    for window, penalty in [('short', 1), ('long', 100)]:
        for n, subset in [(2, '12'), (3, '123')]:
            a, b = get(rows, window, 'mild', '1', n * penalty), get(rows, window, 'mild', subset, penalty)
            weight.append([window, f'ψ1, ρ={n*penalty:g}', f'ψ{subset}, ρ={penalty:g}',
                           pct(a['late']['echo_nrmse']), pct(b['late']['echo_nrmse']),
                           pct(a['late']['driver_nrmse']), pct(b['late']['driver_nrmse'])])
    causal = [[c['window'], c['target'], r['mode'], f"{r['zero_input_drift_over_actual_input_effect']:,.2f}×",
               f"{r['regressed_B_relative_error_to_causal_B']:,.3g}"] for c in audit['causal'] for r in c['modes']]
    attempts = [[r['id'], pct(r['late']['echo_nrmse']), pct(r['late']['driver_nrmse']),
                 pct(r['late']['saturation_fraction']), f"{r['maximum_care_relative_residual']:.3g}"]
                if 'late' in r else [r['id'], 'FAILED', r['failure'], '—', '—']
                for r in sorted(rows.values(), key=lambda x: x['id'])]
    plots = {f'{w}-{t}': response(rows, w, t, p) for w, t, p in
             [('short', 'mild', 1), ('long', 'mild', 100), ('short', 'nonlinear', 100), ('long', 'nonlinear', 100)]}
    sweeps = sweep(rows, 'mild') + sweep(rows, 'nonlinear')
    # Preserve the report's established reading style; no Viewer/UI redesign.
    previous = (ROOT / 'docs/reports/ECHO_MULTIMODE_CONTROL_2026-10-09.html').read_text(encoding='utf-8')
    style = previous.split('<style>', 1)[1].split('</style>', 1)[0]
    report = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Refitting Bψ for every eigenfunction set</title><style>{style}</style></head><body><main>
<p>DynSysLab · echo-only · 11 October 2026 · isolated native C++ experiment</p>
<h1>Refitting B<sub>ψ</sub> for every eigenfunction set</h1>
<p class="note"><strong>The correction changes the answer.</strong> With two latent coordinates and the short window, jointly refitting B for ψ1+ψ2 improves mild-target echo error from 8.50% (ψ1 alone) to <strong>5.13%</strong>, and driver error from 9.58% to <strong>8.63%</strong>, at the same penalty. The earlier inherited-row pair had 14.43% echo and 187.44% driver error. A third mode retains good echo tracking (5.46%) but worsens driver recovery (39.03%). The nonlinear target remains unsuccessful.</p>
<p>Your correction: <q>No, I want to refit B_psi for each set.</q> This report replaces the earlier conclusion for the <em>joint-refit</em> experiment. The <a href="ECHO_MULTIMODE_CONTROL_2026-10-09.html">October 9 report</a> is preserved as the restricted, fixed-row comparator.</p>
<nav><a href="#method">Exact refit</a> · <a href="#results">Control results</a> · <a href="#diagnosis">Why it changes</a> · <a href="#limits">Remaining limitations</a> · <a href="#verification">Verification</a></nav>
<section id="method"><h2>1. Precisely what was refitted</h2>
<p>Each window has three fixed complex eigenfunctions ψ1, ψ2, ψ3 of the same two latent coordinates (ξ₁, ξ₂). For every nonempty subset S = {{1}}, {{2}}, {{3}}, {{1,2}}, {{1,3}}, {{2,3}}, {{1,2,3}}, <strong>every coefficient of every selected B row was refitted</strong> from the calibration records, including cross-mode predictors. The resulting 28 separate maps cover 2 windows × 2 calibration regimes × 7 subsets. Penalty changes reuse the map for that set; penalties do not affect regression.</p>
<div class="eq">zⱼ,ₖ = ψⱼ(ξ₁,ₖ, ξ₂,ₖ), &nbsp; μⱼ = exp(λⱼh), &nbsp; qⱼ,ₖ = zⱼ,ₖ₊₁ − μⱼzⱼ,ₖ<br>
qⱼ,ₖ ≈ [cⱼ₀ + Σ<sub>ℓ∈S</sub>(aⱼℓ Re zℓ,ₖ + bⱼℓ Im zℓ,ₖ)]uₖ, &nbsp; j ∈ S.</div>
<p>For m complex eigenfunctions, the real design matrix has 1+2m columns: u times a constant, then the real and imaginary parts of <em>all</em> selected modes. A simultaneous least-squares solve fits 2m real output columns. It is equivalent to separate output-row regressions sharing the same expanded design—not an additional constraint tying their coefficients together. Columns are normalized by their L² norms; the existing native SVD least-squares routine uses ridge 10⁻⁸. The same single-mode calculation agrees with the existing history-input fitter. This is not a constant B matrix: B is evaluated as an affine function of the selected eigenfunction values, hence ultimately of echo-history coordinates.</p>
<p>The map is stored in the existing conjugate-complete <code>BilinearControlMap</code>: B<sup>c</sup>ⱼ(z) = βⱼ₀ + Σ(Aⱼℓzℓ + Cⱼℓconj(zℓ)). Discrete coefficients are divided by Iⱼ = (exp(λⱼh)−1)/λⱼ for use in the unchanged continuous controller. This preserves the previous frozen-map sampled-to-continuous approximation; it is not claimed to be exact sampled-data LQR.</p>
<h3>Held fixed</h3>
<ul><li><strong>Two latent coordinates in every run.</strong> Long: 8644 echo samples spanning 4.3215 τ (13.018 μs). Short: 101 samples spanning 0.05 τ (0.15062 μs). Sensor step 0.0005 τ; controller hold h=0.01 τ; one τ=3.012383 μs. The existing U matrices are reused, with their first two columns only.</li>
<li><strong>No SVD or eigenfunction refit.</strong> The exact October 9 banks are loaded by reference. They were fitted to eight unforced recordings over 35–80 τ with a degree-five, 21-term polynomial dictionary. Three stable positive-frequency modes were selected using reused free-validation recordings, one each from [4,8), [8,16), [16,24) rad/τ. Their errors and validity limitations are unchanged.</li>
<li><strong>Same input measurements.</strong> Local sine records 1–6 calibrate each regime; 7–8 are withheld complete waveforms. Mild calibration amplitudes: 0.004, 0.008, 0.016 (0.025 withheld). Nonlinear: 0.18, 0.24, 0.36 (0.42 withheld), two phases each. The interval is 12 τ through record end minus 6.02 τ. The target record is score-only, never used in B regression.</li>
<li><strong>Same controller.</strong> wⱼ=ψⱼ/RMS(ψⱼ,target), Q=I/m, R=ρ/0.01; ρ∈{{100,1,0.01}}. Existing continuous CARE recomputed each hold; u=clip[−K(w−wtarget),−0.6,0.6]. No feedforward. Feedback starts at 6 τ and ends at 60 τ. Target normalization uses target echo, not its driver.</li>
<li><strong>No hidden states in identification or feedback.</strong> Only measured echo/history and calibration input enter regression. Only echo/history and the fitted model enter feedback. The simulator alone uses equations of motion; radius and velocity are scoring/plotting truth. The sensor retains the prior intrinsic-echo convention, omitting direct driver scattering.</li></ul></section>
<section id="results"><h2>2. Same-penalty comparison: inherited versus refitted B</h2>
{table(['Window','Subset','ρ','Old echo error','Refit echo error','Old driver error','Refit driver error'], same)}
<p>Errors are relative L²/RMS errors over 30–60 τ; no time/phase alignment is applied. The short ψ1+ψ2 pair has zero saturation, driver cosine similarity 0.9966 and maximum CARE residual 2.45×10⁻¹⁵. Over the final 50–60 τ, echo and driver errors are 4.98% and 8.57%. Its radius-oscillation and velocity errors are 5.19% and 5.12%. Over the entire driven interval 6–60 τ, echo error is 4.94% and driver error 10.85%. This is a useful finite-horizon mild-case result, not proof of indefinite stability or exact driver recovery.</p>
{plots['short-mild']}
<p>The third mode introduces visible unwanted input structure even though the bubble response remains close. Echo tracking and recovery of the known driver are distinct objectives; ranking solely by echo error conceals this.</p>
<h3>Long-window mild response, same ρ=100</h3>{plots['long-mild']}
<p>Long-window echo error improves from 62.37% (one mode) to 60.79% (pair) and 52.15% (triple), but this remains poor tracking. Driver error is best for the pair here (51.09%); the triple is worse (68.94%). Long history does not become a successful controller merely by refitting.</p>
<h3>Weight-matched singles</h3>
{table(['Window','Single','Multi','Single echo','Multi echo','Single driver','Multi driver'], weight)}
<p>Because Q=I/m, a one-mode run at mρ matches the fundamental cost/input ratio of an m-mode run at ρ. The short pair also improves on this weight-matched single, so its benefit is not merely a redistribution of the original fundamental's weight. This does not isolate feedback dimension from B predictor capacity: both change when the selected set changes, intentionally matching your requested experiment. At a <em>fixed</em> set, old-versus-new isolates the changed B fit.</p>
<h3>Best observed echo error per mode count</h3>
{table(['Window','Target','Subset','ρ','Echo error','Driver error','Radius oscillation error','Velocity error'], best)}
<p>These are post-hoc minima within the three-penalty main grid, not independently validated controller choices. Radius error is normalized by target r−1, not equilibrium radius. Values around 100% in the nonlinear cases mean essentially no tracking. The short/nonlinear single at ρ=0.01 also has an elevated CARE residual (3.07×10⁻⁵), so its tiny difference below 100% is not a certified gain.</p>
<details><summary>All subsets and penalties, using standard Viewer plots</summary><p>Subset IDs 1–7 correspond to 1, 2, 3, 12, 13, 23, 123. Log-y values are relative errors (1=100%). Failed runs are omitted, not assigned a fabricated score. Connecting lines guide the eye across discrete sets.</p>{sweeps}</details>
<h3>Nonlinear target: neither window recovers it</h3><p>For comparable completed curves, these plots hold ρ=100 for ψ1, ψ1+ψ2 and ψ1+ψ2+ψ3. They show the known target against near-zero controlled responses. The nonlinear target's known sine has amplitude 0.3 and frequency 0.83 cycles/τ; the mild target uses 0.012 and 0.93 cycles/τ.</p>{plots['short-nonlinear']}{plots['long-nonlinear']}</section>
<section id="diagnosis"><h2>3. What the refit improves—and what it does not</h2>
<p>For short/mild ψ2, the old self-only B row had <strong>102.67% withheld-waveform residual</strong>. In the jointly refitted pair, it is <strong>38.66%</strong>; ψ1's row is 17.51%. The pair's training residuals are 11.58% and 30.16%. Its scaled design is full rank (5/5), condition 1.92. Cross-mode predictors meaningfully improve the fitted forced-response relationship; this was missing from the earlier comparison.</p>
<p>In the short/mild triple, the held-out residuals are 17.97%, 46.68%, and <strong>82.39%</strong>. The third free-trained eigenfunction already has roughly 105% late 4-τ semigroup error. More polynomial features do not guarantee trustworthy extra controlled dynamics. The increased model capacity helps some rows, but the third controlled output is still poorly modeled and its feedback can distort the driver.</p>
<p>For the short/nonlinear triple, the scaled design condition is about 273,017. Held-out residuals are approximately 68%, yet the target-only residuals are orders of magnitude larger. This is severe extrapolation/model failure, not a small tuning error. Full numerical rank alone does not make a regression useful.</p>
<details><summary>Every refitted row: training, withheld-waveform and target-only validation</summary>
<p>Each row reports ‖q−B<sup>d</sup>u‖/‖q‖. Target-only scoring uses the projected target record after history warmup, whereas controller scores use 30–60 τ. Stored summary means are arithmetic means of per-mode relative residuals, not energy-pooled residuals. The files named <code>*-joint-maps.json</code> and <code>*.fit.json</code> contain the actual joint fits; <code>*-maps.json</code> retains scalar comparator diagnostics.</p>
{table(['Window','Regime','Set','Output mode','Train','Withheld','Target only','Rank/columns','Condition'], fits)}</details>
<h3>Does the improved fit now identify the causal input sensitivity?</h3>
<p><strong>No such conclusion is justified.</strong> The reused matched-history probes hold prehistory fixed and vary only the next input. They separate a=znext(0)−μz from the local causal input derivative. The following diagnostics evaluate the refitted <strong>three-mode maps only</strong>; they are not measurements of the pair's own causal-map error.</p>
{table(['Window','Regime','Mode','Zero-input drift / actual input effect','Refitted B relative error vs causal B'], causal)}
<p>The zero-input drift ratios are unchanged because ψ and the history transform are unchanged. On forced histories, the free-trained μz is still not generally the correct autonomous prediction. Regression of the total residual q against current input can absorb that drift and correlations with forcing history. Better observational prediction and useful mild-case feedback do not demonstrate that the fitted B is the physical causal actuation derivative. The pair's target-only scalar inverse still has 24.34% driver error and 21.11% unexplainable normalized lifted residual; this inverse diagnostic is not the feedback law, and its error need not match the closed-loop driver's 8.63%.</p></section>
<section id="limits"><h2>4. Revised conclusion and boundaries</h2>
<ul><li><strong>Use the jointly refitted short-window ψ1+ψ2 case as the next mild-case experimental reference.</strong> It improves both echo tracking and driver recovery over the single fundamental at the same penalty and over the weight-matched single.</li>
<li>Do not assume three eigenfunctions are better than two. Here the third harms driver recovery despite low echo error. Different choices of modes, more robust free validity and different weighting remain separate questions.</li>
<li>The long-window and nonlinear failures remain. Refitting B corrects an important limitation of the earlier test, but does not remove forced-history closure errors, extrapolation, numerical conditioning or lack of independent information in ξ₁,ξ₂.</li>
<li>Extra eigenfunctions remain functions of the same two measurements-derived coordinates. They cannot separate histories that already project to the same ξ₁,ξ₂. This run does not test four latent variables or a new history length.</li>
<li>All feedback laws still vanish at zero lifted tracking error. Without target feedforward, exact sustained tracking of a trajectory requiring nonzero forcing is not generally invariant. This controller limitation is held fixed, not repaired here.</li>
<li>Of 100 attempts (84 main-grid and 16 weight checks), <strong>{audit['completed']} completed and {audit['failed']} failed</strong> with the unchanged CARE solver reporting a non-real Riccati solution. These are numerical failures, not evidence of physical uncontrollability. Completed runs can also be numerically suspect: all maximum CARE residuals and saturation fractions are retained below.</li>
<li>Repeated use of the same validation records and the same two targets makes this a controlled diagnostic comparison, not a broad generalization certificate. No new independent robustness/noise test was performed.</li></ul>
<p>The core pipeline was not modified. No models were promoted into it. A conservative next step would be independent mild-target and perturbation validation of the refitted pair, with the same measurement-only causal diagnostics, before claiming general control capability.</p></section>
<section id="verification"><h2>5. Verification and reproduction</h2>
<p>All fitting, simulations, controller calculations, metrics and the scientific audit ran in native C++. Python only assembles these tables and native Viewer display payloads. The isolated adapter uses the existing least-squares, canonical map, Koopman model, causal history projector and continuous LQR APIs. No new controller or primary pipeline changes were required.</p>
<p>Synthetic dense cross-mode recovery has maximum coefficient error 1.92×10⁻¹⁶, with mode-order permutation checked. Every scalar refit is checked against the existing native scalar history-input fitter. Every refitted map's direct sampled coefficients are checked against its canonical operator-plus-conjugate evaluation. Saved completed control errors were recomputed from recordings; online/offline history projection agrees. The four native gates passed: history Koopman, history input regression, echo-history workflow, and KLQR case-study parity.</p>
<p>A separate archive check confirmed that all 28 subset models retain exactly the same polynomial coefficients, eigenvalues, exponents, centers and scales as their October 9 counterparts. All 16 multi-mode maps contain nonzero cross-mode operators. Single-mode canonical maps agree within 1.82×10⁻¹⁵ relative error; converting with λ directly rather than log(exp(λh))/h introduces floating-point roundoff, not a changed model family. These are algebraically equivalent reruns, not bitwise-identical trajectories. One already-failed, 99%-saturated long/mild ψ2 case at ρ=0.01 shifts from 606.75% to 599.77% echo error under that roundoff. Neither the good mild pair nor the reported conservative comparisons depends on this numerically sensitive case.</p>
<pre>{html.escape(json.dumps(audit['checks'], indent=2))}</pre>
<p>Local evidence: <a href="../../archives/experiments/echo-multimode-joint-20261011-v1/audit.json">complete audit</a>, <a href="../../archives/experiments/echo-multimode-joint-20261011-v1/short-joint-contract.json">refit contract and tests</a>, <a href="../../archives/experiments/echo-multimode-joint-20261011-v1/short-mild-m12.fit.json">pair fit diagnostics</a>. Original bank archives are referenced, not duplicated. Scientific archives and rendered graphics remain local and Git-ignored. No commit or publish was performed.</p>
<pre>powershell -File tools/experiments/Build-EchoMultimode.ps1
# Add build/native-suite-fast/bin to PATH. Use a NEW output directory.
build/echo-multimode/echo_multimode.exe . NEW_OUTPUT long --joint
build/echo-multimode/echo_multimode.exe . NEW_OUTPUT short --joint
build/echo-multimode/echo_multimode.exe . NEW_OUTPUT --audit
python tools/report_echo_multimode_joint.py</pre>
<details><summary>All 100 attempts, including numerical failures</summary>{table(['Attempt','Echo error','Driver error / failure','Saturation','Max CARE residual'], attempts)}</details>
</section></main></body></html>'''
    REPORT.write_text(report, encoding='utf-8')
    print(REPORT)


if __name__ == '__main__':
    main()
