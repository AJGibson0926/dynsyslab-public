"""Presentation and archive checks only; all acquisition/fitting/control is native C++.

Reuses the installed Viewer renderer. Does not fit scientific models or simulate.
"""
from __future__ import annotations
import html
import json
from pathlib import Path
import numpy as np
from report_echo_multimode import read, pct, table, image
from plot_rpe_echo_jet_pipeline_evidence import figure, render, trace

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'archives/experiments/echo-pulse-diagnosis-20261012-v1'
BANK = ROOT / 'archives/experiments/echo-multimode-20261009-v2'
OLD = ROOT / 'archives/experiments/echo-multimode-joint-20261011-v1'
INPUT = ROOT / 'archives/experiments/echo-input-regression-20261006-v1'
REPORT = ROOT / 'docs/reports/ECHO_PULSE_MEMORY_DIAGNOSIS_2026-10-12.html'
BLUE, ORANGE, GREEN, GREY = '#087fa0', '#b75c16', '#25824b', '#58616b'


def payload(path):
    """Load saved numerical arrays from a local NPZ without enabling pickle."""
    with np.load(path, allow_pickle=False) as a:
        return {k: a[k] for k in a.files}


def summarize():
    """Collect pulse, resolution, regression, and control evidence into a presentation audit."""
    acquisition = read(DATA / 'acquire-summary.json')['preparations']
    dense = read(DATA / 'dense-acquisition.json')['records']
    checks = read(DATA / 'resolution-checks/summary-with-immediate.json')['checks']
    responses = [r for w in ('long', 'short')
                 for r in read(DATA / f'{w}-summary.json')['preparations']]
    controls, fits = [], []
    for w in ('long', 'short'):
        for extra in ('', '-dense'):
            doc = read(DATA / f'{w}{extra}-pulse-controls.json')
            controls.extend(doc['controls'])
            fits.extend(dict(f, campaign='dense' if extra else 'sparse') for f in doc['fits'])
            assert doc['joint_fitter_tests']['synthetic_dense_cross_mode_coefficient_error'] < 1e-10
    assert len(acquisition) == 86 and len(responses) == 172
    assert len(fits) == 112 and len(controls) == 40
    assert sum('failure' in r for r in controls) == 3
    for meta in acquisition:
        a = payload(DATA / (meta['id'] + '.pulse.npz'))
        assert not {'states', 'values'}.intersection(a)
        assert a['echo_history'].shape == (8644,)
        assert a['echo_continuations'].shape[1] == 28001
        assert np.isfinite(a['echo_continuations']).all()
        np.testing.assert_allclose(a['echo_continuations'][:, 0], a['echo_history'][-1], atol=1e-12, rtol=0)
    for w in ('long', 'short'):
        for meta in dense:
            a = payload(DATA / (f'{w}-dense-' + meta['source'] + '.samples.npz'))
            assert not {'states', 'values'}.intersection(a)
            assert a['psi'].shape == (meta['points'], 3)
            assert np.isfinite(a['sensitivity']).all()
        # Dense free fits intentionally replicate the sparse baseline.
        for subset in ('1', '2', '3', '12', '13', '23', '123'):
            a = payload(DATA / f'{w}-pulse-free-m{subset}.krom.npz')
            b = payload(DATA / f'{w}-dense-pulse-free-m{subset}.krom.npz')
            np.testing.assert_array_equal(a['sampled_joint_coefficients'], b['sampled_joint_coefficients'])
    summary = dict(schema='dynsyslab.echo-pulse-audit/v1', preparations=len(acquisition),
                   response_analyses=len(responses), dense_points=sum(r['points'] for r in dense),
                   dense_train_points=sum(r['points'] for r in dense if r['split'] == 'train'),
                   map_artifacts=len(fits), distinct_fit_designs=98,
                   controls_completed=37, controls_failed=3,
                   maximum_projection_absolute_error=max(r['projection_max_absolute_error'] for r in responses),
                   maximum_free_replay_relative_error=max(r['free_continuation_relative_error'] for r in acquisition),
                   maximum_fundamental_immediate_refinement_change=max(r['immediate_sensitivity_checks'][0]['refined_relative_change'] for r in checks),
                   checks='timing, paired initial echo, shapes, finite arrays, no regression states, repeated-free coefficient parity')
    (DATA / 'presentation-audit.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return acquisition, dense, checks, responses, controls, fits, summary


def save_figure(name, title, figures):
    """Render a two-column native Viewer figure and return its report image link."""
    render(name=name, label=title, page_label=title, source=DATA / 'presentation-audit.json',
           figures=figures, columns=2, height=1140)
    return image(name, title)


def pulse_plot(preparation, label):
    """Compare long- and short-window impulse responses for a selected preparation."""
    figures = []
    records = {}
    for w in ('long', 'short'):
        a = payload(DATA / f'{w}-{preparation}.response.npz')
        meta = read(DATA / f'{w}-{preparation}.response.json')
        records[w] = (a, meta)
        t, g, prediction = a['time'], a['causal_response'][:, 0], a['predicted_response'][:, 0]
        clean = t >= meta['first_clean_control_time'] - 1e-9
        normalizer = abs(g[1])
        traces = [trace('Measured paired response', t[1:], abs(g[1:]) / normalizer, color=BLUE),
                  trace('From first hold', t[1:], abs(prediction[1:]) / normalizer, color=ORANGE, line_style='dash'),
                  trace('Re-anchored after window clears', t[clean], abs(a['post_clean_prediction'][clean, 0]) / normalizer,
                        color=GREEN, line_style='dash')]
        f = figure('time', f'{w.capitalize()} window: response magnitude', r'\tau since pulse', r'|G_1(\tau)| / |G_1(h)|', traces)
        f['value_scale'] = 'log'; figures.append(f)
    for w in ('long', 'short'):
        a, meta = records[w]
        t, g = a['time'], a['causal_response'][:, 0]
        start = meta['first_clean_control_time']
        ix = (t >= start - 1e-9) & (t <= start + 4.32)
        denom = meta['modes'][0]['peak_sensitivity_magnitude']
        traces = [trace('Measured real response', t[ix], g[ix].real / denom, color=BLUE),
                  trace('Re-anchored exponential', t[ix], a['post_clean_prediction'][ix, 0].real / denom,
                        color=GREEN, line_style='dash', width=2.3)]
        f = figure('time', f'{w.capitalize()}: after window clears', r'\tau since pulse', r'Re(G_1) / peak |G_1|', traces)
        f['value_scale'] = 'linear'; figures.append(f)
    return save_figure(f'echo-pulse-{preparation}-20261012', label, figures)


def coordinate_plot():
    """Plot paired-pulse sensitivity of the echo-history latent coordinates."""
    figures = []
    for w in ('long', 'short'):
        a = payload(DATA / f'{w}-free-train-1-p3.response.npz')
        t, delta = a['time'], a['input_amplitudes'][1]
        response = (a['xi'][1] - a['xi'][2]) / (2 * delta)
        traces = [trace(r'\xi_1 response', t, response[:, 0], color=BLUE),
                  trace(r'\xi_2 response', t, response[:, 1], color=ORANGE)]
        f = figure('time', f'{w.capitalize()}: latent-coordinate pulse response', r'\tau since pulse', r'\Delta\xi / (2\delta)', traces)
        f['value_scale'] = 'linear'; figures.append(f)
    a = payload(DATA / 'short-free-train-1-p3.response.npz')
    t = a['time']; ix = t <= .2
    for key, title, ylabel in [('echo', 'Echo response: immediate onset', r'\Delta e / (2\delta)'),
                                ('xi', 'Short-window latents: first 0.2', r'\Delta\xi / (2\delta)')]:
        g = (a[key][1] - a[key][2]) / (2 * a['input_amplitudes'][1])
        traces = [trace('Echo response', t[ix], g[ix], color=BLUE)] if key == 'echo' else [
            trace(r'\xi_1 response', t[ix], g[ix, 0], color=BLUE), trace(r'\xi_2 response', t[ix], g[ix, 1], color=ORANGE)]
        f = figure('time', title, r'\tau since pulse', ylabel, traces)
        f['value_scale'] = 'linear'; figures.append(f)
    return save_figure('echo-pulse-coordinates-20261012', 'Same late-free pulse: echo and latent responses', figures)


def control_plot(target):
    """Compare sine-trained and pulse-trained control maps on the same saved target."""
    figures = []
    ref = payload(INPUT / f'target-{target}.measurements.npz')
    truth = payload(INPUT / f'target-{target}.truth.npz')
    candidates = [(OLD, f'short-{target}-m12-r1.000000', 'Prior sine fit; pair', GREY, 'dash'),
                  (DATA, f'short-dense-pulse-{target}-m1-{target}-r1.000000', 'Dense pulses; fundamental', ORANGE, 'solid'),
                  (DATA, f'short-dense-pulse-{target}-m12-{target}-r1.000000', 'Dense pulses; pair', BLUE, 'solid')]
    for key, title, ylabel in [('input', 'Control signal', 'Dimensionless pressure u'), ('echo', 'Echo', 'Dimensionless echo'),
                                ('r', 'Radius (evaluation only)', r'r / R_0'), ('v', 'Velocity (evaluation only)', r'd(r/R_0)/d\tau')]:
        ix = np.flatnonzero((ref['time'] >= 54) & (ref['time'] <= 58))
        if key == 'input': ix = ix[::20]
        target_y = ref[key] if key in ('input', 'echo') else truth['states'][:, 0 if key == 'r' else 1]
        traces = [trace('Known sine / target', ref['time'][ix], target_y[ix], color=GREEN, line_style='dash', width=2.5)]
        for folder, stem, label, color, style in candidates:
            if 'failure' in read(folder / f'{stem}.json'):
                continue  # Failed trials have no curve; retain them in the report table.
            a = payload(folder / f'{stem}.measurements.npz')
            assert np.array_equal(a['time'], ref['time'])
            y = a[key] if key in ('input', 'echo') else payload(folder / f'{stem}.truth.npz')['states'][:, 0 if key == 'r' else 1]
            traces.append(trace(label, a['time'][ix], y[ix], color=color, line_style=style))
        f = figure('time', f'Short window / {target}: {title}', r'\tau', ylabel, traces)
        f['value_scale'] = 'linear'; figures.append(f)
    return save_figure(f'echo-pulse-control-{target}-20261012', f'Short window / {target}: same KLQR, newly fitted B', figures)


def main():
    """Assemble the saved local evidence and native Viewer plots into the dated HTML report."""
    acq, dense, checks, responses, controls, fits, audit = summarize()
    pulse = {p: pulse_plot(p, label) for p, label in [
        ('free-train-1-p3', 'Late free preparation: fundamental pulse response'),
        ('target-mild-p0', 'Mild target preparation: forcing removed after pulse'),
        ('target-nonlinear-p0', 'Nonlinear target preparation: forcing removed after pulse')]}
    coords, mild, nonlinear = coordinate_plot(), control_plot('mild'), control_plot('nonlinear')
    bankrows = []
    for w in ('long', 'short'):
        for j, mode in enumerate(read(BANK / f'{w}-bank.json')['selected']):
            bankrows.append([w, f'ψ{j+1}', f"{mode['lambda_re']:.6g} + {mode['lambda_im']:.6g}i",
                             pct(mode['late']['pair']), pct(mode['late']['four_tau']), pct(mode['early']['four_tau'])])
    peakrows = []
    for r in responses:
        if r['id'] not in ('free-train-1-p3', 'target-mild-p0', 'target-nonlinear-p0'): continue
        m = r['modes'][0]
        peakrows.append([r['id'], r['window'], f"{m['peak_over_immediate']:,.1f}×", m['time_to_peak'],
                         pct(m['bands'][3]['pulse_semigroup_relative_error']),
                         pct(m['bands'][3]['post_clean_response_semigroup_error']),
                         pct(m['bands'][4]['post_clean_response_semigroup_error'])])
    groups = []
    for w in ('long', 'short'):
        for regime, late in [('free', False), ('free', True), ('mild', None), ('nonlinear', None)]:
            rows = [r for r in responses if r['window'] == w and r['regime'] == regime
                    and (late is None or (r['source_time'] >= 35) == late)]
            ratios = [r['modes'][0]['peak_over_immediate'] for r in rows]
            post = [r['modes'][0]['bands'][3]['post_clean_response_semigroup_error'] for r in rows]
            groups.append([w, regime + ((' late' if late else ' early') if late is not None else ''), len(rows),
                           f'{np.median(ratios):,.1f}×', f'{min(ratios):,.1f}–{max(ratios):,.1f}×', pct(float(np.median(post)))])
    fitrows, fitappendix = [], []
    for f in fits:
        r = f['training']; scores = f['scores']; test = 'mild' if r in ('free', 'pooled') else r
        row = [f['campaign'], f['window'], r, f['subset'], scores['train-'+r]['samples'],
               pct(scores['train-'+r]['mean_mode_relative_residual']), pct(scores['heldout-'+test]['mean_mode_relative_residual']),
               pct(scores['target-diagnostic-'+test]['mean_mode_relative_residual']), f"{f['rank']}/{f['columns']}", f"{f['condition']:,.3g}"]
        fitappendix.append(row)
        if r in ('mild', 'nonlinear') and f['subset'] == '12': fitrows.append(row)
    compare = []
    for w in ('short', 'long'):
        for target in ('mild', 'nonlinear'):
            for subset in ('1', '12'):
                penalty = 1 if w == 'short' else 100
                old = read(OLD / f'{w}-{target}-m{subset}-r{penalty:.6f}.json')
                new = next(r for r in controls if r['id'] == f'{w}-dense-pulse-{target}-m{subset}-{target}-r{penalty:.6f}')
                compare.append([w, target, subset, penalty, pct(old['late']['echo_nrmse']) if 'late' in old else 'FAILED',
                                pct(new['late']['echo_nrmse']), pct(old['late']['driver_nrmse']) if 'late' in old else 'FAILED',
                                pct(new['late']['driver_nrmse'])])
    controlrows = [[r['id'], pct(r['late']['echo_nrmse']), pct(r['late']['driver_nrmse']),
                    pct(r['final_10tau']['echo_nrmse']), pct(r['late']['saturation_fraction'])]
                   if 'late' in r else [r['id'], 'FAILED', r['failure'], '—', '—'] for r in controls]
    checkrows = [[r['id'], r['window'], f"{r['immediate_sensitivity_checks'][0]['refined_relative_change']:.3g}",
                  f"{r['modes'][0]['double_amplitude_sensitivity_change']:.3g}",
                  pct(r['width_checks'][0]['post_pulse_semigroup_relative_error']),
                  pct(r['width_checks'][1]['post_pulse_semigroup_relative_error'])] for r in checks]
    style = (ROOT / 'docs/reports/ECHO_MULTIMODE_JOINT_BPSI_2026-10-11.html').read_text(encoding='utf-8').split('<style>', 1)[1].split('</style>', 1)[0]
    report = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pulse diagnosis: echo-history memory and Bψ</title><style>{style}</style></head><body><main>
<p>DynSysLab · echo-only · 12 October 2026 · completed isolated native C++ diagnostic</p>
<h1>Pulse diagnosis: the input is visible, but its latent response has memory</h1>
<p class="note"><strong>Your pulse proposal was diagnostically useful.</strong> It exposes a failure beyond B regression: after a brief input ends, the measured eigenfunction response does not initially follow its fitted free-motion exponential. In a late, well-modelled free trajectory, the long-window response grows to <strong>7,287 times its first-hold value</strong>. Extrapolating that first response misses essentially all subsequent response. After the window is entirely post-pulse, a re-anchored exponential predicts it to <strong>0.14–0.20%</strong>. A better immediate B cannot, by itself, supply this missing autonomous evolution.</p>
<p>The short window is less affected, not unaffected. Dense pulse calibration retains useful mild control (two eigenfunctions: <strong>6.68% echo error, 11.25% known-driver error</strong>) but does not outperform the prior sine-trained fit (5.13%, 8.63%). Nonlinear control is still unsuccessful. <strong>The primary pipeline was not changed.</strong></p>
<nav><a href="#setup">Experiment</a> · <a href="#test">Decisive test</a> · <a href="#results">Pulse results</a> · <a href="#fits">B fits</a> · <a href="#control">Control</a> · <a href="#meaning">Diagnosis</a> · <a href="#next">Recommendations</a> · <a href="#verify">Checks and code</a></nav>
<section id="setup"><h2>1. What was tested and what remained fixed</h2>
<p>You proposed replacing continuous sine excitation with localized pings, then watching several history windows to separate natural motion from input response. I tested that directly, using matched continuations from exactly the same preparation and echo prehistory: no pulse, positive/negative pulses, and positive/negative half- and double-amplitude pulses.</p>
{table(['Quantity', 'Long window', 'Short window'], [['History samples', '8,644', '101'], ['History span (τ)', '4.3215', '0.05'], ['History span (µs)', '13.018', '0.15062'], ['Echo sample spacing', '0.0005 τ', '0.0005 τ'], ['Control hold h', '0.01 τ = 20 echo samples', '0.01 τ = 20 echo samples'], ['Latent coordinates used', 'ξ₁, ξ₂', 'ξ₁, ξ₂'], ['Frozen bank', '3 complex polynomial eigenfunctions', '3 complex polynomial eigenfunctions']])}
<p>The plant remains viscous air-in-water RPE: equilibrium radius 10 µm, viscosity 0.001 Pa·s, density 998 kg/m³, ambient pressure 100 kPa, surface tension 0.0725 N/m, κ = 4/3. One τ is 3.01238 µs. Input u is pressure divided by the configured 10,997.9 Pa pressure scale, with the existing negative pressure-input sign. The echo remains the existing intrinsic-echo convention; direct scattering of the driving wave is intentionally excluded.</p>
<p>The two frozen U/S bases and the existing degree-five, 21-term polynomial banks were reused. No SVD, new latent coordinates, new eigenfunctions, or eigenvalue selection was performed here. These banks were originally fitted on eight <strong>unforced</strong> echo records over 35–80 τ; six free records were reused for mode selection/validation. They are not an independent blind validation set. The first, second, and third positive-frequency bands give ψ1, ψ2, ψ3. The table recalls their earlier free-motion validation:</p>
{table(['Window', 'Mode', 'λ (per τ)', 'Late 1-hold error', 'Late 4τ error', 'Early 4τ error'], bankrows)}
<p>This is why the late-free fundamental is the cleanest diagnostic: its free motion is already accurate. The higher modes and early/nonlinear motion are not equally trustworthy.</p>
<h3>Two complementary pulse campaigns</h3>
<ul><li><strong>86 long-response preparations:</strong> 18 free (12 calibration, 6 held out), 22 mild and 22 nonlinear local preparations (12 calibration and 10 held out in each regime), plus 12 phases on each target used only for diagnosis. Every branch is followed for 14 τ, more than three long windows. Early and late free preparations are both included.</li>
<li><strong>12,192 dense immediate-response preparations:</strong> 710 points on each of eight mild calibration records and 814 points on each of eight nonlinear records, spaced by 0.04 τ. Six whole records per regime train B (4,260 mild + 4,884 nonlinear points); two records are held out (1,420 + 1,628 points). Target phases never train B. Adjacent windows overlap, so these are not 12,192 statistically independent observations.</li></ul>
<p>Default pulse width is one hold, 0.01 τ (about 30.1 ns). Amplitude δ is 0.01 for free/mild and 0.1 for nonlinear preparations, with ±δ/2 and ±2δ controls in the long-response campaign. Selected preparations also test widths 0.005 and 0.05 τ. These are finite rectangular pulses, <strong>not unit-area Dirac impulses</strong>; fixed-amplitude width comparisons also change pulse area. Hardware bandwidth was not modelled. The pulse durations are numerical probes, not a claim that an acoustic laboratory can deliver them unchanged.</p>
<p>The simulator and archived r,v are used <em>only to acquire matched branches</em>. Analysis and regression receive echo prehistory, echo continuations and known inputs—not states or equations. Exact replay is an idealized capability; a laboratory would need reproducible preparations, reset/repetition, measurement-noise averaging and actual actuator calibration. At the pulse start, any preceding sine forcing is switched off; after the finite pulse, input is zero for the entire continuation.</p></section>
<section id="test"><h2>2. A test that removes B-regression error from the argument</h2>
<p>At each echo sample the native ring buffer can immediately evaluate the frozen coordinates. With samples ordered oldest to newest, the actual projection is:</p>
<div class="eq">ξⱼ,ₖ = σⱼ⁻¹ Σ<sub>a=0</sub><sup>d−1</sup> U<sub>a,j</sub> e<sub>k−d+1+a</sub>, &nbsp; zⱼ,ₖ = ψⱼ(ξ₁,ₖ, ξ₂,ₖ).</div>
<p>Stored singular-vector and singular-value conventions are reused exactly. There is no new normalization, centering, or SVD in this experiment. Online ring-buffer evaluation and the native batch projection are cross-checked. <strong>After initial warm-up, there is no requirement to wait a whole window to read ξ.</strong> The issue is how the new echo affects the coordinate as it moves through that window.</p>
<p>For identical initial histories, define the measured, centered finite-difference response to a one-hold pulse:</p>
<div class="eq">Gⱼ(ℓh) = [zⱼ<sup>+δ</sup>(ℓh) − zⱼ<sup>−δ</sup>(ℓh)] / (2δ), &nbsp; Gⱼ(0) = 0.</div>
<p>This subtraction cancels the shared zero-input baseline and estimates local causal sensitivity about zero current input. It is not an estimate of the entire response to an arbitrarily large input. Half/double amplitudes test the local linear approximation.</p>
<p>If the assumed current model were adequate,</p>
<div class="eq">zₖ₊₁ = Mzₖ + B<sup>d</sup>(zₖ)uₖ, &nbsp; M = diag(exp(λⱼh)),<br>
G(ℓh) ≈ M<sup>ℓ−1</sup>G(h), &nbsp; ℓ ≥ 1.</div>
<p>After the first pulse hold, u is zero, so <strong>no choice of B can repair failure of the second line</strong> while M and the coordinates stay fixed. This test uses the <em>measured</em> first-hold G, not a fitted B. It therefore bypasses sparse coverage, regression conditioning and B parametrization as explanations for that failure.</p>
<p>A second test re-anchors the exponential after every sample in the buffer is post-pulse. The first such control times are 4.34 τ (long) and 0.06 τ (short). This is a diagnostic reset of the prediction, <strong>not a reset of the simulated trajectory and not an implemented controller</strong>.</p></section>
<section id="results"><h2>3. What the pulses reveal</h2>
<h3>Late free trajectory: a strong separation of memory error from eigenfunction error</h3>
<p>The principal example is free training record 1 at τ = 50.11, well into the eigenfunction’s valid late regime. The long-window fundamental reaches its largest response at 3.43 τ after the pulse, with a peak/first-hold ratio of 7,286.9. Its zero-input one-hold error is only about 0.003–0.004% over the later bands. Nevertheless, propagating the pulse’s first response fails at roughly 100% relative error. After 4.34 τ, the re-anchored complex exponential has 0.138% error over 4.34–8.67 τ and 0.203% over 8.67–14 τ.</p>
{pulse['free-train-1-p3']}
<p>Upper panels are log-y magnitudes, normalized by the first-hold response; the green prediction starts only when the window clears. Lower panels compare <em>real parts</em> for 4.32 τ after that point, so matching magnitudes alone cannot hide phase error. The short-window response peaks at 0.05 τ and is 11.7 times its first-hold response; its re-anchored prediction is useful but drifts more over longer horizons.</p>
{coords}
<p>The echo responds immediately, and both latent coordinates begin responding immediately. Their response continues to change after the input is zero. This is the distinction between <em>online availability</em> and <em>having a closed dynamical model for that online quantity</em>.</p>
<h3>Mild target preparation</h3>
<p>Starting from one phase of the mild sine-driven target gives the same qualitative result. The long-window peak/first-hold ratio is 12,187; after the buffer clears, re-anchored errors are 0.21% and 0.41% in the two later bands. The prior sine is not continued: these are pulse-and-release tests from its history.</p>
{pulse['target-mild-p0']}
<h3>Nonlinear target preparation: clearing the window is not enough</h3>
<p>The nonlinear example is different. Even the short-window zero-input continuation initially violates the fitted exponential strongly; the re-anchored pulse response remains poor. Thus memory is <strong>not the only failure</strong>: the late-trained, two-coordinate polynomial bank also lacks adequate validity over this strongly nonlinear region. Clearing the buffer is not a universal cure.</p>
{pulse['target-nonlinear-p0']}
{table(['Preparation', 'Window', 'Peak / first hold', 'Peak τ', 'From first hold: error 4.34–8.67', 'Re-anchor: error 4.34–8.67', 'Re-anchor: error 8.67–14'], peakrows)}
<h3>Across the preparations, not just one selected trace</h3>
<p>The following are descriptive medians/ranges of the fundamental’s already-computed native metrics. “Late free” means preparation time ≥35 τ. Mild/nonlinear groups combine local and target-diagnostic preparations; they are not confidence intervals or independent-sample population estimates.</p>
{table(['Window', 'Preparation group', 'Count', 'Median peak / first', 'Range', 'Median re-anchor error 4.34–8.67'], groups)}
<p>All quoted exponential-response errors use complex values: √[Σ|G−Ĝ|² / Σ|G|²] on the specified band. Near-zero denominators can make relative errors misleading; the principal late-free response is well resolved, and its first hold was independently checked for numerical convergence.</p></section>
<section id="fits"><h2>4. Refitting B from the pulses</h2>
<p>For every nonempty subset of the three bank modes, all selected B rows were jointly refitted. The existing native affine-in-eigenfunctions/conjugates map, column scaling, SVD least squares and ridge 10⁻⁸ were reused:</p>
<div class="eq">Gⱼ(h) ≈ cⱼ₀ + Σ<sub>ℓ∈S</sub> [aⱼℓ Re(zℓ) + bⱼℓ Im(zℓ)], &nbsp; j ∈ S.</div>
<p>There are 3, 5 or 7 real predictor columns for one, two or three complex modes. Each fit uses either free, mild, nonlinear, or pooled preparations. The dense campaign replaces mild/nonlinear sparse rows; its pooled fit combines both dense regimes plus the 12 free points without balancing the regimes. That weighting is important when interpreting its poor results.</p>
<p>The existing fitter accepts sampled transition residuals. To reuse it without altering its implementation, the experiment supplies a synthetic row with input 1 and next value Mz + G(h). Its residual is exactly G(h), up to a checked roundoff bound. This is a regression adapter, <strong>not a claim that a unit-amplitude pulse was measured</strong>. The fitted sampled coefficient is converted to the existing continuous map by dividing each row by (exp(λⱼh)−1)/λⱼ, exactly as before; this is the pre-existing frozen-map-in-a-hold approximation.</p>
<p>The two campaigns produced 98 distinct fit designs plus 14 deliberately repeated free-baseline fits (112 small map artifacts). Held-out records were never included in their own campaign’s training. The dense and sparse splits are different and explicitly recorded: sparse trains records 1,3,5; dense trains 1–6 and holds out 7–8. Do not interpret their validation scores as identical test sets.</p>
{table(['Campaign', 'Window', 'Training', 'Modes', 'Train points', 'Train error', 'Held-out same-regime error', 'Target-phase error', 'Rank', 'Condition'], fitrows)}
<p>These are errors in <strong>causal first-hold sensitivity</strong>, not the previous report’s full transition residuals, so the two error types must not be equated. For dense mild fits, the pair’s target-phase errors are 3.59%/8.03% for long-window ψ1/ψ2 and <strong>0.746%/6.32%</strong> for short-window ψ1/ψ2. Therefore a fairly accurate local B can coexist with the dramatic pulse-propagation failure. Dense nonlinear fits still have about 94–98% training error; some target relative errors become enormous. More points do not validate that map. This does not, by itself, prove that no more flexible history-dependent map exists.</p>
<details><summary>All B-fit scores, including one-, two- and three-mode sets</summary>
<p>For pooled/free fits, the displayed held-out and target scores refer to the mild regime. Full per-mode and cross-regime scores are retained in the local fit JSONs.</p>
{table(['Campaign', 'Window', 'Training', 'Modes', 'Train points', 'Train error', 'Held-out error', 'Target error', 'Rank', 'Condition'], fitappendix)}</details></section>
<section id="control"><h2>5. Did pulse-fitted B improve actual control?</h2>
<p>There were 40 closed-loop tests: 24 sparse-map and 16 dense-map trials; 37 completed and 3 stopped at native Riccati-solver failures. Control tests use ψ1 and ψ1+ψ2. All seven subsets were <em>fitted and scored</em>, but new three-mode control sweeps were not repeated: the prior experiment already showed the third mode’s poor free validity/driver behaviour, and the decisive structural test does not require it.</p>
<p>The controller is unchanged: normalized complex modes are split into real components, Q = I/m, R = ρ/0.01, a continuous CARE is solved at each 0.01 τ hold, u = −K(w−w<sub>target</sub>), and input is clipped at ±0.6. Feedback begins at 6 τ and runs to 60 τ. Both target and live coordinates use the same frozen causal transform. No target driver, state estimate, equations or future live echo enter feedback. The target driver is read only for post-run scoring. Target echo is used to define the reference and its RMS normalization, as before.</p>
<p>For a controlled comparison, ρ = 1 for short windows and 100 for long windows, matching the earlier displayed settings. These are not freshly optimized penalties. Results below use 30–60 τ, with final-10τ values in the appendix.</p>
{table(['Window', 'Target', 'Modes', 'ρ', 'Prior sine: echo', 'Dense pulse: echo', 'Prior sine: driver', 'Dense pulse: driver'], compare)}
{mild}
<p>The short-window mild pair has 6.68% late echo error and 11.25% driver error; final-10τ errors are 6.44% and 10.70%. It tracks and largely recovers the sine, but does not beat the previous pair. Dense sampling modestly improves the sparse pulse pair’s driver result (12.86% → 11.25%). The long-window local maps give almost zero actuation at the fixed penalty, with errors near 100%; this is not an actuator-saturation failure.</p>
{nonlinear}
<p>The nonlinear target is not recovered. The prior sine-fitted short-window pair at ρ = 1 failed its Riccati solve, so it has no curve in this nonlinear comparison; the table explicitly marks that prior failure. The plots include r and velocity <em>only as evaluation truth</em>; neither is available to the controller or B regression. Pooling incompatible mild/nonlinear data can additionally destabilize the numerics: the three new failed trials are retained as failures, not omitted as missing successes. Successful Riccati solves demonstrate algebraic solvability, not correctness of the lifted plant model.</p>
<details><summary>All 40 control attempts, including failures</summary>{table(['Trial', 'Late echo error', 'Late driver error / failure', 'Final-10τ echo', 'Late saturation'], controlrows)}</details></section>
<section id="meaning"><h2>6. What is actually standing in the way?</h2>
<ol><li><strong>Forced histories leave the free-history relationship on which ψ was trained.</strong> For an unforced deterministic trajectory, a past window may be parameterized by its present physical state. Under forcing, the same current r,v can have different past echo/input histories. A pulse creates a mixed pre-/post-pulse window. A polynomial ψ fitted to free histories need not remain a Koopman eigenfunction on these histories.</li>
<li><strong>The current model assigns all unexplained change to the current input.</strong> In reality the projected pulse response continues evolving differently from Mz even after u = 0. Then a residual fit to continuously forced records can absorb correlations with previous inputs and drifting history into an “effective B.” That can be helpful for one periodic signal without being a true instantaneous causal input map.</li>
<li><strong>Pulse subtraction cleans up B estimation but does not add the missing dynamics.</strong> Fitting the true first-hold response may remove the accidental compensation that made a sine-trained effective map useful. That explains why a better causal calibration does not automatically improve closed-loop tracking.</li>
<li><strong>The nonlinear regime has an additional eigenfunction/domain problem.</strong> Early free errors, post-clean nonlinear errors and poor nonlinear B fits all persist. Two coordinates may also be insufficiently injective, but these tests do not separately quantify injectivity or prove its unique contribution. They establish failure of this compressed representation with its fitted generator, not impossibility of echo-only control.</li>
<li><strong>Pure error feedback is not exact reference tracking for a persistently forced target.</strong> If z equals z<sub>target</sub> exactly, the current law returns zero input. A target requiring a nonzero continuing sine cannot be exactly maintained at zero error by that law alone. This is a separate limitation from history memory; no feedforward or controller redesign was added here.</li></ol>
<p><strong>So your explanation was substantially right, with one important correction:</strong> we can see ξ and ψ online and can measure their immediate response. We are not blind for 8,644 samples. But the future response contains memory that is not represented by the assumed free generator plus instantaneous B. Continuous forcing does not make disentanglement mathematically impossible; the present compressed model and calibration objective do not account for it adequately.</p>
<p>The experiment cannot rule out a more appropriate set of echo-history observables, more latent coordinates, or a different dynamical realization. It also does not claim all sufficiently rich Koopman control frameworks fail. It shows that <strong>holding these ξ, ψ and M fixed and only refitting the present B is insufficient for a faithful pulse-response model</strong>, especially with the long window.</p></section>
<section id="next"><h2>7. Recommended next steps, within your framework</h2>
<ol><li><strong>Keep this pulse test as a model-acceptance gate.</strong> Before control tuning, require useful unforced multi-step accuracy and a correct post-pulse response from the measured first-hold increment. Judge complex error and long-time behaviour, not just one-step residuals.</li>
<li><strong>Use the short window and ψ1/ψ2 in the mild regime as the current practical baseline.</strong> Keep the earlier sine-calibrated pair as the better demonstrated tracking controller. Label its B as an effective fitted map, not a validated causal state-Jacobian. Pulse maps are a diagnostic comparator, not a replacement that has earned deployment.</li>
<li><strong>Do not spend the next effort merely adding more sine waves or denser immediate B samples.</strong> The dense campaign has already removed sparse coverage as the sole explanation in the mild regime. Gain sweeps may improve tracking but cannot make G(ℓh) = M<sup>ℓ−1</sup>G(h) true when it is false.</li>
<li><strong>For stronger control, fix representation/closure first.</strong> Within the existing experimental policy, a next separately authorized test would repeat this same gate with four latent coordinates and eigenfunctions trained only on unforced echoes spanning the intended nonlinear domain. More coordinates may retain missing information, but they do not automatically make the currently assumed generator correct on mixed forced histories.</li>
<li><strong>If four coordinates still fail the pulse gate, explicitly decide whether memory belongs in the model.</strong> Possibilities include a state containing relevant history/input-history information or a pulse/delayed-response realization. These are future modelling choices, not changes made here. Waiting for a clean buffer and calling the resulting multi-window effect “instantaneous B” would mislabel the physics and would not fix the current controller.</li>
<li><strong>Only after the dynamical model passes, revisit persistent-reference tracking and laboratory robustness.</strong> Then consider the appropriate reference/feedforward treatment, measurement noise, finite actuator bandwidth and the intentionally omitted direct scattering. No claim of laboratory readiness follows from these noiseless numerical tests.</li></ol></section>
<section id="verify"><h2>8. Numerical checks, boundaries and reproducibility</h2>
{table(['Preparation', 'Window', 'First-hold change on half integrator step', 'Whole-response change at double amplitude', '0.005τ pulse: post-pulse prediction error', '0.05τ pulse: post-pulse prediction error'], checkrows)}
<p>The last two columns extrapolate from the first control-grid point at or after each pulse ends; for a 0.005τ pulse this point is 0.01τ. Longer pulses help the short-window late/mild mismatch somewhat but do not eliminate it. They do not repair the long-window response, and nonlinear failures remain. These are finite-width checks, not a zero-width convergence proof.</p>
<p>The largest immediate fundamental sensitivity change under halving the integration maximum step in the six explicit resolution cases is {audit['maximum_fundamental_immediate_refinement_change']:.3g}. The largest online/batch coordinate discrepancy over the main campaign is {audit['maximum_projection_absolute_error']:.3g}; the free replay relative discrepancy is at most {audit['maximum_free_replay_relative_error']:.3g}. Pulse signs, widths, sample counts and shared initial echoes were checked. Surrogate cancellation, conjugate-map conversion, single-mode compatibility and a synthetic joint cross-mode recovery/permutation test also pass. The four native history/Koopman/KLQR parity CTest targets pass; the installed native doctor passes. Existing unrelated repository audit debt was left untouched.</p>
<p>All new acquisition, projection analysis, fitting and control are C++. Python only verifies archive structure, assembles already-computed metrics and formats native Viewer plots. No scientific algorithm was moved into Python. New code lives beside the core under <code>tools/experiments/echo_pulse*</code>. Only a main guard and a forwarding entry point were added to the earlier isolated multimode experiment, so its existing control implementation can be linked unchanged.</p>
{table(['Source', 'Responsibility'], [['echo_pulse_acquire.cpp', 'Only new module permitted to access archived truth and the native simulator; emits echo-only measurement branches'], ['echo_pulse_analyze.cpp', 'Frozen online/batch projection, paired response, native exponential-error metrics and width/amplitude checks'], ['echo_pulse_dense.cpp', 'Dense echo-only immediate sensitivity; references existing prehistories rather than duplicating them'], ['echo_pulse_fit.cpp', 'Reuse of existing joint B fitter, canonical map storage, and unchanged KLQR trial'], ['echo_pulse.cpp / echo_pulse.hpp', 'Small command dispatcher and experiment contract'], ['report_echo_pulse_diagnosis.py', 'This report, archive checks and native Viewer presentation']])}
<p>Local results: <code>archives/experiments/echo-pulse-diagnosis-20261012-v1</code>. Inputs remain immutable; large scientific data, plots and Viewer payloads remain local/ignored. Failed trial JSONs are retained. This report does not commit or publish anything.</p>
<pre><code>powershell -ExecutionPolicy Bypass -File tools/experiments/Build-EchoPulse.ps1
# Add build/native-suite-fast/bin to this process PATH.
# Use a fresh output directory for a full new acquisition.
build/echo-pulse/echo_pulse.exe . OUTPUT acquire
build/echo-pulse/echo_pulse.exe . OUTPUT long
build/echo-pulse/echo_pulse.exe . OUTPUT short
build/echo-pulse/echo_pulse.exe . OUTPUT fit-long
build/echo-pulse/echo_pulse.exe . OUTPUT fit-short
build/echo-pulse/echo_pulse.exe . OUTPUT control-long
build/echo-pulse/echo_pulse.exe . OUTPUT control-short
build/echo-pulse/echo_pulse.exe . OUTPUT dense-acquire
build/echo-pulse/echo_pulse.exe . OUTPUT dense-fit-long
build/echo-pulse/echo_pulse.exe . OUTPUT dense-fit-short
build/echo-pulse/echo_pulse.exe . OUTPUT dense-control-long
build/echo-pulse/echo_pulse.exe . OUTPUT dense-control-short
build/echo-pulse/echo_pulse.exe . OUTPUT verify
python tools/report_echo_pulse_diagnosis.py</code></pre>
<p>The presentation helper reads the named October 12 output folder. Original comparison: <a href="ECHO_MULTIMODE_JOINT_BPSI_2026-10-11.html">October 11 subset-specific B report</a>. This report’s plots are actual native Viewer renders; click any to inspect the full-size PNG.</p>
<p class="note"><strong>Bottom line:</strong> pulses successfully isolated the problem. There is a genuine echo-history dynamical-closure failure, not just an inability to collect cleaner B-training data. Mild short-window control remains useful; the present fixed representation does not support trustworthy nonlinear control merely by refitting B.</p>
</section></main></body></html>'''
    REPORT.write_text(report, encoding='utf-8')
    print(json.dumps(audit, indent=2))
    print(REPORT)


if __name__ == '__main__':
    main()
