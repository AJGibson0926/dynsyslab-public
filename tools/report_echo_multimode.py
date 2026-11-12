"""Presentation only: native C++ metrics/archives, rendered with native Viewer."""
from __future__ import annotations
import html
import json
from pathlib import Path
import numpy as np
from plot_rpe_echo_jet_pipeline_evidence import figure, render, trace

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'archives/experiments/echo-multimode-20261009-v2'
INPUTS = ROOT / 'archives/experiments/echo-input-regression-20261006-v1'
REPORT = ROOT / 'docs/reports/ECHO_MULTIMODE_CONTROL_2026-10-09.html'
SUBSETS = ('1', '2', '3', '12', '13', '23', '123')
COLORS = ('#007f8b', '#b15a16', '#7857a4')

def read(path):
    """Read the saved JSON evidence used to assemble this local report."""
    return json.loads(Path(path).read_text(encoding='utf-8'))

def pct(value):
    """Format a fractional error as a percentage with two decimal places."""
    return f'{100*value:,.2f}%'

def table(headers, rows):
    """Build an HTML table while escaping all supplied headings and cell text."""
    return '<div class="scroll"><table><thead><tr>'+''.join(f'<th>{html.escape(str(x))}</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join(f'<td>{html.escape(str(x))}</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'

def image(name, alt):
    """Link and embed a local Viewer PNG with escaped alternative text."""
    return f'<a href="../../viewer_outputs/{name}.png"><img src="../../viewer_outputs/{name}.png" alt="{html.escape(alt)}"></a>'

def get(rows, window, target, subset, penalty):
    """Find a saved controller result by window, target, subset, and penalty."""
    return rows[f'{window}-{target}-m{subset}-r{penalty:.6f}']

def comparison(rows, target):
    """Render error comparisons across all tested eigenfunction subsets and penalties."""
    figures=[]
    for key, ylabel in [('echo_nrmse','Echo relative RMS error'), ('driver_nrmse','Driver relative RMS error')]:
        for window in ('long','short'):
            traces=[]
            for penalty,color in zip((100,1,.01),COLORS):
                points=[(i+1,get(rows,window,target,s,penalty)) for i,s in enumerate(SUBSETS)]
                points=[(i,r['late'][key]) for i,r in points if 'late' in r]
                traces.append(trace(f'R multiplier {penalty:g}',np.array([p[0] for p in points]),
                    np.array([p[1] for p in points]),color=color,marker=True,marker_size=6))
            f=figure('time',f'{window.capitalize()} history / {target}', 'Subset ID (key in report)',ylabel,traces)
            f['value_scale']='log';figures.append(f)
    name=f'echo-multimode-{target}-all-subsets-20261009'
    render(name=name,label=f'Multi-eigenfunction control / {target}',page_label='2 latent coordinates; subset 1,2,3,12,13,23,123; late errors 30-60',
           source=DATA/'audit.json',figures=figures,columns=2,height=1100)
    return image(name,f'All seven eigenfunction subsets: {target}')

def details(rows,window,target):
    """Plot the post-hoc best observed echo tracking at each tested eigenfunction count."""
    chosen=[]
    for n in (1,2,3):
        eligible=[r for r in rows.values() if r['id'].startswith(f'{window}-{target}-') and r['modes']==n
                  and r['R_multiplier'] in (100,1,.01) and 'late' in r]
        chosen.append(min(eligible,key=lambda r:r['late']['echo_nrmse']))
    reference=np.load(INPUTS/f'target-{target}.measurements.npz')
    truth=np.load(INPUTS/f'target-{target}.truth.npz')
    figures=[]
    for key,label in [('input','u = p / p0'),('echo','Echo'),('r','r / R0'),('v','dr / dτ')]:
        keep=(reference['time']>=54)&(reference['time']<=58)
        # Display-only decimation: retain the original 0.0005 grid for echo/state,
        # and interval-start samples for zero-order-held input.
        ix=np.flatnonzero(keep)
        if key=='input':ix=ix[::20]
        y=reference[key] if key in ('input','echo') else truth['states'][:,0 if key=='r' else 1]
        traces=[trace('Known driver / target',reference['time'][ix],y[ix],color='#58616b',line_style='dash',width=2.4)]
        for row,color in zip(chosen,COLORS):
            a=np.load(DATA/(row['id']+'.measurements.npz'))
            v=a[key] if key in ('input','echo') else np.load(DATA/(row['id']+'.truth.npz'))['states'][:,0 if key=='r' else 1]
            tag=row['id'].split('-m')[-1].split('-r')[0]
            traces.append(trace(f"{row['modes']} mode(s): {tag}; R={row['R_multiplier']:g}",a['time'][ix],v[ix],color=color))
        f=figure('time',f'{window.capitalize()} / {target}: {label}',r'\tau',label,traces);f['value_scale']='linear';figures.append(f)
    name=f'echo-multimode-{window}-{target}-responses-20261009'
    render(name=name,label=f'{window} / {target} / best observed echo tracking per count',
        page_label='Post-hoc best of tested cases; not a tuned or independently validated selection',source=DATA/'audit.json',figures=figures,columns=2,height=1100)
    return image(name,f'{window} {target}: driver, echo, radius, velocity'),chosen

def main():
    """Assemble the saved local evidence and native Viewer plots into the dated HTML report."""
    audit=read(DATA/'audit.json');rows={r['id']:r for r in audit['results']}
    banks={w:read(DATA/f'{w}-bank.json') for w in ('long','short')}
    eigen=[]
    for w,b in banks.items():
        for j,r in enumerate(b['selected'],1):
            eigen.append([w,f'ψ{j}',f"{r['lambda_re']:.6f} + {r['lambda_im']:.6f}i",
                pct(r['late']['pair']),pct(r['late']['one_tau']),pct(r['late']['four_tau']),pct(r['early']['four_tau'])])
    best=[];details_html=[]
    for t in ('mild','nonlinear'):
        for w in ('long','short'):
            plot,chosen=details(rows,w,t);details_html.append(f'<h3>{w.capitalize()} history · {t}</h3>'+plot)
            for r in chosen:
                best.append([w,t,r['id'].split('-m')[-1].split('-r')[0],r['R_multiplier'],pct(r['late']['echo_nrmse']),
                    pct(r['late']['driver_nrmse']),pct(r['truth_late']['radius_oscillation_nrmse']),pct(r['truth_late']['velocity_nrmse'])])
    matched=[]
    for w,penalty in [('short',1),('long',100)]:
        for n,subset in [(2,'12'),(3,'123')]:
            a=get(rows,w,'mild','1',penalty*n);b=get(rows,w,'mild',subset,penalty)
            matched.append([w,f'ψ1 only, R={penalty*n:g}',f'ψ{subset}, R={penalty:g}',
                pct(a['late']['echo_nrmse']),pct(b['late']['echo_nrmse']),pct(a['late']['driver_nrmse']),pct(b['late']['driver_nrmse'])])
    maps=[]
    for w in ('long','short'):
        for r in read(DATA/f'{w}-maps.json'):
            maps.append([w,r['target'],r['mode'],pct(r['train_residual']),pct(r['heldout_waveform_residual']),pct(r['target_residual']),f"{r['condition']:.2f}"])
    causal=[]
    for case in audit['causal']:
        for r in case['modes']:
            causal.append([case['window'],case['target'],r['mode'],f"{r['zero_input_drift_over_actual_input_effect']:,.2f}×",
                           f"{r['regressed_B_relative_error_to_causal_B']:,.3g}"])
    inverse=[]
    for w in ('long','short'):
        for t in ('mild','nonlinear'):
            for s in ('1','12','123'):
                r=read(DATA/f'{w}-{t}-m{s}.inverse.json')
                inverse.append([w,t,s,pct(r['driver_nrmse']),pct(r['unachievable_residual_fraction']),pct(r['known_driver_model_residual'])])
    all_rows=[]
    for r in sorted(rows.values(),key=lambda r:r['id']):
        if 'failure' in r:all_rows.append([r['id'],'FAILED',r['failure'],'—','—'])
        else:all_rows.append([r['id'],pct(r['late']['echo_nrmse']),pct(r['late']['driver_nrmse']),pct(r['late']['saturation_fraction']),f"{r['maximum_care_relative_residual']:.2g}"])
    plots=comparison(rows,'mild')+comparison(rows,'nonlinear')
    report=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Echo-only control: one, two and three eigenfunctions</title>
<style>body{{font:17px/1.65 system-ui,sans-serif;color:#253344;background:#f3f6f8;margin:0}}main{{max-width:1280px;margin:auto;padding:35px 25px}}h1{{font-size:36px;line-height:1.25}}section{{background:white;border:1px solid #d9e2e9;border-radius:10px;margin:24px 0;padding:25px}}h2{{margin-top:0}}table{{border-collapse:collapse;font-size:14px;width:100%}}td,th{{padding:9px;border-bottom:1px solid #dce4ea;text-align:left}}th{{background:#e8f2f5}}.scroll{{overflow:auto}}.note{{border-left:5px solid #b77515;padding:16px 20px;background:#fff8e9}}.eq{{background:#eef5f7;padding:15px;font:19px/1.6 Georgia,serif}}img{{width:100%;height:auto}}a{{color:#066778}}code,pre{{font:14px/1.6 ui-monospace,monospace;white-space:pre-wrap;overflow-wrap:anywhere}}summary{{cursor:pointer}}@media print{{body{{background:white}}section{{border:0}}}}</style></head><body><main>
<p>DynSysLab · echo-only · 9 October 2026 · isolated native C++ experiment</p>
<h1>Do more eigenfunctions improve echo-history control?</h1>
<p class="note"><strong>Not in this controlled comparison.</strong> The strongest result remains one fundamental eigenfunction with the short window on the mild target: 8.50% echo error and 9.58% driver error. Adding modes does not reliably improve either window and does not rescue the nonlinear target. The reasons are partly identifiable: cost-weight redistribution, poor forced-response models for the added modes, free-to-forced history mismatch, and—in aggressive cases—numerical failure or saturation. This is not a proof that every multi-eigenfunction controller must fail.</p>
<nav><a href="#design">Design</a> · <a href="#free">Eigenfunctions</a> · <a href="#results">Results</a> · <a href="#weights">Weight check</a> · <a href="#models">Model diagnosis</a> · <a href="#responses">Responses</a> · <a href="#limits">Limits / next step</a> · <a href="#reproduce">Reproduce</a></nav>

<section id="design"><h2>1. Exactly what changed</h2>
<p><strong>Two latent coordinates in every run.</strong> The long and short history bases are the previously archived U matrices; no new Hankel SVD was performed. The long window has 8644 echo samples spanning 4.3215 τ (13.018 μs); the short window has 101 spanning 0.05 τ (0.15062 μs). Sensor spacing is 0.0005 τ, control hold 0.01 τ, with τ = 3.012383 μs. The frozen bases store four columns, but only ξ₁ and ξ₂ enter every eigenfunction.</p>
<p>For each window, one degree-five EDMD operator was fitted to the same eight <em>unforced</em> training records over 35–80 τ. Its dictionary contains 21 terms. Three stable, positive-frequency modes were selected, one from each predeclared band [4,8), [8,16), [16,24) rad/τ. Selection minimized the sum of 1- and 4-τ free-trajectory errors over 35–80 τ on six existing validation records. Those records were also used in previous work, so these are reused selection/validation data, not a fresh blind test. All 21 candidates were preserved this time.</p>
<p>We call the selected fundamental, second-frequency and third-frequency candidates ψ1, ψ2 and ψ3. These labels are local to each window; they do not denote latent coordinates. The seven tested subsets were ψ1, ψ2, ψ3, ψ1+ψ2, ψ1+ψ3, ψ2+ψ3, and ψ1+ψ2+ψ3. Negative-frequency conjugates were not double-counted.</p>
<p>This is a <strong>new matched mode-bank comparison</strong>, not an exact rerun of the earlier degree-one single-mode baseline. Both windows now use the same polynomial degree and training interval. The old models and primary pipeline remain unchanged.</p>
<h3>Same model family and controller</h3>
<div class="eq">zⱼ = ψⱼ(ξ₁, ξ₂), &nbsp; μⱼ = exp(λⱼh)<br>zⱼ,ₖ₊₁ − μⱼzⱼ,ₖ ≈ [cⱼ₀ + cⱼ₁ Re(zⱼ,ₖ) + cⱼ₂ Im(zⱼ,ₖ)]uₖ.</div>
<p>Each row of B was fitted using the existing native history-input fitter, on six local sine calibration recordings (records 1–6); two complete waveforms (7–8) were withheld. The target waveform itself was not used to regress B. Mild calibration amplitudes were 0.004, 0.008, 0.016, with 0.025 withheld; nonlinear calibration amplitudes were 0.18, 0.24, 0.36, with 0.42 withheld. Each amplitude has two phases. The mild target is amplitude 0.012 at 0.93 cycles/τ; nonlinear is 0.3 at 0.83 cycles/τ. B sees only measured echo, its history and recorded input. Eigenfunctions see only unforced training data.</p>
<p>The same independently fitted B rows are retained across subsets. They are assembled into the standard conjugate-complete bilinear control-map class with diagonal predictor blocks. <strong>This isolates adding feedback eigenfunctions; it does not test a richer cross-mode B regression.</strong> Each added mode supplies another output row, not new predictor terms to an existing row.</p>
<div class="eq">wⱼ = zⱼ / RMS(zⱼ,target), &nbsp; Q = I/m, &nbsp; R = ρ/0.1²,<br>uₖ = clip[−K(wₖ − wtarget,ₖ), −0.6, 0.6].</div>
<p>The existing continuous-time CARE solver is called each control hold; sampled B is converted using Bᶜⱼ = Bᵈⱼ / [(μⱼ−1)/λⱼ]. This is the prior frozen-map conversion approximation, not a claim of exact sampled-data LQR. The target is projected from its echo; the known target driver is excluded from feedback and used only for scoring. There is no feedforward term. Runs start at r=1, v=0, turn feedback on at 6 τ and end at 60 τ.</p>
<p>The main factorial sweep comprises 2 windows × 2 targets × 7 subsets × 3 penalties (ρ=100,1,0.01) = <strong>84 attempts</strong>. Sixteen supplementary single-mode attempts isolate weight redistribution. Total: <strong>{audit['completed']} completed, {audit['failed']} failed</strong>. Errors below use 30–60 τ unless noted; final 50–60 τ metrics are also saved. This is finite-horizon late-time evidence, not proof of indefinite stability.</p></section>

<section id="free"><h2>2. How accurate are the eigenfunctions before control?</h2>
{table(['Window','Mode','λ (per τ)','Late one-step','Late 1 τ','Late 4 τ','Early 4 τ (5–15)'],eigen)}
<p>Errors are means of per-trajectory relative L² errors against exp(λΔt)ψ(t), computed by the existing native semigroup validator. Late means 35–80 τ. Long ψ1 is very accurate late, but this does not establish its validity on forced histories. Short ψ3 already fails the late long-horizon test. Every selected function is substantially less accurate during early nonlinear free evolution. Adding these functions is not adding three equally trustworthy coordinates.</p></section>

<section id="results"><h2>3. Closed-loop results</h2>
<p>The table selects the lowest observed late echo error for each mode count from the <em>main</em> grid. This is an optimistic, post-hoc summary—not controller selection with independent validation. R settings and actual mode subsets are shown so improved tuning is not confused with dimensionality.</p>
{table(['Window','Target','Modes included','ρ','Echo error','Driver error','Radius oscillation error','Velocity error'],best)}
<p>Radius error is normalized by the target oscillation rtarget−1, not the equilibrium radius. Roughly 100% echo and driver error in the nonlinear cases is approximately the no-control baseline, not successful tracking. Small differences around 100% should not be interpreted as a meaningful improvement.</p>
<h3>All tested subsets and penalties</h3><p>Subset IDs: 1 → ψ1; 2 → ψ2; 3 → ψ3; 4 → ψ1+ψ2; 5 → ψ1+ψ3; 6 → ψ2+ψ3; 7 → ψ1+ψ2+ψ3. Vertical axes are logarithmic, showing relative error (1 = 100%). Failed runs are absent from the curves and listed in the complete table below. Connecting lines guide the eye across discrete subsets.</p>{plots}
<p>The short mild pair ψ1+ψ2 at ρ=1 has 14.43% echo error but <strong>187.44% driver error</strong>, versus 8.50% and 9.58% for ψ1 alone. Thus a moderately plausible response does not mean the sinusoidal input has been rediscovered. At ρ=0.01, short ψ1 alone has only 10.90% echo error but 5784% driver error and about 48% saturation: high-frequency control activity can hide behind a comparatively smooth response.</p></section>

<section id="weights"><h2>4. Is it merely weaker weighting of ψ1?</h2>
<p>With Q=I/m, adding modes divides the fundamental's cost weight by m. A one-mode controller with penalty mρ has the same fundamental cost/input ratio as an m-mode controller with penalty ρ. The following supplementary cases distinguish dilution from what the extra modes do.</p>
{table(['Window','Weight-matched single','Multiple','Single echo','Multi echo','Single driver','Multi driver'],matched)}
<p><strong>Long mild:</strong> ψ1 at ρ=200 and ψ1+ψ2 at ρ=100 both give approximately 73.8% echo error. Most of this apparent worsening relative to ψ1 at ρ=100 is just cost redistribution; it is not convincing evidence that ψ2 independently destroys the response.</p>
<p><strong>Short mild:</strong> the weight-matched single gives 11.93% echo and 12.89% driver error, while the pair gives 14.43% and 187.44%. Here the extra mode is doing something genuinely harmful to driver recovery beyond cost dilution. For the triple at ρ=1, even the weight-matched single is much better. The evidence does not support a monotonic “more modes is better” rule.</p></section>

<section id="models"><h2>5. What is wrong with the models?</h2>
<h3>Forced-response validation</h3>
{table(['Window','Target','Mode','Calibration residual','Withheld-waveform residual','Target residual','Scaled design condition'],maps)}
<p>These residuals are ‖q−B̂ᵈu‖/‖q‖, where q=znext−μz. A value near 100% means that the B fit explains little of that residual. Higher-frequency rows are poor even for the mild case. Extremely large nonlinear target errors show serious failure outside the fitted free-history behavior; a small eigenfunction free one-step error is not a certificate of controlled-model accuracy.</p>
<h3>Reusing the matched-history causal probes</h3>
<p>At 12 previously recorded phases per target, the echo prehistory is identical and only the next held input varies. No state or equation is supplied to this analysis. The zero-input residual a=znext(0)−μz and local input sensitivity are separately measurable in these simulator-generated echo probes. These are diagnostic interventions along the targets, not independent validation or regression training.</p>
{table(['Window','Target','Mode','‖a‖ / ‖actual input effect‖','Relative error B̂ vs causal B'],causal)}
<p>The original structural problem survives adding modes: free-trained history eigenfunctions have nonzero autonomous residual on forced histories, while regression attributes the total residual to current input. The newly added rows also suffer this mismatch. A larger feedback vector cannot make those incorrect dynamics identities true.</p>
<h3>Can one scalar input satisfy all the target mode increments?</h3>
<p>Using target echoes only, we solve the scalar least-squares inverse u*=Re(B̂* q)/‖B̂‖² after the same per-mode normalization. The remaining perpendicular residual measures inconsistency with the fitted input direction, not a formal physical uncontrollability certificate.</p>
{table(['Window','Target','Subset','Inverse driver error','Residual not explained by any scalar u','Known-driver model residual'],inverse)}
<p>For short/mild, adding ψ2 need not produce an enormous inverse-input error, but it raises the unexplained lifted residual. In closed loop the modes compete through one actuator and their incorrect dynamics/map estimates; this can inject unwanted frequencies. The present study diagnoses incompatibility but does not uniquely apportion it among projection noninjectivity, omitted history/input dependence, polynomial extrapolation and regression-family error.</p></section>

<section id="responses"><h2>6. Actual responses, not just scores</h2><p>Native Viewer tiles show driver, echo, radius and velocity over 54–58 τ. Each mode count uses its best observed echo-error case from the main grid; differing penalties are explicit in legends. Simulator radius and velocity are evaluation-only ground truth and were never supplied to regression or feedback.</p>{''.join(details_html)}</section>

<section id="limits"><h2>7. What this does—and does not—establish</h2>
<ul><li>For these two-coordinate representations, these candidates and these maps, adding eigenfunctions gives no reliable improvement and can materially worsen input recovery.</li>
<li>Long-history late free ψ1 is much more accurate than its control performance suggests. The failure is not simply “all eigenfunctions are bad”; the forced-history extension and B interpretation matter.</li>
<li>Harmonic-like functions are additional features of the same ξ₁,ξ₂, not additional measurements. They cannot separate two echo histories that already project to identical ξ. A second-frequency function alone can also introduce phase ambiguity.</li>
<li>Only three frequency bands, one fixed degree-five dictionary and one late training interval were tested. Other modes, real decay modes, more global training, different weights or a full cross-mode B fit could behave differently.</li>
<li>The experiment deliberately held each B row fixed across subsets. It does not answer whether allowing Bⱼ to depend on <em>all</em> selected ψ values improves approximation. That is a separate model-capacity factor, not yet tested.</li>
<li>Seven aggressive short-window cases failed in the unchanged CARE implementation (“Riccati solution is not real”). Coordinate normalization did not remove all such failures. Those runs do not establish physical instability; they are recorded as numerical failures, not assigned fabricated tracking scores.</li>
<li>Three completed short/nonlinear runs at ρ=0.01 (ψ1, ψ1+ψ2, ψ1+ψ3) also reached elevated CARE residuals, up to 1.42×10⁻⁴. Their near-100% tracking errors are retained transparently, but these are not numerically certified controller solutions. No favorable nonlinear conclusion rests on them.</li>
<li>The first unnormalized-cost sweep is retained as v1. The main reported v2 uses the algebraically equivalent normalized formulation. Some saturated nonlinear closed-loop runs are sensitive to floating-point differences. Good mild runs agree closely; numerical equivalence is separately tested.</li>
<li>No feedforward was added. Because u=−K error vanishes at zero tracking error, exact sustained tracking of a nonzero-input target is not generally an invariant solution of this controller. This limitation was held fixed across mode counts.</li></ul>
<p class="note"><strong>Recommendation:</strong> do not promote the multi-mode banks into the primary pipeline. Keep short-window ψ1 as the mild-case reference. Before testing four latent coordinates, the clean next two-coordinate comparison is the same modes with a jointly regressed, conjugate-complete B map versus the fixed-row maps here, retaining causal and held-out validation. That would separate “more feedback coordinates” from “more capacity to model actuation.” If clean forced-response validation still fails, adding more controller modes is not the appropriate fix. No such extension was implemented in this run.</p></section>

<section id="reproduce"><h2>8. Verification and reproducibility</h2>
<p>All scientific computation, fitting, integration, errors and audits were compiled C++; Python only assembles tables and existing Viewer display payloads. Existing primary source modules were not edited. New code is confined to tools/experiments/echo_multimode* and its standalone build helper. Source measurements and prior archives were not rewritten. Local scientific data are Git-ignored; no commit or publish was performed.</p>
<pre>{html.escape(json.dumps(audit['checks'],indent=2))}</pre>
<p>The duplication test confirms that repeating an identical mode twice or three times while averaging its cost does not artificially strengthen feedback. The scaling test confirms equivalence of normalized coordinates and inverse-square cost weighting on well-conditioned test cases. Every saved completed trial's echo/driver/saturation metrics were recomputed from its recordings. Online/offline projection and canonical map assembly were checked during each target/window experiment.</p>
<p>All four focused existing native test gates passed: history Koopman, history input regression, echo-history workflow, and KLQR case-study parity. The native doctor passed. The pre-existing architecture audit still reports unrelated debt; this experiment did not attempt to repair it.</p>
<p>Evidence: <a href="../../archives/experiments/echo-multimode-20261009-v2/audit.json">complete machine-readable audit</a>; <a href="../../archives/experiments/echo-multimode-20261009-v2/long-bank.json">long mode bank</a>; <a href="../../archives/experiments/echo-multimode-20261009-v2/short-bank.json">short mode bank</a>. Each subset has a standard polynomial KROM and saved canonical bilinear map. The custom map name is loaded with the existing loadApproximateControlMap(path, "control_fitted", modes, 1) API. Basis provenance points to the prior archive; it is not duplicated.</p>
<pre>powershell -File tools/experiments/Build-EchoMultimode.ps1
# Put build/native-suite-fast/bin on PATH. Use a NEW output directory.
build/echo-multimode/echo_multimode.exe . NEW_OUTPUT long
build/echo-multimode/echo_multimode.exe . NEW_OUTPUT short
build/echo-multimode/echo_multimode.exe . NEW_OUTPUT --audit
python tools/report_echo_multimode.py  # presentation of the recorded v2 experiment</pre>
<details><summary>All 100 attempts, including failures and saturation</summary>{table(['Attempt','Late echo error','Late driver error / failure','Saturation','Maximum CARE residual'],all_rows)}</details>
</section></main></body></html>'''
    REPORT.write_text(report,encoding='utf-8')
    print(REPORT)

if __name__=='__main__':main()
