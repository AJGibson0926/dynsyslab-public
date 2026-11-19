"""Run a DynSysLab pipeline and/or build its per-stage HTML report.

    python tools/pipeline_report.py run JOB.json [--param N=V]... [--set P=V]... [--until S]
    python tools/pipeline_report.py report MANIFEST.pipeline.json [...]

A thin client: computation stays in the native CLI; plots come from the native
Viewer CLI (DynSysLabViewerCli render). The report has a dependency graph of
the whole job, then one section per stage with that stage bolded in the graph,
a method table, measured results and the stage's Viewer plots. Output is one
self-contained HTML file (images embedded) under viewer_outputs/pipeline-reports/.
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "build" / "native-suite-fast" / "bin"
REF = re.compile(r"\$\{([A-Za-z][A-Za-z0-9_-]*)\.(?:archive|manifest|archives\[\d+\]|output:[A-Za-z0-9_]+)\}")
VIEWS = {  # preferred Viewer routes per stage role, first available wins
    "training": ["time", "phase"], "identification": ["control", "time"],
    "hdmd": ["koopman_spectrum"], "coordinates": ["time", "observable_time"],
    "sindy": ["koopman_spectrum", "krom_dashboard"], "krom": ["krom_dashboard", "koopman_spectrum"],
    "controller": ["control"], "input_map": ["control_map_error", "koopman_spectrum"],
}


# ---------------------------------------------------------------- classification
def coordinate_text(settings: dict) -> str:
    """Short name of a coordinate change, chains joined by arrows."""
    kind = settings.get("kind", "identity")
    if kind == "chain":
        return "(r, v) -> " + " -> ".join(coordinate_text(s) for s in settings.get("steps", []))
    if kind == "identity":
        return "state (r, v)"
    if kind == "observable":
        return {"echo_jet2": "(e, e')"}.get(settings.get("observable"), str(settings.get("observable")))
    if kind == "history":
        rank = settings.get("coordinate_rank", settings.get("hankel_rank", 4))
        source = settings.get("measurement") or settings.get("observable") or "measurement"
        return f"xi_1..xi_{rank} of {'e' if source == 'far_field_acoustic_echo' else source}"
    if kind == "channels":
        return ", ".join(str(c) for c in settings.get("channels", []))
    return kind


def role(job: dict) -> tuple[str, str, str]:
    """(role, layer label, short description) for a stage job."""
    kind = job.get("kind")
    if kind == "solve":
        ctrl = (job.get("controller") or {}).get("kind", "none")
        if ctrl in ("klqr", "kmpc"):
            rep = job["controller"]["settings"].get("control_representation", "exact")
            return "controller", "controller", f"{ctrl.upper()} ({rep} B) + target"
        if ctrl != "none":
            return "identification", "solver", f"forced run ({ctrl})"
        return "training", "solver", "training data"
    methods = [m["name"] for m in job.get("methods", [])]
    names = {"hdmd": "H-DMD", "sindy": "SINDy", "krom_extract": "K-ROM",
             "history_basis": "Hankel-SVD basis", "history_projection": "xi projection",
             "history_control_regression": "input-map regression"}
    text = " + ".join(names.get(m, m) for m in methods)
    if "krom_extract" in methods:
        return "krom", "extracted data", text
    if "history_control_regression" in methods:
        return "input_map", "identification", text
    if "coordinate_transform" in methods:
        settings = next(m["settings"] for m in job["methods"] if m["name"] == "coordinate_transform")
        return "coordinates", "coordinate change", coordinate_text(settings)
    if set(methods) & {"history_basis", "history_projection"}:
        return "coordinates", "coordinate change", text
    if "sindy" in methods:
        s = next(m["settings"] for m in job["methods"] if m["name"] == "sindy")
        sources = json.dumps(job.get("sources", []))
        coords = ("xi (history coordinates)" if "history_projection" in sources or "history-projection" in sources
                  else "coordinates" if "coordinate-transform" in sources
                  else s.get("coordinate_observable", "state"))
        return "sindy", "post-processing", f"{text} in {coords}"
    return "hdmd" if "hdmd" in methods else "post", "post-processing", text


def edges(stages: list[dict]) -> list[tuple[str, str]]:
    out = []
    for stage in stages:
        for dep in sorted(set(REF.findall(json.dumps(stage["job"])))):
            out.append((dep, stage["id"]))
    return out


# ---------------------------------------------------------------- dependency graph
def graph_svg(stages: list[dict], links, active: str | None, status: dict[str, str]) -> str:
    ids = [s["id"] for s in stages]
    depth = {i: 0 for i in ids}
    for _ in ids:  # longest-path layering (stages are already topologically ordered)
        for a, b in links:
            depth[b] = max(depth[b], depth[a] + 1)
    layers: dict[int, list[str]] = {}
    for i in ids:
        layers.setdefault(depth[i], []).append(i)
    w, h, gx, gy = 190, 58, 46, 22
    skips = [(a, b) for a, b in links if depth[b] - depth[a] > 1]  # routed over the top
    lane, top = 9, 20 + 9 * len(skips)
    pos = {}
    for d, members in layers.items():
        for k, i in enumerate(members):
            pos[i] = (20 + d * (w + gx), top + k * (h + gy))
    width = 40 + (max(layers) + 1) * (w + gx) - gx
    height = top + 20 + max(len(m) for m in layers.values()) * (h + gy) - gy
    info = {s["id"]: role(s["job"]) for s in stages}
    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" role="img" '
             f'style="max-width:{width}px"><defs><marker id="a" viewBox="0 0 10 10" refX="9" refY="5" '
             'markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0L10,5L0,10z" '
             'fill="var(--muted)"/></marker></defs>']
    edge = 'fill="none" stroke="var(--muted)" stroke-width="1.4" marker-end="url(#a)"'
    for a, b in links:
        (x1, y1), (x2, y2) = pos[a], pos[b]
        if (a, b) not in skips:
            parts.append(f'<line x1="{x1 + w}" y1="{y1 + h / 2}" x2="{x2 - 2}" y2="{y2 + h / 2}" {edge}/>')
            continue
        # Skip-layer edge: leave above the adjacent-edge port, rise in the gap column
        # to its own lane above every box, cross, and descend in the gap before the target.
        k = skips.index((a, b))
        ly, off = 14 + k * lane, 6 + 4 * k
        sy, ty = y1 + h / 2 - 14, y2 + h / 2 - 14
        gx1, gx2 = x1 + w + gx / 2 - off, x2 - gx / 2 + off
        parts.append(f'<path d="M{x1 + w},{sy} H{gx1} V{ly} H{gx2} V{ty} H{x2 - 2}" {edge}/>')
    for i in ids:
        x, y = pos[i]
        _, layer, text = info[i]
        on = i == active
        state = status.get(i, "")
        stroke = "var(--bad)" if state == "failed" else ("var(--accent)" if on else "var(--line)")
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="var(--card)" '
                     f'stroke="{stroke}" stroke-width="{3.2 if on else 1.2}"/>')
        weight = "700" if on else "400"
        parts.append(f'<text x="{x + 10}" y="{y + 19}" font-size="11.5" fill="var(--muted)" '
                     f'font-weight="{weight}">({html.escape(layer)}) {html.escape(i)}'
                     f'{" — " + state if state in ("failed", "reused", "not run") else ""}</text>')
        label = text if len(text) <= 30 else text[:29] + "…"
        parts.append(f'<text x="{x + 10}" y="{y + 40}" font-size="13" fill="var(--ink)" '
                     f'font-weight="{weight}">{html.escape(label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------- stage results
def load(path: str):
    return np.load(path, allow_pickle=False)


def meta(z) -> dict:
    try:
        return json.loads(str(z["metadata_json"]))
    except Exception:
        return {}


def fmt(x, digits=4) -> str:
    if isinstance(x, complex) or (hasattr(x, "imag") and np.iscomplexobj(x)):
        return f"{x.real:.{digits}g} {'+' if x.imag >= 0 else '−'} {abs(x.imag):.{digits}g}i"
    return f"{x:.{digits}g}" if isinstance(x, (float, np.floating)) else str(x)


def results(role_name: str, job: dict, archives: list[str]) -> list[tuple[str, str]]:
    if not archives:
        return []
    z = load(archives[0])
    rows: list[tuple[str, str]] = []
    if role_name in ("training", "identification", "controller"):
        t, x = z["time"], z["values"]
        rows += [("samples", f"{len(t)} over τ = {t[0]:g}…{t[-1]:g} (Δτ = {t[1] - t[0]:.3g})"),
                 ("r range", f"{x[:, 0].min():.4f} … {x[:, 0].max():.4f}"),
                 ("dr/dτ range", f"{x[:, 1].min():.3f} … {x[:, 1].max():.3f}")]
        if "control_values" in z.files:
            u = z["control_values"].astype(float).ravel()
            rows.append(("control u (p₀ units)", f"{u.min():.3f} … {u.max():.3f}  (|u|max = {abs(u).max() * 100:.0f} kPa)"))
        if role_name == "controller" and "control_target_values" in z.files:
            ref = z["control_target_values"]
            if ref.shape == x.shape:
                e = x[:, 0] - ref[:, 0]
                for lo, hi in [(0, 10), (10, t[-1])]:
                    m = (t >= lo) & (t <= hi)
                    rows.append((f"r tracking error τ {lo:g}–{hi:g}",
                                 f"RMS {np.sqrt(np.mean(e[m] ** 2)):.5f}, max {np.abs(e[m]).max():.5f}"))
                rows.append(("target r range", f"{ref[:, 0].min():.4f} … {ref[:, 0].max():.4f}"))
    elif role_name == "hdmd":
        md, lam = meta(z), z["continuous_eigenvalues"]
        rows += [("Hankel", f"{md.get('matrix_rows')} × {md.get('matrix_columns')}, Δτ = {md.get('delta_t'):.3g}, "
                  f"sample stride {md.get('sample_stride')}"),
                 ("rank retained", f"{md.get('rank')} (full spectrum: {md.get('full_spectrum_retained')})")]
        if "eigenfunction_values" in z.files and "reconstruction_matrix" in z.files:
            w = np.abs(z["reconstruction_matrix"][0] * z["eigenfunction_values"][0])
            order = np.argsort(-w)[:8]
            rows.append(("spectral measure |c_j φ_j(0)|, top 8",
                         "<br>".join(f"λ = {fmt(lam[j])}  →  {w[j]:.4g}" for j in order)))
            rows.append(("weight share of top 20 modes", f"{w[np.argsort(-w)[:20]].sum() / w.sum():.1%}"))
    elif role_name in ("sindy", "krom"):
        lam = z["continuous_eigenvalues"] if "continuous_eigenvalues" in z.files else z.get("spectrum_rates")
        md = meta(z)
        if lam is not None:
            rows.append(("modes", f"{len(lam)}"))
            rows.append(("leading eigenvalues λ", ", ".join(fmt(v) for v in lam[: min(6, len(lam))])))
        if "training_residual_rms" in z.files:
            r = z["training_residual_rms"]
            rows.append(("generator residual RMS (leading 6 | median)",
                         ", ".join(f"{v:.4g}" for v in r[:6]) + f" | {np.median(r):.4g}"))
        for key in ("coordinate_provider", "gamma_method", "training_samples", "polynomial_degree", "selected_modes"):
            if key in md:
                rows.append((key.replace("_", " "), str(md[key])))
    elif role_name == "coordinates":
        md = meta(z)
        if "singular_values" in z.files:
            s = z["singular_values"]
            rows += [("Hankel window", f"{z['left_singular_vectors'].shape[0]} delays"),
                     ("singular values σ₁…", ", ".join(f"{v:.4g}" for v in s[:8])),
                     ("cumulative energy", ", ".join(f"{v:.3f}" for v in z["cumulative_energy"][:8]))]
        elif "values" in z.files:
            q = z["values"]
            rows += [("ξ samples", f"{q.shape[0]} × {q.shape[1]} (after {md.get('warmup_samples', '?')} warm-up samples)"),
                     ("ξ ranges", ", ".join(f"[{a:.3f}, {b:.3f}]" for a, b in zip(q.min(0), q.max(0))))]
    elif role_name == "input_map":
        md = meta(z)
        rows += [("regression pairs", str(md.get("regression_pairs"))),
                 ("relative residual", f"{md.get('relative_residual', float('nan')):.4f}")]
        if "control_constant_matrix" in z.files:
            rows.append(("constant B (per mode)", ", ".join(fmt(v) for v in z["control_constant_matrix"].ravel()[:4])))
    return rows


def method_rows(job: dict) -> list[tuple[str, str]]:
    def compact(value) -> str:
        text = json.dumps(value)
        return text if len(text) <= 160 else text[:157] + "…"
    if job.get("kind") == "solve":
        rows = [("system", compact(job.get("system", {}).get("name"))), ("time", compact(job.get("time"))),
                ("initial condition", compact(job.get("initial_conditions", {}).get("values"))),
                ("solver", compact(job.get("solver")))]
        if job.get("controller"):
            rows += [("controller", compact(job["controller"].get("kind"))),
                     ("controller settings", compact({k: v for k, v in job["controller"].get("settings", {}).items()
                                                       if k != "krom"})),
                     ("target", compact(job["controller"].get("target")))]
        return rows
    return [(m["name"], compact({k: v for k, v in m.get("settings", {}).items()
                                 if not str(v).startswith(("C:", "archives"))})) for m in job.get("methods", [])]


# ---------------------------------------------------------------- Viewer plots
def viewer_png(archive: str, role_name: str, out: Path) -> tuple[str | None, str]:
    cli = str(BIN / "DynSysLabViewerCli.exe")
    try:
        panels = subprocess.run([cli, "panels", "--source", archive, "--archives", "archives"],
                                cwd=ROOT, capture_output=True, text=True, timeout=300, stdin=subprocess.DEVNULL)
        routes = {p["route"] for p in json.loads(panels.stdout).get("panels", [])}
    except Exception as error:  # Viewer unavailable or archive not in its catalog
        return None, f"Viewer could not list views: {error}"
    route = next((r for r in VIEWS.get(role_name, []) if r in routes), None)
    command = [cli, "render", "--source", archive, "--archives", "archives", "--output", str(out),
               "--width", "1500", "--height", "1000"] + (["--view", route] if route else [])
    if role_name == "sindy":  # show the first oscillatory eigenfunction (what the KROM rule selects)
        z = load(archive)
        lam = z["continuous_eigenvalues"] if "continuous_eigenvalues" in z.files else z["spectrum_rates"]
        mode = next((k for k, v in enumerate(lam) if abs(v.imag) > 1e-9), 0)
        command += ["--modes", str(mode)]
        route = f"{route}, mode {mode}"
    run = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=600, stdin=subprocess.DEVNULL)
    if run.returncode or not out.exists():
        return None, f"Viewer render failed: {(run.stderr or run.stdout)[-300:]}"
    return base64.b64encode(out.read_bytes()).decode(), route or "default"


# ---------------------------------------------------------------- report
CSS = """
:root{--bg:#f6f7f9;--card:#fff;--ink:#16202b;--muted:#5b6876;--line:#c9d1da;--accent:#0f766e;--bad:#b42318}
@media (prefers-color-scheme:dark){:root{--bg:#11161c;--card:#19212a;--ink:#e6ebf0;--muted:#9aa7b4;--line:#33404d;--accent:#2dd4bf;--bad:#f97066}}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,Segoe UI,sans-serif}
main{max-width:1180px;margin:0 auto;padding:24px 16px}h1{font-size:24px;margin:0 0 4px}h2{margin:36px 0 8px;font-size:19px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;margin:12px 0;overflow-x:auto}
table{border-collapse:collapse;width:100%}td{border-top:1px solid var(--line);padding:6px 8px;vertical-align:top}
td:first-child{color:var(--muted);width:30%}code{font-size:12.5px}.bad{color:var(--bad);font-weight:600}
img{max-width:100%;border:1px solid var(--line);border-radius:6px}.muted{color:var(--muted)}
"""


def build(manifest_path: Path) -> Path:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    stages = manifest["expanded_job"]["stages"]
    records = {s["id"]: s for s in manifest["stages"]}
    failed_stage = None
    if manifest.get("status") == "failed":
        m = re.match(r"stage '([^']+)'", manifest.get("error", ""))
        failed_stage = m.group(1) if m else None
    status = {}
    for s in stages:
        rec = records.get(s["id"])
        status[s["id"]] = ("reused" if rec and rec.get("reused") else "ran") if rec else (
            "failed" if s["id"] == failed_stage else "not run")
    links = edges(stages)
    out_dir = ROOT / "viewer_outputs" / "pipeline-reports" / manifest_path.stem.replace(".pipeline", "")
    out_dir.mkdir(parents=True, exist_ok=True)
    title = manifest.get("label", manifest_path.stem)
    body = [f"<h1>{html.escape(title)}</h1><p class='muted'>Pipeline manifest "
            f"<code>{html.escape(manifest_path.name)}</code> · status <b>{manifest.get('status', 'complete')}</b>"
            + (f" · extends {', '.join(html.escape(e) for e in manifest.get('extends_chain', []))}"
               if manifest.get("extends_chain") else "") + "</p>"]
    if manifest.get("error"):
        body.append(f"<p class='bad'>{html.escape(manifest['error'])}</p>")
    if manifest.get("expanded_job", {}).get("command_inputs"):
        body.append(f"<p class='muted'>Command inputs: <code>{html.escape(json.dumps(manifest['expanded_job']['command_inputs']))}</code></p>")
    body.append("<div class='card'>" + graph_svg(stages, links, None, status) + "</div>")
    for s in stages:
        sid, job = s["id"], (records.get(s["id"]) or {}).get("job", s["job"])
        rname, layer, text = role(job)
        archives = (records.get(sid) or {}).get("archives", [])
        body.append(f"<h2>({html.escape(layer)}) {html.escape(sid)}: {html.escape(text)}</h2>")
        body.append("<div class='card'>" + graph_svg(stages, links, sid, status) + "</div>")
        if status[sid] == "failed":
            body.append(f"<p class='bad'>Failed: {html.escape(manifest.get('error', ''))}</p>")
        elif status[sid] == "not run":
            body.append("<p class='muted'>Not run (an earlier stage failed or --until stopped before it).</p>")
        body.append("<div class='card'><b>Method</b><table>" + "".join(
            f"<tr><td>{html.escape(k)}</td><td><code>{html.escape(v)}</code></td></tr>" for k, v in method_rows(job))
            + "</table></div>")
        if archives:
            rows = results(rname, job, archives)
            body.append("<div class='card'><b>Results</b>" + (" <span class='muted'>(reused from an identical earlier run)</span>"
                        if status[sid] == "reused" else "") + "<table>" + "".join(
                f"<tr><td>{html.escape(k)}</td><td>{v}</td></tr>" for k, v in rows)
                + f"<tr><td>archive</td><td><code>{html.escape(Path(archives[0]).name)}</code></td></tr></table></div>")
            image, note = viewer_png(archives[0], rname, out_dir / f"{sid}.png")
            body.append("<div class='card'><b>Viewer</b> <span class='muted'>(" + html.escape(note) + ")</span><br>"
                        + (f"<img alt='Viewer plot for {html.escape(sid)}' src='data:image/png;base64,{image}'>"
                           if image else "") + "</div>")
    page = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
            f"content='width=device-width,initial-scale=1'><title>{html.escape(title)}</title><style>{CSS}</style>"
            f"</head><body><main>{''.join(body)}</main></body></html>")
    path = out_dir / "report.html"
    path.write_text(page, encoding="utf-8")
    return path


def run(job: str, extra: list[str]) -> Path:
    started = time.time()
    cli = str(BIN / "dynsyslab-cli.exe")
    proc = subprocess.run([cli, "run", job, "--archive-root", "archives", *extra], cwd=ROOT,
                          capture_output=True, text=True, stdin=subprocess.DEVNULL)
    sys.stderr.write(proc.stderr[-2000:])
    manifests = [p for p in (ROOT / "archives" / "campaigns" / "pipelines").glob("*.pipeline.json")
                 if p.stat().st_mtime >= started - 1]
    if not manifests:
        raise SystemExit(f"pipeline produced no manifest (exit {proc.returncode})")
    return max(manifests, key=lambda p: p.stat().st_mtime)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run"); r.add_argument("job")
    p = sub.add_parser("report"); p.add_argument("manifests", nargs="+")
    args, extra = parser.parse_known_args()
    manifests = [run(args.job, extra)] if args.command == "run" else [Path(m) for m in args.manifests]
    for manifest in manifests:
        print(build(manifest.resolve()))


if __name__ == "__main__":
    main()
