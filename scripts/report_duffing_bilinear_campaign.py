"""Build retained Markdown and native Viewer reports for a Duffing KROM campaign."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from viewer import ArchiveSource, ControlErrorSource, PlotRequest, TilePlan, Viewer


DEFAULT_CAMPAIGN = ROOT / "archives" / "experiments" / "duffing_bilinear_krom_search_20260829"
DEFAULT_MARKDOWN = ROOT / "docs" / "case-studies" / "DUFFING_BILINEAR_KROM_SEARCH_2026-08-29.md"
DEFAULT_VIEWER_ROOT = ROOT / "viewer_outputs" / "duffing-bilinear-krom-search"
DEFAULT_SUMMARY_VIEWER = ROOT / "viewer_outputs" / "duffing-bilinear-krom-search-final.html"


def _fmt(value: Any, default: str = "-") -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return f"{number:.6g}" if np.isfinite(number) else default


def _case_cell(case: dict[str, Any] | None) -> str:
    if not case:
        return "missing"
    status = str(case.get("status", "unknown"))
    metrics = case.get("metrics") or {}
    if status == "complete":
        return f"complete; release={_fmt(metrics.get('distance_to_target_at_release'))}"
    return status


def _best_completed_bilinear(attempts: list[dict[str, Any]]) -> tuple[int, str, float] | None:
    rows = []
    for attempt in attempts:
        for name, case in attempt.get("cases", {}).items():
            if not name.startswith("bilinear") or case.get("status") != "complete":
                continue
            error = float((case.get("metrics") or {}).get("distance_to_target_at_release", np.inf))
            if np.isfinite(error):
                rows.append((int(attempt["attempt_number"]), name, error))
    return min(rows, key=lambda item: item[2]) if rows else None


def _write_attempt_viewers(attempts: list[dict[str, Any]], viewer_root: Path) -> list[Path]:
    viewer_root.mkdir(parents=True, exist_ok=True)
    outputs = []
    for attempt in attempts:
        number = int(attempt["attempt_number"])
        archive = attempt.get("artifacts", {}).get("control_error_archive")
        if not archive:
            continue
        viewer = Viewer(
            tile_plan=TilePlan(f"duffing_bilinear_search_attempt_{number:02d}", columns=3),
            renderer="plotly",
        )
        viewer.add(PlotRequest(
            "control_error",
            ControlErrorSource(archive),
            control_error_measure="relative",
            control_error_log_scale=True,
        ))
        destination = viewer_root / f"attempt-{number:02d}-control-map-errors.html"
        outputs.append(Path(viewer.write_html(viewer.render(), destination)).resolve())
    return outputs


def _write_summary_viewer(attempt: dict[str, Any], destination: Path) -> Path:
    viewer = Viewer(
        tile_plan=TilePlan("duffing_bilinear_search_final", columns=3),
        renderer="plotly",
    )
    viewer.add(PlotRequest(
        "control_error",
        ControlErrorSource(attempt["artifacts"]["control_error_archive"]),
        control_error_measure="relative",
        control_error_log_scale=True,
    ))
    exact = attempt.get("cases", {}).get("exact", {})
    if exact.get("status") == "complete" and exact.get("archive_path"):
        source = ArchiveSource(exact["archive_path"], label="Exact B_psi reference")
        viewer.add(PlotRequest(
            "phase", source, title="Final-candidate exact reference: phase portrait",
            max_display_samples=5000,
        ))
        viewer.add(PlotRequest(
            "time", source, title="Final-candidate exact reference",
            time_series_layout="separate", max_display_samples=5000,
        ))
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(viewer.write_html(viewer.render(), destination)).resolve()


def _markdown(report: dict[str, Any], viewer_outputs: list[Path], summary_viewer: Path) -> str:
    attempts = list(report["attempts"])
    exact_success = [
        int(item["attempt_number"])
        for item in attempts
        if "exact" in item.get("successful_cases", [])
    ]
    best = _best_completed_bilinear(attempts)
    lines = [
        "# Duffing bilinear KROM control search",
        "",
        "## Bottom line",
        "",
        f"The bounded campaign completed **{len(attempts)} durable attempts** and found **no successful bilinear control approximation**. "
        "One earlier pilot was interrupted before it entered the durable ledger; counting it conservatively gives 25 physical extraction-regression-control sequences, exactly the user-specified ceiling.",
        "",
        "This is a negative result about the retained approximation families, not about KLQR or the selected KROMs in general. "
        f"Exact analytic `B_psi` passed the unchanged transfer/release gate in attempts {', '.join(map(str, exact_success))}. "
        "The same KROMs failed when their input maps were replaced by a constant matrix, a controlled-dynamics bilinear fit, a global pointwise bilinear map, a corridor-local map, a conjugate-complete real-linear map, or the path-conditioned candidates.",
        "",
        f"The final native Viewer dashboard is [{summary_viewer.name}]({summary_viewer.as_posix()}). "
        f"All {len(viewer_outputs)} per-attempt error-map dashboards were generated under `{viewer_outputs[0].parent if viewer_outputs else 'viewer_outputs'}`.",
        "",
        "## Frozen experiment contract",
        "",
        "- Existing Duffing degree-4 SINDy KROM and existing 100-member random-forced controlled dataset; neither was regenerated or expanded.",
        "- Initial state `(0.5, 0)`, constant target `(-0.25, 0)`, control active on `[1, 10]`, release at 10, final time 20.",
        "- Dimensionless Duffing dynamics `x_dot=v`, `v_dot=-0.5 v + x - 4 x^3 + u`.",
        "- KLQR law `u=-K(Psi(x)-Psi(target))`, `Q=12500 I`, `R=0.01`, no control cap, continuous RHS control evaluation.",
        "- SciPy BDF, `rtol=1e-8`, `atol=1e-10`, `max_step=0.01`; screening output interval `0.001`.",
        "- Success gate: release target error at most 0.02 and final left-equilibrium error at most 0.10.",
        "",
        "## Tooling and efficiency",
        "",
        "The campaign is available through `dynsyslab-duffing-bilinear-search` and the public `rom.control_campaign` API. "
        "It enforces a hard maximum of 25 attempts, writes an atomic resume ledger after every attempt, isolates every controller solve in a timed process, archives each extracted independent KROM and fitted control representation, and saves every full phase-space error field.",
        "",
        "The Korda-Mezic constant-matrix multistep design was algebraically reduced from a nested convolution to an O(horizon) recurrence. "
        "Representative fresh constant fits now take roughly 0.4 seconds; pointwise bilinear fits take about 1-4 seconds. "
        "The long wall time comes from stiff or unstable closed-loop integration, not matrix regression.",
        "",
        "## Attempt ledger",
        "",
        "| # | retained KROM indices | exact | constant | bilinear dynamics | direct bilinear map | successful cases |",
        "|---:|---|---|---|---|---|---|",
    ]
    for attempt in attempts:
        cases = attempt.get("cases", {})
        indices = ",".join(map(str, attempt["candidate"]["retained_indices"]))
        successes = ", ".join(attempt.get("successful_cases", [])) or "none"
        lines.append(
            f"| {attempt['attempt_number']} | `{indices}` | {_case_cell(cases.get('exact'))} | "
            f"{_case_cell(cases.get('constant'))} | {_case_cell(cases.get('bilinear-dynamics'))} | "
            f"{_case_cell(cases.get('bilinear-map'))} | {successes} |"
        )
    best_text = "No bilinear case completed with finite metrics."
    if best is not None:
        best_text = (
            f"The smallest release error among completed bilinear solves was `{best[2]:.6g}` "
            f"for attempt {best[0]} (`{best[1]}`), still far outside the 0.02 success gate."
        )
    lines.extend([
        "",
        "## Evidence and diagnosis",
        "",
        best_text,
        "",
        "1. **Exact controllability is a necessary independent gate.** Several subsets fail even with exact `B_psi`; those KROMs cannot be rescued by a better regression. Modes 6 and 9 and the pairs `(6,8)`, `(6,9)`, and `(8,9)` are exact-controllable for this target.",
        "2. **Global affine closure is absent.** Across exact-controllable KROMs, global bilinear-map median relative errors remain approximately order one. This means `B_psi(x)` is not in the affine span of the retained eigenfunctions over the training rectangle.",
        "3. **Adding a conjugate term is insufficient.** The candidate `M Psi + N conjugate(Psi) + b` correctly spans a general real-linear map on split coordinates, but its global median errors also remain near one and its feedback solves fail.",
        "4. **Naive locality can make the global map worse.** The rectangular target-corridor dataset is nonuniform and highly correlated; local fits often have worse conditioning and larger out-of-corridor errors.",
        "5. **Time-weighted path fitting hid the engagement defect.** For mode 6 the reference-path median relative error fell to about 0.059, but the error at engagement was about 2.69 because most samples represented the long target dwell. Both the path and path-tube controllers timed out.",
        "6. **Arc-length collocation improved pathwise error but not closed-loop robustness.** For `(8,9)`, seven- and nine-node conjugate-complete arc fits reduced reference-path maximum relative error to roughly 0.18, yet all closed-loop solves timed out. Once the approximate controller leaves the exact reference path, its errors and CARE conditioning amplify each other.",
        "7. **Regression error alone is not a sufficient selection score.** Useful pre-screens must combine exact-control success, target displacement in lift space, controllability/CARE conditioning, engagement-point error, path-tube error, and sensitivity of the gain to perturbations in `B_psi`.",
        "",
        "## Selection criteria learned",
        "",
        "A future candidate should be admitted only if it passes these gates in order: (a) exact `B_psi` control succeeds; (b) local controllability rank is complete with a usable smallest singular value; (c) the fitted map is accurate at the initial state and target, not merely in aggregate; (d) path-tube error stays small; (e) CARE solutions and gains vary smoothly under bounded map perturbations; and only then (f) full nonlinear integration is run.",
        "",
        "The polynomial closure score, eigenvalue frequency spread, and global median map error were individually weak predictors. Mode 8 demonstrates why: it fails alone under exact control but succeeds when paired with mode 9, so complementary target geometry matters. Conversely, adding modes can reduce controllability conditioning enough to make exact control impractical.",
        "",
        "## Considerations for later work (not performed)",
        "",
        "- Expand the controlled dataset with better state-space coverage and persistent excitation, especially around the transfer tube and engagement state, while preserving a held-out validation set.",
        "- Build a larger or control-aware parent KROM whose eigenfunction span is closer to derivative closure under the Duffing input vector field. The derivative of a degree-4 polynomial eigenfunction is degree 3, but arbitrary small subsets of learned eigenfunctions need not span those derivatives.",
        "- Seek eigenfunctions with different eigenvalues/frequencies only after exact-control and conditioning gates; frequency diversity by itself did not predict success here.",
        "- Consider regularized local bases, orthonormalized lift coordinates, and gain-sensitivity-aware regression objectives. Do not infer that simply adding modes will help: larger exact KLQR systems became much more expensive and sometimes ill-conditioned.",
        "- Preserve exact `B_psi` as the objectively correct default. Constant and bilinear maps remain optional approximations with explicit validation metadata.",
        "",
        "## Reproduction",
        "",
        "```text",
        "dynsyslab-duffing-bilinear-search --dry-run --max-attempts 24",
        "dynsyslab-duffing-bilinear-search --resume --max-attempts 24 --case-timeout-seconds 60",
        "dynsyslab-duffing-bilinear-report",
        "```",
        "",
        f"Durable campaign record: `{Path(DEFAULT_CAMPAIGN / 'campaign-report.json').resolve()}`",
        "",
        "Every attempt directory contains its independent KROM archive references, fitted matrix arrays, controller records, and `control-map-errors.npz`. Failed solves are retained as explicit statuses; they are never converted into missing or successful data.",
    ])
    return "\n".join(lines)


def build(
    campaign_root: Path,
    markdown_path: Path,
    viewer_root: Path,
    summary_viewer_path: Path,
) -> dict[str, Any]:
    """Build Viewer artifacts and Markdown from one retained search campaign."""
    report_path = campaign_root / "campaign-report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    attempts = list(report["attempts"])
    if not attempts:
        raise ValueError("Campaign report contains no attempts.")
    viewers = _write_attempt_viewers(attempts, viewer_root)
    summary = _write_summary_viewer(attempts[-1], summary_viewer_path)
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.write_text(_markdown(report, viewers, summary) + "\n", encoding="utf-8")
    return {
        "campaign_report": str(report_path.resolve()),
        "attempt_count": len(attempts),
        "markdown_report": str(markdown_path.resolve()),
        "summary_viewer": str(summary),
        "attempt_viewer_count": len(viewers),
        "attempt_viewer_root": str(viewer_root.resolve()),
    }


def main() -> int:
    """Create the standard report bundle for a retained Duffing search."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, default=DEFAULT_CAMPAIGN)
    parser.add_argument("--markdown", type=Path, default=DEFAULT_MARKDOWN)
    parser.add_argument("--viewer-root", type=Path, default=DEFAULT_VIEWER_ROOT)
    parser.add_argument("--summary-viewer", type=Path, default=DEFAULT_SUMMARY_VIEWER)
    args = parser.parse_args()
    result = build(
        args.campaign.expanduser().resolve(),
        args.markdown.expanduser().resolve(),
        args.viewer_root.expanduser().resolve(),
        args.summary_viewer.expanduser().resolve(),
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
