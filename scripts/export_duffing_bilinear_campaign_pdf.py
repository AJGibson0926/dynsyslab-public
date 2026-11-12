"""Export the bounded Duffing bilinear-KROM campaign as a formal static PDF."""
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

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

from rom import load_control_approximation_error_fields


DEFAULT_REPORT = ROOT / "archives" / "experiments" / "duffing_bilinear_krom_search_20260829" / "campaign-report.json"
DEFAULT_OUTPUT = ROOT / "output" / "pdf" / "duffing-bilinear-krom-control-search.pdf"
DEFAULT_TEMP = ROOT / "tmp" / "pdfs" / "duffing-bilinear-krom-control-search"

BLUE = "#2468a2"
RED = "#c9362b"
GOLD = "#d69428"
GREEN = "#27805d"
PURPLE = "#7357a5"


def _load_solution(path: str | Path) -> dict[str, np.ndarray]:
    with np.load(Path(path), allow_pickle=False) as data:
        time = np.asarray(data["time"], dtype=float)
        values = np.asarray(data["values"], dtype=float)
        control = (
            np.asarray(data["control_values"], dtype=float)
            if "control_values" in data.files else np.empty((time.size, 0))
        )
        target = (
            np.asarray(data["control_target_values"], dtype=float)
            if "control_target_values" in data.files else np.tile((-0.25, 0.0), (time.size, 1))
        )
    return {"time": time, "state": values, "control": control, "target": target}


def _display_indices(arrays: list[np.ndarray], limit: int = 6000) -> np.ndarray:
    count = int(np.asarray(arrays[0]).shape[0])
    if count <= limit:
        return np.arange(count)
    selected = set(np.linspace(0, count - 1, limit, dtype=int).tolist())
    for values in arrays:
        flat = np.asarray(values, dtype=float).reshape(count, -1)
        for column in range(flat.shape[1]):
            finite = np.flatnonzero(np.isfinite(flat[:, column]))
            if finite.size:
                selected.add(int(finite[np.argmin(flat[finite, column])]))
                selected.add(int(finite[np.argmax(flat[finite, column])]))
    return np.asarray(sorted(selected), dtype=int)


def _solution_figure(path: Path, attempt: dict[str, Any]) -> None:
    exact = attempt["cases"]["exact"]
    data = _load_solution(exact["archive_path"])
    indices = _display_indices([data["time"], data["state"], data["control"]])
    time = data["time"][indices]
    state = data["state"][indices]
    target = data["target"][indices]
    control = data["control"][indices, 0]
    release = int(np.argmin(np.abs(data["time"] - 10.0)))
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 7.55), constrained_layout=True)
    fig.suptitle(
        "Final candidate exact B_psi reference - retained indices (8, 9)\n"
        f"Release error {exact['metrics']['distance_to_target_at_release']:.6g}; "
        f"final left-equilibrium error {exact['metrics']['final_distance_to_left_equilibrium']:.6g}",
        fontsize=15, fontweight="bold",
    )
    phase, position, velocity, driver = axes.ravel()
    phase.plot(state[:, 0], state[:, 1], color=BLUE, lw=1.3, label="solution")
    phase.plot(target[:, 0], target[:, 1], color=RED, ls="--", lw=1.7, label="target")
    phase.scatter(*data["state"][0], color=GOLD, edgecolor="black", s=42, label="initial")
    phase.scatter(*data["state"][release], color=GREEN, edgecolor="black", s=42, label="release")
    phase.set(xlabel="x", ylabel="x_dot", title="Phase portrait")
    phase.set_box_aspect(1)
    position.plot(time, state[:, 0], color=BLUE, lw=1.2, label="solution")
    position.plot(time, target[:, 0], color=RED, ls="--", lw=1.6, label="target")
    position.set(xlabel="dimensionless time", ylabel="x", title="Position")
    velocity.plot(time, state[:, 1], color=PURPLE, lw=1.2, label="solution")
    velocity.plot(time, target[:, 1], color=RED, ls="--", lw=1.6, label="target")
    velocity.set(xlabel="dimensionless time", ylabel="x_dot", title="Velocity")
    driver.plot(time, control, color=GREEN, lw=1.0)
    driver.set(xlabel="dimensionless time", ylabel="dimensionless u", title="Control history")
    active = (data["time"] >= 1.0) & (data["time"] <= 10.0)
    values = data["control"][active, 0]
    low, high = np.percentile(values, (1.0, 99.0))
    if np.ptp(data["control"][:, 0]) / max(high - low, np.finfo(float).eps) > 12.0:
        inset = driver.inset_axes([0.53, 0.51, 0.44, 0.42])
        inset.plot(time, control, color=GREEN, lw=0.8)
        inset.set(xlim=(1.0, 10.0), ylim=(low - 0.08 * (high - low), high + 0.08 * (high - low)))
        inset.set_title("active detail (1-99%)", fontsize=7)
        inset.tick_params(labelsize=6)
        inset.grid(alpha=0.2)
    for axis in axes.ravel():
        axis.grid(alpha=0.22)
        if axis is not driver:
            axis.legend(fontsize=8, loc="best")
    for axis in (position, velocity, driver):
        axis.axvspan(1.0, 10.0, color="#dce6ef", alpha=0.35, zorder=-5)
        axis.axvline(10.0, color="#555555", lw=0.8, ls=":")
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def _campaign_figure(path: Path, report: dict[str, Any]) -> None:
    attempts = report["attempts"]
    numbers = np.asarray([item["attempt_number"] for item in attempts])
    exact_error = np.asarray([
        (item.get("cases", {}).get("exact", {}).get("metrics") or {}).get("distance_to_target_at_release", np.nan)
        for item in attempts
    ], dtype=float)
    map_error = np.asarray([
        (item.get("diagnostics", {}).get("error_fields", {}).get("bilinear-map", {}).get("relative") or {}).get("median", np.nan)
        for item in attempts
    ], dtype=float)
    statuses = {}
    for item in attempts:
        for name, case in item.get("cases", {}).items():
            if name.startswith("bilinear"):
                statuses[case.get("status", "unknown")] = statuses.get(case.get("status", "unknown"), 0) + 1
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 7.55), constrained_layout=True)
    fig.suptitle("Bounded Duffing bilinear-KROM search - campaign diagnostics", fontsize=16, fontweight="bold")
    axes[0, 0].semilogy(numbers, exact_error, "o-", color=BLUE, ms=4)
    axes[0, 0].axhline(0.02, color=RED, ls="--", lw=1.2, label="success tolerance")
    axes[0, 0].set(xlabel="attempt", ylabel="release target error", title="Exact B_psi controllability gate")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(numbers, map_error, "o-", color=PURPLE, ms=4)
    axes[0, 1].axhline(1.0, color=RED, ls="--", lw=1.1, label="order-one error")
    axes[0, 1].set(xlabel="attempt", ylabel="median relative Frobenius error", title="Global direct bilinear-map error")
    axes[0, 1].legend(fontsize=8)
    labels = list(statuses)
    axes[1, 0].bar(labels, [statuses[key] for key in labels], color=[GREEN, GOLD, RED, PURPLE][:len(labels)])
    axes[1, 0].set(ylabel="case count", title="Bilinear closed-loop outcomes")
    axes[1, 0].tick_params(axis="x", rotation=20)
    final = attempts[-1]["diagnostics"]["timings"]
    names = ["constant", "bilinear dynamics", "direct map", "path tube", "error fields"]
    values = [
        final.get("constant_regression_seconds", np.nan),
        final.get("bilinear_dynamics_regression_seconds", np.nan),
        final.get("bilinear_map_regression_seconds", np.nan),
        final.get("bilinear_path_tube_map_regression_seconds", np.nan),
        final.get("error_field_seconds", np.nan),
    ]
    axes[1, 1].barh(names, values, color=BLUE)
    axes[1, 1].set(xlabel="seconds", title="Final-attempt identification and map timing")
    for axis in axes.ravel():
        axis.grid(alpha=0.2)
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def _heatmap_pages(temp: Path, error_archive: str | Path) -> list[Path]:
    fields, _ = load_control_approximation_error_fields(error_archive)
    items = list(fields.items())
    outputs = []
    for page_index, start in enumerate(range(0, len(items), 4), start=1):
        group = items[start:start + 4]
        fig, axes = plt.subplots(2, 2, figsize=(12.2, 7.55), constrained_layout=True)
        fig.suptitle(
            f"Final-candidate control-map errors - group {page_index}",
            fontsize=16, fontweight="bold",
        )
        for axis in axes.ravel():
            axis.set_visible(False)
        for axis, (name, field) in zip(axes.ravel(), group):
            axis.set_visible(True)
            relative = np.asarray(field.relative_error, dtype=float)
            positive = relative[relative > 0.0]
            floor = max(float(np.min(positive)) if positive.size else 1e-16, 1e-16)
            image = axis.imshow(
                np.log10(np.maximum(relative, floor)), origin="lower", aspect="equal",
                cmap="viridis", extent=(field.x_values[0], field.x_values[-1], field.y_values[0], field.y_values[-1]),
            )
            axis.set(xlabel="x", ylabel="x_dot", title=name.replace("-", " "))
            colorbar = fig.colorbar(image, ax=axis, shrink=0.82)
            colorbar.set_label("log10 relative Frobenius error", fontsize=8)
            axis.text(
                0.02, 0.02,
                f"median={np.median(relative):.3g}\np95={np.percentile(relative, 95):.3g}",
                transform=axis.transAxes, color="white", fontsize=8, va="bottom",
                bbox={"facecolor": "black", "alpha": 0.45, "pad": 3, "edgecolor": "none"},
            )
        output = temp / f"final-errors-{page_index}.png"
        fig.savefig(output, dpi=180, facecolor="white")
        plt.close(fig)
        outputs.append(output)
    return outputs


def _status(case: dict[str, Any] | None) -> str:
    if not case:
        return "missing"
    if case.get("status") != "complete":
        return str(case.get("status"))
    error = (case.get("metrics") or {}).get("distance_to_target_at_release")
    return f"complete ({float(error):.3g})" if error is not None else "complete"


def export(report_path: Path, output_path: Path, temp_path: Path) -> Path:
    """Render one retained Duffing campaign result as the standard PDF."""
    if plt is None:
        raise ModuleNotFoundError(
            "PDF reporting dependencies are unavailable; install dynsyslab[reporting]."
        )
    try:
        from postprocessing.formal_report import FormalLandscapePDF
    except ModuleNotFoundError as error:
        raise ModuleNotFoundError(
            "PDF reporting dependencies are unavailable; install dynsyslab[reporting]."
        ) from error
    report = json.loads(report_path.read_text(encoding="utf-8"))
    attempts = report["attempts"]
    temp_path.mkdir(parents=True, exist_ok=True)
    campaign_png = temp_path / "campaign-diagnostics.png"
    solution_png = temp_path / "final-exact-solution.png"
    _campaign_figure(campaign_png, report)
    _solution_figure(solution_png, attempts[-1])
    heatmaps = _heatmap_pages(temp_path, attempts[-1]["artifacts"]["control_error_archive"])

    document = FormalLandscapePDF(
        output_path,
        title="Duffing bilinear KROM control search",
        footer_label="DynSysLab - bounded Duffing bilinear KROM control search",
    )
    pdf = document.canvas
    pdf.setFillColor("#17212b")
    pdf.rect(0, 0, document.width, document.height, fill=1, stroke=0)
    pdf.setFillColor("#f4f7f9")
    pdf.setFont("Helvetica-Bold", 25)
    pdf.drawString(42, document.height - 58, "Duffing bilinear KROM control search")
    pdf.setFont("Helvetica", 11)
    pdf.setFillColor("#c9d5df")
    pdf.drawString(44, document.height - 81, "Formal retained report - 24 durable attempts, hard-capped search")
    pdf.setFont("Helvetica-Bold", 17)
    pdf.setFillColor("#e7edf2")
    pdf.drawString(44, document.height - 132, "Outcome: no successful bilinear control approximation")
    document.paragraph(
        "Exact analytic <b>B_psi</b> repeatedly reaches the transfer target and releases into the left well. "
        "Every fitted constant or bilinear approximation fails the unchanged closed-loop gate. "
        "The result isolates a control-map closure and robustness problem, not a regression-speed problem or a failure of KLQR with the exact map.",
        46, document.height - 225, document.width - 92, 86, light=True, font_size=11, leading=16,
    )
    cover_rows = [
        ["Contract", "Value"],
        ["State / target", "(0.5, 0) -> (-0.25, 0); release to (-0.5, 0)"],
        ["Control", "u=-K(Psi-Psi_target); Q=12500 I; R=0.01; no cap"],
        ["Timing", "active 1 to 10; final 20; dimensionless"],
        ["Gate", "release <=0.02; final left-equilibrium <=0.10"],
        ["Attempt accounting", "24 durable + 1 interrupted pilot = conservative ceiling 25"],
    ]
    document.table(cover_rows, x=45, y=92, column_widths=(150, 510), row_heights=28, font_size=9)
    document.finish_page()

    document.image_page(campaign_png, "Campaign-level evidence")
    document.image_page(solution_png, "Exact-control reference for the final candidate")

    for start in (0, 12):
        document.title(f"Attempt ledger {start + 1}-{min(start + 12, len(attempts))}")
        rows = [["#", "indices", "exact", "constant", "bilinear dyn.", "direct map", "bilinear success"]]
        for attempt in attempts[start:start + 12]:
            cases = attempt["cases"]
            rows.append([
                attempt["attempt_number"],
                ",".join(map(str, attempt["candidate"]["retained_indices"])),
                _status(cases.get("exact")),
                _status(cases.get("constant")),
                _status(cases.get("bilinear-dynamics")),
                _status(cases.get("bilinear-map")),
                "yes" if attempt.get("bilinear_success") else "no",
            ])
        document.table(rows, x=36, y=56, column_widths=(30, 55, 104, 104, 112, 112, 82), row_heights=35, font_size=7.3)
        document.finish_page()

    for image in heatmaps:
        document.image_page(image, "Final-candidate retained error maps")

    text_pages = [
        (
            "Diagnosis: what the campaign establishes",
            "<b>1. Exact control and approximation quality are separate gates.</b><br/>"
            "Modes 6 and 9 and selected two-mode pairs pass exact control. Other subsets fail before regression quality matters.<br/><br/>"
            "<b>2. Global bilinear closure is absent.</b><br/>"
            "Median relative B_psi errors stay near one over the training rectangle. The retained small eigenfunction subsets do not span their input-direction derivatives.<br/><br/>"
            "<b>3. Conjugate completion is necessary in principle but insufficient here.</b><br/>"
            "M Psi + N conjugate(Psi) + b represents a general real-linear split-coordinate map, yet its fitted errors remain too large.<br/><br/>"
            "<b>4. Closed-loop sensitivity dominates.</b><br/>"
            "Arc-length fits reduce error along the exact path, but a small early deviation changes B_psi, the CARE gain, and the state path together. The approximate trajectory then exits the accurately fitted region and the feedback solve becomes stiff or unstable.",
        ),
        (
            "Selection criteria learned",
            "Future searches should gate candidates in this order:<br/><br/>"
            "1. Exact B_psi control succeeds for the frozen target and tuning.<br/>"
            "2. Local controllability rank is complete and its smallest singular value is usable.<br/>"
            "3. Map error is small at engagement and at the target, not merely in median.<br/>"
            "4. Error remains small in a tube around the expected closed-loop path.<br/>"
            "5. CARE gains remain finite and smooth under bounded B_psi perturbations.<br/>"
            "6. Only then run the full nonlinear controlled solve.<br/><br/>"
            "Frequency spread, polynomial closure score, and global median error were weak predictors when used alone. Complementary target geometry matters: mode 8 fails alone but succeeds exactly when paired with mode 9.",
        ),
        (
            "Later options considered but not performed",
            "<b>Controlled data.</b> Improve phase-space coverage and persistent excitation, especially around engagement and the transfer tube, with a held-out validation set.<br/><br/>"
            "<b>Parent KROM.</b> Consider a larger or control-aware eigenfunction span whose derivatives under the input vector field are closer to closure. Do not assume more modes automatically help; larger exact KLQR systems were often less well conditioned.<br/><br/>"
            "<b>Eigenvalues and frequencies.</b> Explore additional spectral content only after exact-control and conditioning gates. Frequency diversity alone did not predict success.<br/><br/>"
            "<b>Regression objectives.</b> Consider orthonormalized lift coordinates, regularized local bases, gain-sensitivity-aware loss, and explicit endpoint/path-tube validation. Preserve exact B_psi as the objectively correct default; approximations remain optional and metadata-tagged.",
        ),
        (
            "Reproduction and retained evidence",
            "The lightweight workflow is available through:<br/><br/>"
            "<font name='Courier'>dynsyslab-duffing-bilinear-search --dry-run --max-attempts 24</font><br/>"
            "<font name='Courier'>dynsyslab-duffing-bilinear-search --resume --max-attempts 24 --case-timeout-seconds 60</font><br/>"
            "<font name='Courier'>dynsyslab-duffing-bilinear-report</font><br/>"
            "<font name='Courier'>dynsyslab-duffing-bilinear-pdf</font><br/><br/>"
            "Every attempt retains an independent extracted KROM, fitted control representations, matrices, full phase-space error fields, case records, and completed solution archives where available. The atomic campaign ledger supports resume without recounting completed work. All 24 per-attempt error maps are also rendered through Viewer as native square heat maps.",
        ),
    ]
    for title, body in text_pages:
        document.title(title)
        document.paragraph_top(
            body, 48, document.height - 86, document.width - 96,
            document.height - 145, font_size=11, leading=17,
        )
        document.finish_page()
    return document.save()


def main() -> int:
    """Export a bounded bilinear-search result to the repeatable PDF format."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--temp", type=Path, default=DEFAULT_TEMP)
    args = parser.parse_args()
    output = export(
        args.report.expanduser().resolve(),
        args.output.expanduser().resolve(),
        args.temp.expanduser().resolve(),
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
