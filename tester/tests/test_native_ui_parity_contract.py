from __future__ import annotations

import json
import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs" / "architecture" / "NATIVE_UI_PARITY_CONTRACT.json"


def _check_every_retained_capability_has_a_native_ui_surface() -> None:
    ledger = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert ledger["schema"] == "dynsyslab.native-ui-parity/v1"
    assert ledger["release_gate"] is True
    assert ledger["runtime_python_allowed"] is False

    capabilities = ledger["capabilities"]
    identifiers = [entry["id"] for entry in capabilities]
    assert len(identifiers) == len(set(identifiers))
    assert {entry["app"] for entry in capabilities} == {
        "DynSysLab.exe",
        "DynSysLabViewer.exe",
    }

    for entry in capabilities:
        assert entry["status"] == "implemented", entry["id"]
        service = ROOT / entry["native_service"]
        assert service.is_file(), (entry["id"], service)
        assert service.suffix in {".cpp", ".hpp"}, entry["id"]
        assert "python" not in service.as_posix().lower(), entry["id"]

        evidence = entry["ui_evidence"]
        assert evidence, entry["id"]
        for item in evidence:
            source = ROOT / item["path"]
            assert source.is_file(), (entry["id"], source)
            assert source.suffix in {".cpp", ".hpp"}, entry["id"]
            text = source.read_text(encoding="utf-8")
            for token in item["tokens"]:
                assert token in text, (entry["id"], source, token)

        tests = entry["test_evidence"]
        assert tests, entry["id"]
        for test in tests:
            assert (ROOT / test).is_file(), (entry["id"], test)


def _check_core_human_workflow_groups_are_covered() -> None:
    ledger = json.loads(CONTRACT.read_text(encoding="utf-8"))
    identifiers = {entry["id"] for entry in ledger["capabilities"]}
    for prefix in ("solver.", "operations.", "postprocess.", "viewer."):
        assert any(identifier.startswith(prefix) for identifier in identifiers)

    required = {
        "solver.controllers_targets_and_kroms",
        "solver.control_evaluation_policy",
        "solver.krom_training_domain_expansion",
        "operations.durable_queue",
        "operations.archive_audit_and_deduplication",
        "operations.legacy_tagged_csv_import",
        "operations.krom_control_candidate_search",
        "postprocess.hdmd_scalar_observable",
        "postprocess.sindy_eigenfunctions",
        "postprocess.boundary_eigenfunctions",
        "postprocess.eigenvalue_optimization",
        "postprocess.krom_extraction",
        "postprocess.koopman_validation",
        "postprocess.control_map_regression",
        "postprocess.observability_and_klqe",
        "viewer.solution_and_control_views",
        "viewer.batch_cluster_and_campaign_aggregates",
        "viewer.complete_view_selection",
        "viewer.automatic_detail_insets",
        "viewer.koopman_spectrum_and_eigenfunctions",
        "viewer.koopman_validation",
        "viewer.eigenvalue_optimization",
        "viewer.control_map_error_heatmaps",
        "viewer.krom_control_campaigns",
        "viewer.clipboard_and_disk_export",
    }
    assert required <= identifiers


def _check_compiled_cli_scientific_capabilities_are_not_cli_only() -> None:
    ledger = json.loads(CONTRACT.read_text(encoding="utf-8"))
    source = (ROOT / "compute_workbench" / "native" / "src" / "Cli.cpp").read_text(
        encoding="utf-8"
    )
    capabilities = source.split("void capabilities()", 1)[1].split("void systems()", 1)[0]
    controller_text, postprocess_and_workflows = capabilities.split("postprocessing", 1)
    controller_text = controller_text.split("controllers", 1)[1]
    postprocess_text, workflow_text = postprocess_and_workflows.split("workflows", 1)
    workflow_text = workflow_text.split("desktop_ui", 1)[0]
    quoted = re.compile(r'\\"([a-z][a-z0-9_]*)\\"')
    native_controllers = set(quoted.findall(controller_text))
    native_postprocessing = set(quoted.findall(postprocess_text))
    native_workflows = set(quoted.findall(workflow_text))

    controller_ui = set()
    postprocess_ui = set()
    workflow_ui = set()
    for entry in ledger["capabilities"]:
        covered = set(entry.get("covers", []))
        if entry["id"].startswith("solver."):
            controller_ui |= covered
        if entry["id"].startswith("postprocess."):
            postprocess_ui |= covered
        # Scientific commands such as Fourier have their human surface in Viewer.
        if entry["id"].startswith("viewer."):
            postprocess_ui |= covered & native_postprocessing
        if entry["id"].startswith("operations."):
            workflow_ui |= covered
    assert native_controllers == controller_ui
    assert native_postprocessing == postprocess_ui
    assert native_workflows == workflow_ui


class NativeUiParityContractTests(unittest.TestCase):
    """Release gates for complete, Python-free retained-workflow UI parity."""

    def test_every_retained_capability_has_a_native_ui_surface(self) -> None:
        _check_every_retained_capability_has_a_native_ui_surface()

    def test_core_human_workflow_groups_are_covered(self) -> None:
        _check_core_human_workflow_groups_are_covered()

    def test_compiled_cli_scientific_capabilities_are_not_cli_only(self) -> None:
        _check_compiled_cli_scientific_capabilities_are_not_cli_only()
