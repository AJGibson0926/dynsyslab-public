"""Regression tests for the neutral event channel and audit front door."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from infrastructure import source_events
from quality.architecture import audit_project, compare_baseline, render_markdown


PROJECT_ROOT = Path(__file__).resolve().parents[2]


class ArchitectureQualityTests(unittest.TestCase):
    def test_document_links_distinguish_declared_local_assets_from_broken_paths(self):
        from quality.architecture import _documentation_metrics
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "docs").mkdir()
            (root / "docs" / "README.md").write_text(
                "[Local figure](../viewer_outputs/known.png)\n"
                "[Unknown](../viewer_outputs/missing.png)\n", encoding="utf-8")
            (root / "docs" / "LOCAL_ASSET_POINTERS.json").write_text(json.dumps({
                "schema": "dynsyslab.local-asset-pointers/v1", "assets": [
                    {"path": "viewer_outputs/known.png", "sha256": "a" * 64, "bytes": 20}]}),
                encoding="utf-8")
            result = _documentation_metrics(root)
            self.assertEqual([row["target"] for row in result["local_asset_links"]],
                             ["../viewer_outputs/known.png"])
            self.assertEqual([row["target"] for row in result["broken_links"]],
                             ["../viewer_outputs/missing.png"])

    def test_audit_skips_side_research(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "pyproject.toml").write_text('[project]\nname = "x"\n', encoding="utf-8")
            (root / "docs").mkdir()
            (root / "docs" / "README.md").write_text("# Docs\n", encoding="utf-8")
            (root / "pkg").mkdir()
            (root / "pkg" / "kept.py").write_text('"""Kept."""\n', encoding="utf-8")
            (root / "research" / "study" / "_vendor").mkdir(parents=True)
            (root / "research" / "study" / "_vendor" / "huge.py").write_text(
                '"""Vendored."""\n' + "x = 1\n" * 1200, encoding="utf-8")
            python = audit_project(root)["python"]
            self.assertFalse(python["over_1000_lines"])
        from project_paths import is_side_research, source_files
        self.assertTrue(is_side_research(PROJECT_ROOT / "research" / "anything.py"))
        self.assertFalse(is_side_research(PROJECT_ROOT / "quality" / "architecture.py"))
        self.assertFalse(any(is_side_research(path) for path in source_files()))

    def test_audit_is_versioned_and_covers_python_native_and_docs(self):
        result = audit_project(PROJECT_ROOT)
        self.assertEqual(result["schema"], "dynsyslab.architecture-audit/v1")
        self.assertGreater(result["python"]["files"], 100)
        self.assertGreater(result["native"]["files"], 20)
        self.assertEqual(result["documentation"]["broken_links"], [])
        self.assertEqual(result["python"]["package_cycles"], [])
        self.assertGreaterEqual(result["python"]["public_docstring_coverage"], 0.70)
        self.assertIn("compute_workbench",
                      result["python"]["package_public_docstrings"])
        self.assertIn("# DynSysLab architecture audit", render_markdown(result))

    def test_audit_cli_is_non_executing_json(self):
        completed = subprocess.run(
            (sys.executable, "-m", "quality", "--format", "json"),
            cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=30,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(json.loads(completed.stdout)["schema"],
                         "dynsyslab.architecture-audit/v1")

    def test_checked_baseline_enforces_the_clean_architecture_floor(self):
        result = audit_project(PROJECT_ROOT)
        baseline = json.loads((PROJECT_ROOT / "quality" /
                               "architecture_baseline.json").read_text(encoding="utf-8"))
        comparison = compare_baseline(result, baseline)
        self.assertTrue(comparison["passes"], comparison)
        self.assertTrue(result["passes_strict"])
        self.assertEqual(comparison["regressions"]["hard_size_paths"], [])
        self.assertEqual(comparison["regressions"]["boundary_violations"], [])

    def test_refactored_modules_remain_below_preferred_ceiling(self):
        paths = [PROJECT_ROOT / "viewer" / "core.py",
                 PROJECT_ROOT / "viewer" / "contracts.py",
                 PROJECT_ROOT / "viewer" / "records.py",
                 PROJECT_ROOT / "job_manager" / "manager.py",
                 PROJECT_ROOT / "job_manager" / "_runtime.py",
                 PROJECT_ROOT / "scripts" / "search_duffing_bilinear_krom.py",
                 PROJECT_ROOT / "rom" / "koopman.py",
                 PROJECT_ROOT / "rom" / "control_campaign.py",
                 PROJECT_ROOT / "rom" / "control_candidate_screening.py",
                 PROJECT_ROOT / "rom" / "candidates" /
                 "control_identification.py",
                 PROJECT_ROOT / "rom" / "candidates" /
                 "control_identification_contracts.py",
                 PROJECT_ROOT / "postprocessing" / "koopman" /
                 "koopman_sindy.py",
                 PROJECT_ROOT / "postprocessing" / "koopman" /
                 "sindy_nullspace.py",
                 PROJECT_ROOT / "postprocessing" / "koopman" /
                 "workflows.py",
                 PROJECT_ROOT / "postprocessing" / "koopman" /
                 "duffing_pipeline.py",
                 PROJECT_ROOT / "estimation" / "koopman.py",
                 PROJECT_ROOT / "estimation" / "klqe_storage.py",
                 PROJECT_ROOT / "infrastructure" / "controlled_data.py",
                 PROJECT_ROOT / "infrastructure" / "persistence.py",
                 PROJECT_ROOT / "operations" / "rpe_krom_screening.py"]
        paths.extend((PROJECT_ROOT / "viewer" / "rendering").glob("*.py"))
        paths.extend((PROJECT_ROOT / "job_manager").glob("*_services.py"))
        paths.extend((PROJECT_ROOT / "operations" /
                      "duffing_bilinear").glob("*.py"))
        paths.extend((PROJECT_ROOT / "rom").glob("krom_*.py"))
        oversized = {
            path.relative_to(PROJECT_ROOT).as_posix():
            len(path.read_text(encoding="utf-8").splitlines())
            for path in paths
            if len(path.read_text(encoding="utf-8").splitlines()) > 500
        }
        self.assertEqual(oversized, {})

    def test_quality_workflow_covers_supported_python_wheel_and_native_apps(self):
        workflow = (PROJECT_ROOT / ".github" / "workflows" /
                    "quality.yml").read_text(encoding="utf-8")
        self.assertIn('python-version: ["3.10", "3.12"]', workflow)
        self.assertIn("python -m build", workflow)
        self.assertIn("dynsyslab-wheel-smoke", workflow)
        self.assertIn("cmake --preset windows-fast-release", workflow)
        self.assertIn("windows-fast-release-gates", workflow)
        self.assertIn("dynsyslab_package_smoke", workflow)

    def test_native_release_gate_builds_every_registered_test_executable(self):
        from project_paths import source_files
        cmake_files = source_files("CMakeLists.txt")
        sources = "\n".join(path.read_text(encoding="utf-8") for path in cmake_files)
        registered = set(re.findall(
            r"add_test\s*\(\s*NAME\s+(dynsyslab_[a-z0-9_]+)", sources))
        executables = set(re.findall(
            r"add_executable\s*\(\s*(dynsyslab_[a-z0-9_]+)", sources))
        built_tests = registered & executables
        root_cmake = (PROJECT_ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
        gate = root_cmake.split("if(DYNSYSLAB_BUILD_TESTS)", 1)[1].split("if(WIN32)", 1)[0]
        # Exercise the real recursive collector with the registered target names,
        # without configuring numerical dependencies or invoking a compiler.
        with TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            child = source / "nested"
            child.mkdir(parents=True)
            (child / "CMakeLists.txt").write_text(
                "\n".join(f"add_custom_target({name})" for name in sorted(built_tests)),
                encoding="utf-8")
            (source / "CMakeLists.txt").write_text(
                "cmake_minimum_required(VERSION 3.24)\nproject(GateProbe NONE)\n"
                "add_subdirectory(nested)\nif(TRUE)\n" + gate +
                '\nget_target_property(deps dynsyslab_test_binaries MANUALLY_ADDED_DEPENDENCIES)\n'
                'file(WRITE "${CMAKE_BINARY_DIR}/dependencies.txt" "${deps}")\n',
                encoding="utf-8")
            build = Path(temporary) / "build"
            result = subprocess.run(["cmake", "-S", str(source), "-B", str(build)],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            dependencies = set((build / "dependencies.txt").read_text().split(";"))
            self.assertFalse(built_tests - dependencies, built_tests - dependencies)

    def test_eigen_is_private_to_explicit_native_numerical_units(self):
        native = PROJECT_ROOT / "native"
        public_headers = list((native / "include").rglob("*.hpp"))
        leaking_headers = [
            path.relative_to(PROJECT_ROOT).as_posix()
            for path in public_headers
            if re.search(r"#\s*include\s*[<\"](?:unsupported/)?Eigen/",
                         path.read_text(encoding="utf-8"))
        ]
        self.assertEqual(leaking_headers, [])
        allowed = {
            "native/src/EigenPch.hpp",
            # Matrix-free Hankel kernels; Eigen stays inside these private units.
            "native/src/HistoryKoopman.cpp",
            "native/src/HistoryCorrelation.hpp",
            # Randomized Hankel operator (dynsyslab_observation); Eigen-free header.
            "native/src/HankelOperator.cpp",
            # Private Eigen helpers of the linear-systems kernels in dynsyslab_numerics.
            "native/src/EigenNumericsLinearSystems.hpp",
            # Dense conic solver and bilinear K-ROM units; their public headers are Eigen-free.
            "native/src/DenseSdp.cpp",
            "native/src/BilinearKrom.cpp",
            "native/src/CertifiedDesign.cpp",
        }
        allowed.update(
            path.relative_to(PROJECT_ROOT).as_posix()
            for path in (native / "src").glob("EigenNumerics*.cpp")
        )
        exposed = {
            path.relative_to(PROJECT_ROOT).as_posix()
            for suffix in ("*.cpp", "*.hpp")
            for path in (native / "src").glob(suffix)
            if re.search(r"#\s*include\s*[<\"](?:unsupported/)?Eigen/",
                         path.read_text(encoding="utf-8"))
        }
        self.assertEqual(exposed, allowed)
        cmake = (native / "CMakeLists.txt").read_text(encoding="utf-8")
        self.assertNotRegex(cmake,
            r"target_link_libraries\([^)]*PUBLIC[^)]*Eigen3::Eigen")
        self.assertEqual(set(re.findall(r"dynsyslab_eigen_module\((\w+)\)", cmake)), {
            "dynsyslab_numerics", "dynsyslab_observation", "dynsyslab_history_koopman",
            "dynsyslab_control_certificates"})
        self.assertIn("set(DYNSYSLAB_EIGEN_COMPILE_JOBS 3 CACHE STRING", cmake)

    def test_native_sources_use_narrow_numerical_interfaces(self):
        native_sources = list((PROJECT_ROOT / "native" / "src").glob("*.cpp"))
        aggregate_users = [
            path.relative_to(PROJECT_ROOT).as_posix()
            for path in native_sources
            if "dynsyslab/numerics/LinearAlgebra.hpp" in
            path.read_text(encoding="utf-8")
        ]
        self.assertEqual(aggregate_users, [])

    def test_distribution_manifest_is_client_only(self):
        manifest = (PROJECT_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
        self.assertIn("include docs/START_HERE.md", manifest)
        self.assertIn("prune archives", manifest)
        self.assertNotIn("include archives", manifest)

    def test_distribution_contains_only_the_optional_native_client(self):
        project = (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        self.assertIn('packages = ["dynsyslab_client"]', project)
        self.assertIn("include-package-data = false", project)
        self.assertIn("dependencies = []", project)
        for package in ("archiver*", "compute_workbench*", "estimation*",
                        "job_manager*", "postprocessing*", "scripts*"):
            self.assertNotIn(f'"{package}"', project)

    def test_optional_python_client_does_not_import_scientific_runtimes(self):
        script = (
            "import json,sys; import dynsyslab_client; "
            "from dynsyslab_client.cli import diagnostics; "
            "assert 'numpy' not in sys.modules; assert 'scipy' not in sys.modules; "
            "print(json.dumps(diagnostics()))"
        )
        with TemporaryDirectory() as directory:
            completed = subprocess.run(
                (sys.executable, "-c", script), cwd=directory,
                env={**dict(os.environ), "PYTHONPATH": str(PROJECT_ROOT)},
                capture_output=True, text=True, timeout=30, check=False,
            )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(completed.stdout)
        self.assertEqual(result["schema"],
                         "dynsyslab.python-client-diagnostics/v1")
        self.assertFalse(result["python_required_by_native"])

    def test_source_events_are_viewer_independent_and_passive(self):
        source = (PROJECT_ROOT / "infrastructure" /
                  "source_events.py").read_text(encoding="utf-8")
        self.assertNotIn("from viewer", source)
        self.assertNotIn("import viewer", source)
        with TemporaryDirectory() as directory, patch.object(
                source_events, "ROOT", Path(directory)):
            observed = Path(directory) / "solution.npz"
            source_events.register_sources("session", [{
                "name": "test", "path": str(observed.resolve()),
            }])
            self.assertEqual(source_events.notify_source_update(
                observed, source_type="solution", reason="test"), 1)
            events, cursor = source_events.poll_events("session")
            self.assertEqual(cursor, 1)
            self.assertEqual(events[0]["reason"], "test")


if __name__ == "__main__":
    unittest.main()
