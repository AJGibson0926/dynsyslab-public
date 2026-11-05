"""Regression tests for the safe, read-only MCP harness surface."""

from __future__ import annotations

import unittest

from codex_harness.server import project_structure, read_project_document


class CodexHarnessTests(unittest.TestCase):
    def test_structure_lists_source_and_hides_archive_payloads(self) -> None:
        packages = project_structure()["source_packages"]
        self.assertIn("solver", packages)
        self.assertNotIn("archives", packages)

    def test_document_reader_rejects_arbitrary_paths(self) -> None:
        self.assertIn("error", read_project_document("../config/antenna_policy.json"))

    def test_document_reader_allows_public_architecture_document(self) -> None:
        result = read_project_document("ARCHITECTURE.md")
        self.assertEqual(result["document"], "ARCHITECTURE.md")
        self.assertTrue(result["content"])
