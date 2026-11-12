from __future__ import annotations
import tempfile, unittest
from pathlib import Path
import numpy as np
from archiver.gateway import ArchiveAccessError, ArchiveGateway

class ArchiveGatewayContractTests(unittest.TestCase):
    def test_index_refresh_requires_management_and_preserves_selected_root(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            with patch("archives.build_index", return_value={"solutions": 3}) as refresh:
                for capability in ("read", "write"):
                    with self.assertRaises(ArchiveAccessError):
                        ArchiveGateway(capability, root).refresh_index()
                refresh.assert_not_called()
                self.assertEqual(ArchiveGateway("manage", root).refresh_index(), {"solutions": 3})
                refresh.assert_called_once_with(reconcile=True, archive_root=root)

    def test_read_gateway_cannot_create_a_session(self):
        with tempfile.TemporaryDirectory() as directory:
            gateway=ArchiveGateway("read",Path(directory))
            with self.assertRaises(ArchiveAccessError): gateway.new_session(object())
    def test_manage_gateway_can_discover_empty_sections(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(ArchiveGateway("manage",Path(directory)).iter_paths("solutions"),[])
    def test_read_gateway_cannot_update_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ArchiveAccessError): ArchiveGateway("read",Path(directory)).update_metadata(Path(directory)/"x.npz",{})
    def test_open_read_rejects_missing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                with ArchiveGateway("read",Path(directory)).open_read(Path(directory)/"missing.npz"): pass

    def test_archive_root_is_resolved_absolutely(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertTrue(ArchiveGateway("read",Path(directory)).archive_root.is_absolute())

    def test_write_gateway_can_discover_an_empty_section(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(ArchiveGateway("write",Path(directory)).iter_paths("history"),[])

    def test_read_gateway_cannot_create_hierarchy_archive(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ArchiveAccessError): ArchiveGateway("read",Path(directory)).save_hierarchy(object())

    def test_read_gateway_cannot_extend_a_solution(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ArchiveAccessError): ArchiveGateway("read",Path(directory)).extend_solution("a.npz","b.npz")

    def test_open_read_exposes_existing_arrays_without_write_permission(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"sample.npz"; np.savez(path,values=np.array([1.,2.]))
            with ArchiveGateway("read",Path(directory)).open_read(path) as archive:
                self.assertTrue(np.array_equal(archive["values"],[1.,2.]))

    def test_metadata_and_solution_discovery_do_not_materialize_solution_arrays(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory).resolve(); target=root/"solutions"/"finished.npz"; target.parent.mkdir()
            np.savez(target,metadata_json=np.array('{"label":"finished","status":"complete"}'),values=np.arange(100))
            partial=root/"solutions"/"active.partial.npz"
            np.savez(partial,metadata_json=np.array('{"label":"active"}'),values=np.arange(2))
            gateway=ArchiveGateway("read",root)
            self.assertEqual(gateway.metadata(target)["label"],"finished")
            self.assertEqual(gateway.solution_paths(),[target])
            self.assertEqual(gateway.partial_solution_paths(),[partial])

    def test_source_discovery_excludes_generated_build_output(self):
        from project_paths import PROJECT_ROOT, is_generated, source_files
        self.assertTrue(is_generated(Path(__file__).resolve().parents[1]/"build"/"lib"/"x.py"))
        self.assertFalse(any("build" in path.relative_to(PROJECT_ROOT).parts for path in source_files()))

    def test_read_gateway_reads_small_manifest_without_write_capability(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"campaign.json"; path.write_text('{"campaign_id":"demo"}',encoding="utf-8")
            self.assertEqual(ArchiveGateway("read",Path(directory)).read_json(path)["campaign_id"],"demo")
