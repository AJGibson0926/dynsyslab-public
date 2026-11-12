from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest

import numpy as np

from dynsyslib.legacy_dynbatch_csv import (
    export_legacy_dynbatch_csv,
    read_legacy_dynbatch_csv,
)


class LegacyDynBatchCSVTests(unittest.TestCase):
    def _write_fixture(self, root: Path) -> Path:
        path = root / "legacy.csv"
        rows = [
            ["Var1", "Var2_1", "Var2_2", "Var2_3"],
            ["dynbatch", "2", "", ""],
        ]
        for index, initial in ((1, (0.2, 0.9)), (2, (-0.4, -0.1))):
            rows.extend(
                [
                    ["dynsys", str(index), "Duffing", "0.1", "2", "2"],
                    ["time", "0", "0.1", "0.2"],
                    ["ctrl", "random", "uniform", ""],
                    ["u", "0", str(index), str(index + 0.5)],
                    ["dynclst", "1", "", ""],
                    ["inst", "1", "", ""],
                    ["time", "0", "0.1", "0.2"],
                    ["sol", str(initial[0]), str(initial[0] + 0.1), str(initial[0] + 0.2)],
                    ["", str(initial[1]), str(initial[1] - 0.1), str(initial[1] - 0.2)],
                ]
            )
        with path.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerows(rows)
        return path

    def test_decodes_hierarchy_and_removes_control_initializer_for_regression(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = self._write_fixture(Path(temporary))
            decoded = read_legacy_dynbatch_csv(source)
            self.assertEqual(decoded.declared_system_count, 2)
            self.assertEqual(len(decoded.systems), 2)
            self.assertEqual(decoded.trajectory_count, 2)
            np.testing.assert_allclose(decoded.systems[0].trajectories[0].initial_condition, (0.2, 0.9))
            np.testing.assert_allclose(decoded.systems[0].control_history.values[:, 0], (0.0, 1.0, 1.5))
            controlled = decoded.controlled_dataset(input_scale=0.5)
            np.testing.assert_allclose(controlled.trajectories[0].inputs[:, 0], (0.5, 0.75))
            self.assertEqual(controlled.transition_count, 4)
            pooled = decoded.pooled_time_series()
            self.assertEqual(pooled.values.shape, (3, 2, 2))

    def test_exports_numeric_archive_and_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoded = read_legacy_dynbatch_csv(self._write_fixture(root))
            archive, manifest = export_legacy_dynbatch_csv(
                decoded, root / "decoded", include_controlled_dataset=True
            )
            self.assertTrue(archive.is_file())
            self.assertTrue(manifest.is_file())
            with np.load(archive) as values:
                self.assertEqual(values["controlled_states"].shape, (2, 3, 2))
                self.assertEqual(values["controlled_interval_inputs"].shape, (2, 2, 1))

    def test_rejects_unattached_continuation_rows(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = self._write_fixture(Path(temporary))
            with source.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.reader(stream))
            rows.insert(2, ["", "1", "2", "3"])
            with source.open("w", newline="", encoding="utf-8") as stream:
                csv.writer(stream).writerows(rows)
            with self.assertRaisesRegex(ValueError, "Unexpected continuation"):
                read_legacy_dynbatch_csv(source)


if __name__ == "__main__":
    unittest.main()
