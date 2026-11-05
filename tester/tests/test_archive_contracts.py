from __future__ import annotations
import tempfile
import unittest
import sqlite3
from contextlib import closing
from pathlib import Path
import numpy as np
from archives.store import _atomic_npz, _uniform_time_steps
from archives.store import extend_solution, load_time_series
from time_series import TimeSeries
from archives.schema import ArchiveSchemaError, require_format, require_keys
from archives.index import index_archive_object

class ArchiveContractTests(unittest.TestCase):
    def test_atomic_npz_publish_is_readable(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"solution.npz"
            _atomic_npz(path,time=np.array([0.,1.]),values=np.array([[1.],[2.]]))
            with np.load(path,allow_pickle=False) as record:
                self.assertEqual(record["values"].shape,(2,1))
            self.assertFalse(any(Path(directory).glob("*.partial")))

    def test_uniform_grid_indicator(self):
        self.assertTrue(_uniform_time_steps(np.array([0.,.1,.2])))
        self.assertFalse(_uniform_time_steps(np.array([0.,.1,.25])))

    def test_versioned_record_contract_rejects_missing_or_wrong_formats(self):
        require_format({"format":"dynsyslab.job/v1"},"dynsyslab.job/")
        with self.assertRaises(ArchiveSchemaError): require_format({},"dynsyslab.job/")
        with self.assertRaises(ArchiveSchemaError): require_keys({"one":1},"one","two")
        with self.assertRaises(ArchiveSchemaError): require_format({"format":"dynsyslab.job/v99"},"dynsyslab.job/")

    def test_extension_keeps_a_strictly_increasing_combined_time_series(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"solution.npz"
            _atomic_npz(path,time=np.array([0.,1.]),values=np.array([[0.],[1.]]),metadata_json="{}")
            extend_solution(path,TimeSeries(np.array([1.,2.]),np.array([[1.],[2.]])))
            combined=load_time_series(path)
            self.assertTrue(np.array_equal(combined.time,np.array([0.,1.,2.])))
            self.assertTrue(np.all(np.diff(combined.time)>0))

    def test_extension_rejects_a_gap_in_state_dimension(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"solution.npz"; _atomic_npz(path,time=np.array([0.,1.]),values=np.array([[0.],[1.]]),metadata_json="{}")
            with self.assertRaises(ValueError): extend_solution(path,TimeSeries(np.array([1.,2.]),np.array([[1.,2.],[2.,3.]])))

    def test_uniform_indicator_accepts_single_sample(self):
        self.assertTrue(_uniform_time_steps(np.array([0.])))

    def test_archive_publication_invalidates_catalog_only_on_change(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); path=root/"derived"/"record.npz"
            _atomic_npz(path,values=np.array([1.]),metadata_json='{"label":"one"}')
            self.assertTrue(index_archive_object(path,archive_root=root))
            database=root/"management"/"archive_index.sqlite"
            with closing(sqlite3.connect(database)) as connection:
                generation=connection.execute(
                    "SELECT generation FROM archive_index_state").fetchone()[0]
                connection.execute("INSERT INTO archive_catalog_snapshots VALUES(?,?,?,?)",
                                   ("test",generation,1,"{}"))
                connection.commit()
            self.assertTrue(index_archive_object(path,archive_root=root))
            with closing(sqlite3.connect(database)) as connection:
                self.assertEqual(connection.execute(
                    "SELECT generation FROM archive_index_state").fetchone()[0],generation)
                self.assertEqual(connection.execute(
                    "SELECT count(*) FROM archive_catalog_snapshots").fetchone()[0],1)
            _atomic_npz(path,values=np.array([2.]),metadata_json='{"label":"two"}')
            self.assertTrue(index_archive_object(path,archive_root=root))
            with closing(sqlite3.connect(database)) as connection:
                self.assertEqual(connection.execute(
                    "SELECT generation FROM archive_index_state").fetchone()[0],generation+1)
                self.assertEqual(connection.execute(
                    "SELECT count(*) FROM archive_catalog_snapshots").fetchone()[0],0)
