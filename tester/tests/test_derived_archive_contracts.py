from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import zipfile

import numpy as np

from archiver.gateway import ArchiveGateway
from postprocessing.core import (
    DERIVED_DEFAULT_COMPRESSION,
    Trajectory,
    _archive_compression,
    _write_derived_npz,
    discover_derived,
    load_derived_array,
    recompress_derived,
    save_derived,
)
from postprocessing.koopman.pipeline import _storage_arrays
from postprocessing.koopman.dictionary import PolynomialDictionary
from time_series import TimeSeries


class DerivedArchiveContractTests(unittest.TestCase):
    def test_immutable_publication_preserves_prior_results_and_is_discoverable(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source_path=root/"source.npz"
            source=Trajectory(TimeSeries(np.array((0.,1.)),np.array(((1.,),(2.,)))),{"system":"test"},source_path)
            with patch("postprocessing.core.POSTPROCESSING_ARCHIVE",ArchiveGateway("write",root)):
                first=save_derived(source,"unit-test",{"rank":2},values=np.array((1.,2.)))
                second=save_derived(source,"unit-test",{"rank":2},values=np.array((3.,4.)))
                records=discover_derived(analysis="unit-test",source_archive=source_path)
            self.assertNotEqual(first,second)
            self.assertTrue(first.exists() and second.exists())
            self.assertEqual([item.path for item in records],[second.resolve(),first.resolve()])
            self.assertEqual(records[0].analysis_signature_sha256,records[1].analysis_signature_sha256)
            self.assertNotEqual(records[0].publication_id,records[1].publication_id)
            with np.load(first,allow_pickle=False) as data: metadata=json.loads(str(data["metadata_json"]))
            self.assertEqual(metadata["format"],"dynsyslab.derived/v1")
            self.assertTrue(metadata["publication"]["immutable"])

    def test_replacement_requires_explicit_publication_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source_path=root/"source.npz"
            source=Trajectory(TimeSeries(np.array((0.,1.)),np.array(((1.,),(2.,)))),{"system":"test"},source_path)
            with patch("postprocessing.core.POSTPROCESSING_ARCHIVE",ArchiveGateway("write",root)):
                first=save_derived(source,"unit-test",publication="replace",values=np.array((1.,)))
                second=save_derived(source,"unit-test",publication="replace",values=np.array((2.,)))
            self.assertEqual(first,second)
            with np.load(second,allow_pickle=False) as data:
                np.testing.assert_array_equal(data["values"],np.array((2.,)))
                self.assertFalse(json.loads(str(data["metadata_json"]))["publication"]["immutable"])

    def test_default_derived_codec_is_fast_deflate_and_metadata_is_truthful(self):
        self.assertEqual(DERIVED_DEFAULT_COMPRESSION,"deflate")
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source_path=root/"source.npz"
            source=Trajectory(TimeSeries(np.array((0.,1.)),np.array(((1.,),(2.,)))),{"system":"test"},source_path)
            with patch("postprocessing.core.POSTPROCESSING_ARCHIVE",ArchiveGateway("write",root)):
                path=save_derived(source,"unit-test",values=np.arange(20,dtype=float))
            self.assertEqual(_archive_compression(path),"deflate")
            with np.load(path,allow_pickle=False) as data:
                metadata=json.loads(str(data["metadata_json"]))
            self.assertEqual(metadata["storage_compression"],"npz-deflate")

    def test_recompression_updates_codec_provenance_without_changing_values(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"derived.npz"; values=np.linspace(0.,1.,100)
            _write_derived_npz(path,{"values":values,"metadata_json":json.dumps({"storage_compression":"npz-deflate"})})
            recompress_derived(path,algorithm="lzma")
            self.assertEqual(_archive_compression(path),"lzma")
            with np.load(path,allow_pickle=False) as data:
                np.testing.assert_array_equal(data["values"],values)
                self.assertEqual(json.loads(str(data["metadata_json"]))["storage_compression"],"npz-lzma")

    def test_explicit_array_alias_resolves_without_duplicate_member(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"derived.npz"; modes=np.arange(12).reshape(4,3)
            metadata={"analysis_metadata":{"storage":{"array_aliases":{"eigenfunction_values":"modes"}}}}
            _write_derived_npz(path,{"modes":modes,"metadata_json":json.dumps(metadata)})
            with zipfile.ZipFile(path) as archive: self.assertNotIn("eigenfunction_values.npy",archive.namelist())
            np.testing.assert_array_equal(load_derived_array(path,"eigenfunction_values"),modes)

    def test_virtual_polynomial_eigenfunctions_reconstruct_only_requested_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"coefficient-first.npz"
            states=np.array(((1.,2.),(2.,4.),(3.,5.)))
            center=np.array((1.,2.)); scale=np.array((2.,3.)); multiplier=np.array((.5,2.))
            transformed=(states-center)/scale*multiplier
            features=PolynomialDictionary(degree=2).transform(transformed)
            coefficients=(np.arange(18).reshape(3,6)+1j*np.arange(18,36).reshape(3,6))/10.
            specification={"encoding":"polynomial-dictionary-evaluation/v1",
                "states":"training_states","coefficients":"coefficients",
                "dictionary":{"encoding":"polynomial-dictionary/v1","class":"PolynomialDictionary",
                              "degree":2,"include_constant":True},
                "state_preprocessing":{"center":center.tolist(),"scale":scale.tolist(),
                                       "multiplier":multiplier.tolist()},
                "shape":[3,3],"dtype":"complex128"}
            metadata={"analysis_metadata":{"storage":{"virtual_arrays":{
                "eigenfunction_values":specification}}}}
            _write_derived_npz(path,{"training_states":states,"coefficients":coefficients,
                                     "metadata_json":json.dumps(metadata)})
            expected=features@coefficients.T
            np.testing.assert_allclose(load_derived_array(path,"eigenfunction_values"),expected)
            np.testing.assert_allclose(load_derived_array(path,"eigenfunction_values",columns=[2,0]),
                                       expected[:,[2,0]])
            with self.assertRaises(IndexError):
                load_derived_array(path,"eigenfunction_values",columns=[3])

    def test_direct_and_aliased_arrays_support_column_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"columns.npz"; values=np.arange(20).reshape(5,4)
            metadata={"analysis_metadata":{"storage":{"array_aliases":{"alias":"values"}}}}
            _write_derived_npz(path,{"values":values,"metadata_json":json.dumps(metadata)})
            np.testing.assert_array_equal(load_derived_array(path,"values",columns=[3,1]),values[:,[3,1]])
            np.testing.assert_array_equal(load_derived_array(path,"alias",columns=[0,2]),values[:,[0,2]])

    def test_hankel_identity_storage_uses_an_explicit_alias(self):
        values=np.arange(12,dtype=float).reshape(4,3)
        source=SimpleNamespace(values=np.ones((4,2)))
        result=SimpleNamespace(eigenvalues=np.arange(3),operator=np.eye(3),modes=values,
                               eigenfunction_values=values,spectrum=SimpleNamespace(eigenfunctions=[]))
        arrays,aliases=_storage_arrays(source,result)
        self.assertNotIn("eigenfunction_values",arrays)
        self.assertEqual(aliases,{"eigenfunction_values":"modes"})


if __name__=="__main__":
    unittest.main()
