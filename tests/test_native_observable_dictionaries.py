"""Independent formula, CLI and native workflow checks for observable dictionaries."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from scipy import special
from dynsyslab_client.workflow import NativeWorkflow, NativeWorkflowError

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("DYNSYSLAB_NATIVE_BIN", ROOT / "build/native-suite-fast/bin"))
LIB = BIN / ("dynsyslab_workflow_c_api.dll" if os.name == "nt" else "libdynsyslab_workflow_c_api.so")


@unittest.skipUnless(LIB.is_file(), "Build the native workflow first")
class ObservableDictionaryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.api = NativeWorkflow(LIB)
        cls.x = np.random.default_rng(42).uniform(-1, 1, (96, 2))

    def dictionary(self, kind, **settings):
        return self.api.koopman_spectral({"action": "dictionary",
            "data": {"dimension": 2, "x": self.x.ravel().tolist()},
            "settings": {"dictionary": {"kind": kind, **settings}}})

    def test_product_formulas_and_frozen_evaluation(self):
        for kind, function in (("chebyshev", special.eval_chebyt), ("laguerre", special.eval_laguerre)):
            with self.subTest(kind=kind):
                result = self.dictionary(kind, degree=5, index_set="tensor_product", scale=[2.5, 2.5])
                d = result["dictionary"]
                indices = np.array(d["exponents"]).reshape(-1, 2)
                expected = np.prod(function(indices[None, :, :], self.x[:, None, :]/2.5), axis=2)
                np.testing.assert_allclose(np.array(result["values"]).reshape(expected.shape), expected, atol=2e-14)
                repeat = self.api.koopman_spectral({"action": "dictionary", "dictionary": d, "points": self.x.ravel().tolist()})
                self.assertEqual(repeat["values"], result["values"])
        h = self.api.koopman_spectral({"action": "dictionary",
            "dictionary": {"kind": "hermite", "dimension": 1, "exponents": [0, 1, 2, 20, 100]},
            "points": [-4., -.3, 0., 5.]})
        orders = np.array([0, 1, 2, 20, 100])
        points = np.array([-4., -.3, 0., 5.])[:, None]
        expected = special.eval_hermite(orders, points)*np.exp(-points**2/2-
            (orders*np.log(2)+special.gammaln(orders+1))/2)/np.pi**.25
        np.testing.assert_allclose(np.array(h["values"]).reshape(expected.shape), expected, atol=2e-14)

    def test_radial_formulas(self):
        for kind in ("gaussian", "laplacian", "lorentzian", "matern"):
            with self.subTest(kind=kind):
                result = self.dictionary(kind, centers=7, bandwidth=.8, smoothness=2., amplitude=2., center_selection="kmeans")
                nodes = np.array(result["dictionary"]["nodes"]).reshape(-1, 2)
                z = np.linalg.norm(self.x[:, None, :]-nodes[None, :, :], axis=2)/.8
                if kind == "gaussian": expected = 2*np.exp(-z*z/2)
                elif kind == "laplacian": expected = 2*np.exp(-z)
                elif kind == "lorentzian": expected = 2/(1+z*z)
                else:
                    expected = np.full_like(z, 2.)
                    positive = z > 0
                    expected[positive] = z[positive]**2*special.kv(2, z[positive])
                np.testing.assert_allclose(np.array(result["values"]).reshape(z.shape), expected, atol=1e-13)

    def test_fit_and_archive_for_each_family(self):
        with tempfile.TemporaryDirectory() as folder:
            for kind in ("chebyshev", "laguerre", "hermite", "fourier_hermite", "laplacian", "lorentzian", "matern"):
                with self.subTest(kind=kind):
                    path = str(Path(folder) / (kind+".npz"))
                    result = self.api.koopman_spectral({"action": "fit",
                        "data": {"dimension": 2, "x": self.x.ravel().tolist(), "y": (.9*self.x).ravel().tolist()},
                        "settings": {"method": "resdmd", "rank": 5, "dictionary": {
                            "kind": kind, "degree": 2, "centers": 8, "angular_order": 1, "hermite_order": 2}},
                        "output_archive": path})
                    loaded = self.api.koopman_spectral({"action": "inspect", "archive": path})
                    self.assertEqual(result["model_identity"], loaded["model_identity"])
                    evaluated = self.api.koopman_spectral({"action": "evaluate", "archive": path, "points": [.1, .2]})
                    self.assertTrue(np.isfinite(evaluated["values"]).all())

    def test_complex_observable_api_and_component_flow(self):
        for system in ("duffing", "nonlinear_pendulum"):
            query = {"system": {"name": system}, "observable": "phase_sech", "state": [.7, -.3]}
            result = self.api.observable(query)
            self.assertEqual(result["value_type"], "complex")
            self.assertEqual(np.array(result["jacobian"]).shape, (2, 2))
            components = self.api.observable({**query, "observable": "phase_sech_components"})
            np.testing.assert_allclose(result["values"][0], components["values"])
            rate = self.api.observable({**query, "observable": "phase_sech_components", "flow_order": 1})
            self.assertTrue(np.isfinite(rate["values"]).all())
        with self.assertRaises(NativeWorkflowError):
            self.dictionary("fourier_hermite", normalize=True)
        with self.assertRaises(NativeWorkflowError):
            self.dictionary("chebyshev", normalize=True, scale=[2, 2])

    def test_cli_and_pipeline_dictionary(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            query = root / "query.json"
            query.write_text(json.dumps({"action": "dictionary", "data": {"dimension": 2, "x": [0, 0, 1, 1]},
                "settings": {"dictionary": {"kind": "laguerre", "degree": 12}}}))
            cli = BIN / ("dynsyslab-cli.exe" if os.name == "nt" else "dynsyslab-cli")
            result = subprocess.run([str(cli), "koopman-spectral", str(query)], check=True, capture_output=True, text=True)
            self.assertEqual(json.loads(result.stdout)["features"], 91)
            source = root / "series.npz"
            t = np.arange(150)*.05
            np.savez(source, time=t, values=np.column_stack((np.cos(t), np.sin(t))),
                metadata_json=np.array(json.dumps({"coordinate_system": "native"})))
            job = {"schema": "dynsyslab.compute-job/v1", "kind": "postprocess",
                "label": "dictionary test", "sources": [str(source)], "methods": [
                    {"name": "resdmd", "settings": {"rank": 5, "dictionary": {"kind": "chebyshev", "degree": 2}}}]}
            self.api.plan(job)
            self.api.run(job, root / "archives")
            fitted = []
            for p in (root / "archives").rglob("*.npz"):
                with np.load(p) as archive:
                    if "spectral_krom_contract_json" in archive.files: fitted.append(p)
            self.assertEqual(len(fitted), 1)
            model = self.api.koopman_spectral({"action": "inspect", "archive": str(fitted[0])})
            self.assertEqual(model["model"]["dictionary"]["kind"], "chebyshev")

    def test_named_observables_taps_and_derivatives_in_hankel(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "pendulum.npz"
            t = np.arange(80)*.02
            states = np.column_stack((.7*np.cos(t), -.7*np.sin(t)))
            np.savez(source, time=t, values=states, metadata_json=np.array(json.dumps({
                "system": "nonlinear_pendulum", "archive_coordinate_system": "dimensional"})))
            job = {"schema": "dynsyslab.compute-job/v1", "kind": "postprocess",
                "label": "named observable test", "sources": [str(source)], "methods": [
                    {"name": "hankel_matrix", "settings": {"observables": [
                        {"observable": "phase_sech_real", "taps": [0, 2]},
                        {"observable": "phase_sech_imag", "derivative_order": 1}],
                        "window": {"delays": 3, "columns": 20}}}]}
            self.api.plan(job)
            self.api.run(job, root / "archives")
            matrices = []
            for path in (root / "archives").rglob("*.npz"):
                with np.load(path) as archive:
                    if "hankel_matrix" in archive.files: matrices.append(archive["hankel_matrix"].copy())
            self.assertEqual(len(matrices), 1)
            matrix = matrices[0]
            self.assertEqual(matrix.shape, (9, 20))
            for block in range(3):
                for col in (0, 7, 19):
                    current = block+col+2
                    q = {"system": {"name": "nonlinear_pendulum"}, "state": states[current].tolist()}
                    real = self.api.observable({**q, "observable": "phase_sech_real"})["values"][0]
                    delayed = self.api.observable({**q, "state": states[current-2].tolist(),
                        "observable": "phase_sech_real"})["values"][0]
                    rate = self.api.observable({**q, "observable": "phase_sech_imag", "flow_order": 1})["values"][0]
                    np.testing.assert_allclose(matrix[3*block:3*block+3, col], [real, delayed, rate], atol=2e-13)


if __name__ == "__main__":
    unittest.main()
