"""Scientific and boundary regression for the compiled spectral Koopman API."""
from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from dynsyslab_client.workflow import NativeWorkflow, NativeWorkflowError

ROOT = Path(__file__).resolve().parents[1]
BIN = Path(os.environ.get("DYNSYSLAB_NATIVE_BIN", ROOT / "build/native-suite-fast/bin"))
LIB = BIN / ("dynsyslab_workflow_c_api.dll" if os.name == "nt" else "libdynsyslab_workflow_c_api.so")


@unittest.skipUnless(LIB.is_file(), "Build the native workflow first")
class SpectralKoopmanTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.api = NativeWorkflow(LIB)
        theta = np.arange(128) * 2 * np.pi / 128
        cls.x = np.column_stack((np.cos(theta), np.sin(theta)))
        cls.y = np.column_stack((np.cos(theta + .31), np.sin(theta + .31)))

    def fit(self, method="resdmd", **extra):
        settings = {"method": method, "dictionary": {"kind": "linear", "constant": False},
                    "sample_interval": .1, "assume_measure_preserving": True}
        settings.update(extra.pop("settings", {}))
        return self.api.koopman_spectral({"action": "fit", "data": {
            "dimension": 2, "x": self.x.ravel().tolist(), "y": self.y.ravel().tolist()},
            "settings": settings, **extra})

    def test_wrapper_rollout_and_archive_reference(self):
        with tempfile.TemporaryDirectory() as folder:
            path = str(Path(folder) / "rotation.npz")
            fit = self.fit(output_archive=path)
            result = self.api.koopman_spectral({"action": "rollout", "archive": path,
                "initial": [1, 0], "steps": 50})
            expected = np.column_stack((np.cos(.31*np.arange(51)), np.sin(.31*np.arange(51))))
            np.testing.assert_allclose(np.array(result["values"]).reshape(-1, 2), expected, atol=2e-12)
            loaded = self.api.koopman_spectral({"action": "inspect", "archive": path})
            self.assertEqual(loaded["model_identity"], fit["model_identity"])
            reference = {"artifact": path, "model_identity": fit["model_identity"], "component": "coordinates"}
            value = self.api.koopman_spectral({"action": "observable", "reference": reference, "points": [.2, .7]})
            np.testing.assert_allclose(np.array(value["values"])[:, 0], [.2, .7], atol=1e-12)
            reference["model_identity"] = "wrong-model"
            with self.assertRaises(NativeWorkflowError):
                self.api.koopman_spectral({"action": "observable", "reference": reference, "points": [.2, .7]})
            query = Path(folder) / "query.json"
            query.write_text(json.dumps({"action": "inspect", "archive": path, "include_model": False}))
            cli = BIN / ("dynsyslab-cli.exe" if os.name == "nt" else "dynsyslab-cli")
            result = subprocess.run([str(cli), "koopman-spectral", str(query)], check=True, capture_output=True, text=True)
            self.assertEqual(json.loads(result.stdout)["model_identity"], fit["model_identity"])
            info = subprocess.run([str(cli), "krom-info", path], check=True, capture_output=True, text=True)
            self.assertEqual(json.loads(info.stdout)["model_identity"], fit["model_identity"])

    def test_independent_weighted_edmd_normal_equations(self):
        rng = np.random.default_rng(14)
        x = rng.normal(size=(90, 2))
        y = np.column_stack((.9*x[:, 0]+.1*x[:, 1]**2, .7*x[:, 1]))
        w = np.linspace(.2, 2, len(x))
        result = self.api.koopman_spectral({"action": "fit", "data": {"dimension": 2,
            "x": x.ravel().tolist(), "y": y.ravel().tolist(), "weights": w.tolist()},
            "settings": {"method": "resdmd", "dictionary": {"kind": "linear", "constant": False}}})
        model = result["model"]
        basis = np.array(model["basis"]).reshape(2, 2)
        operator = basis @ np.array(model["action"]).reshape(2, 2) @ np.linalg.inv(basis)
        direct = np.linalg.solve(x.T @ (w[:, None]*x), x.T @ (w[:, None]*y))
        np.testing.assert_allclose(operator, direct, atol=1e-12)
        c = np.array(model["eigenfunctions"])
        c = (c[:, 0]+1j*c[:, 1]).reshape(2, 2)
        lam = np.array(model["multipliers"])
        lam = lam[:, 0]+1j*lam[:, 1]
        residual = np.sqrt(np.sum(w[:, None]*abs(y@c-(x@c)*lam)**2, axis=0)/
                           np.sum(w[:, None]*abs(x@c)**2, axis=0))
        np.testing.assert_allclose(result["spectral"]["residuals"]["eigenpair"], residual, atol=1e-12)

    def test_shared_hankel_postprocessor(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "trajectory.npz"
            theta = .31*np.arange(160)
            values = np.column_stack((np.cos(theta), np.sin(theta)))
            np.savez(source, time=.1*np.arange(len(values)), values=values,
                     metadata_json=np.array(json.dumps({"coordinate_system": "native"})))
            job = {"schema": "dynsyslab.compute-job/v1", "kind": "postprocess",
                "label": "spectral regression", "sources": [str(source)], "methods": [
                    {"name": "hankel_matrix", "settings": {"observables": [0, 1],
                        "window": {"delays": 3, "columns": 100}}},
                    {"name": "resdmd", "settings": {"use_previous_hankel": True,
                        "dictionary": {"kind": "linear", "constant": False}, "rank": 2}}]}
            self.api.plan(job)
            result = self.api.run(job, root / "archives")
            paths = list((root / "archives").rglob("*.npz"))
            fitted = [p for p in paths if "spectral_krom_contract_json" in np.load(p).files]
            self.assertEqual(len(fitted), 1, result)
            model = self.api.koopman_spectral({"action": "inspect", "archive": str(fitted[0])})
            self.assertEqual(model["rank"], 2)
            self.assertLess(max(model["spectral"]["residuals"]["eigenpair"]), 1e-10)
            self.assertAlmostEqual(model["model"]["sample_interval"], .1)

    @unittest.skipUnless((ROOT / "build/colbrook-reference-parity.mat").is_file(), "Run pinned MATLAB reference script")
    def test_pinned_matlab_mp_and_rigged_parity(self):
        from scipy.io import loadmat
        reference = loadmat(ROOT / "build/colbrook-reference-parity.mat")
        angles = reference["TH"].ravel().tolist()
        result = self.fit("mpedmd", settings={"rigged_measures": [
            {"name": f"cos{order}", "coefficients": [1, 0], "angles": angles, "bandwidth": .08, "order": order}
            for order in (1, 4)]})
        model = result["model"]
        basis = np.array(model["basis"]).reshape(2, 2)
        action = np.array(model["action"]).reshape(2, 2)
        np.testing.assert_allclose(basis@action@np.linalg.inv(basis), reference["mpK"], atol=2e-12)
        for index, order in enumerate((1, 4)):
            measure = model["spectral"]["measures"][index]
            np.testing.assert_allclose(measure["density"], reference[f"xi{order}"].ravel(), atol=2e-10)
            coeff = np.array(model["spectral"]["wave_packets"][index]["coefficients"])
            coeff = (coeff[:, 0]+1j*coeff[:, 1]).reshape(-1, 2)
            np.testing.assert_allclose(coeff.T, reference[f"g{order}"], atol=2e-10)

if __name__ == "__main__":
    unittest.main()
