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

    def test_cyclic_state_pair_physical_matrix_orientation(self):
        from dynsyslab_client.cyclic_dynamics import fit_cyclic_state_pair_models
        x = np.array([[.03, .01], [-.02, .04], [.01, -.02], [-.04, -.01]])
        matrices = [np.array([[1., .5], [0., .8]]), np.array([[.6, 0.], [-.2, 1.]])]
        pairs = [dict(dimension=2, samples=len(x), x=x.ravel().tolist(),
                      y=(x @ a.T).ravel().tolist()) for a in matrices]
        fitted = fit_cyclic_state_pair_models(self.api, pairs, validation_pairs=pairs)
        np.testing.assert_allclose(np.array(fitted['values']).reshape(2, 2, 2), matrices, atol=2e-13)
        self.assertEqual(len(fitted['validation']), 2)
        bad = [dict(dimension=2, samples=2, x=[1, 0, -1, 0], y=[1, 0, -1, 0])]
        with self.assertRaises((ValueError, NativeWorkflowError)):
            fit_cyclic_state_pair_models(self.api, bad)

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
            self.assertEqual(json.loads(info.stdout)["compression"]["identity"], fit["model_identity"])

    def test_independent_weighted_edmd_normal_equations(self):
        rng = np.random.default_rng(14)
        x = rng.normal(size=(90, 2))
        y = np.column_stack((.9*x[:, 0]+.1*x[:, 1]**2, .7*x[:, 1]))
        w = np.linspace(.2, 2, len(x))
        result = self.api.koopman_spectral({"action": "fit", "data": {"dimension": 2,
            "x": x.ravel().tolist(), "y": y.ravel().tolist(), "weights": w.tolist()},
            "settings": {"method": "resdmd", "dictionary": {"kind": "linear", "constant": False}}})
        model = result["compression"]
        basis = np.array(model["basis"]).reshape(2, 2)
        action = np.array(model["action"]).reshape(2, 2)
        operator = basis @ action @ np.linalg.inv(basis)
        direct = np.linalg.solve(x.T @ (w[:, None]*x), x.T @ (w[:, None]*y))
        np.testing.assert_allclose(operator, direct, atol=1e-12)
        # The eigenpairs follow from the operator; each pair's residual is independent of its scale.
        lam, vectors = np.linalg.eig(action)
        c = basis @ vectors
        residual = np.sqrt(np.sum(w[:, None]*abs(y@c-(x@c)*lam)**2, axis=0)/
                           np.sum(w[:, None]*abs(x@c)**2, axis=0))
        np.testing.assert_allclose(sorted(result["spectral"]["residuals"]["eigenpair"]), sorted(residual), atol=1e-12)

    def test_batch_free_rollout_with_independent_rotation_reference(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        model = self.fit('mpedmd', output_archive=str(Path(folder.name) / "rotation.npz"))['archive']
        points = np.array([[1., 0.], [.3, -.7], [0., 0.]])
        for steps in (0, 31):
            result = self.api.koopman_spectral(dict(action='rollout_batch', krom=model,
                initial_points=points.ravel().tolist(), steps=steps))
            theta = .31*np.arange(steps+1)
            expected = np.stack([np.column_stack((x*np.cos(theta)-y*np.sin(theta),
                x*np.sin(theta)+y*np.cos(theta))) for x, y in points])
            np.testing.assert_allclose(np.array(result['values']).reshape(expected.shape), expected, atol=2e-12)
            self.assertEqual(result['layout'], 'trajectory/time/output')
            self.assertEqual(result['trajectories'], 3)
            self.assertEqual(result['samples'], steps+1)
        for points_bad, steps_bad in (([], 1), ([1], 1), ([1, 0], 10000000)):
            with self.assertRaises(NativeWorkflowError):
                self.api.koopman_spectral(dict(action='rollout_batch', krom=model,
                    initial_points=points_bad, steps=steps_bad))

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
            fitted = [p for p in paths if "compression_action" in np.load(p).files]
            self.assertEqual(len(fitted), 1, result)
            model = self.api.koopman_spectral({"action": "inspect", "krom": str(fitted[0])})
            self.assertEqual(model["rank"], 2)
            self.assertLess(max(model["spectral"]["residuals"]["eigenpair"]), 1e-10)
            self.assertAlmostEqual(model["compression"]["sample_interval"], .1)

if __name__ == "__main__":
    unittest.main()
