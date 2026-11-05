from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

from dynsyslib import load_system
from postprocessing.linalg_backends import (
    LinearAlgebraBackend,
    available_linear_algebra_backends,
    select_linear_algebra_backend,
    torch_cuda_backend,
)
from postprocessing.koopman.dmd import hankel_dmd
from provenance import serializable


def _nearest_set_distance(left: np.ndarray,right: np.ndarray) -> float:
    return max(float(np.min(np.abs(value-right))) for value in left)


class LinearAlgebraBackendTests(unittest.TestCase):
    def _problem(self):
        time=np.arange(30,dtype=float)*.01
        states=np.column_stack((np.cos(2*time)+.1*np.sin(5*time),np.sin(2*time)))
        system=load_system("rpe")
        metadata={"uniform_time_steps":True,"system":system.name,
                  "system_parameters":serializable(system)}
        return time,states,metadata

    def test_auto_is_the_admitted_scipy_reference(self):
        selected=select_linear_algebra_backend("auto")
        self.assertEqual(selected.name,"scipy-lapack")
        self.assertTrue(selected.auto_admitted)
        self.assertTrue(selected.full_spectrum)
        self.assertIn("gesvd",selected.svd_drivers)

    def test_capability_inventory_keeps_cuda_explicit_only(self):
        backends=available_linear_algebra_backends()
        self.assertEqual(set(backends),{"scipy","torch-cuda"})
        self.assertTrue(backends["scipy"].available)
        self.assertFalse(backends["torch-cuda"].auto_admitted)

    def test_explicit_unavailable_cuda_request_never_falls_back(self):
        unavailable=LinearAlgebraBackend("torch-cuda",False,True,"cuda",True,("gesvd",),False,
                                         reason="test runtime has no CUDA")
        with patch("postprocessing.linalg_backends.torch_cuda_backend",return_value=unavailable):
            with self.assertRaisesRegex(RuntimeError,"no CUDA"):
                select_linear_algebra_backend("torch-cuda")

    def test_invalid_backend_name_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"linear_algebra_backend"):
            select_linear_algebra_backend("fortran")

    def test_scipy_dispatch_preserves_spectrum_and_records_runtime(self):
        time,states,metadata=self._problem()
        auto=hankel_dmd(states,time=time,matrix_size=5,sampling_window=20,observable=0,
                        matrix_layout="window-snapshots",svd_driver="gesvd",svd_tolerance=0.,
                        linear_algebra_backend="auto",metadata=metadata)
        explicit=hankel_dmd(states,time=time,matrix_size=5,sampling_window=20,observable=0,
                            matrix_layout="window-snapshots",svd_driver="gesvd",svd_tolerance=0.,
                            linear_algebra_backend="scipy",metadata=metadata)
        self.assertLess(_nearest_set_distance(auto.eigenvalues,explicit.eigenvalues),1e-11)
        backend=explicit.metadata["linear_algebra_backend"]
        self.assertEqual(backend["requested"],"scipy")
        self.assertEqual(backend["name"],"scipy-lapack")
        self.assertIsInstance(backend["libraries"],tuple)

    def test_torch_cuda_full_kernel_matches_reference_when_available(self):
        if not torch_cuda_backend().available:
            self.skipTest("PyTorch CUDA backend is not available on this machine")
        time,states,metadata=self._problem()
        reference=hankel_dmd(states,time=time,matrix_size=5,sampling_window=20,observable=0,
                             matrix_layout="window-snapshots",svd_driver="gesvd",svd_tolerance=0.,
                             linear_algebra_backend="scipy",metadata=metadata)
        candidate=hankel_dmd(states,time=time,matrix_size=5,sampling_window=20,observable=0,
                             matrix_layout="window-snapshots",svd_driver="gesvd",svd_tolerance=0.,
                             linear_algebra_backend="torch-cuda",metadata=metadata)
        self.assertLess(_nearest_set_distance(reference.eigenvalues,candidate.eigenvalues),1e-8)
        self.assertEqual(candidate.metadata["linear_algebra_backend"]["name"],"torch-cuda")


if __name__=="__main__":
    unittest.main()
