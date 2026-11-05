from __future__ import annotations

import inspect
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import call, patch

import numpy as np
from scipy.integrate import solve_ivp

from dynsyslib import load_system
from postprocessing.koopman.dmd import (
    HDMD_DEFAULT_MATRIX_SIZE,
    HDMD_DEFAULT_SAMPLING_WINDOW,
    hankel_dmd,
    hankel_dmd_ensemble,
)
from postprocessing.koopman.pipeline import hankel_dmd_archive, hankel_dmd_processor
from postprocessing.core import Trajectory
from postprocessing.processor import processor
from postprocessing.method_handlers import METHOD_HANDLERS, SUPPORTED_METHODS
from provenance import serializable
from time_series import TimeSeries


class HDMDDefaultTests(unittest.TestCase):
    def test_processor_method_registry_is_complete_and_callable(self):
        expected={"dmd","hdmd","edmd","sindy","pca","pod","bagged_pca",
                  "robust_pca","svd","hankel_pod","havok"}
        self.assertEqual(set(SUPPORTED_METHODS),expected)
        self.assertEqual(set(METHOD_HANDLERS),expected)
        for handler in METHOD_HANDLERS.values():
            self.assertTrue(callable(handler.run))
            self.assertTrue(callable(handler.publish))

    def test_reference_hankel_geometry_is_the_public_default(self):
        self.assertEqual(HDMD_DEFAULT_MATRIX_SIZE,2000)
        self.assertEqual(HDMD_DEFAULT_SAMPLING_WINDOW,20000)
        for function in (hankel_dmd,hankel_dmd_archive,hankel_dmd_processor):
            parameters=inspect.signature(function).parameters
            self.assertEqual(parameters["matrix_size"].default,2000)
            self.assertEqual(parameters["sampling_window"].default,20000)

    def test_sampling_window_limits_hankel_rows(self):
        time=np.arange(20,dtype=float)*.1
        states=np.column_stack((np.cos(time),np.sin(time)))
        system=load_system("rpe")
        result=hankel_dmd(states,time=time,matrix_size=4,sampling_window=5,rank=2,
                          metadata={"uniform_time_steps":True,"system":system.name,"system_parameters":serializable(system)})
        embedding=result.metadata["delay_embedding"]
        self.assertEqual(embedding["matrix_size"],4)
        self.assertEqual(embedding["sampling_window"],5)
        self.assertEqual(embedding["hankel_rows"],5)
        self.assertEqual(embedding["hankel_columns"],8)
        self.assertTrue(embedding["used_full_sampling_window"])

    def test_trajectory_pipeline_encodes_full_scalar_shifted_hankel_policy(self):
        workflow=processor.preset("trajectory-hdmd-sindy",sources=[])
        hdmd,sindy=workflow.methods
        self.assertEqual(hdmd.name,"hdmd")
        self.assertEqual(hdmd.settings["matrix_size"],2000)
        self.assertEqual(hdmd.settings["sampling_window"],20000)
        self.assertEqual(hdmd.settings["observable"],0)
        self.assertEqual(hdmd.settings["matrix_layout"],"window-snapshots")
        self.assertEqual(hdmd.settings["mode_sort"],"terminal-amplitude-energy")
        self.assertIsNone(hdmd.settings["rank"])
        self.assertEqual(hdmd.settings["svd_tolerance"],1e-10)
        self.assertEqual(hdmd.settings["svd_driver"],"gesvd")
        self.assertIsNone(hdmd.settings["svd_threads"])
        self.assertIsNone(hdmd.settings["eigendecomposition_threads"])
        self.assertEqual(hdmd.settings["linear_algebra_backend"],"auto")
        self.assertEqual(hdmd.settings["analysis_coordinates"],"integration")
        self.assertEqual(sindy.settings["spectrum"],"previous")
        self.assertEqual(sindy.settings["polynomial_degree"],4)
        self.assertEqual(sindy.settings["gamma_method"],"system_rhs")
        self.assertEqual(sindy.settings["null_space_backend"],"gram-eigh")
        self.assertEqual(sindy.settings["analysis_coordinates"],"integration")
        self.assertEqual(sindy.settings["archive_storage"],"coefficient-first")
        self.assertFalse(sindy.settings["compute_training_values"])
        self.assertNotIn("archive_storage",hdmd.settings)

    def test_window_snapshot_layout_matches_direct_shifted_hankel_pair(self):
        time=np.arange(12,dtype=float)*.1
        radius=np.sin(time)+.2*np.cos(3*time)
        states=np.column_stack((radius,np.arange(time.size,dtype=float)))
        system=load_system("rpe")
        result=hankel_dmd(states,time=time,matrix_size=4,sampling_window=6,observable=0,
                          matrix_layout="window-snapshots",svd_tolerance=0.,
                          metadata={"uniform_time_steps":True,"system":system.name,"system_parameters":serializable(system)})
        windows=np.lib.stride_tricks.sliding_window_view(radius,4)
        X,Y=windows[:6],windows[1:7]
        U,s,Vh=np.linalg.svd(X,full_matrices=False)
        operator=(U.conj().T@Y@Vh.conj().T)@np.diag(1/s)
        self.assertEqual(result.metadata["delay_embedding"]["hankel_rows"],6)
        self.assertEqual(result.metadata["delay_embedding"]["hankel_columns"],4)
        self.assertEqual(result.metadata["preprocessing"]["observable"],{"kind":"state_index","index":0})
        np.testing.assert_allclose(result.operator,operator,rtol=1e-10,atol=1e-10)
        self.assertEqual(result.metadata["delay_embedding"]["interpolation_state_anchor"],"window-end")
        np.testing.assert_array_equal(
            result.spectrum.eigenfunctions[0].state_callable.samples,
            states[3:9],
        )

    def test_ensemble_window_snapshot_interpolation_uses_each_window_end_state(self):
        time=np.arange(12,dtype=float)*.1
        first=np.column_stack((np.sin(time),np.arange(time.size,dtype=float)))
        second=np.column_stack((10.+np.sin(time),100.+np.arange(time.size,dtype=float)))
        pooled=TimeSeries(time,np.stack((first,second),axis=1),pooled_members=2)
        system=load_system("rpe")
        result=hankel_dmd_ensemble(pooled,matrix_size=4,sampling_window=6,rank=2,
            observable=0,matrix_layout="window-snapshots",
            metadata={"uniform_time_steps":True,"system":system.name,
                      "system_parameters":serializable(system)})
        expected=np.vstack((first[3:9],second[3:9]))
        self.assertEqual(result.metadata["delay_embedding"]["interpolation_state_anchor"],"window-end")
        np.testing.assert_array_equal(
            result.spectrum.eigenfunctions[0].state_callable.samples,expected)

    def test_terminal_amplitude_energy_sort_preserves_full_spectrum(self):
        time=np.arange(16,dtype=float)*.1
        radius=np.sin(time)+.2*np.cos(3*time)+.04*np.sin(7*time)
        states=np.column_stack((radius,np.zeros_like(radius)))
        system=load_system("rpe")
        metadata={"uniform_time_steps":True,"system":system.name,"system_parameters":serializable(system)}
        unsorted=hankel_dmd(states,time=time,matrix_size=5,sampling_window=9,observable=0,
                            matrix_layout="window-snapshots",rank=4,mode_sort="none",metadata=metadata)
        sorted_result=hankel_dmd(states,time=time,matrix_size=5,sampling_window=9,observable=0,
                                 matrix_layout="window-snapshots",rank=4,
                                 mode_sort="terminal-amplitude-energy",metadata=metadata)
        np.testing.assert_allclose(np.sort_complex(sorted_result.eigenvalues),
                                   np.sort_complex(unsorted.eigenvalues))
        self.assertEqual(sorted_result.eigenvalues.size,unsorted.eigenvalues.size)
        self.assertTrue(sorted_result.metadata["mode_sort"]["full_spectrum_retained"])
        order=sorted_result.metadata["mode_sort"]["original_mode_indices"]
        self.assertEqual(sorted(order),list(range(unsorted.eigenvalues.size)))
        energy=sorted_result.metadata["mode_sort"]["sorted_amplitude_magnitudes"]
        self.assertTrue(all(left>=right for left,right in zip(energy,energy[1:])))

    def test_ensemble_hdmd_uses_shared_svd_policy_without_cross_boundary_pairs(self):
        time=np.arange(12,dtype=float)*.1
        first=np.column_stack((np.cos(time),np.sin(time)))
        second=np.column_stack((10.+np.cos(2*time),-5.+np.sin(2*time)))
        pooled=TimeSeries(time,np.stack((first,second),axis=1),pooled_members=2)
        system=load_system("rpe")
        result=hankel_dmd_ensemble(pooled,matrix_size=3,sampling_window=6,rank=2,
            svd_driver="gesvd",metadata={"uniform_time_steps":True,"system":system.name,
                                          "system_parameters":serializable(system)})
        blocks=[]
        for values in (first,second):
            used=values[:8]
            blocks.append(np.hstack([used[index:6+index] for index in range(3)]))
        X=np.hstack([block[:-1].T for block in blocks]); Y=np.hstack([block[1:].T for block in blocks])
        from scipy.linalg import svd
        U,s,Vh=svd(X,full_matrices=False,lapack_driver="gesvd",check_finite=False)
        U,s,Vh=U[:,:2],s[:2],Vh[:2]
        expected=(U.conj().T@Y)@Vh.conj().T/s[None,:]
        # Singular-vector signs are free between LAPACK calls, so the reduced
        # operators can differ by a diagonal similarity transform while their
        # spectra remain identical.
        np.testing.assert_allclose(np.sort_complex(result.eigenvalues),
                                   np.sort_complex(np.linalg.eigvals(expected)),rtol=2e-13,atol=2e-13)
        self.assertEqual(result.metadata["svd_driver"],"gesvd")
        self.assertEqual(result.metadata["delay_embedding"]["paired_snapshots"],10)
        self.assertEqual(result.metadata["delay_embedding"]["trajectory_boundary_transitions"],0)
        self.assertEqual(result.eigenfunction_values.shape[0],12)

    def test_explicit_svd_thread_limit_is_applied_and_recorded(self):
        time=np.arange(12,dtype=float)*.1
        states=np.column_stack((np.cos(time),np.sin(time)))
        system=load_system("rpe")
        with patch("postprocessing.linalg_backends.threadpool_limits") as limit:
            result=hankel_dmd(states,time=time,matrix_size=3,sampling_window=6,rank=2,
                svd_threads=1,eigendecomposition_threads=2,
                metadata={"uniform_time_steps":True,"system":system.name,
                                        "system_parameters":serializable(system)})
        self.assertEqual(limit.call_args_list,[
            call(limits=1,user_api="blas"),
            call(limits=2,user_api="blas"),
        ])
        self.assertEqual(result.metadata["svd_threads"],1)
        self.assertEqual(result.metadata["eigendecomposition_threads"],2)

    def test_processor_reports_algorithm_archive_and_total_timings(self):
        time=np.arange(14,dtype=float)*.1
        states=np.column_stack((np.cos(time),np.sin(time)))
        system=load_system("rpe")
        source=Trajectory(TimeSeries(time,states),{"uniform_time_steps":True,"system":system.name,
            "system_parameters":serializable(system)},None)
        run=processor([source],archive=False).add_method("hdmd",matrix_size=4,sampling_window=6,rank=2).run()[0]
        timing=run.results[0].timing
        self.assertEqual(set(timing),{"algorithm_seconds","archive_seconds","total_seconds"})
        self.assertGreaterEqual(timing["total_seconds"],timing["algorithm_seconds"])
        self.assertGreaterEqual(timing["total_seconds"],timing["archive_seconds"])
        self.assertEqual(run.results[0].result.metadata["processor_timing"]["clock"],"time.perf_counter")

    def test_sindy_archive_records_exact_preceding_spectrum_archive(self):
        time=np.arange(14,dtype=float)*.1
        states=np.column_stack((np.cos(time),np.sin(time)))
        system=load_system("rpe")
        source=Trajectory(TimeSeries(time,states),{"uniform_time_steps":True,"system":system.name,
            "system_parameters":serializable(system)},Path("source.npz"))
        hdmd_spectrum=SimpleNamespace()
        hdmd_result=SimpleNamespace(metadata={},spectrum=hdmd_spectrum)
        sindy_result=SimpleNamespace(metadata={},spectrum=SimpleNamespace())
        archives=[Path("hdmd-derived.npz"),Path("sindy-derived.npz")]
        with patch("postprocessing.method_handlers.koopman.hankel_dmd",return_value=hdmd_result), \
             patch("postprocessing.method_handlers.koopman.koopman_with_sindy",return_value=sindy_result) as sindy, \
             patch.object(processor,"_archive",side_effect=archives):
            run=processor.preset(sources=[source]).run()[0]
        self.assertIs(sindy.call_args.kwargs["spectrum"],hdmd_spectrum)
        self.assertFalse(sindy.call_args.kwargs["compute_training_values"])
        self.assertEqual(run.results[1].result.metadata["eigenvalue_source_archive"],"hdmd-derived.npz")

    def test_full_preset_preserves_dimensionless_krom_coordinate_capabilities(self):
        system=load_system("rpe"); tau=np.linspace(0.,.2,101)
        integrated=solve_ivp(system.dimensionless_rhs,(tau[0],tau[-1]),np.asarray((1.1,0.)),t_eval=tau,rtol=1e-10,atol=1e-12)
        self.assertTrue(integrated.success)
        physical_time=np.asarray(system.to_dimensional_time(tau))
        physical_states=np.asarray([system.to_dimensional_state(row) for row in integrated.y.T])
        source=Trajectory(TimeSeries(physical_time,physical_states),{"uniform_time_steps":True,"delta_t":float(np.diff(physical_time)[0]),
            "integration_form":"dimensionless","nondimensionalization":serializable(system.dimensionless_parameters),
            "system":system.name,"system_parameters":serializable(system)},None)
        run=processor.preset(sources=[source],archive=False,matrix_size=8,sampling_window=40,rank=4,
            sindy_settings={"polynomial_degree":2,"max_iterations":3}).run()[0]
        hdmd,sindy=(item.result for item in run.results)
        self.assertEqual(hdmd.krom.coordinate_system,"dimensionless")
        self.assertEqual(hdmd.krom.state_representation,"delay-observable")
        self.assertFalse(hdmd.krom.capabilities()["state_reconstruction"])
        self.assertTrue(hdmd.krom.capabilities()["exact_control_coupling"])
        self.assertFalse(hdmd.krom.capabilities()["analytic_control_coupling"])
        self.assertTrue(all(function.representation_kind=="scattered-interpolation"
                            for function in hdmd.krom.eigenfunctions))
        self.assertTrue(all(function.gradient_backend=="finite_difference"
                            for function in hdmd.krom.eigenfunctions))
        self.assertIsNotNone(hdmd.krom.observable_reconstruction)
        self.assertEqual(hdmd.krom.C(output="complex").shape[1],hdmd.krom.lifted_dimension)
        self.assertEqual(sindy.krom.coordinate_system,"dimensionless")
        self.assertEqual(sindy.krom.state_representation,"eigenfunction-only")
        self.assertIsNone(sindy.training_eigenfunction_values)
        self.assertTrue(all(function.values is None for function in sindy.spectrum.eigenfunctions))
        self.assertTrue(sindy.krom.capabilities()["analytic_control_coupling"])
        self.assertTrue(all(left is right for left,right in zip(hdmd.spectrum.eigenvalues,sindy.spectrum.eigenvalues)))
        physical_state=physical_states[10]; model_state=integrated.y.T[10]; scale=system.dimensionless_parameters.time_scale_s
        np.testing.assert_allclose(sindy.krom.psi(physical_state,state_coordinates="physical"),sindy.krom.psi(model_state,state_coordinates="model"))
        np.testing.assert_allclose(sindy.krom.A(time_coordinates="physical"),sindy.krom.A(time_coordinates="model")/scale)
        hdmd_input=hdmd.krom.B_psi(model_state,output="complex",state_coordinates="model",
                                   time_coordinates="model")
        self.assertEqual(hdmd_input.shape[0],hdmd.krom.lifted_dimension)
        self.assertTrue(np.all(np.isfinite(hdmd_input)))
        np.testing.assert_allclose(sindy.krom.B_psi(physical_state,output="complex",state_coordinates="physical",time_coordinates="physical"),
                                   sindy.krom.B_psi(model_state,output="complex",state_coordinates="model",time_coordinates="model")/scale)


if __name__=="__main__":
    unittest.main()
