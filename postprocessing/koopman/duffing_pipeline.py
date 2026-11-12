"""Composable Duffing H-DMD-to-SINDy orchestration.

The workflow policy and data generators remain in :mod:`workflows`; this
module owns the comparatively high-level artifact assembly step so each file
stays within DynSysLab's preferred source-size boundary.
"""
from __future__ import annotations


def duffing_hdmd_sindy(sop=None):
    """Run the retained scalar-H-DMD/full-state-SINDy operating procedure."""
    from dynsyslib import Observable
    from postprocessing.core import Trajectory
    from rom import fit_observable_reconstruction, transfer_observable_reconstruction

    from .composable import KoopmanPipelineSpecification, run_composable_koopman_pipeline
    from .dmd import hankel_dmd_ensemble
    from .workflows import DuffingKoopmanSOP, generate_duffing_free_training

    policy = sop or DuffingKoopmanSOP()
    series, system, metadata = generate_duffing_free_training(policy)
    scalar = Observable.component(2, policy.hdmd_scalar_state_index, name="x")
    hdmd = hankel_dmd_ensemble(
        series,
        matrix_size=policy.hdmd_matrix_size,
        sampling_window=policy.hdmd_sampling_window,
        rank=policy.hdmd_rank,
        observable=scalar,
        matrix_layout=policy.hdmd_matrix_layout,
        mode_sort=policy.hdmd_mode_sort,
        svd_tolerance=1e-10,
        metadata=metadata,
        max_memory_mb=2048.0,
    )
    full_state = Observable.identity(2, name="Duffing state (x, xdot)")
    trajectories = [
        Trajectory(
            member,
            {
                **metadata,
                "pooled_input": True,
                "pooled_member_index": index,
                "pooled_member_count": series.ensemble_size,
            },
            None,
        )
        for index, member in enumerate(series.members())
    ]
    pipeline = run_composable_koopman_pipeline(
        KoopmanPipelineSpecification("hdmd", "sindy"),
        eigenvalue_options={"result": hdmd},
        eigenfunction_options={
            "trajectories": trajectories,
            "system": system,
            "coordinate_observable": full_state,
            "polynomial_degree": policy.sindy_polynomial_degree,
            "threshold": policy.sindy_threshold,
            "normalize": policy.sindy_normalize,
            "null_space_backend": policy.sindy_null_space_backend,
            "gamma_method": "system_rhs",
            "compute_training_errors": True,
            "compute_training_values": True,
            "max_memory_mb": 2048.0,
        },
    )
    sindy = pipeline.eigenfunctions.artifact
    common = series.values[: min(20, series.sample_count), :, :].reshape(-1, 2)
    transferred = transfer_observable_reconstruction(hdmd.krom, sindy.krom, common)
    if transferred.validated:
        reconstruction = transferred
    else:
        all_states = series.values.reshape(-1, 2)
        reconstruction = fit_observable_reconstruction(
            sindy.krom, all_states, scalar, ridge=1e-12
        )
        reconstruction.provenance["hdmd_native_C_transfer_diagnostic"] = (
            transferred.summary()
        )
    sindy.krom = sindy.krom.with_observable_reconstruction(reconstruction)
    sindy.spectrum.metadata["hdmd_native_observable_reconstruction"] = (
        hdmd.krom.observable_reconstruction.summary()
    )
    sindy.spectrum.metadata["hdmd_C_transfer_diagnostic"] = transferred.summary()
    sindy.spectrum.metadata["observable_reconstruction"] = reconstruction.summary()
    sindy.metadata["krom"] = sindy.krom.metadata
    return {
        "series": series,
        "system": system,
        "metadata": metadata,
        "hdmd": hdmd,
        "sindy": sindy,
        "pipeline": pipeline,
        "sop": policy,
    }


__all__ = ["duffing_hdmd_sindy"]
