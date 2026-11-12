"""Read-only planning and compact result summaries for post-processing workflows."""
from __future__ import annotations
from math import comb
from pathlib import Path
from typing import Any
import numpy as np

from .core import Trajectory, load_solution
from .resources import estimate_dense_bytes
from .workflow import WorkflowContext


def plan_trajectory_hdmd_sindy(source: str | Path | Trajectory, hdmd: dict[str,Any], sindy: dict[str,Any]) -> dict[str,Any]:
    """Inspect one source and predict standard pipeline geometry without computing it."""
    trajectory=source if isinstance(source,Trajectory) else load_solution(source)
    prepared=WorkflowContext().prepare(trajectory,hdmd)
    analysis_trajectory=prepared.source; configured=prepared.settings
    samples=int(analysis_trajectory.time.size); state_dimension=int(np.asarray(analysis_trajectory.values).reshape(samples,-1).shape[1])
    stride=int(configured["sample_stride"]); analysis_samples=(samples-1)//stride+1
    columns=int(configured["matrix_size"]); requested_rows=int(configured["sampling_window"])
    realized_rows=min(requested_rows,max(0,analysis_samples-columns))
    full_samples=columns+requested_rows; shape=(realized_rows,columns)
    estimate=estimate_dense_bytes(shape,shape,itemsize=8,workspace_factor=4.) if realized_rows and columns else 0
    rank_limit=min(shape) if all(shape) else 0
    requested_rank=configured.get("rank"); expected_upper=rank_limit if requested_rank is None else min(rank_limit,int(requested_rank))
    degree=int(sindy["polynomial_degree"]); feature_count=comb(state_dimension+degree,degree)
    issues=[]
    if realized_rows<2: issues.append("The source does not contain enough analysis samples for two shifted Hankel rows.")
    if not trajectory.metadata.get("uniform_time_steps",False): issues.append("The source is not declared uniformly sampled; explicit Koopman resampling is required.")
    budget=float(configured["max_memory_mb"])
    if estimate>budget*1024**2: issues.append(f"The H-DMD preflight estimate exceeds max_memory_mb={budget:g}.")
    return {"format":"dynsyslab.postprocessing-plan/v1","valid":not issues,"issues":issues,
            "source":{"archive":None if trajectory.source is None else str(trajectory.source),"samples":samples,
                      "state_dimension":state_dimension,"uniform_time_steps":bool(trajectory.metadata.get("uniform_time_steps",False)),
                      "archive_delta_t":trajectory.metadata.get("delta_t"),"analysis_delta_t":analysis_trajectory.metadata.get("delta_t"),
                      "analysis_coordinates":prepared.coordinates,"resampled":prepared.resampled},
            "hdmd":{"observable":configured["observable"],"matrix_layout":configured["matrix_layout"],
                    "requested_shape":[requested_rows,columns],"realized_shape":list(shape),"full_window_source_samples":full_samples,
                    "uses_full_window":realized_rows==requested_rows,"rank_policy":"full-above-absolute-tolerance" if requested_rank is None else "explicit-rank",
                    "retained_rank_upper_bound":expected_upper,"svd_tolerance":configured["svd_tolerance"],"svd_driver":configured["svd_driver"],
                    "svd_threads":configured.get("svd_threads"),
                    "eigendecomposition_threads":configured.get("eigendecomposition_threads"),
                    "linear_algebra_backend":configured.get("linear_algebra_backend","auto"),
                    "estimated_dense_working_set_mib":estimate/1024**2,"max_memory_mb":budget},
            "sindy":{"polynomial_degree":degree,"feature_count":feature_count,"null_space_backend":sindy["null_space_backend"],
                     "reuse_conjugate_fits":bool(sindy["reuse_conjugate_fits"]),"generator_eigenvalue_domain":"continuous",
                     "display_eigenvalue_domain":"discrete","maximum_input_eigenvalues":expected_upper,"max_memory_mb":float(sindy["max_memory_mb"]),
                     "archive_storage":sindy.get("archive_storage","coefficient-first"),
                     "algorithm_materializes_training_values":bool(sindy.get("compute_training_values",True)),
                     "publication_materializes_training_values":sindy.get("archive_storage","coefficient-first")=="full",
                     "estimated_training_values_mib":samples*expected_upper*np.dtype(complex).itemsize/1024**2,
                     "estimated_full_values_mib":samples*expected_upper*np.dtype(complex).itemsize/1024**2,
                     "estimated_coefficient_first_mib":(samples*state_dimension*8+expected_upper*feature_count*np.dtype(complex).itemsize)/1024**2}}


def summarize_processing_runs(runs) -> dict[str,Any]:
    """Return stable JSON output without serializing numerical arrays."""
    records=[]
    for run in runs:
        methods=[]
        for item in run.results:
            spectrum=getattr(item.result,"spectrum",None); krom=getattr(item.result,"krom",None)
            methods.append({"method":item.method,"archive":None if item.archive is None else str(item.archive),"timing":dict(item.timing),
                            "spectrum_eigenvalues":None if spectrum is None else len(spectrum.eigenvalues),
                            "krom_retained_eigenvalues":None if krom is None else len(krom.eigenvalues)})
        records.append({"source":None if run.source.source is None else str(run.source.source),
                        "workflow_policy":getattr(run,"workflow_policy",None),"methods":methods})
    return {"format":"dynsyslab.postprocessing-run/v1","runs":records}
