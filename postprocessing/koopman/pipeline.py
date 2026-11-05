"""Archive-to-Koopman convenience calls and solver callback factories."""
from __future__ import annotations
from pathlib import Path
import numpy as np
from ..core import integration_coordinates, load_solution, save_derived
from .dmd import HDMD_DEFAULT_MATRIX_SIZE, HDMD_DEFAULT_SAMPLING_WINDOW, dmd, hankel_dmd
from .edmd import edmd
from .dictionary import ObservableDictionary, PolynomialDictionary
from .resampling import KoopmanResamplingOptions, resample_for_koopman
from .validation import ResamplingDecisionRequired, require_uniform_time_steps
from .koopman_sindy import KoopmanSINDyResult, koopman_with_sindy
from .storage import sindy_archive_payload

def _storage_arrays(source, result):
    """Build canonical arrays and explicit aliases without duplicating payloads."""
    coordinates=None if result.eigenfunction_values is None else np.asarray(result.eigenfunction_values)
    modes=None if result.modes is None else np.asarray(result.modes)
    arrays={"eigenvalues":result.eigenvalues,"operator":result.operator,
            "modes":[] if modes is None else modes,
            "initial_modal_amplitudes":[] if coordinates is None or not coordinates.size else coordinates[0],
            "eigenfunction_state_samples":np.asarray(source.values)}
    modal_amplitudes=getattr(result,"modal_amplitudes",None)
    krom=getattr(result,"krom",None)
    if modal_amplitudes is not None:
        arrays["terminal_modal_amplitudes"]=np.asarray(modal_amplitudes)
    if krom is not None and krom.observable_reconstruction is not None:
        reconstruction=krom.observable_reconstruction
        arrays["observable_reconstruction_coefficients"]=reconstruction.coefficients
        arrays["observable_reconstruction_intercept"]=reconstruction.intercept
    aliases={}
    if coordinates is None: arrays["eigenfunction_values"]=[]
    elif result.eigenfunction_values is result.modes:
        aliases["eigenfunction_values"]="modes"
    else: arrays["eigenfunction_values"]=coordinates
    coefficients=[function.coefficients for function in result.spectrum.eigenfunctions]
    if coefficients and all(item is not None for item in coefficients): arrays["eigenfunction_coefficients"]=np.column_stack(coefficients)
    return arrays,aliases

def _save(source, result, *, save_full_data: bool=False):
    """Persist the complete algorithm spectrum plus its nonredundant KROM view."""
    arrays,aliases=_storage_arrays(source,result)
    coordinates=None if result.eigenfunction_values is None else np.asarray(result.eigenfunction_values)
    metadata={**result.metadata,"spectrum":result.spectrum.summary(),"krom":None if result.krom is None else result.krom.metadata,
              "storage":{"format":"complete-spectrum","full_coordinate_series_saved":coordinates is not None,"operator_saved":True,
                         "terminal_modal_amplitudes_saved":getattr(result,"modal_amplitudes",None) is not None,
                         "observable_reconstruction_saved":getattr(result,"krom",None) is not None and result.krom.observable_reconstruction is not None,
                         "array_aliases":aliases,"deduplicated_array_count":len(aliases)}}
    return save_derived(source,"koopman-"+result.method,analysis_metadata=metadata,**arrays)

def _fixed_lag_source(source, resampling: KoopmanResamplingOptions | None):
    try: require_uniform_time_steps(source.time,source.metadata)
    except ResamplingDecisionRequired:
        if resampling is None: raise
        return resample_for_koopman(source,resampling)
    return source

def dmd_archive(path: str | Path, rank: int|None=None, *, svd_tolerance: float|None=None,
                svd_driver: str="gesdd",svd_threads: int|None=None,
                eigendecomposition_threads: int|None=None,linear_algebra_backend: str="auto",
                resampling: KoopmanResamplingOptions | None=None, save_full_data: bool=False):
    source=_fixed_lag_source(load_solution(path),resampling); return _save(source,dmd(source.values,rank=rank,
        svd_tolerance=svd_tolerance,svd_driver=svd_driver,svd_threads=svd_threads,
        eigendecomposition_threads=eigendecomposition_threads,time=source.time,metadata=source.metadata,
        source_archive=str(source.source),linear_algebra_backend=linear_algebra_backend),save_full_data=save_full_data)

def hankel_dmd_archive(path: str | Path, matrix_size: int=HDMD_DEFAULT_MATRIX_SIZE,
                       rank: int|None=None, *,
                       sampling_window: int=HDMD_DEFAULT_SAMPLING_WINDOW,
                       delay_stride: int=1, sample_stride: int=1, observable=None,
                       observable_series=None, observable_name: str|None=None,
                       matrix_layout: str="delay-state", mode_sort: str="none",
                       svd_tolerance: float|None=None, svd_driver: str="gesdd",
                       svd_threads: int|None=None,
                       eigendecomposition_threads: int|None=None,
                       linear_algebra_backend: str="auto",
                       analysis_coordinates: str="integration",
                       resampling: KoopmanResamplingOptions | None=None,
                       save_full_data: bool=False):
    canonical=load_solution(path)
    if analysis_coordinates=="integration": source=integration_coordinates(canonical)
    elif analysis_coordinates=="archive": source=canonical
    else: raise ValueError("analysis_coordinates must be 'integration' or 'archive'.")
    source=_fixed_lag_source(source,resampling)
    result=hankel_dmd(source.values,matrix_size=matrix_size,rank=rank,
        sampling_window=sampling_window,delay_stride=delay_stride,
        sample_stride=sample_stride,observable=observable,
        observable_series=observable_series,observable_name=observable_name,
        matrix_layout=matrix_layout,mode_sort=mode_sort,
        svd_tolerance=svd_tolerance,svd_driver=svd_driver,
        svd_threads=svd_threads,eigendecomposition_threads=eigendecomposition_threads,
        linear_algebra_backend=linear_algebra_backend,time=source.time,
        coordinate_system="dimensional" if analysis_coordinates=="archive" else None,
        metadata=source.metadata,source_archive=str(canonical.source))
    result.metadata["analysis_coordinates"]=analysis_coordinates
    return _save(source,result,save_full_data=save_full_data)

def edmd_archive(path: str | Path, dictionary: ObservableDictionary|None=None, regularization: float=0., *, resampling: KoopmanResamplingOptions | None=None, save_full_data: bool=False):
    source=_fixed_lag_source(load_solution(path),resampling); return _save(source,edmd(source.values,dictionary or PolynomialDictionary(),regularization=regularization,time=source.time,metadata=source.metadata,source_archive=str(source.source)),save_full_data=save_full_data)

def dmd_processor(rank: int|None=None, **settings):
    def processor(path: str): return dmd_archive(path,rank,**settings)
    processor.__name__="koopman_dmd"; return processor

def hankel_dmd_processor(matrix_size: int=HDMD_DEFAULT_MATRIX_SIZE, rank: int|None=None, *, sampling_window: int=HDMD_DEFAULT_SAMPLING_WINDOW, delay_stride: int=1, sample_stride: int=1, observable=None, observable_series=None, observable_name: str|None=None, matrix_layout: str="delay-state", mode_sort: str="none", svd_tolerance: float|None=None, svd_driver: str="gesdd", svd_threads: int|None=None, eigendecomposition_threads: int|None=None, linear_algebra_backend: str="auto", analysis_coordinates: str="integration", save_full_data: bool=False):
    """Return a solver-compatible callback that archives a Hankel-DMD result."""
    def processor(path: str): return hankel_dmd_archive(path,matrix_size=matrix_size,rank=rank,sampling_window=sampling_window,delay_stride=delay_stride,sample_stride=sample_stride,observable=observable,observable_series=observable_series,observable_name=observable_name,matrix_layout=matrix_layout,mode_sort=mode_sort,svd_tolerance=svd_tolerance,svd_driver=svd_driver,svd_threads=svd_threads,eigendecomposition_threads=eigendecomposition_threads,linear_algebra_backend=linear_algebra_backend,analysis_coordinates=analysis_coordinates,save_full_data=save_full_data)
    processor.__name__="koopman_hankel_dmd"; return processor

def koopman_with_sindy_archive(trajectories, *, archive_storage: str="full", save_full_data: bool=False, **kwargs):
    """Archive wrapper; writes only the derived spectrum, never new trajectories."""
    kwargs.setdefault("compute_training_values",archive_storage=="full")
    result=koopman_with_sindy(trajectories=trajectories,**kwargs)
    primary=result.metadata.get("primary_source")
    if primary is None: raise ValueError("Domain-only results have no trajectory archive; use koopman_with_sindy directly.")
    source=load_solution(primary)
    result.metadata["archive_storage_policy"]=archive_storage
    arrays,metadata=sindy_archive_payload(result,storage_policy=archive_storage)
    return save_derived(source,"koopman-with-sindy",analysis_metadata=metadata,**arrays)
