"""Numerically explicit standard and Hankel DMD implementations."""
from __future__ import annotations
from typing import Any, Literal
import numpy as np
from time import perf_counter
from time_series import TimeSeries
from .result import KoopmanResult
from .validation import require_uniform_time_steps
from .spectrum import discrete_spectrum
from postprocessing.resources import require_memory_budget
from postprocessing.linalg_backends import LinearAlgebraBackendPreference, run_dmd_kernel

# Reference H-DMD geometry retained from the established DynSysLab workflow:
# 2,000 fixed delay samples (columns per observable) and 20,000 sliding
# windows (rows).  Shorter trajectories retain the fixed delay requirement
# and use every available sliding row.
HDMD_DEFAULT_MATRIX_SIZE=2000
HDMD_DEFAULT_SAMPLING_WINDOW=20000
HDMD_DEFAULT_SVD_TOLERANCE=1e-10
HankelDMDLayout=Literal["delay-state","window-snapshots"]
HDMDModeSort=Literal["none","terminal-amplitude-energy"]
SVDDriver=Literal["gesdd","gesvd"]

def _sort_by_terminal_amplitude_energy(result: KoopmanResult, terminal_snapshot: np.ndarray) -> float:
    """Sort projected H-DMD modes by terminal modal-amplitude magnitude."""
    from scipy.linalg import pinv
    started=perf_counter()
    modes=np.asarray(result.modes)
    amplitudes=pinv(modes,check_finite=False)@np.asarray(terminal_snapshot)
    magnitudes=np.abs(amplitudes)
    order=np.argsort(-magnitudes,kind="stable")
    result.eigenvalues=np.asarray(result.eigenvalues)[order]
    result.modes=modes[:,order]
    result.modal_amplitudes=np.asarray(amplitudes)[order]
    elapsed=perf_counter()-started
    result.metadata["mode_sort"]={
        "method":"terminal-amplitude-energy",
        "definition":"abs(pinv(projected_modes) @ terminal_hankel_snapshot)",
        "descending":True,
        "full_spectrum_retained":True,
        "terminal_snapshot_column":-1,
        "original_mode_indices":order.astype(int).tolist(),
        "sorted_amplitude_magnitudes":magnitudes[order].astype(float).tolist(),
    }
    return elapsed

def _dmd_pair(X: np.ndarray, Y: np.ndarray, rank: int|None=None, exact: bool=True, *,
              svd_tolerance: float|None=None, svd_driver: SVDDriver="gesdd", compute_coordinates: bool=True,
              svd_threads: int|None=None,
              eigendecomposition_threads: int|None=None,
              linear_algebra_backend: LinearAlgebraBackendPreference="auto",
              max_memory_mb: float|None=1024.0) -> KoopmanResult:
    """DMD for matrices whose columns are paired snapshots."""
    started=perf_counter(); X=np.asarray(X,dtype=float); Y=np.asarray(Y,dtype=float)
    if X.ndim!=2 or X.shape!=Y.shape or min(X.shape)<2:
        raise ValueError("DMD snapshot matrices must be equal two-dimensional arrays with at least two rows and columns.")
    if svd_tolerance is not None and (not np.isfinite(svd_tolerance) or svd_tolerance<0):
        raise ValueError("svd_tolerance must be a finite nonnegative absolute singular-value threshold.")
    if svd_driver not in {"gesdd","gesvd"}: raise ValueError("svd_driver must be 'gesdd' or 'gesvd'.")
    if svd_threads is not None and (int(svd_threads)!=svd_threads or svd_threads<1):
        raise ValueError("svd_threads must be a positive integer or None.")
    if eigendecomposition_threads is not None and (int(eigendecomposition_threads)!=eigendecomposition_threads or eigendecomposition_threads<1):
        raise ValueError("eigendecomposition_threads must be a positive integer or None.")
    require_memory_budget("DMD",X.shape,Y.shape,max_memory_mb=max_memory_mb)
    prepared=perf_counter()
    kernel=run_dmd_kernel(X,Y,rank=rank,exact=exact,svd_tolerance=svd_tolerance,
        svd_driver=svd_driver,compute_coordinates=compute_coordinates,svd_threads=svd_threads,
        eigendecomposition_threads=eigendecomposition_threads,linear_algebra_backend=linear_algebra_backend)
    completed=perf_counter(); timing=dict(kernel.timing); timing["input_preflight_seconds"]=prepared-started
    timing["total_seconds"]=completed-started
    return KoopmanResult("dmd",kernel.eigenvalues,kernel.modes,kernel.coordinates,kernel.operator,
                         {"rank":kernel.retained_rank,"exact":exact,"svd_tolerance":svd_tolerance,"svd_driver":svd_driver,
                          "svd_threads":None if svd_threads is None else int(svd_threads),
                          "eigendecomposition_threads":None if eigendecomposition_threads is None else int(eigendecomposition_threads),
                          "linear_algebra_backend":kernel.backend,
                          "svd_rank_policy":"explicit-rank" if rank is not None else ("absolute-tolerance" if svd_tolerance is not None else "numerical-rank"),
                          "kernel_timing":timing})

def _dmd(states: np.ndarray, rank: int|None=None, exact: bool=True, *,
         svd_tolerance: float|None=None, svd_driver: SVDDriver="gesdd",
         svd_threads: int|None=None,
         eigendecomposition_threads: int|None=None,
         linear_algebra_backend: LinearAlgebraBackendPreference="auto",
         max_memory_mb: float|None=1024.0) -> KoopmanResult:
    x=np.asarray(states,dtype=float)
    if x.ndim==1: x=x[:,None]
    if x.ndim!=2 or x.shape[0]<3: raise ValueError("DMD requires at least three time-ordered state samples.")
    X,Y=x[:-1].T,x[1:].T
    result=_dmd_pair(X,Y,rank=rank,exact=exact,svd_tolerance=svd_tolerance,
                     svd_driver=svd_driver,svd_threads=svd_threads,
                     eigendecomposition_threads=eigendecomposition_threads,
                     linear_algebra_backend=linear_algebra_backend,max_memory_mb=max_memory_mb)
    result.eigenfunction_values=(np.linalg.pinv(result.modes)@x.T).T
    return result

def _series_arrays(states: np.ndarray | TimeSeries, time: np.ndarray | None):
    if isinstance(states,TimeSeries):
        if time is not None: raise ValueError("Pass either a TimeSeries or separate states and time, not both.")
        if states.is_pooled: raise ValueError("This DMD call received a pooled TimeSeries. Use series.members() to analyze each trajectory independently, preserving trajectory boundaries.")
        return states.values,states.time
    if time is None: raise ValueError("Array DMD input requires an explicit time array.")
    return states,time

def _hankel_settings(matrix_size: int, delay_stride: int, sampling_window: int) -> tuple[int,int,int]:
    """Validate fixed-delay columns, stride, and requested sliding rows."""
    if int(matrix_size)!=matrix_size or matrix_size<2: raise ValueError("matrix_size must be an integer of at least 2.")
    if int(delay_stride)!=delay_stride or delay_stride<1: raise ValueError("delay_stride must be a positive integer number of analysis samples.")
    if int(sampling_window)!=sampling_window or sampling_window<3: raise ValueError("sampling_window must be an integer of at least 3 Hankel rows.")
    return int(matrix_size),int(delay_stride),int(sampling_window)

def _subsample_hankel_source(states: np.ndarray, time: np.ndarray, *, sample_stride: int) -> tuple[np.ndarray,np.ndarray,dict]:
    """Return a lower-rate grid-aligned trajectory without changing its source archive."""
    x=np.asarray(states,dtype=float); t=np.asarray(time,dtype=float).reshape(-1)
    native_dt=float(t[1]-t[0]); native_rate=1./native_dt
    if int(sample_stride)!=sample_stride or sample_stride<1: raise ValueError("sample_stride must be a positive integer.")
    stride=int(sample_stride); return x[::stride],t[::stride],{"sample_stride":stride,"native_sample_rate":native_rate,"analysis_sample_rate":native_rate/stride,"interpolated":False}

def _observable_coordinate_system(metadata: dict|None, explicit: str|None) -> str|None:
    if explicit is not None: return explicit
    details=(metadata or {}).get("postprocessing_preprocessing") or {}
    return details.get("coordinate_system")

def _observable_system(observable, metadata: dict|None, supplied: Any|None):
    if supplied is not None or not isinstance(observable,str): return supplied
    if not metadata: return None
    from dynsyslib import resolve_archived_system
    try: return resolve_archived_system(metadata.get("system"),metadata.get("system_parameters"))
    except Exception as error:
        raise ValueError("Named H-DMD observable could not resolve its archived dynamical system.") from error

def _coerce_scalar_observable(observable, state_dimension: int, *, metadata: dict|None,
                              system: Any|None, coordinate_system: str|None):
    from dynsyslib.observables import coerce_observable
    selected_coordinates=_observable_coordinate_system(metadata,coordinate_system)
    if isinstance(observable,str) and selected_coordinates is None:
        raise ValueError(
            "Named H-DMD observables require an explicit coordinate_system or "
            "prepared analysis-coordinate metadata."
        )
    item=coerce_observable(observable,state_dimension,
        system=_observable_system(observable,metadata,system),
        coordinate_system=selected_coordinates)
    if item.output_dimension!=1: raise ValueError("Window-snapshot H-DMD requires one scalar Observable.")
    if item.value_type!="real": raise ValueError("H-DMD currently requires a real-valued scalar Observable.")
    return item

def _sampled_scalar_values(source, time: np.ndarray, *, sample_stride: int,
                           ensemble_members: int=0) -> np.ndarray:
    """Validate a precomputed scalar signal on the exact trajectory grid."""
    if isinstance(source,TimeSeries):
        if not np.array_equal(np.asarray(source.time),np.asarray(time)):
            raise ValueError("Observable TimeSeries must use the exact H-DMD source time grid.")
        if ensemble_members and (not source.is_pooled or source.ensemble_size!=ensemble_members):
            raise ValueError("Pooled observable TimeSeries must match the state ensemble size.")
        if not ensemble_members and source.is_pooled:
            raise ValueError("An individual H-DMD trajectory cannot use a pooled observable TimeSeries.")
        values=np.asarray(source.values)
    else: values=np.asarray(source)
    if np.iscomplexobj(values): raise ValueError("H-DMD currently requires a real-valued sampled observable.")
    if ensemble_members:
        if values.ndim==3 and values.shape[2]==1: values=values[:,:,0]
        if values.shape!=(np.asarray(time).size,ensemble_members):
            raise ValueError("Sampled ensemble observable must have shape (samples, members) or (samples, members, 1).")
    else:
        if values.ndim==2 and values.shape[1]==1: values=values[:,0]
        if values.shape!=(np.asarray(time).size,):
            raise ValueError("Sampled scalar observable must have shape (samples,) or (samples, 1).")
    values=np.asarray(values,dtype=float)
    if np.any(~np.isfinite(values)): raise ValueError("Sampled scalar observable contains non-finite values.")
    return values[::int(sample_stride)]

def _sampled_observable_summary(name: str|None, state_dimension: int, *, pooled: bool) -> dict:
    label=str(name or "provided scalar time series")
    return {"encoding":"sampled-scalar-observable/v1","name":label,
        "state_dimension":int(state_dimension),"output_dimension":1,
        "analytic_jacobian":False,"value_type":"real","descriptor":{
            "kind":"sampled-time-series","label":label,"pooled":bool(pooled),
            "state_callable":False,"time_grid_verified":True}}

def dmd(states: np.ndarray | TimeSeries, rank: int|None=None, exact: bool=True, *,
        svd_tolerance: float|None=None, svd_driver: SVDDriver="gesdd",
        svd_threads: int|None=None, eigendecomposition_threads: int|None=None,
        linear_algebra_backend: LinearAlgebraBackendPreference="auto", time: np.ndarray|None=None,
        metadata: dict|None=None, source_archive: str|None=None,
        max_memory_mb: float|None=1024.0) -> KoopmanResult:
    """Standard DMD on a verified fixed-lag time series."""
    states,time=_series_arrays(states,time); validation=require_uniform_time_steps(time,metadata)
    if np.asarray(states).shape[0]!=np.asarray(time).size: raise ValueError("DMD states and time arrays must have equal sample counts.")
    result=_dmd(states,rank=rank,exact=exact,svd_tolerance=svd_tolerance,svd_driver=svd_driver,
                svd_threads=svd_threads,eigendecomposition_threads=eigendecomposition_threads,
                linear_algebra_backend=linear_algebra_backend,max_memory_mb=max_memory_mb); result.metadata.update(validation)
    result.metadata["preprocessing"]={"state_dimension":int(np.asarray(states).reshape(len(time),-1).shape[1]),"sample_count":int(len(time)),"analysis_delta_t":validation.get("delta_t"),"source_metadata_flag":None if metadata is None else metadata.get("uniform_time_steps")}
    result.metadata["delay_embedding"]={"type":"none","delays":1,"window_samples":1,"window_duration":0.}
    result.spectrum=discrete_spectrum(result.eigenvalues,result.eigenfunction_values,result.modes,time=time,source_metadata=metadata,source_archive=source_archive,
                                      method="dmd",representation="DMD trajectory eigenfunction interpolation",method_metadata=result.metadata,
                                      interpolation_states=np.asarray(states,dtype=float))
    from rom import KROM
    result.krom=KROM(result.spectrum)
    result.metadata["krom"]=result.krom.metadata
    result.metadata["spectrum"]=result.spectrum.summary()
    return result

def hankel_dmd(states: np.ndarray | TimeSeries, matrix_size: int=HDMD_DEFAULT_MATRIX_SIZE, rank: int|None=None, *, sampling_window: int=HDMD_DEFAULT_SAMPLING_WINDOW, delay_stride: int=1, sample_stride: int=1, observable=None, observable_series: np.ndarray|TimeSeries|None=None, observable_name: str|None=None, system: Any|None=None, coordinate_system: str|None=None, matrix_layout: HankelDMDLayout="delay-state", mode_sort: HDMDModeSort="none", svd_tolerance: float|None=None, svd_driver: SVDDriver="gesdd", svd_threads: int|None=None, eigendecomposition_threads: int|None=None, linear_algebra_backend: LinearAlgebraBackendPreference="auto", time: np.ndarray|None=None, metadata: dict|None=None, source_archive: str|None=None, max_memory_mb: float|None=1024.0) -> KoopmanResult:
    """Fit H-DMD from all states or one explicit scalar observable.

    ``observable`` accepts an :class:`Observable`, state index, or system
    catalog name. ``observable_series`` accepts an already evaluated scalar
    signal on the exact source time grid. The two forms are mutually exclusive.
    """
    states,time=_series_arrays(states,time); x=np.asarray(states,dtype=float)
    if x.ndim==1: x=x[:,None]
    validation=require_uniform_time_steps(time,metadata)
    if x.shape[0]!=np.asarray(time).size: raise ValueError("Hankel DMD states and time arrays must have equal sample counts.")
    source_time=np.asarray(time,dtype=float)
    if observable is not None and observable_series is not None:
        raise ValueError("Pass either observable or observable_series, not both.")
    sampled_signal=None if observable_series is None else _sampled_scalar_values(
        observable_series,source_time,sample_stride=sample_stride)
    x,time,sampling=_subsample_hankel_source(x,time,sample_stride=sample_stride)
    source_sample_count=int(x.shape[0]); source_state_dimension=int(x.shape[1]); dt=float(time[1]-time[0]); delays,delay_stride,requested_rows=_hankel_settings(matrix_size,delay_stride,sampling_window); last=(delays-1)*delay_stride
    if matrix_layout not in {"delay-state","window-snapshots"}: raise ValueError("matrix_layout must be 'delay-state' or 'window-snapshots'.")
    if mode_sort not in {"none","terminal-amplitude-energy"}: raise ValueError("mode_sort must be 'none' or 'terminal-amplitude-energy'.")
    if mode_sort!="none" and matrix_layout!="window-snapshots": raise ValueError("H-DMD terminal-amplitude energy sorting requires matrix_layout='window-snapshots'.")
    observable_object=None; observable_contract=None
    if matrix_layout=="window-snapshots":
        embedding_started=perf_counter()
        if observable is None and sampled_signal is None:
            raise ValueError("matrix_layout='window-snapshots' requires one scalar Observable, name, state index, or sampled signal.")
        if delay_stride!=1: raise ValueError("matrix_layout='window-snapshots' currently requires delay_stride=1.")
        if sampled_signal is None:
            observable_object=_coerce_scalar_observable(observable,x.shape[1],metadata=metadata,
                system=system,coordinate_system=coordinate_system)
            signal=observable_object.evaluate(x)[:,0]
            observable_contract=observable_object.summary()
        else:
            signal=sampled_signal
            observable_contract=_sampled_observable_summary(observable_name,x.shape[1],pooled=False)
        embedded_samples=min(requested_rows,signal.size-delays)
        if embedded_samples<2: raise ValueError("H-DMD matrix size and sampling window leave fewer than two shifted window rows.")
        used_samples=embedded_samples+delays
        signal=signal[:used_samples]; analysis_states=x[:used_samples]; analysis_time=np.asarray(time)[:used_samples]
        require_memory_budget("Hankel-DMD",(embedded_samples,delays),(embedded_samples,delays),max_memory_mb=max_memory_mb,workspace_factor=4.0)
        windows=np.lib.stride_tricks.sliding_window_view(signal,delays)
        X,Y=np.asarray(windows[:embedded_samples]),np.asarray(windows[1:embedded_samples+1])
        embedding_seconds=perf_counter()-embedding_started
        result=_dmd_pair(X,Y,rank=rank,exact=False,svd_tolerance=svd_tolerance,svd_driver=svd_driver,compute_coordinates=False,svd_threads=svd_threads,eigendecomposition_threads=eigendecomposition_threads,linear_algebra_backend=linear_algebra_backend,max_memory_mb=max_memory_mb)
        energy_sort_seconds=0.
        if mode_sort=="terminal-amplitude-energy": energy_sort_seconds=_sort_by_terminal_amplitude_energy(result,X[:,-1])
        else:
            from scipy.linalg import pinv
            amplitude_started=perf_counter()
            result.modal_amplitudes=pinv(np.asarray(result.modes),check_finite=False)@X[:,-1]
            energy_sort_seconds=perf_counter()-amplitude_started
            result.metadata["mode_sort"]={"method":"none","full_spectrum_retained":True,
                "terminal_modal_amplitudes_computed":True,
                "definition":"pinv(projected_modes) @ terminal_hankel_snapshot"}
        result.method="hankel-dmd"; result.eigenfunction_values=result.modes
        result.metadata["delays"]=delays; result.metadata.update({**validation,"delta_t":dt})
        legacy_observable={"kind":"state_index","index":int(observable)} if isinstance(observable,(int,np.integer)) and not isinstance(observable,(bool,np.bool_)) else observable_contract
        result.metadata["preprocessing"]={"state_dimension":source_state_dimension,"source_sample_count":source_sample_count,"analysis_sample_count":int(used_samples),"analysis_delta_t":dt,"source_metadata_flag":None if metadata is None else metadata.get("uniform_time_steps"),"subsampling":sampling,"observable":legacy_observable,"observable_contract":observable_contract}
        if metadata is not None: result.metadata["preprocessing"].update(metadata.get("postprocessing_preprocessing") or {})
        result.metadata["delay_embedding"]={"type":"hankel","matrix_layout":matrix_layout,"matrix_size":int(delays),"sampling_window":requested_rows,"time_delay":dt,"delay_stride_samples":1,"embedded_state_dimension":int(X.shape[0]),"embedded_samples":int(X.shape[1]),"hankel_rows":int(X.shape[0]),"hankel_columns":int(X.shape[1]),"window_samples":int(delays),"window_duration":float((delays-1)*dt),"sampling_window_duration":float((embedded_samples-1)*dt),"used_full_sampling_window":bool(embedded_samples==requested_rows),"shifted_pair":True,"interpolation_state_anchor":"window-end"}
        result.metadata["kernel_timing"]["energy_sort_seconds"]=energy_sort_seconds
        result.metadata["kernel_timing"]["total_seconds"]+=energy_sort_seconds
        result.metadata["pipeline_timing"]={"matrix_construction_seconds":embedding_seconds,**result.metadata["kernel_timing"]}
        interpolation_slice=slice(delays-1,delays-1+embedded_samples)
        result.spectrum=discrete_spectrum(result.eigenvalues,result.eigenfunction_values,result.modes,time=analysis_time[interpolation_slice],source_metadata=metadata,source_archive=source_archive,
                                          method="hankel-dmd",representation="Hankel-DMD trajectory eigenfunction interpolation",method_metadata=result.metadata,
                                          interpolation_states=analysis_states[interpolation_slice])
        from rom import KROM
        result.krom=KROM(result.spectrum)
        if result.modal_amplitudes is not None and observable_object is not None:
            from rom.observation import ObservableReconstruction
            coefficients=np.asarray(result.modal_amplitudes,dtype=complex)[list(result.krom._retained_indices)][None,:]
            result.krom.observable_reconstruction=ObservableReconstruction(
                observable_object,coefficients,np.zeros(1),
                {"passed":False,"reason":"native H-DMD terminal-snapshot modes are not a held-out fit"},
                {"method":"H-DMD terminal modal amplitudes","definition":"pinv(projected_modes) @ terminal_hankel_snapshot","mode_sort":"terminal-amplitude-energy"})
            result.krom._validate_observable_reconstruction()
            result.metadata["observable_reconstruction"]=result.krom.observable_reconstruction.summary()
        result.metadata["krom"]=result.krom.metadata; result.metadata["spectrum"]=result.spectrum.summary()
        return result
    if observable is not None or observable_series is not None:
        raise ValueError("Scalar observable selection requires matrix_layout='window-snapshots'.")
    if x.shape[0]<=last+2: raise ValueError("H-DMD matrix size and time delay leave fewer than three embedded samples.")
    available_rows=x.shape[0]-last; embedded_samples=min(requested_rows,available_rows); used_samples=embedded_samples+last
    x=x[:used_samples]; time=np.asarray(time)[:used_samples]
    require_memory_budget("Hankel-DMD",(embedded_samples,x.shape[1]*delays),max_memory_mb=max_memory_mb,workspace_factor=4.0)
    hankel=np.hstack([x[index*delay_stride:x.shape[0]-last+index*delay_stride] for index in range(delays)])
    result=_dmd(hankel,rank=rank,exact=True,svd_tolerance=svd_tolerance,
                svd_driver=svd_driver,svd_threads=svd_threads,eigendecomposition_threads=eigendecomposition_threads,linear_algebra_backend=linear_algebra_backend,max_memory_mb=max_memory_mb); result.method="hankel-dmd"; result.metadata["delays"]=delays; result.metadata.update({**validation,"delta_t":dt})
    result.metadata["preprocessing"]={"state_dimension":source_state_dimension,"source_sample_count":source_sample_count,"analysis_sample_count":int(x.shape[0]),"analysis_delta_t":dt,"source_metadata_flag":None if metadata is None else metadata.get("uniform_time_steps"),"subsampling":sampling,"observable":{"kind":"all_states"}}
    result.metadata["delay_embedding"]={"type":"hankel","matrix_layout":matrix_layout,"matrix_size":int(delays),"sampling_window":requested_rows,"time_delay":float(delay_stride*dt),"delay_stride_samples":int(delay_stride),"embedded_state_dimension":int(hankel.shape[1]),"embedded_samples":int(hankel.shape[0]),"hankel_rows":int(hankel.shape[0]),"hankel_columns":int(hankel.shape[1]),"window_samples":int(last+1),"window_duration":float(last*dt),"sampling_window_duration":float((embedded_samples-1)*dt),"used_full_sampling_window":bool(embedded_samples==requested_rows),"shifted_pair":False}
    result.spectrum=discrete_spectrum(result.eigenvalues,result.eigenfunction_values,result.modes,time=time,source_metadata=metadata,source_archive=source_archive,
                                      method="hankel-dmd",representation="Hankel-DMD trajectory eigenfunction interpolation",method_metadata=result.metadata,
                                      interpolation_states=x[:result.eigenfunction_values.shape[0]])
    from rom import KROM
    result.krom=KROM(result.spectrum)
    result.metadata["krom"]=result.krom.metadata
    result.metadata["spectrum"]=result.spectrum.summary()
    return result

def hankel_dmd_ensemble(series: TimeSeries, matrix_size: int=HDMD_DEFAULT_MATRIX_SIZE, rank: int|None=None, *, sampling_window: int=HDMD_DEFAULT_SAMPLING_WINDOW, delay_stride: int=1, sample_stride: int=1, observable=None, observable_series: np.ndarray|TimeSeries|None=None, observable_name: str|None=None, system: Any|None=None, coordinate_system: str|None=None, matrix_layout: HankelDMDLayout="delay-state", mode_sort: HDMDModeSort="none", svd_tolerance: float|None=None, svd_driver: SVDDriver="gesdd", svd_threads: int|None=None, eigendecomposition_threads: int|None=None, linear_algebra_backend: LinearAlgebraBackendPreference="auto", metadata: dict|None=None, source_archive: str|None=None, max_memory_mb: float|None=1024.0) -> KoopmanResult:
    """Hankel-DMD for pooled trajectories with no cross-boundary pairs."""
    if not series.is_pooled:
        return hankel_dmd(series,matrix_size=matrix_size,rank=rank,
            sampling_window=sampling_window,delay_stride=delay_stride,
            sample_stride=sample_stride,observable=observable,observable_series=observable_series,
            observable_name=observable_name,system=system,coordinate_system=coordinate_system,matrix_layout=matrix_layout,
            mode_sort=mode_sort,svd_tolerance=svd_tolerance,
            svd_driver=svd_driver,svd_threads=svd_threads,
            eigendecomposition_threads=eigendecomposition_threads,linear_algebra_backend=linear_algebra_backend,metadata=metadata,
            source_archive=source_archive,max_memory_mb=max_memory_mb)
    validation=require_uniform_time_steps(series.time,metadata)
    delays,delay_stride,requested_rows=_hankel_settings(matrix_size,delay_stride,sampling_window)
    if int(sample_stride)!=sample_stride or sample_stride<1:
        raise ValueError("sample_stride must be a positive integer.")
    stride=int(sample_stride); analysis_time=np.asarray(series.time)[::stride]
    if analysis_time.size<2:
        raise ValueError("sample_stride leaves fewer than two ensemble time samples.")
    if observable is not None and observable_series is not None:
        raise ValueError("Pass either observable or observable_series, not both.")
    dt=float(analysis_time[1]-analysis_time[0]); last=(delays-1)*delay_stride
    members=series.members(); blocks=[]; interpolation=[]
    if matrix_layout=="window-snapshots":
        if delay_stride!=1: raise ValueError("Ensemble window-snapshots currently requires delay_stride=1.")
        if observable is None and observable_series is None:
            raise ValueError("Ensemble window-snapshots requires one scalar Observable, name, state index, or sampled signal.")
        state_dimension=int(members[0].values.shape[1])
        sampled_signals=None if observable_series is None else _sampled_scalar_values(
            observable_series,series.time,sample_stride=stride,ensemble_members=series.ensemble_size)
        observable_object=None if sampled_signals is not None else _coerce_scalar_observable(
            observable,state_dimension,metadata=metadata,system=system,coordinate_system=coordinate_system)
        observable_contract=(observable_object.summary() if observable_object is not None else
            _sampled_observable_summary(observable_name,state_dimension,pooled=True))
        embedding_started=perf_counter(); x_blocks=[]; y_blocks=[]
        for member_index,item in enumerate(members):
            state=np.asarray(item.values,dtype=float)[::stride]
            signal=(observable_object.evaluate(state)[:,0] if observable_object is not None
                    else sampled_signals[:,member_index])
            embedded=min(requested_rows,signal.size-delays)
            if embedded<2: raise ValueError("Every ensemble member must leave at least two shifted H-DMD window rows.")
            used=embedded+delays; windows=np.lib.stride_tricks.sliding_window_view(signal[:used],delays)
            x_blocks.append(np.asarray(windows[:embedded])); y_blocks.append(np.asarray(windows[1:embedded+1]))
            # Each modal row belongs to the physical state at the right edge
            # of its delay window, because the DMD snapshots advance along
            # Hankel columns and the terminal column is used for amplitudes.
            interpolation.append(state[delays-1:delays-1+embedded])
        X=np.vstack(x_blocks); Y=np.vstack(y_blocks)
        require_memory_budget("ensemble window-snapshot H-DMD",X.shape,Y.shape,max_memory_mb=max_memory_mb,workspace_factor=5.0)
        embedding_seconds=perf_counter()-embedding_started
        result=_dmd_pair(X,Y,rank=rank,exact=False,svd_tolerance=svd_tolerance,svd_driver=svd_driver,
            compute_coordinates=False,svd_threads=svd_threads,eigendecomposition_threads=eigendecomposition_threads,
            linear_algebra_backend=linear_algebra_backend,max_memory_mb=max_memory_mb)
        if mode_sort=="terminal-amplitude-energy": sort_seconds=_sort_by_terminal_amplitude_energy(result,X[:,-1])
        elif mode_sort=="none":
            from scipy.linalg import pinv
            sort_started=perf_counter(); result.modal_amplitudes=pinv(np.asarray(result.modes),check_finite=False)@X[:,-1]
            sort_seconds=perf_counter()-sort_started
            result.metadata["mode_sort"]={"method":"none","full_spectrum_retained":True,
                "terminal_modal_amplitudes_computed":True,"definition":"pinv(projected_modes) @ terminal_hankel_snapshot"}
        else: raise ValueError("mode_sort must be 'none' or 'terminal-amplitude-energy'.")
        result.method="hankel-dmd-ensemble"; result.eigenfunction_values=result.modes
        result.metadata.update({"delays":delays,"ensemble_members":series.ensemble_size,**validation,"delta_t":dt})
        result.metadata["preprocessing"]={"state_dimension":int(members[0].values.shape[1]),
            "source_sample_count_per_trajectory":int(members[0].sample_count),
            "analysis_sample_count_per_trajectory":int(x_blocks[0].shape[0]+delays),
            "trajectory_count":series.ensemble_size,"analysis_delta_t":dt,
            "trajectory_boundaries_preserved":True,"subsampling":{"sample_stride":stride,
            "native_sample_rate":1./float(series.time[1]-series.time[0]),"analysis_sample_rate":1./dt,
            "interpolated":False},"observable":observable_contract,
            "observable_contract":observable_contract}
        result.metadata["delay_embedding"]={"type":"hankel","matrix_layout":"window-snapshots",
            "matrix_size":delays,"sampling_window":requested_rows,"hankel_rows":int(X.shape[0]),
            "hankel_columns":int(X.shape[1]),"rows_per_trajectory":int(x_blocks[0].shape[0]),
            "window_samples":delays,"window_duration":float((delays-1)*dt),
            "trajectory_boundary_transitions":0,"shifted_pair":True,
            "interpolation_state_anchor":"window-end"}
        result.metadata["kernel_timing"]["energy_sort_seconds"]=sort_seconds
        result.metadata["kernel_timing"]["total_seconds"]+=sort_seconds
        result.metadata["pipeline_timing"]={"matrix_and_pair_construction_seconds":embedding_seconds,**result.metadata["kernel_timing"]}
        states=np.vstack(interpolation)
        result.spectrum=discrete_spectrum(result.eigenvalues,result.eigenfunction_values,result.modes,time=analysis_time,
            source_metadata=metadata,source_archive=source_archive,method="hankel-dmd-ensemble",
            representation="ensemble scalar-observable window H-DMD interpolation",method_metadata=result.metadata,
            interpolation_states=states)
        from rom import KROM
        result.krom=KROM(result.spectrum)
        if result.modal_amplitudes is not None and observable_object is not None:
            from rom.observation import ObservableReconstruction
            coefficients=np.asarray(result.modal_amplitudes,dtype=complex)[list(result.krom._retained_indices)][None,:]
            result.krom.observable_reconstruction=ObservableReconstruction(observable_object,coefficients,np.zeros(1),
                {"passed":False,"reason":"native H-DMD terminal-snapshot modes are not a held-out fit"},
                {"method":"H-DMD terminal modal amplitudes","definition":"pinv(projected_modes) @ terminal_hankel_snapshot","mode_sort":mode_sort})
            result.krom._validate_observable_reconstruction()
            result.metadata["observable_reconstruction"]=result.krom.observable_reconstruction.summary()
        result.metadata["krom"]=result.krom.metadata; result.metadata["spectrum"]=result.spectrum.summary()
        return result
    if observable is not None or observable_series is not None:
        raise ValueError("Scalar observable selection requires matrix_layout='window-snapshots'.")
    if matrix_layout!="delay-state": raise ValueError("matrix_layout must be 'delay-state' or 'window-snapshots'.")
    embedding_started=perf_counter()
    for item in members:
        x=np.asarray(item.values,dtype=float)[::stride]
        if x.ndim!=2 or x.shape[0]<=last+2:
            raise ValueError("Every pooled member must leave at least three embedded samples after delay and sample strides.")
        rows=min(requested_rows,x.shape[0]-last); used=rows+last; x=x[:used]
        block=np.hstack([x[index*delay_stride:x.shape[0]-last+index*delay_stride] for index in range(delays)])
        blocks.append(block); interpolation.append(x[:block.shape[0]])
    X=np.hstack([block[:-1].T for block in blocks]); Y=np.hstack([block[1:].T for block in blocks])
    require_memory_budget("ensemble Hankel-DMD",X.shape,Y.shape,max_memory_mb=max_memory_mb,workspace_factor=5.0)
    embedding_seconds=perf_counter()-embedding_started
    result=_dmd_pair(X,Y,rank=rank,exact=True,svd_tolerance=svd_tolerance,
                     svd_driver=svd_driver,compute_coordinates=False,svd_threads=svd_threads,
                     eigendecomposition_threads=eigendecomposition_threads,linear_algebra_backend=linear_algebra_backend,max_memory_mb=max_memory_mb)
    coordinate_started=perf_counter()
    all_embedded=np.hstack([block.T for block in blocks])
    result.eigenfunction_values=(np.linalg.pinv(result.modes)@all_embedded).T
    coordinate_seconds=perf_counter()-coordinate_started
    result.method="hankel-dmd-ensemble"
    result.metadata.update({"delays":delays,"ensemble_members":series.ensemble_size,**validation,"delta_t":dt})
    result.metadata["kernel_timing"]["coordinate_recovery_seconds"]=coordinate_seconds
    result.metadata["kernel_timing"]["total_seconds"]+=coordinate_seconds
    sampling={"sample_stride":stride,"native_sample_rate":1./float(series.time[1]-series.time[0]),
              "analysis_sample_rate":1./dt,"interpolated":False}
    result.metadata["preprocessing"]={"state_dimension":int(members[0].values.shape[1]),
        "source_sample_count_per_trajectory":int(members[0].sample_count),
        "analysis_sample_count_per_trajectory":int(blocks[0].shape[0]+last),
        "trajectory_count":series.ensemble_size,"analysis_delta_t":dt,
        "trajectory_boundaries_preserved":True,"subsampling":sampling}
    paired_snapshots=int(sum(block.shape[0]-1 for block in blocks))
    result.metadata["delay_embedding"]={"type":"hankel","matrix_layout":"delay-state",
        "matrix_size":int(delays),"sampling_window":requested_rows,
        "delay_stride_samples":int(delay_stride),"time_delay":float(delay_stride*dt),
        "embedded_state_dimension":int(blocks[0].shape[1]),
        "embedded_samples_per_trajectory":int(blocks[0].shape[0]),
        "hankel_rows_per_trajectory":int(blocks[0].shape[0]),
        "hankel_columns":int(blocks[0].shape[1]),"window_samples":int(last+1),
        "window_duration":float(last*dt),"used_full_sampling_window":bool(blocks[0].shape[0]==requested_rows),
        "shifted_pair":True,"paired_snapshots":paired_snapshots,
            "trajectory_boundary_transitions":0,"interpolation_state_anchor":"window-start"}
    result.metadata["pipeline_timing"]={"matrix_and_pair_construction_seconds":embedding_seconds,
                                         **result.metadata["kernel_timing"]}
    interpolation_states=np.vstack(interpolation)
    result.spectrum=discrete_spectrum(result.eigenvalues,result.eigenfunction_values,result.modes,time=analysis_time,source_metadata=metadata,source_archive=source_archive,method="hankel-dmd-ensemble",representation="ensemble Hankel-DMD trajectory eigenfunction interpolation",method_metadata=result.metadata,interpolation_states=interpolation_states)
    from rom import KROM
    result.krom=KROM(result.spectrum)
    result.metadata["krom"]=result.krom.metadata
    result.metadata["spectrum"]=result.spectrum.summary(); return result
