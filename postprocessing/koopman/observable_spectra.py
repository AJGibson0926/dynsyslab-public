"""Multi-observable H-DMD spectrum construction and comparison."""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Literal

import numpy as np

from postprocessing.core import Trajectory, integration_coordinates, load_solution, save_derived
from time_series import TimeSeries

from .dmd import (HDMD_DEFAULT_MATRIX_SIZE, HDMD_DEFAULT_SAMPLING_WINDOW,
                  HDMD_DEFAULT_SVD_TOLERANCE, hankel_dmd, hankel_dmd_ensemble)


SpectrumDomain = Literal["discrete", "continuous"]


@dataclass
class ObservableSpectrumEntry:
    """One observable's H-DMD spectrum with its exact input contract."""

    name: str
    observable: dict[str, Any]
    discrete_eigenvalues: np.ndarray
    continuous_eigenvalues: np.ndarray
    result: Any | None = field(default=None, repr=False)

    def values(self, domain: SpectrumDomain) -> np.ndarray:
        if domain == "discrete": return np.asarray(self.discrete_eigenvalues, dtype=complex)
        if domain == "continuous": return np.asarray(self.continuous_eigenvalues, dtype=complex)
        raise ValueError("domain must be 'discrete' or 'continuous'.")


def _symmetric_nearest_rms(left: np.ndarray, right: np.ndarray) -> float:
    left=np.asarray(left,dtype=complex).reshape(-1); right=np.asarray(right,dtype=complex).reshape(-1)
    if not left.size or not right.size: return float("nan")
    distances=np.abs(left[:,None]-right[None,:])**2
    return float(np.sqrt(.5*(np.mean(np.min(distances,axis=1))+
                              np.mean(np.min(distances,axis=0)))))


@dataclass
class HDMDObservableComparison:
    """Spectra obtained with identical H-DMD settings and different observables."""

    entries: tuple[ObservableSpectrumEntry, ...]
    metadata: dict[str, Any] = field(default_factory=dict)
    archive: Path | None = None

    def __post_init__(self) -> None:
        if len(self.entries)<2: raise ValueError("Observable spectrum comparison requires at least two entries.")
        names=[entry.name for entry in self.entries]
        if len(set(names))!=len(names): raise ValueError("Observable spectrum comparison names must be unique.")

    def pairwise_distances(self, domain: SpectrumDomain="discrete") -> np.ndarray:
        values=[entry.values(domain) for entry in self.entries]
        output=np.zeros((len(values),len(values)),dtype=float)
        for row in range(len(values)):
            for column in range(row+1,len(values)):
                output[row,column]=output[column,row]=_symmetric_nearest_rms(values[row],values[column])
        return output

    def summary(self) -> dict[str, Any]:
        return {"encoding":"hdmd-observable-spectrum-comparison/v1",
            "observable_names":[entry.name for entry in self.entries],
            "spectrum_sizes":[int(entry.discrete_eigenvalues.size) for entry in self.entries],
            "observable_contracts":[entry.observable for entry in self.entries],
            "pairwise_discrete_nearest_rms":self.pairwise_distances("discrete").tolist(),
            "pairwise_continuous_nearest_rms":self.pairwise_distances("continuous").tolist(),
            "metadata":dict(self.metadata),
            "archive":None if self.archive is None else str(self.archive)}


def _selection_items(observables) -> list[tuple[str,Any]]:
    if isinstance(observables,Mapping): items=[(str(name),value) for name,value in observables.items()]
    elif isinstance(observables,Sequence) and not isinstance(observables,(str,bytes)):
        items=[]
        for index,value in enumerate(observables):
            name=(value if isinstance(value,str) else getattr(value,"name",None))
            items.append((str(name or f"observable_{index}"),value))
    else: raise TypeError("observables must be a mapping or a sequence of scalar observable selections.")
    if len(items)<2: raise ValueError("Supply at least two observables to compare H-DMD spectra.")
    return items


def compare_hankel_observables(
    states: np.ndarray | TimeSeries,
    observables,
    *,
    time: np.ndarray | None=None,
    system: Any | None=None,
    coordinate_system: str | None=None,
    metadata: dict[str,Any] | None=None,
    source_archive: str | None=None,
    matrix_size: int=HDMD_DEFAULT_MATRIX_SIZE,
    sampling_window: int=HDMD_DEFAULT_SAMPLING_WINDOW,
    rank: int | None=None,
    sample_stride: int=1,
    svd_tolerance: float | None=HDMD_DEFAULT_SVD_TOLERANCE,
    svd_driver: str="gesvd",
    svd_threads: int | None=None,
    eigendecomposition_threads: int | None=None,
    linear_algebra_backend: str="auto",
    max_memory_mb: float | None=8192.,
) -> HDMDObservableComparison:
    """Run one standard scalar H-DMD fit per observable on identical states.

    Mapping values may be named catalog entries, :class:`Observable` objects,
    state indices, or already evaluated scalar ``TimeSeries``/arrays. Mapping
    keys are the stable comparison labels and sampled-signal names.
    """
    call=hankel_dmd_ensemble if isinstance(states,TimeSeries) and states.is_pooled else hankel_dmd
    common={"matrix_size":matrix_size,"sampling_window":sampling_window,"rank":rank,
        "sample_stride":sample_stride,"matrix_layout":"window-snapshots",
        "mode_sort":"terminal-amplitude-energy","svd_tolerance":svd_tolerance,
        "svd_driver":svd_driver,"svd_threads":svd_threads,
        "eigendecomposition_threads":eigendecomposition_threads,
        "linear_algebra_backend":linear_algebra_backend,"metadata":metadata,
        "source_archive":source_archive,"max_memory_mb":max_memory_mb,"system":system}
    if coordinate_system is not None: common["coordinate_system"]=coordinate_system
    if call is hankel_dmd: common["time"]=time
    entries=[]
    for name,selection in _selection_items(observables):
        sampled=isinstance(selection,(TimeSeries,np.ndarray,list,tuple)) and not isinstance(selection,(str,bytes))
        options={**common,("observable_series" if sampled else "observable"):selection}
        if sampled: options["observable_name"]=name
        result=call(states,**options)
        eigenvalues=tuple(result.spectrum.eigenvalues)
        contract=dict(result.metadata["preprocessing"]["observable_contract"])
        entries.append(ObservableSpectrumEntry(name,contract,
            np.asarray([value.discrete_time for value in eigenvalues],dtype=complex),
            np.asarray([value.continuous_time for value in eigenvalues],dtype=complex),result))
    return HDMDObservableComparison(tuple(entries),{
        "matrix_size":int(matrix_size),"sampling_window":int(sampling_window),
        "rank":None if rank is None else int(rank),"sample_stride":int(sample_stride),
        "svd_tolerance":svd_tolerance,"svd_driver":svd_driver,
        "matrix_layout":"window-snapshots","mode_sort":"terminal-amplitude-energy",
        "coordinate_system":coordinate_system or
            ((metadata or {}).get("postprocessing_preprocessing") or {}).get("coordinate_system"),
        "same_state_training_data":True,"source_archive":source_archive})


def _comparison_arrays(comparison: HDMDObservableComparison) -> dict[str,np.ndarray]:
    counts=np.asarray([entry.discrete_eigenvalues.size for entry in comparison.entries],dtype=int)
    width=int(np.max(counts)); shape=(len(comparison.entries),width)
    discrete=np.full(shape,np.nan+1j*np.nan,dtype=complex)
    continuous=np.full(shape,np.nan+1j*np.nan,dtype=complex)
    for index,entry in enumerate(comparison.entries):
        discrete[index,:counts[index]]=entry.discrete_eigenvalues
        continuous[index,:counts[index]]=entry.continuous_eigenvalues
    return {"observable_names":np.asarray([entry.name for entry in comparison.entries]),
        "observable_contracts_json":np.asarray([json.dumps(entry.observable,sort_keys=True)
                                                 for entry in comparison.entries]),
        "mode_counts":counts,"discrete_eigenvalues":discrete,
        "continuous_eigenvalues":continuous,
        "pairwise_discrete_nearest_rms":comparison.pairwise_distances("discrete"),
        "pairwise_continuous_nearest_rms":comparison.pairwise_distances("continuous")}


def save_observable_spectrum_comparison(source: Trajectory,
                                        comparison: HDMDObservableComparison) -> Path:
    path=save_derived(source,"koopman-observable-spectrum-comparison",
        analysis_metadata={"comparison":comparison.summary()},**_comparison_arrays(comparison))
    comparison.archive=path
    return path


def load_observable_spectrum_comparison(path: str|Path) -> HDMDObservableComparison:
    target=Path(path).expanduser().resolve()
    with np.load(target,allow_pickle=False) as data:
        names=np.asarray(data["observable_names"]).astype(str)
        contracts=[json.loads(str(value)) for value in data["observable_contracts_json"]]
        counts=np.asarray(data["mode_counts"],dtype=int)
        discrete=np.asarray(data["discrete_eigenvalues"],dtype=complex)
        continuous=np.asarray(data["continuous_eigenvalues"],dtype=complex)
        outer=json.loads(str(data["metadata_json"]))
    entries=tuple(ObservableSpectrumEntry(str(name),contract,
        discrete[index,:counts[index]],continuous[index,:counts[index]])
        for index,(name,contract) in enumerate(zip(names,contracts)))
    analysis=outer.get("analysis_metadata") or {}
    metadata=dict((analysis.get("comparison") or {}).get("metadata") or {})
    metadata["source_archive"]=outer.get("source_archive")
    return HDMDObservableComparison(entries,metadata,target)


def compare_hankel_observables_archive(path: str|Path, observables, *,
                                       archive: bool=True,
                                       analysis_coordinates: str="integration", **settings):
    source=load_solution(path)
    if analysis_coordinates=="integration": analysis=integration_coordinates(source)
    elif analysis_coordinates=="archive": analysis=source
    else: raise ValueError("analysis_coordinates must be 'integration' or 'archive'.")
    comparison=compare_hankel_observables(analysis.time_series,observables,
        metadata=analysis.metadata,source_archive=str(source.source),**settings)
    comparison.metadata["analysis_coordinates"]=analysis_coordinates
    if archive: save_observable_spectrum_comparison(source,comparison)
    return comparison


def observable_spectrum_comparison_report(source: HDMDObservableComparison|str|Path,
                                          *, domain: SpectrumDomain="discrete"):
    """Compatibility gateway to the Viewer-owned presentation adapter."""
    comparison=(load_observable_spectrum_comparison(source)
                if isinstance(source,(str,Path)) else source)
    from importlib import import_module
    adapter=import_module("viewer.observable_spectra")
    return adapter.observable_spectrum_comparison_report(comparison,domain=domain)


__all__=["HDMDObservableComparison","ObservableSpectrumEntry","SpectrumDomain",
    "compare_hankel_observables","compare_hankel_observables_archive",
    "load_observable_spectrum_comparison","observable_spectrum_comparison_report",
    "save_observable_spectrum_comparison"]
