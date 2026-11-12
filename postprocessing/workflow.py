"""Typed source resolution and per-trajectory post-processing workflow state."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Literal

from provenance import serializable_system

from .core import Trajectory, integration_coordinates, load_solution


AnalysisCoordinates = Literal["archive", "integration"]


def _source_metadata(system: Any) -> dict[str, Any]:
    return {
        "system":getattr(system,"name",type(system).__name__),
        "system_parameters":serializable_system(system),
        "uniform_time_steps":True,
    }


class SourceResolver:
    """Flatten supported archive, result, hierarchy, and campaign sources."""

    def resolve(self, sources: Iterable[Any]) -> list[Trajectory]:
        resolved: list[Trajectory]=[]
        active: set[int]=set()
        for source in sources:
            self._visit(source,resolved,active)
        if not resolved:
            raise ValueError("processor has no completed trajectory sources.")
        return resolved

    def _visit(self,item: Any,resolved: list[Trajectory],active: set[int],system: Any=None) -> None:
        if isinstance(item,Trajectory):
            resolved.append(item); return
        if isinstance(item,(str,Path)):
            resolved.append(load_solution(item)); return
        archive_path=getattr(item,"archive_path",None)
        if archive_path is not None:
            resolved.append(load_solution(archive_path)); return

        identity=id(item)
        if identity in active:
            raise ValueError(f"Cyclic post-processing source graph at {type(item).__name__}.")
        active.add(identity)
        try:
            if (hasattr(item,"time_series") and item.time_series is not None
                    and not hasattr(item,"trajectories")):
                if system is None:
                    raise ValueError("In-memory series require their parent dynamical system for provenance.")
                resolved.append(Trajectory(item.time_series,_source_metadata(system),None))
            elif hasattr(item,"trajectories"):
                for child in item.trajectories:
                    self._visit(child,resolved,active,system)
            elif hasattr(item,"clusters"):
                parent_system=getattr(item,"dynamical_system",None)
                if parent_system is None:
                    raise ValueError("Cluster sources require their parent dynamical system for provenance.")
                for child in item.clusters:
                    self._visit(child,resolved,active,parent_system)
            elif hasattr(item,"systems"):
                for child in item.systems:
                    self._visit(child,resolved,active)
            else:
                raise TypeError(f"Unsupported post-processing source: {type(item).__name__}.")
        finally:
            active.remove(identity)


@dataclass(frozen=True)
class PreparedAnalysis:
    """One method's validated analysis source and algorithm-only settings."""

    source: Trajectory
    settings: dict[str,Any]
    coordinates: AnalysisCoordinates
    resampled: bool


@dataclass
class WorkflowContext:
    """Typed shared state for one ordered, single-trajectory workflow."""

    spectrum: Any=None
    spectrum_archive: Path|None=None
    _analysis_cache: dict[AnalysisCoordinates,Trajectory]=field(default_factory=dict,repr=False)

    def prepare(self,source: Trajectory,settings: dict[str,Any]) -> PreparedAnalysis:
        """Validate workflow controls and return settings owned by the method."""
        algorithm_settings=dict(settings)
        requested=algorithm_settings.pop("analysis_coordinates",None)
        if requested is None:
            coordinates: AnalysisCoordinates="archive"
        elif requested in {"archive","integration"}:
            coordinates=requested
        else:
            raise ValueError("analysis_coordinates must be 'archive', 'integration', or omitted.")

        resampling=algorithm_settings.pop("resampling",None)
        if resampling is None:
            analysis=self._analysis_cache.get(coordinates)
            if analysis is None:
                analysis=source if coordinates=="archive" else integration_coordinates(source)
                self._analysis_cache[coordinates]=analysis
            return PreparedAnalysis(analysis,algorithm_settings,coordinates,False)

        from .resampling import ResamplingOptions, resample
        if not isinstance(resampling,ResamplingOptions):
            raise TypeError("resampling must be a postprocessing.ResamplingOptions instance or None.")
        analysis=resample(source,resampling)
        if coordinates=="integration":
            analysis=integration_coordinates(analysis)
        return PreparedAnalysis(analysis,algorithm_settings,coordinates,True)

    def observe(self,result: Any,archive: Path|None) -> None:
        """Advance spectrum and lineage together when a method emits a spectrum."""
        spectrum=getattr(result,"spectrum",None)
        if spectrum is not None:
            self.spectrum=spectrum
            self.spectrum_archive=archive


__all__=["AnalysisCoordinates","PreparedAnalysis","SourceResolver","WorkflowContext"]
