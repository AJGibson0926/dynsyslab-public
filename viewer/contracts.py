"""Versioned in-process contracts for read-only Viewer composition."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, Sequence
import numpy as np

Level=Literal["trajectory","cluster","system","batch"]
PlotKind=Literal["time","fourier","phase","three_d","observable_time","observable_heatmap","koopman_spectrum","eigenfunction","krom","pod_mode","bifurcation","control_error","observability","klqe"]
EigenfunctionDisplay=Literal["domain_color","real","imaginary"]
EigenfunctionView=Literal["domain","phase_portrait"]
TimeSeriesLayout=Literal["separate","combined"]
Renderer=Literal["svg","plotly"]
ControlErrorMeasure=Literal["relative","absolute","exact_norm"]
ObservableComponent=Literal["auto","real","imaginary","magnitude","phase"]
_KROM_MODE_COLORS=("#e41a1c","#ff7f00","#4daf4a","#377eb8","#984ea3","#a65628","#f781bf","#17becf")
_DEFAULT_TRACE_COLORS=("#3478c8","#e45756","#54a24b","#f2cf5b","#b279a2","#ff9da6","#9d755d","#bab0ac")

@dataclass(frozen=True)
class PlotStyle:
    """Uniform defaults, with scalar-or-per-trajectory style inputs."""
    colors: str | Sequence[str] = _DEFAULT_TRACE_COLORS
    opacity: float | Sequence[float] = 1.0
    line_width: float | Sequence[float] = 1.5
    line_style: str | Sequence[str] = "solid"
    marker_size: float = 3.0

@dataclass(frozen=True)
class TilePlan:
    """Layout directives for grouping generated figures."""
    name: str="single"
    columns: int=1
    rows: int | None=None
    @classmethod
    def preset(cls,name: str) -> "TilePlan":
        presets={"single":(1,1),"vertical":(1,None),"horizontal":(None,1),"grid_2x2":(2,2),"time_spectrum":(2,1),"time_phase":(2,1),"dashboard_2x2":(2,2)}
        if name not in presets: raise ValueError(f"Unknown tile preset {name!r}; choose from {sorted(presets)}.")
        columns,rows=presets[name]
        return cls(name,columns or 1,rows)

@dataclass(frozen=True)
class ArchiveSource:
    """An existing solution archive or hierarchy manifest to view."""
    path: str | Path
    level: Level="trajectory"
    system_index: int | None=None
    cluster_index: int | None=None
    trajectory_index: int | None=None
    label: str | None=None
    # Set explicitly for a solver checkpoint/partial record. ``None`` asks
    # Viewer to infer incompleteness from the archived status metadata.
    incomplete: bool | None=None

@dataclass(frozen=True)
class KoopmanSource:
    """An existing derived archive or in-memory KoopmanSpectrum to view."""
    source: str | Path | Any
    label: str | None=None
    # Older sparse-eigenfunction archives may reuse, rather than duplicate,
    # eigenvalues.  Supply the exact archived source spectrum when its durable
    # spectrum-id link cannot be resolved from the archive index.
    eigenvalue_source: str | Path | None=None

@dataclass(frozen=True)
class PODSource:
    """An in-memory PROM or SVD result containing a PROM."""
    source: Any
    label: str | None=None

@dataclass(frozen=True)
class ControlErrorSource:
    """Typed control-error field(s) or a saved control-error archive."""
    source: str | Path | Any
    label: str | None=None

@dataclass(frozen=True)
class ObservabilitySource:
    """An in-memory or saved :class:`rom.ObservabilityAnalysis`."""
    source: str | Path | Any
    label: str | None=None

@dataclass(frozen=True)
class KLQESource:
    """An in-memory or saved :class:`estimation.KLQEResult`."""
    source: str | Path | Any
    label: str | None=None

@dataclass(frozen=True)
class BifurcationSource:
    """A JSON manifest mapping parameter values to existing solution archives."""
    manifest_path: str | Path
    parameter_key: str="parameter"
    tail_fraction: float=.25

@dataclass(frozen=True)
class CampaignSource:
    """An existing Campaigner manifest of any supported campaign type.

    Viewer classifies it as bifurcation, Lyapunov, trajectory-ensemble,
    parameter-sweep, or generic, then delegates only disposable display-data
    preparation to its internal compute service. It can never trigger a job,
    extension, solve, or primary-archive postprocessing operation.
    """
    campaign_path: str | Path
    tail_fraction: float=.25

@dataclass(frozen=True)
class _ViewerArtifact:
    """Private bridge from the internal compute service back to Viewer."""
    path: str | Path
    label: str | None=None

@dataclass(frozen=True)
class ViewerDataSource:
    """A named source registered for the current viewing session."""
    name: str
    data_type: str
    source: ArchiveSource | KoopmanSource | PODSource | ControlErrorSource | ObservabilitySource | KLQESource | BifurcationSource | CampaignSource
    campaign_type: str | None=None

@dataclass(frozen=True)
class PlotRequest:
    """Describe one read-only plot without performing scientific computation."""

    kind: PlotKind
    source: ArchiveSource | KoopmanSource | PODSource | ControlErrorSource | ObservabilitySource | KLQESource | BifurcationSource | CampaignSource | None=None
    source_name: str | None=None
    state_components: tuple[int,...]=(0,1)
    # Stable catalog key or an explicit dynsyslib Observable. Observable
    # displays are read-only evaluations of the supplied trajectory/system.
    observable: Any | None=None
    observable_component: ObservableComponent="auto"
    eigenfunction_index: int=0
    # Dashboard selection is a Viewer concern only; it never changes the KROM
    # spectrum.  By default, display the first five nonredundant KROM modes.
    # Set explicit spectrum indices, or set the limit to None to request every
    # retained mode (potentially a very large report).
    krom_mode_indices: tuple[int,...] | None=None
    krom_mode_limit: int | None=5
    pod_mode_index: int=0
    eigenfunction_display: EigenfunctionDisplay="domain_color"
    # ``None`` selects the dimension-aware default: domain in 2-D and phase
    # portrait in 3-D.
    eigenfunction_view: EigenfunctionView | None=None
    domain: tuple[tuple[float,float],tuple[float,float]] | None=None
    domain_density: int=100
    state_samples: np.ndarray | None=None
    max_display_samples: int | None=None
    # Separate state-component panels are the legible default. ``combined``
    # retains the compact multi-trace display for callers that prefer it.
    time_series_layout: TimeSeriesLayout="separate"
    color_scheme: str="domain"  # domain (phase/magnitude) or jet (phase mapped through jet)
    show_colorbar: bool=False
    style: PlotStyle=field(default_factory=PlotStyle)
    spectrum_scale: str="linear"  # linear (default) or log_power
    control_error_measure: ControlErrorMeasure="relative"
    control_error_log_scale: bool=True
    # None selects automatic behavior: unfinished campaigns refresh; finished
    # and ordinary archived sources remain static unless explicitly enabled.
    refreshable: bool | None=None
    title: str | None=None

@dataclass
class ViewerFigure:
    """Renderer-neutral traces and metadata for one Viewer panel."""

    kind: PlotKind
    title: str
    traces: list[dict[str,Any]]
    metadata: dict[str,Any]=field(default_factory=dict)

@dataclass
class ViewerReport:
    """An ordered set of Viewer panels with layout and warning metadata."""

    figures: list[ViewerFigure]
    tile_plan: TilePlan
    warnings: list[str]=field(default_factory=list)
