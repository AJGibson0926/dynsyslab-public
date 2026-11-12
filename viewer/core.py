"""Read-only archive display with an internal disposable-compute delegate.

Viewer never imports ``solver``, ``job_manager``, or ``campaigner`` and never
changes the primary archive. Requested display derivations are delegated to
the internal ``viewer_compute`` service, which writes only its temporary
secondary archive. Viewer then reads those artifacts and renders graphics.
"""
from __future__ import annotations
from dataclasses import replace
import json
from pathlib import Path
from typing import Any
import numpy as np
from archiver.gateway import VIEWER_ARCHIVE
from time_series import TimeSeries
from ._compute import ViewerCompute
from .settings import get_viewer_settings
from . import sources as _sources
from . import routing as _routing
from .session import ViewingSession

from .contracts import (
    _ViewerArtifact, ArchiveSource, BifurcationSource, CampaignSource,
    ControlErrorMeasure, ControlErrorSource, EigenfunctionView, KLQESource,
    KoopmanSource, ObservableComponent, ObservabilitySource, PlotKind, PlotRequest,
    PlotStyle, PODSource, Renderer, TilePlan, TimeSeriesLayout,
    ViewerDataSource, ViewerFigure, ViewerReport,)

from .rendering.helpers import (
    _automatic_time_detail_inset, _control_display, _display_samples,
    _interpolation_kind, _krom_colors, _krom_display_indices,
    _limit_display_traces, _sequence, _trace_display_indices,)

from .records import (
    _axis_context, _display_analysis_series, _display_series,
    _load_archive_source, _manifest, _solution_record,
    _trajectory_paths, _viewer_compute_record,)

from .rendering.analysis import AnalysisRenderMixin
from .rendering.campaign import CampaignRenderMixin
from .rendering.koopman import KoopmanRenderMixin
from .rendering.serialization import (SerializationMixin, _axis_ticks,
    _jet_rgb, _svg_colorbar, _svg_figure, _trace_color,
    _write_interactive_html)

class Viewer(AnalysisRenderMixin, KoopmanRenderMixin, CampaignRenderMixin, SerializationMixin):
    """Configure sources and requests first; call render only after configuration."""
    _viewer_compute_record = staticmethod(_viewer_compute_record)
    def __init__(self, *, tile_plan: TilePlan | None=None, default_style: PlotStyle | None=None,
                 renderer: Renderer="svg"):
        self.tile_plan=tile_plan or TilePlan.preset("single")
        self.default_style=default_style or PlotStyle()
        if renderer not in {"svg","plotly"}: raise ValueError("renderer must be 'svg' or 'plotly'.")
        self.renderer=renderer
        self._requests: list[PlotRequest]=[]
        self._sources: dict[str,ViewerDataSource]={}
        self._refreshable_paths: set[str]=set()
        # This is a lifecycle delegation only: Viewer never derives data or
        # writes primary archives.  The compute service owns its temporary
        # archive and clears it for every fresh viewing session.
        self._session=ViewingSession()
        self.viewer_compute_session_id=self._session.identifier
        self._compute=self._session.compute

    def add(self, request: PlotRequest) -> "Viewer":
        self._requests.append(request); return self

    @staticmethod
    def _source_type(source: Any) -> str:
        value=_routing.classify(source)
        if value!="unknown": return value
        raise TypeError(f"Unsupported Viewer data source {type(source).__name__}.")

    def attach_source(self, source: ArchiveSource | KoopmanSource | PODSource | ControlErrorSource | ObservabilitySource | KLQESource | BifurcationSource | CampaignSource,
                      *, name: str | None=None) -> "Viewer":
        """Register a named source for automatic plot-request resolution."""
        data_type=self._source_type(source); proposed=name or getattr(source,"label",None) or f"{data_type}-{len(self._sources)+1}"
        if proposed in self._sources: raise ValueError(f"Viewer already has a data source named {proposed!r}.")
        campaign_type=None
        if isinstance(source,CampaignSource):
            payload=VIEWER_ARCHIVE.read_json(source.campaign_path); campaign_type=ViewerCompute.campaign_type(payload)
        self._sources[proposed]=ViewerDataSource(proposed,data_type,source,campaign_type)
        self._sync_subscriptions()
        return self

    def _sync_subscriptions(self) -> None:
        """Publish passive source interests to the secondary session registry."""
        self._session.synchronize(self._sources.values())

    def poll_refresh_events(self) -> list[dict[str,Any]]:
        """Read new passive source-update events without changing any archive."""
        return self._session.poll()

    def refresh_if_updated(self) -> ViewerReport | None:
        """Re-render only when a subscribed source has reported an update."""
        events=self.poll_refresh_events()
        return self.render() if any(event.get("source_path") in self._refreshable_paths for event in events) else None

    def attach_sources(self, *sources: ArchiveSource | KoopmanSource | PODSource | ControlErrorSource | ObservabilitySource | KLQESource | BifurcationSource | CampaignSource) -> "Viewer":
        """Attach several sources using their labels or generated session names."""
        for source in sources: self.attach_source(source)
        return self

    @property
    def data_sources(self) -> tuple[ViewerDataSource,...]:
        """Read-only registry of current session sources and their data types."""
        return tuple(self._sources.values())

    def configure(self, *, tile_plan: TilePlan | None=None, default_style: PlotStyle | None=None,
                  renderer: Renderer | None=None) -> "Viewer":
        if tile_plan is not None: self.tile_plan=tile_plan
        if default_style is not None: self.default_style=default_style
        if renderer is not None:
            if renderer not in {"svg","plotly"}: raise ValueError("renderer must be 'svg' or 'plotly'.")
            self.renderer=renderer
        return self

    def render(self) -> ViewerReport:
        if any(request.kind=="krom" for request in self._requests) and self.tile_plan.name=="single":
            self.tile_plan=TilePlan("krom_dashboard",columns=4)
        figures=[]; warnings=[]; self._refreshable_paths=set()
        for requested in self._requests:
            try:
                resolved=self._resolve_request_source(requested)
                refreshable=self._is_refreshable(resolved)
                if refreshable:
                    path=self._source_path(resolved.source)
                    if path is not None: self._refreshable_paths.add(path)
                for request in self._materialize_requests(resolved):
                    rendered=self._render_request(request)
                    for figure in rendered:
                        figure.metadata["refreshable"]=refreshable
                        if figure.kind in {"time","observable_time","klqe"} and "detail_inset" not in figure.metadata:
                            inset=_automatic_time_detail_inset(figure.traces)
                            if inset is not None: figure.metadata["detail_inset"]=inset
                    figures.extend(rendered)
            except (ValueError,KeyError,FileNotFoundError,NotImplementedError,RuntimeError) as error:
                warnings.append(f"{requested.kind}: {error}")
        if self.tile_plan.name=="time_phase":
            for figure in figures:
                if figure.kind=="phase" or figure.metadata.get("time_series_layout")=="combined" or figure.metadata.get("control_history"):
                    figure.metadata.setdefault("grid_span",self.tile_plan.columns)
        return ViewerReport(figures,self.tile_plan,warnings)

    @staticmethod
    def _source_path(source: Any) -> str | None:
        if isinstance(source,ArchiveSource): return str(Path(source.path).expanduser().resolve())
        if isinstance(source,CampaignSource): return str(Path(source.campaign_path).expanduser().resolve())
        if isinstance(source,BifurcationSource): return str(Path(source.manifest_path).expanduser().resolve())
        if isinstance(source,KoopmanSource) and isinstance(source.source,(str,Path)): return str(Path(source.source).expanduser().resolve())
        if isinstance(source,ControlErrorSource) and isinstance(source.source,(str,Path)): return str(Path(source.source).expanduser().resolve())
        if isinstance(source,(ObservabilitySource,KLQESource)) and isinstance(source.source,(str,Path)): return str(Path(source.source).expanduser().resolve())
        return None

    def _is_refreshable(self, request: PlotRequest) -> bool:
        if request.refreshable is not None: return bool(request.refreshable)
        if isinstance(request.source,ArchiveSource):
            if request.source.incomplete is not None: return bool(request.source.incomplete)
            path=Path(request.source.path).expanduser()
            try:
                metadata=_manifest(path)
                status=str(metadata.get("status",metadata.get("solve_status","complete"))).lower()
                if status in {"failed","failure","error"}: return False
                return status not in {"complete","completed","solve completed"}
            except Exception: return "history" in str(path).lower() or "partial" in str(path).lower()
        if not isinstance(request.source,CampaignSource): return False
        path=Path(request.source.campaign_path).expanduser()
        if not path.is_file(): return False
        campaign=VIEWER_ARCHIVE.read_json(path); jobs=VIEWER_ARCHIVE.coordination_state().get("jobs",{})
        return any(jobs.get(job_id,{}).get("status") not in {"completed","failed"} for job_id in campaign.get("jobs",[]))

    @staticmethod
    def _required_source_type(kind: PlotKind) -> tuple[str,...]:
        return _routing.required(kind)

    def _resolve_request_source(self, request: PlotRequest) -> PlotRequest:
        if request.source is not None:
            if request.source_name is not None: raise ValueError("Specify either source or source_name, not both.")
            return request
        if request.source_name is not None:
            if request.source_name not in self._sources: raise ValueError(f"No Viewer source named {request.source_name!r}; available: {sorted(self._sources)}.")
            return replace(request,source=self._sources[request.source_name].source)
        expected=self._required_source_type(request.kind)
        matches=_routing.matching(list(self._sources.values()),request.kind)
        if len(matches)==1: return replace(request,source=matches[0].source)
        if not matches: raise ValueError(f"No {', '.join(expected)} source is attached for a {request.kind!r} plot.")
        names=[entry.name for entry in matches]
        raise ValueError(f"Multiple compatible sources are attached ({names}); specify PlotRequest(source_name=...).")

    def _materialize_requests(self, request: PlotRequest) -> list[PlotRequest]:
        """Delegate needed display derivations to the session-local service."""
        if isinstance(request.source,_ViewerArtifact): return [request]
        if isinstance(request.source,CampaignSource):
            self._compute.register_campaign_view(request.source.campaign_path)
            payload=VIEWER_ARCHIVE.read_json(request.source.campaign_path); campaign_type=ViewerCompute.campaign_type(payload)
            if request.kind=="bifurcation":
                if campaign_type!="bifurcation": raise ValueError(f"Campaign type {campaign_type!r} does not supply bifurcation data.")
                artifact=self._compute.bifurcation_from_campaign(request.source.campaign_path,tail_fraction=request.source.tail_fraction,label="campaign-bifurcation")
            elif request.kind=="fourier": artifact=self._compute.campaign_fourier(request.source.campaign_path,components=request.state_components,label=f"campaign-{campaign_type}-fourier")
            elif request.kind in {"time","phase","three_d"}: artifact=self._compute.campaign_trajectories(request.source.campaign_path,label=f"campaign-{campaign_type}-trajectories")
            else: raise ValueError(f"Campaign type {campaign_type!r} cannot supply a {request.kind!r} display directly.")
            return [replace(request,source=_ViewerArtifact(artifact))]
        if request.kind=="fourier":
            if not isinstance(request.source,ArchiveSource): raise ValueError("Fourier plots require an ArchiveSource.")
            entries=_load_archive_source(request.source)
            return [replace(request,source=_ViewerArtifact(self._compute.fourier(metadata["source_archive"],components=request.state_components,label=f"fourier-{index}"))) for index,(_,_,metadata) in enumerate(entries)]
        if request.kind=="bifurcation":
            if isinstance(request.source,CampaignSource):
                artifact=self._compute.bifurcation_from_campaign(request.source.campaign_path,tail_fraction=request.source.tail_fraction,label="campaign-bifurcation")
            elif isinstance(request.source,BifurcationSource):
                artifact=self._compute.bifurcation_from_manifest(request.source.manifest_path,parameter_key=request.source.parameter_key,tail_fraction=request.source.tail_fraction,label="manifest-bifurcation")
            else: raise ValueError("Bifurcation plots require CampaignSource or BifurcationSource.")
            return [replace(request,source=_ViewerArtifact(artifact))]
        if request.kind=="eigenfunction" and isinstance(request.source,KoopmanSource):
            spectrum=request.source.source
            if isinstance(spectrum,(str,Path)): raise ValueError("Eigenfunction computation requires an in-memory spectrum with eigenfunctions.")
            if request.state_samples is None: raise ValueError("Eigenfunction computation requires state_samples.")
            system=getattr(spectrum,"dynamical_system",None); dimension=int(getattr(system,"state_dimension",np.asarray(request.state_samples).shape[1]))
            view=request.eigenfunction_view or ("phase_portrait" if dimension==3 else "domain")
            artifact=self._compute.koopman_eigenfunction(spectrum,request.eigenfunction_index,state_samples=request.state_samples,view=view,domain=request.domain,density=request.domain_density,label=f"koopman-eigenfunction-{request.eigenfunction_index}")
            return [replace(request,source=_ViewerArtifact(artifact))]
        if request.kind=="pod_mode" and isinstance(request.source,PODSource):
            prom=getattr(request.source.source,"prom",request.source.source)
            if request.state_samples is None: raise ValueError("POD-mode computation requires state_samples.")
            dimension=np.asarray(request.state_samples).shape[1]; view=request.eigenfunction_view or ("phase_portrait" if dimension==3 else "domain")
            artifact=self._compute.pod_mode(prom,request.pod_mode_index,state_samples=request.state_samples,view=view,domain=request.domain,density=request.domain_density,label=f"pod-mode-{request.pod_mode_index}")
            return [replace(request,source=_ViewerArtifact(artifact))]
        return [request]

    def _render_request(self, request: PlotRequest) -> list[ViewerFigure]:
        if isinstance(request.source,_ViewerArtifact):
            _,metadata=_viewer_compute_record(request.source)
            if metadata.get("kind")=="campaign_trajectories": return [self._campaign_trajectory_figure(request)]
            if metadata.get("kind")=="campaign_fourier": return [self._campaign_fourier_figure(request)]
        if request.kind in {"time","phase","three_d"}: return self._state_figures(request)
        if request.kind in {"observable_time","observable_heatmap"}: return self._observable_figures(request)
        if request.kind=="fourier":
            return [self._viewer_compute_figure(request,"fourier")]
        if request.kind in {"koopman_spectrum","eigenfunction"}: return [self._koopman_figure(request)]
        if request.kind=="krom": return self._krom_figures(request)
        if request.kind=="control_error": return self._control_error_figures(request)
        if request.kind=="observability": return self._observability_figures(request)
        if request.kind=="klqe": return self._klqe_figures(request)
        # POD amplitudes are real by construction.  An imaginary-component
        # request is omitted instead of producing a blank panel.
        if request.kind=="pod_mode":
            if request.eigenfunction_display=="imaginary": return []
            return [self._viewer_compute_figure(request,"pod_mode")] if isinstance(request.source,_ViewerArtifact) else [self._pod_mode_figure(request)]
        if request.kind=="bifurcation": return [self._viewer_compute_figure(request,"bifurcation")]
        raise ValueError(f"Unsupported plot kind {request.kind!r}.")

    @staticmethod
    def _observable_display_values(values: np.ndarray, component: ObservableComponent) -> tuple[np.ndarray,str]:
        data=np.asarray(values).reshape(-1)
        if component not in {"auto","real","imaginary","magnitude","phase"}:
            raise ValueError("observable_component must be auto, real, imaginary, magnitude, or phase.")
        if not np.iscomplexobj(data):
            if component in {"imaginary","phase"}: return np.zeros(data.size),component
            return np.asarray(data,dtype=float),("value" if component=="auto" else component)
        selected="magnitude" if component=="auto" else component
        converted={"real":data.real,"imaginary":data.imag,"magnitude":np.abs(data),"phase":np.angle(data)}.get(selected)
        if converted is None: raise ValueError("Complex observables require real, imaginary, magnitude, phase, or auto display.")
        return np.asarray(converted,dtype=float),selected

    @staticmethod
    def _system_observable(metadata: dict[str,Any], requested: Any):
        from dynsyslib import Observable, resolve_archived_system
        system=resolve_archived_system(str(metadata.get("system")),dict(metadata.get("system_parameters") or {}))
        coordinate=str(metadata.get("integration_form") or "native")
        observable=requested if isinstance(requested,Observable) else system.observable(str(requested),coordinate_system=coordinate)
        return system,observable,coordinate

    def _observable_figures(self, request: PlotRequest) -> list[ViewerFigure]:
        """Evaluate a scalar model observable along trajectories or over a 2-D domain."""
        if not isinstance(request.source,ArchiveSource):
            raise ValueError("Observable displays require an ArchiveSource.")
        if request.observable is None:
            raise ValueError("Observable displays require PlotRequest(observable=...).")
        entries=_load_archive_source(request.source)
        if not entries: raise ValueError("The observable source contains no trajectories.")
        groups: dict[int,list[tuple[str,TimeSeries,dict[str,Any]]]]={}
        for entry in entries: groups.setdefault(entry[2].get("system_index",0),[]).append(entry)
        figures=[]
        for group_index,group in groups.items():
            system,observable,coordinate=self._system_observable(group[0][2],request.observable)
            if observable.output_dimension!=1:
                raise ValueError(f"Observable {observable.name!r} is vector-valued; scalar Viewer displays require output_dimension=1.")
            label=str(observable.descriptor.get("label") or observable.name)
            units=str(observable.descriptor.get("units") or "")
            suffix=f" ({units})" if units and units!="dimensionless" else ""
            if request.kind=="observable_time":
                style=request.style or self.default_style; colors=_sequence(style.colors,len(group)); opacity=_sequence(style.opacity,len(group)); widths=_sequence(style.line_width,len(group)); traces=[]
                displayed_component=None; context=None
                for index,(trajectory_label,series,metadata) in enumerate(group):
                    _,_,current=_display_series(series,metadata); context=context or current
                    values,displayed_component=self._observable_display_values(observable.evaluate(series.values)[:,0],request.observable_component)
                    traces.append({"label":f"{trajectory_label} · {label}","x":series.time*current["time_scale"],"y":values,
                                   "color":colors[index],"opacity":opacity[index],"width":widths[index]})
                decorated=(f"|{label}|" if displayed_component=="magnitude" else
                           f"arg({label})" if displayed_component=="phase" else
                           f"{displayed_component}({label})" if displayed_component in {"real","imaginary"} else label)
                shown,sampling=_limit_display_traces(traces,request.max_display_samples)
                figures.append(ViewerFigure("observable_time",request.title or f"observable: {decorated}",shown,
                    {"system_index":group_index,"system":group[0][2].get("system"),"observable":observable.summary(),
                     "observable_component":displayed_component,"axes":{"x":context["time_label"],"y":decorated+suffix},**sampling}))
                continue
            if int(getattr(system,"state_dimension",0))!=2:
                raise ValueError("Observable heat maps are available only for two-dimensional systems.")
            all_states=np.vstack([series.values for _,series,_ in group])
            if request.domain is None:
                lower=np.nanmin(all_states,axis=0); upper=np.nanmax(all_states,axis=0); span=np.maximum(upper-lower,1e-9)
                domain=((float(lower[0]-.05*span[0]),float(upper[0]+.05*span[0])),
                        (float(lower[1]-.05*span[1]),float(upper[1]+.05*span[1])))
            else: domain=request.domain
            density=int(request.domain_density)
            if density<2: raise ValueError("Observable heat-map domain_density must be at least 2.")
            x_values=np.linspace(*domain[0],density); y_values=np.linspace(*domain[1],density)
            xx,yy=np.meshgrid(x_values,y_values); states=np.column_stack((xx.ravel(),yy.ravel()))
            values,displayed_component=self._observable_display_values(observable.evaluate(states)[:,0],request.observable_component)
            labels=system.state_component_labels(coordinate)
            color_label=(f"|{label}|" if displayed_component=="magnitude" else
                         f"arg({label})" if displayed_component=="phase" else
                         f"{displayed_component}({label})" if displayed_component in {"real","imaginary"} else label)+suffix
            figures.append(ViewerFigure("observable_heatmap",request.title or f"observable heat map: {label}",[
                {"x":x_values,"y":y_values,"z":values.reshape(density,density),"heatmap":True,"colorbar_title":color_label}],
                {"system_index":group_index,"system":group[0][2].get("system"),"observable":observable.summary(),
                 "observable_component":displayed_component,"coordinate_system":coordinate,
                 "axes":{"x":labels[0],"y":labels[1]},"colorbar_title":color_label,"height":500,"square":True}))
        return figures

    def _state_figures(self, request: PlotRequest) -> list[ViewerFigure]:
        if not isinstance(request.source,ArchiveSource): raise ValueError("State, Fourier, phase, and 3D plots require ArchiveSource.")
        entries=_load_archive_source(request.source)
        # A batch deliberately becomes one independent plot group per System.
        groups: dict[int,list[tuple[str,TimeSeries,dict[str,Any]]]]={}
        for entry in entries: groups.setdefault(entry[2].get("system_index",0),[]).append(entry)
        figures=[]
        for group_index,group in groups.items():
            style=request.style or self.default_style; colors=_sequence(style.colors,len(group)); opacity=_sequence(style.opacity,len(group)); widths=_sequence(style.line_width,len(group)); traces=[]
            time_colors=_sequence(style.colors,len(group)*len(request.state_components)) if request.kind=="time" else []
            component_traces={component:[] for component in request.state_components} if request.kind=="time" else {}
            control_traces: dict[int,list[dict[str,Any]]]={}
            context=None
            for index,(label,series,metadata) in enumerate(group):
                display_time,display_values,current=_display_series(series,metadata); context=context or current
                target_series=_sources.control_target_series(metadata)
                control_series=_sources.control_series(metadata) if request.kind=="time" else None
                target_time,target_values=(None,None)
                if target_series is not None:
                    target_time,target_values,_=_display_series(target_series,metadata)
                if control_series is not None:
                    for control_component in range(control_series.values.shape[1]):
                        control_name,control_scale=_control_display(metadata,control_component)
                        control_traces.setdefault(control_component,[]).append({"label":f"{label} · {control_name}","x":control_series.time*current["time_scale"],"y":control_series.values[:,control_component]*control_scale,"color":colors[index],"opacity":opacity[index],"width":widths[index],"role":"control_history"})
                if request.kind=="time":
                    for component_index,component in enumerate(request.state_components):
                        state_label=current["state_labels"][component] if component<len(current["state_labels"]) else f"x{component}"
                        trace={"label":f"{label} · {state_label}","x":display_time,"y":display_values[:,component],"color":time_colors[index*len(request.state_components)+component_index],"opacity":opacity[index],"width":widths[index]}
                        traces.append(trace); component_traces[component].append(trace)
                        if target_values is not None and component<target_values.shape[1]:
                            target_trace={"label":f"{label} · control target · {state_label}","x":target_time,"y":target_values[:,component],"color":"#d62728","opacity":1.,"width":2.,"line_style":"dash","role":"control_target"}
                            traces.append(target_trace); component_traces[component].append(target_trace)
                elif request.kind=="phase":
                    if len(request.state_components)!=2: raise ValueError("A phase portrait needs exactly two state components.")
                    a,b=request.state_components; traces.append({"label":label,"x":display_values[:,a],"y":display_values[:,b],"color":colors[index],"opacity":opacity[index],"width":widths[index]})
                    if target_values is not None and max(a,b)<target_values.shape[1]:
                        traces.append({"label":f"{label} · control target","x":target_values[:,a],"y":target_values[:,b],"color":"#d62728","opacity":1.,"width":2.,"line_style":"dash","role":"control_target"})
                else:
                    if len(request.state_components)!=3: raise ValueError("A 3D plot needs exactly three state components.")
                    a,b,c=request.state_components; traces.append({"label":label,"x":display_values[:,a],"y":display_values[:,b],"z":display_values[:,c],"color":colors[index],"opacity":opacity[index],"width":widths[index]})
            title=request.title or f"{request.kind}: system {group_index}"
            if request.kind=="time":
                if request.time_series_layout not in {"separate","combined"}: raise ValueError("time_series_layout must be 'separate' or 'combined'.")
                if request.time_series_layout=="combined":
                    axes={"x":context["time_label"],"y":" / ".join(context["state_labels"][item] if item<len(context["state_labels"]) else f"x{item}" for item in request.state_components)}
                    display_traces,sampling=_limit_display_traces(traces,request.max_display_samples)
                    figures.append(ViewerFigure(request.kind,title,display_traces,{"system_index":group_index,"system":group[0][2].get("system"),"axes":axes,"time_series_layout":"combined","control_target_overlay":any(trace.get("role")=="control_target" for trace in traces),**sampling}))
                else:
                    for component in request.state_components:
                        state_label=context["state_labels"][component] if component<len(context["state_labels"]) else f"x{component}"
                        component_title=f"{title} · {state_label}" if request.title else f"time: system {group_index} · {state_label}"
                        display_traces,sampling=_limit_display_traces(component_traces[component],request.max_display_samples)
                        figures.append(ViewerFigure(request.kind,component_title,display_traces,{"system_index":group_index,"system":group[0][2].get("system"),"state_component":component,"axes":{"x":context["time_label"],"y":state_label},"time_series_layout":"separate","control_target_overlay":any(trace.get("role")=="control_target" for trace in component_traces[component]),**sampling}))
                for control_component,control_component_traces in sorted(control_traces.items()):
                    control_name,_=_control_display(group[0][2],control_component)
                    display_traces,sampling=_limit_display_traces(control_component_traces,request.max_display_samples)
                    figures.append(ViewerFigure("time",f"control: system {group_index} · {control_name}",display_traces,{"system_index":group_index,"system":group[0][2].get("system"),"control_component":control_component,"axes":{"x":context["time_label"],"y":control_name},"control_history":True,**sampling}))
                continue
            elif request.kind=="phase": axes={"x":context["state_labels"][request.state_components[0]] if request.state_components[0]<len(context["state_labels"]) else f"x{request.state_components[0]}","y":context["state_labels"][request.state_components[1]] if request.state_components[1]<len(context["state_labels"]) else f"x{request.state_components[1]}"}
            else: axes={}
            figure_metadata={"system_index":group_index,"system":group[0][2].get("system"),"axes":axes,"control_target_overlay":any(trace.get("role")=="control_target" for trace in traces)}
            if request.kind=="phase": figure_metadata.update({"height":520,"square":True})
            display_traces,sampling=_limit_display_traces(traces,request.max_display_samples)
            figure_metadata.update(sampling)
            figures.append(ViewerFigure(request.kind,title,display_traces,figure_metadata))
        return figures

    def _krom_figures(self, request: PlotRequest) -> list[ViewerFigure]:
        """Render the standard archived-KROM spectrum and retained-mode dashboard."""
        if not isinstance(request.source,KoopmanSource) or not isinstance(request.source.source,(str,Path)):
            raise ValueError("The standard KROM dashboard requires an archived KoopmanSource.")
        path=Path(request.source.source).expanduser().resolve()
        with VIEWER_ARCHIVE.open_read(path) as data:
            required={"training_states","metadata_json"}
            missing=required.difference(data.files)
            if missing: raise ValueError(f"Archived KROM dashboard data are incomplete; missing {sorted(missing)}.")
            states=np.asarray(data["training_states"],dtype=float)
            metadata=json.loads(str(data["metadata_json"]))
        analysis=metadata.get("analysis_metadata") or {}; krom=analysis.get("krom")
        if not krom: raise ValueError("Derived archive has no KROM metadata.")
        retained=tuple(int(index) for index in krom.get("retained_indices",()))
        if not retained: raise ValueError("Archived KROM has no retained eigenfunction indices.")
        displayed=_krom_display_indices(request,retained)
        try: values=np.asarray(_sources.derived_array(path,"eigenfunction_values",columns=displayed),dtype=complex)
        except KeyError as error: raise ValueError("Archived KROM has neither stored nor reconstructable eigenfunction values.") from error
        eigenvalue_path=request.source.eigenvalue_source or analysis.get("eigenvalue_source_archive") or path
        try: eigenvalues=self._archived_eigenvalues(Path(eigenvalue_path))
        except FileNotFoundError as error:
            raise FileNotFoundError(f"{error} Supply KoopmanSource(..., eigenvalue_source=<exact spectrum archive>).") from error
        if values.ndim!=2 or states.ndim!=2 or states.shape[0]!=values.shape[0]:
            raise ValueError("Archived KROM states and eigenfunction values must be aligned two-dimensional arrays.")
        if values.shape[1]!=len(displayed) or max(retained)>=eigenvalues.size:
            raise ValueError("Archived KROM eigenvalues, eigenfunctions, and retained indices are inconsistent.")
        source_archive=metadata.get("source_archive")
        if not source_archive: raise ValueError("Archived KROM has no primary trajectory source for its time axis.")
        series,source_metadata=_solution_record(Path(source_archive))
        if series.time.size!=states.shape[0]: raise ValueError("Archived KROM training states do not align with the primary trajectory time axis.")
        display_time,display_states,context=_display_analysis_series(TimeSeries(series.time,states),source_metadata,analysis)
        colors=_krom_colors(request.style,len(displayed))
        return self._krom_dashboard_figures(eigenvalues,displayed,display_states,display_time,values,context,colors,
            request.title,retained_indices=retained)

    @staticmethod
    def _archived_eigenvalues(path: Path, _visited_paths: set[Path] | None=None, _visited_spectrum_ids: set[str] | None=None) -> np.ndarray:
        """Load local eigenvalues or resolve a nonredundant source-spectrum link."""
        target=Path(path).resolve()
        paths=set() if _visited_paths is None else _visited_paths
        identifiers=set() if _visited_spectrum_ids is None else _visited_spectrum_ids
        if target in paths: raise ValueError(f"Cyclic Koopman eigenvalue archive link detected at {target}.")
        paths.add(target)
        with VIEWER_ARCHIVE.open_read(target) as data:
            if "eigenvalues" in data.files: return np.asarray(data["eigenvalues"],dtype=complex)
            metadata=json.loads(str(data["metadata_json"]))
        analysis=metadata.get("analysis_metadata") or {}
        exact_source=analysis.get("eigenvalue_source_archive")
        if exact_source: return Viewer._archived_eigenvalues(Path(exact_source),paths,identifiers)
        source_id=(analysis.get("eigenvalue_storage") or {}).get("source_spectrum_id")
        if not source_id: raise KeyError("Derived Koopman archive contains neither eigenvalues nor an eigenvalue-source link.")
        if source_id in identifiers: raise ValueError(f"Cyclic Koopman eigenvalue spectrum link detected for spectrum '{source_id}'.")
        identifiers.add(source_id)
        for candidate in VIEWER_ARCHIVE.iter_paths("derived"):
            with VIEWER_ARCHIVE.open_read(candidate) as data:
                if "metadata_json" not in data.files: continue
                candidate_metadata=json.loads(str(data["metadata_json"]))
            spectrum=((candidate_metadata.get("analysis_metadata") or {}).get("spectrum") or {})
            if (spectrum.get("metadata") or {}).get("spectrum_id")==source_id:
                return Viewer._archived_eigenvalues(candidate,paths,identifiers)
        raise FileNotFoundError(f"Cannot resolve eigenvalue source spectrum '{source_id}' for {target}.")
