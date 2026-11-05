"""KoopmanRender methods extracted from Viewer."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence
import numpy as np

from ..contracts import (_ViewerArtifact, KoopmanSource, PlotRequest, PODSource,
                         ViewerFigure)
from ..settings import get_interpolation_sampling
from .helpers import (_display_samples, _interpolation_kind,
                      _krom_colors, _sequence)

class KoopmanRenderMixin:
    """Render Koopman spectra, KROM dashboards, and sampled mode fields."""

    def _koopman_figure(self, request: PlotRequest) -> ViewerFigure:
        if isinstance(request.source,_ViewerArtifact): return self._viewer_compute_figure(request,"koopman_eigenfunction")
        if not isinstance(request.source,KoopmanSource): raise ValueError("Koopman plots require KoopmanSource.")
        source=request.source.source
        if request.kind=="koopman_spectrum":
            if isinstance(source,(str,Path)):
                values=self._archived_eigenvalues(Path(request.source.eigenvalue_source or source))
            else: values=np.asarray([item.discrete_time for item in source.eigenvalues],dtype=complex)
            circle=np.exp(1j*np.linspace(0,2*np.pi,361))
            marker_color=_sequence(request.style.colors,1)[0]
            return ViewerFigure("koopman_spectrum",request.title or "Koopman spectrum",[{"label":"unit circle","x":circle.real,"y":circle.imag,"color":"#777777","width":1},{"label":"eigenvalues","x":values.real,"y":values.imag,"color":marker_color,"marker":True}],{"domain":"discrete-time","eigenvalue_symbol":"mu","axes":{"x":"Re(μ) (unitless)","y":"Im(μ) (unitless)"}})
        if not hasattr(source,"eigenfunctions"): raise ValueError("Eigenfunction display requires an in-memory KoopmanSpectrum.")
        function=source.eigenfunctions[request.eigenfunction_index]
        component=request.eigenfunction_display
        if component not in {"domain_color","real","imaginary"}:
            raise ValueError("eigenfunction_display must be 'domain_color', 'real', or 'imaginary'.")
        system=getattr(source,"dynamical_system",None)
        state_dimension=int(getattr(system,"state_dimension",0) or 0)
        samples=None if request.state_samples is None else np.asarray(request.state_samples,dtype=float)
        if not state_dimension and samples is not None and samples.ndim==2: state_dimension=samples.shape[1]
        if state_dimension==3:
            eigenfunction_view=request.eigenfunction_view or "phase_portrait"
            if eigenfunction_view!="phase_portrait":
                raise ValueError("A full-domain view is not available for 3D Koopman eigenfunctions; use eigenfunction_view='phase_portrait'.")
            if samples is None or samples.ndim!=2 or samples.shape[1]!=3:
                raise ValueError("A 3D eigenfunction display requires source-data state_samples with shape (sample_count, 3).")
            if function.values is None: raise ValueError("Viewer will not evaluate ψ(x). This eigenfunction needs archived sampled values for display.")
            values=np.asarray(function.values,dtype=complex).reshape(-1)
            if values.size!=samples.shape[0]: raise ValueError("3D source states and eigenfunction values must have matching sample counts.")
            display_values=values if component=="domain_color" else (values.real if component=="real" else values.imag)
            return ViewerFigure("three_d",request.title or f"Koopman eigenfunction {request.eigenfunction_index}",[{
                "label":"eigenfunction source trajectories","x":samples[:,0],"y":samples[:,1],"z":samples[:,2],"value":display_values,
                "encoding":"complex-domain-color" if component=="domain_color" else component,"marker":True}],{"system":getattr(system,"name",None),"state_dimension":3,"representation":function.representation,"component":component,"eigenfunction_view":eigenfunction_view,"colored_by":"koopman_eigenfunction"})
        if state_dimension!=2: raise ValueError(f"Eigenfunction display supports associated 2D or 3D systems, not state dimension {state_dimension or 'unknown'}.")
        eigenfunction_view=request.eigenfunction_view or "domain"
        if eigenfunction_view not in {"domain","phase_portrait"}:
            raise ValueError("eigenfunction_view must be 'domain' or 'phase_portrait'.")
        if eigenfunction_view=="phase_portrait":
            if samples is None or samples.ndim!=2 or samples.shape[1]!=2:
                raise ValueError("A Koopman eigenfunction phase portrait requires 2D source-data state_samples.")
            if function.values is None: raise ValueError("Viewer will not evaluate ψ(x). This eigenfunction needs archived sampled values for display.")
            values=np.asarray(function.values,dtype=complex).reshape(-1)
            if values.size!=samples.shape[0]: raise ValueError("Archived eigenfunction values must match supplied phase-portrait samples.")
            display_values=values if component=="domain_color" else (values.real if component=="real" else values.imag)
            return ViewerFigure("phase",request.title or f"Koopman eigenfunction phase portrait {request.eigenfunction_index}",[{"label":"eigenfunction phase portrait","x":samples[:,0],"y":samples[:,1],"value":display_values,"encoding":"complex-domain-color" if component=="domain_color" else component,"colored_line":True,"color_scheme":request.color_scheme,"show_colorbar":request.show_colorbar}],{"system":getattr(system,"name",None),"state_dimension":2,"representation":function.representation,"component":component,"eigenfunction_view":"phase_portrait","axes":{"x":"x1","y":"x2"}})
        raise ValueError("Viewer does not evaluate or interpolate eigenfunctions across a domain. Archive a domain-grid result in postprocessing first.")
        if request.domain is None: raise ValueError("A 2D eigenfunction display requires a domain.")
        (x0,x1),(y0,y1)=request.domain; density=request.domain_density
        if density<2: raise ValueError("domain_density must be at least 2.")
        xx,yy=np.meshgrid(np.linspace(x0,x1,density),np.linspace(y0,y1,density)); states=np.column_stack((xx.ravel(),yy.ravel()))
        if function.has_state_callable:
            values=np.asarray([function.psi(row) for row in states]).reshape(xx.shape)
        else:
            if request.state_samples is None or function.values is None:
                raise ValueError("This trajectory-valued eigenfunction needs matching 2D state_samples to interpolate psi(x) over a domain.")
            samples=np.asarray(request.state_samples,dtype=float); samples_values=np.asarray(function.values,dtype=complex).reshape(-1)
            if samples.ndim!=2 or samples.shape[1]!=2 or samples.shape[0]!=samples_values.size:
                raise ValueError("state_samples must have shape (sample_count, 2) and match the eigenfunction value count.")
            from scipy.interpolate import LinearNDInterpolator
            real=LinearNDInterpolator(samples,samples_values.real,fill_value=np.nan)(states)
            imaginary=LinearNDInterpolator(samples,samples_values.imag,fill_value=np.nan)(states)
            values=(real+1j*imaginary).reshape(xx.shape)
        display_values=values if component=="domain_color" else (values.real if component=="real" else values.imag)
        interpolation_kind=_interpolation_kind(getattr(source,"metadata",None),"hdmd")
        sampling=get_interpolation_sampling(interpolation_kind)
        return ViewerFigure("eigenfunction",request.title or f"Koopman eigenfunction {request.eigenfunction_index}",[{"x":xx,"y":yy,"value":display_values,"encoding":"complex-domain-color" if component=="domain_color" else component}],{"state_dimension":2,"representation":function.representation,"component":component,"interpolation_sampling":{"kind":sampling.kind,"recommended_domain_density":sampling.recommended_domain_density,"requested_domain_density":density,"notes":sampling.notes}})

    @staticmethod
    def _krom_dashboard_figures(eigenvalues: np.ndarray, retained: tuple[int,...], states: np.ndarray,
                                time: np.ndarray, values: np.ndarray, context: dict[str,Any],
                                colors: Sequence[str], title: str | None=None, *,
                                retained_indices: tuple[int,...] | None=None) -> list[ViewerFigure]:
        """Build the canonical one-spectrum plus five-panels-per-mode KROM view."""
        eigenvalues=np.asarray(eigenvalues,dtype=complex).reshape(-1); states=np.asarray(states,dtype=float)
        time=np.asarray(time,dtype=float).reshape(-1); values=np.asarray(values,dtype=complex)
        circle=np.exp(1j*np.linspace(0,2*np.pi,361)); selected=list(retained)
        model_retained=tuple(retained if retained_indices is None else retained_indices)
        if not model_retained:
            raise ValueError("A KROM dashboard requires at least one retained eigenvalue.")
        if any(index<0 or index>=eigenvalues.size for index in model_retained):
            raise ValueError("KROM retained indices exceed the eigenvalue spectrum.")
        if any(index not in model_retained for index in selected):
            raise ValueError("Displayed KROM modes must be members of the retained model spectrum.")
        if values.ndim!=2 or values.shape[0]!=states.shape[0]:
            raise ValueError("KROM eigenfunction values must align with the displayed states.")
        if values.shape[1]==eigenvalues.size: value_columns={index:index for index in selected}
        elif values.shape[1]==len(selected): value_columns={index:position for position,index in enumerate(selected)}
        else: raise ValueError("KROM eigenfunction values must contain either the full spectrum or the retained subset.")
        figures=[ViewerFigure("koopman_spectrum",title or "KROM discrete spectrum (unitless)",[
            {"label":"unit circle","x":circle.real,"y":circle.imag,"color":"#777777","width":1.},
            {"label":"full spectrum","x":eigenvalues.real,"y":eigenvalues.imag,"color":"#b7bdc8","marker":True,"marker_size":5.},
            {"label":"KROM retained spectrum","x":eigenvalues[list(model_retained)].real,"y":eigenvalues[list(model_retained)].imag,
             "color":"#596273","marker":True,"marker_size":6.},
            *[{"label":f"ψ{index}","x":np.asarray([eigenvalues[index].real]),"y":np.asarray([eigenvalues[index].imag]),
               "color":str(color),"marker":True,"marker_size":10.} for index,color in zip(selected,colors)]
        ],{"domain":"discrete-time","axes":{"x":"Re(μ) (unitless)","y":"Im(μ) (unitless)"},"height":500,"grid_span":4,
           "full_spectrum_count":int(eigenvalues.size),"retained_mode_count":len(model_retained),
           "displayed_mode_indices":selected,"omitted_dashboard_mode_count":len(model_retained)-len(selected),
           "display_selection_is_model_truncation":False})]
        state_labels=context.get("state_labels") or ["x0","x1"]
        if states.shape[1]!=2: raise ValueError("The canonical KROM phase dashboard currently requires a two-dimensional physical state.")
        for index,color in zip(selected,colors):
            value_column=value_columns[index]
            eigenvalue=eigenvalues[index]; title_color=str(color)
            mode_title=f"ψ{index} — μ = {eigenvalue.real:.4f} {eigenvalue.imag:+.4f}i (unitless)"
            figures.append(ViewerFigure("phase",mode_title,[{"label":f"ψ{index}: hue = phase, brightness = magnitude",
                "x":states[:,0],"y":states[:,1],"value":values[:,value_column],"encoding":"complex-domain-color","marker":True,"marker_size":3.4}],
                {"title_color":title_color,"component":"domain_color","axes":{"x":state_labels[0],"y":state_labels[1]},"height":430,"grid_span":4}))
            components=(("Phase",np.angle(values[:,value_column]),"arg(ψ) (rad)"),("Magnitude",np.abs(values[:,value_column]),"|ψ|"),
                        ("Real part",values[:,value_column].real,"Re(ψ)"),("Imaginary part",values[:,value_column].imag,"Im(ψ)"))
            for component,component_values,ylabel in components:
                figures.append(ViewerFigure("time",f"ψ{index} · {component}",[{"label":f"ψ{index} {component.lower()}",
                    "x":time,"y":component_values,"color":title_color,"width":1.7}],
                    {"title_color":title_color,"axes":{"x":context.get("time_label","t"),"y":ylabel},"height":280,"grid_span":2}))
        return figures

    def _pod_mode_figure(self, request: PlotRequest) -> ViewerFigure:
        raise ValueError("Viewer will not evaluate or interpolate POD modes. Archive sampled POD display data in postprocessing first.")
        if not isinstance(request.source,PODSource): raise ValueError("POD-mode plots require PODSource.")
        prom=getattr(request.source.source,"prom",request.source.source)
        if not hasattr(prom,"modes"): raise ValueError("PODSource must contain a PROM or an SVD result with .prom.")
        mode=prom.modes[request.pod_mode_index]; component=request.eigenfunction_display
        if component not in {"domain_color","real","imaginary"}: raise ValueError("eigenfunction_display must be 'domain_color', 'real', or 'imaginary'.")
        samples=None if request.state_samples is None else np.asarray(request.state_samples,dtype=float)
        if samples is None or samples.ndim!=2: raise ValueError("A POD-mode display requires source-data state_samples with shape (sample_count, state_dimension).")
        interpolation_kind=_interpolation_kind(prom.metadata,"pod")
        samples,sampling_metadata=_display_samples(samples,interpolation_kind,request.max_display_samples)
        values=np.asarray([mode.phi(row) for row in samples],dtype=float)
        dimension=samples.shape[1]
        if dimension==3:
            displayed=values.astype(complex) if component=="domain_color" else (values if component=="real" else np.zeros_like(values))
            return ViewerFigure("three_d",request.title or f"POD mode {request.pod_mode_index}",[{"label":"POD-mode source trajectory","x":samples[:,0],"y":samples[:,1],"z":samples[:,2],"value":displayed,"encoding":"complex-domain-color" if component=="domain_color" else component,"marker":True}],{"system":prom.system.system_name,"state_dimension":3,"component":component,"colored_by":"pod_mode","singular_value":mode.singular_value,"interpolation_sampling":sampling_metadata})
        if dimension!=2: raise ValueError(f"POD-mode display supports 2D or 3D source states, not dimension {dimension}.")
        if request.domain is None: raise ValueError("A 2D POD-mode display requires a domain.")
        (x0,x1),(y0,y1)=request.domain; density=request.domain_density; xx,yy=np.meshgrid(np.linspace(x0,x1,density),np.linspace(y0,y1,density)); states=np.column_stack((xx.ravel(),yy.ravel()))
        scalar=np.asarray([mode.phi(row) for row in states]).reshape(xx.shape)
        displayed=scalar.astype(complex) if component=="domain_color" else (scalar if component=="real" else np.zeros_like(scalar))
        sampling=get_interpolation_sampling(interpolation_kind)
        return ViewerFigure("eigenfunction",request.title or f"POD mode {request.pod_mode_index}",[{"x":xx,"y":yy,"value":displayed,"encoding":"complex-domain-color" if component=="domain_color" else component}],{"system":prom.system.system_name,"state_dimension":2,"component":component,"representation":"POD trajectory-amplitude interpolation","singular_value":mode.singular_value,"interpolation_sampling":{"kind":sampling.kind,"recommended_domain_density":sampling.recommended_domain_density,"requested_domain_density":density,"notes":sampling.notes}})

    def _viewer_compute_figure(self, request: PlotRequest, expected_kind: str) -> ViewerFigure:
        if not isinstance(request.source,_ViewerArtifact): raise RuntimeError("Internal Viewer artifact was not materialized.")
        arrays,metadata=self._viewer_compute_record(request.source)
        if metadata.get("viewing_session_id")!=self.viewer_compute_session_id:
            raise ValueError("ViewerCompute artifact belongs to another viewing session; regenerate it after this Viewer was initialized.")
        if metadata.get("kind")!=expected_kind: raise RuntimeError(f"Internal Viewer artifact contains {metadata.get('kind')!r}, not {expected_kind!r}.")
        if expected_kind=="fourier":
            frequency=np.asarray(arrays["frequency"],dtype=float); power=np.asarray(arrays["power"],dtype=float); components=np.asarray(arrays.get("components",np.arange(power.shape[1])),dtype=int)
            colors=_sequence(request.style.colors,power.shape[1]); widths=_sequence(request.style.line_width,power.shape[1])
            traces=[{"label":f"x{component}","x":frequency,"y":power[:,index],"color":colors[index],"width":widths[index]} for index,component in enumerate(components)]
            return ViewerFigure("fourier",request.title or "Archived Fourier spectrum",traces,{"axes":{"x":"frequency","y":"spectral power"},"viewer_compute":metadata})
        if expected_kind=="bifurcation":
            campaign_state=self._campaign_view_state(metadata)
            return ViewerFigure("bifurcation",request.title or "Archived bifurcation data",[{"label":"post-transient maxima","x":np.asarray(arrays["pressure_kpa"]),"y":np.asarray(arrays["radius_ratio"]),"color":request.style.colors,"marker":True}],{"axes":{"x":"acoustic pressure amplitude (kPa)","y":"radius maximum (R/R₀)"},"viewer_compute":metadata,"campaign_run_state":campaign_state})
        states=np.asarray(arrays["states"],dtype=float); values=np.asarray(arrays["values"])
        if expected_kind=="pod_mode": values=np.asarray(values,dtype=float)
        else: values=np.asarray(values,dtype=complex)
        component=request.eigenfunction_display
        displayed=values if component=="domain_color" else (values.real if component=="real" else values.imag)
        if metadata.get("view")=="domain":
            (x0,x1),(y0,y1)=metadata["domain"]; density=int(metadata["domain_density"]); xx,yy=np.meshgrid(np.linspace(x0,x1,density),np.linspace(y0,y1,density))
            return ViewerFigure("eigenfunction",request.title or "Archived Koopman eigenfunction",[{"x":xx,"y":yy,"value":displayed.reshape(xx.shape),"encoding":"complex-domain-color" if component=="domain_color" else component}],{"viewer_compute":metadata,"state_dimension":2,"component":component})
        if states.shape[1]==3:
            return ViewerFigure("three_d",request.title or "Archived Koopman eigenfunction",[{"x":states[:,0],"y":states[:,1],"z":states[:,2],"value":displayed,"encoding":"complex-domain-color" if component=="domain_color" else component,"marker":True}],{"viewer_compute":metadata,"state_dimension":3,"component":component})
        return ViewerFigure("phase",request.title or "Archived Koopman eigenfunction",[{"x":states[:,0],"y":states[:,1],"value":displayed,"encoding":"complex-domain-color" if component=="domain_color" else component,"colored_line":True,"color_scheme":request.color_scheme,"show_colorbar":request.show_colorbar}],{"viewer_compute":metadata,"state_dimension":2,"component":component})
