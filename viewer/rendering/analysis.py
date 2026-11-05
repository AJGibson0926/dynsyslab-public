"""AnalysisRender methods extracted from Viewer."""
from __future__ import annotations

from pathlib import Path
import numpy as np

from ..contracts import (ControlErrorSource, KLQESource, ObservableComponent,
                         ObservabilitySource, PlotRequest, ViewerFigure)
from .helpers import _sequence

class AnalysisRenderMixin:
    """Render observability, estimation, and control-error analysis products."""

    def _observability_figures(self, request: PlotRequest) -> list[ViewerFigure]:
        """Render strengths, split-coordinate directions, and mode participation."""
        if not isinstance(request.source,ObservabilitySource):
            raise ValueError("Observability plots require ObservabilitySource.")
        source=request.source.source
        if isinstance(source,(str,Path)):
            from rom import load_observability_analysis
            analysis=load_observability_analysis(source)
        else:
            from rom import ObservabilityAnalysis
            if not isinstance(source,ObservabilityAnalysis):
                raise TypeError("ObservabilitySource must contain an ObservabilityAnalysis or saved artifact.")
            analysis=source
        indices=np.arange(1,analysis.singular_values.size+1)
        floor=max(float(np.max(analysis.singular_values))*np.finfo(float).eps,np.finfo(float).tiny)
        color=request.style.colors if isinstance(request.style.colors,str) else request.style.colors[0]
        figures=[ViewerFigure("observability",request.title or "KROM observability strengths",[
            {"label":"log10 singular strength","x":indices,
             "y":np.log10(np.maximum(analysis.singular_values,floor)),
             "color":color,"marker":True}],
            {"axes":{"x":"observability direction index","y":"log10 singular strength"},
             "observability":analysis.summary()})]
        coordinates=np.arange(analysis.state_dimension)
        figures.append(ViewerFigure("observability","Observability directions (columns)",[
            {"x":indices,"y":coordinates,"z":analysis.directions,"heatmap":True,
             "colorbar_title":"signed split-coordinate coefficient"}],
            {"axes":{"x":"direction index","y":"split lifted coordinate"},
             "colorbar_title":"signed split-coordinate coefficient","directions_are_columns":True,
             "observability":analysis.summary(),"height":500}))
        if analysis.mode_participation is not None:
            groups=np.arange(analysis.mode_participation.shape[0])
            figures.append(ViewerFigure("observability","Koopman-mode participation by observable direction",[
                {"x":indices,"y":groups,"z":analysis.mode_participation,"heatmap":True,
                 "colorbar_title":"fractional squared participation"}],
                {"axes":{"x":"direction index","y":"nonredundant Koopman mode"},
                 "colorbar_title":"fractional squared participation","height":500}))
        return figures

    def _klqe_figures(self, request: PlotRequest) -> list[ViewerFigure]:
        """Render archived state estimates, measured outputs, and innovations."""
        if not isinstance(request.source,KLQESource):
            raise ValueError("KLQE plots require KLQESource.")
        source=request.source.source
        if isinstance(source,(str,Path)):
            from estimation import load_klqe_result
            result=load_klqe_result(source)
        else:
            from estimation import KLQEResult
            if not isinstance(source,KLQEResult):
                raise TypeError("KLQESource must contain a KLQEResult or saved artifact.")
            result=source
        colors=_sequence(request.style.colors,max(result.state_estimate.values.shape[1],
                                                   result.measurement.values.shape[1]))
        time=result.state_estimate.time; figures=[]
        for component in range(result.state_estimate.values.shape[1]):
            traces=[{"label":f"estimated x{component}","x":time,
                     "y":result.state_estimate.values[:,component],"color":colors[component]}]
            if result.true_state is not None:
                traces.append({"label":f"true x{component}","x":time,
                    "y":result.true_state.values[:,component],"color":"#222222","line_style":"dash"})
            figures.append(ViewerFigure("klqe",f"KLQE state component x{component}",traces,
                {"axes":{"x":"time","y":f"x{component}"},"role":"state_estimate",
                 "klqe":result.summary()}))
        for component in range(result.measurement.values.shape[1]):
            figures.append(ViewerFigure("klqe",f"KLQE measurement y{component}",[
                {"label":"measured","x":time,"y":result.measurement.values[:,component],"color":colors[component]},
                {"label":"predicted","x":time,"y":result.measurement_estimate.values[:,component],
                 "color":"#d62728","line_style":"dash"}],
                {"axes":{"x":"time","y":f"y{component}"},"role":"measurement"}))
            figures.append(ViewerFigure("klqe",f"KLQE innovation y{component}",[
                {"label":"innovation","x":time,"y":result.innovation.values[:,component],"color":colors[component]}],
                {"axes":{"x":"time","y":f"y{component} - yhat{component}"},"role":"innovation"}))
        if result.state_estimate.values.shape[1]>=2:
            traces=[{"label":"estimated","x":result.state_estimate.values[:,0],
                     "y":result.state_estimate.values[:,1],"color":colors[0]}]
            if result.true_state is not None:
                traces.append({"label":"true","x":result.true_state.values[:,0],
                    "y":result.true_state.values[:,1],"color":"#222222","line_style":"dash"})
            figures.append(ViewerFigure("phase","KLQE state phase portrait",traces,
                {"axes":{"x":"x0","y":"x1"},"role":"state_phase","square":True,"height":520}))
        if result.state_error is not None:
            norm=np.linalg.norm(result.state_error.values,axis=1)
            figures.append(ViewerFigure("klqe","KLQE state-estimation error norm",[
                {"label":"||x - xhat||2","x":time,"y":norm,"color":"#d62728"}],
                {"axes":{"x":"time","y":"Euclidean state error"},"role":"error_norm"}))
        poles=((result.metadata.get("estimator") or {}).get("observer_poles") or [])
        if poles:
            values=np.asarray([complex(item["real"],item["imag"]) for item in poles])
            figures.append(ViewerFigure("klqe","KLQE observer poles",[
                {"label":"observer poles","x":values.real,"y":values.imag,"marker":True,"color":"#984ea3"}],
                {"axes":{"x":"real part","y":"imaginary part"},"role":"observer_poles"}))
        return figures

    def _control_error_figures(self, request: PlotRequest) -> list[ViewerFigure]:
        """Render typed control-map matrix-norm errors as native heat maps."""
        if not isinstance(request.source,ControlErrorSource):
            raise ValueError("Control-error heat maps require ControlErrorSource.")
        source=request.source.source; archive_metadata={}
        if isinstance(source,(str,Path)):
            from rom import (load_control_approximation_error_fields,
                             load_control_approximation_evaluations)
            try:
                fields,archive_metadata=load_control_approximation_error_fields(source)
            except ValueError as field_error:
                try:
                    evaluations,archive_metadata=load_control_approximation_evaluations(source)
                except ValueError:
                    raise field_error
                fields={item.label:item.error_field for item in evaluations}
        else:
            from rom import ControlApproximationErrorField, ControlApproximationEvaluation
            if isinstance(source,ControlApproximationErrorField): fields={source.approximation:source}
            elif isinstance(source,ControlApproximationEvaluation): fields={source.label:source.error_field}
            elif isinstance(source,(tuple,list)) and source and all(isinstance(item,ControlApproximationEvaluation) for item in source):
                fields={item.label:item.error_field for item in source}
            elif isinstance(source,dict) and source and all(isinstance(item,ControlApproximationErrorField) for item in source.values()): fields=dict(source)
            else: raise TypeError("ControlErrorSource must contain error fields/evaluations or a saved error archive.")
        measure=request.control_error_measure
        if measure not in {"relative","absolute","exact_norm"}:
            raise ValueError("control_error_measure must be 'relative', 'absolute', or 'exact_norm'.")
        attribute={"relative":"relative_error","absolute":"absolute_error","exact_norm":"exact_norm"}[measure]
        label={"relative":"relative Frobenius error","absolute":"absolute Frobenius error","exact_norm":"exact ||B_psi||_F"}[measure]
        figures=[]
        for name,field in fields.items():
            raw=np.asarray(getattr(field,attribute),dtype=float); values=raw
            if request.control_error_log_scale:
                positive=raw[raw>0.]; floor=max(float(np.min(positive)) if positive.size else np.finfo(float).tiny,1e-16)
                values=np.log10(np.maximum(raw,floor)); color_label=f"log10 {label}"
            else: color_label=label
            axes=((archive_metadata.get("metadata") or {}).get("axes") or field.metadata.get("axes") or
                  {"x":"state x1","y":"state x2"})
            figures.append(ViewerFigure("control_error",request.title or f"{name}: {label}",[
                {"x":field.x_values,"y":field.y_values,"z":values,"heatmap":True,"colorbar_title":color_label}],
                {"axes":axes,"colorbar_title":color_label,"measure":measure,
                 "log_scale":bool(request.control_error_log_scale),"approximation":field.approximation,
                 "field_metadata":field.metadata,"archive_metadata":archive_metadata,"height":500,"square":True}))
        return figures
