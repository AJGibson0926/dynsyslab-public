"""Archive-only Viewer public API with lazy implementation loading."""
from __future__ import annotations

from importlib import import_module


_LAZY_EXPORTS={
    name:(".contracts",name) for name in (
        "ArchiveSource","BifurcationSource","CampaignSource","KoopmanSource","ControlErrorSource","ObservabilitySource","KLQESource",
        "PODSource","PlotRequest","PlotStyle","ViewerDataSource","TilePlan",
        "ViewerFigure","ViewerReport","Renderer","EigenfunctionView",
        "TimeSeriesLayout","ControlErrorMeasure","ObservableComponent",
    )
}
_LAZY_EXPORTS["Viewer"]=(".core","Viewer")
_LAZY_EXPORTS.update({name:(".settings",name) for name in (
    "InterpolationSamplingSettings","VIEWER_INTERPOLATION_SAMPLING",
    "VIEWER_SYSTEM_SETTINGS","ViewerSystemSettings","get_interpolation_sampling",
    "get_viewer_settings",
)})
_LAZY_EXPORTS.update({name:(".sources",name) for name in (
    "load_archive_source","axis_context","display_series",
)})
_LAZY_EXPORTS.update({
    "classify_source":(".routing","classify"),
    "required_source_types":(".routing","required"),
})


def __getattr__(name: str):
    target=_LAZY_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value=getattr(import_module(target[0],__name__),target[1])
    globals()[name]=value
    return value


def __dir__():
    return sorted(set(globals())|set(_LAZY_EXPORTS))


__all__=["Viewer","ArchiveSource","KoopmanSource","PODSource","ControlErrorSource","ObservabilitySource","KLQESource","BifurcationSource","CampaignSource","ViewerDataSource","PlotRequest","PlotStyle","TilePlan","ViewerFigure","ViewerReport","Renderer","EigenfunctionView","TimeSeriesLayout","ControlErrorMeasure","ObservableComponent","ViewerSystemSettings","InterpolationSamplingSettings","VIEWER_SYSTEM_SETTINGS","VIEWER_INTERPOLATION_SAMPLING","get_viewer_settings","get_interpolation_sampling","load_archive_source","axis_context","display_series","classify_source","required_source_types"]
