"""Viewer adapter for scalar-observable spectrum comparisons."""
from __future__ import annotations

from pathlib import Path
import numpy as np

from .contracts import TilePlan, ViewerFigure, ViewerReport


def observable_spectrum_comparison_report(comparison, *, domain: str="discrete"):
    """Render a typed observable comparison as a normal Viewer report."""
    if isinstance(comparison, (str, Path)):
        from importlib import import_module
        comparison = import_module(
            "postprocessing.koopman.observable_spectra"
        ).load_observable_spectrum_comparison(comparison)
    colors=("#3478c8","#e45756","#54a24b","#f2cf5b","#b279a2","#ff9da6","#9d755d")
    traces=[]
    if domain=="discrete":
        angles=np.linspace(0.,2*np.pi,361); circle=np.exp(1j*angles)
        traces.append({"label":"unit circle","x":circle.real,"y":circle.imag,
                       "color":"#777777","width":1.})
        axes={"x":"Re(μ) (unitless)","y":"Im(μ) (unitless)"}
    elif domain=="continuous":
        axes={"x":"Re(λ)","y":"Im(λ)"}
    else:
        raise ValueError("domain must be 'discrete' or 'continuous'.")
    for index,entry in enumerate(comparison.entries):
        values=entry.values(domain)
        traces.append({"label":entry.name,"x":values.real,"y":values.imag,
                       "color":colors[index%len(colors)],"marker":True,
                       "marker_size":7.})
    figure=ViewerFigure("koopman_spectrum",
        f"H-DMD spectra by scalar observable ({domain} time)",traces,
        {"domain":f"{domain}-time","axes":axes,"square":True,
         "observable_spectrum_comparison":comparison.summary()})
    return ViewerReport([figure],TilePlan("observable-spectrum-comparison",1))


__all__=["observable_spectrum_comparison_report"]
