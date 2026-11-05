"""HTML, Plotly, and SVG serialization for already-computed figures."""
from __future__ import annotations

import colorsys
from html import escape
import hashlib
import json
from pathlib import Path
import numpy as np

from archiver.gateway import VIEWER_ARCHIVE
from ..contracts import Renderer, ViewerFigure, ViewerReport
from .helpers import _automatic_time_detail_inset

class SerializationMixin:
    """Write renderer-neutral reports as self-contained SVG or Plotly HTML."""

    def write_html(self, report: ViewerReport, path: str | Path, *, renderer: Renderer | None=None,
                   interactive: bool | None=None) -> Path:
        """Write with the selected renderer; ``interactive`` remains a legacy alias."""
        path=Path(path).resolve()
        try:
            path.relative_to(VIEWER_ARCHIVE.archive_root)
        except ValueError:
            pass
        else:
            raise ValueError("Viewer graphics may not be written inside the solution archive.")
        path.parent.mkdir(parents=True,exist_ok=True)
        selected=("plotly" if interactive else "svg") if interactive is not None else (renderer or self.renderer)
        if selected not in {"svg","plotly"}: raise ValueError("renderer must be 'svg' or 'plotly'.")
        if selected=="plotly":
            return _write_interactive_html(report,path)
        figures="".join(f"<div style='grid-column:span {int(item.metadata.get('grid_span',1))}'>{_svg_figure(item)}</div>" for item in report.figures)
        warnings="".join(f"<li>{escape(item)}</li>" for item in report.warnings)
        columns=report.tile_plan.columns
        path.write_text(f"<!doctype html><html><head><meta charset='utf-8'><style>.grid{{display:grid;grid-template-columns:repeat({columns},minmax(0,1fr));gap:1rem}}svg{{width:100%;height:auto}}@media(max-width:700px){{.grid{{grid-template-columns:1fr}}}}</style></head><body><main><h1>dynsyslab Viewer</h1><div class='grid'>{figures}</div>{'<ul>'+warnings+'</ul>' if warnings else ''}</main></body></html>",encoding="utf-8")
        return path

def _write_interactive_html(report: ViewerReport, path: Path) -> Path:
    """Write client-side interactive figures; Plotly is loaded only by the output page."""
    figures=[]; shared_arrays=[]; array_indices={}
    def encoded_array(value):
        array=np.ascontiguousarray(np.asarray(value))
        if array.size<16: return array.tolist()
        digest=hashlib.blake2b(array.view(np.uint8),digest_size=16).hexdigest()
        key=(array.dtype.str,array.shape,digest)
        index=array_indices.get(key)
        if index is None:
            index=len(shared_arrays); array_indices[key]=index; shared_arrays.append(array.tolist())
        return {"__dynsyslab_array__":index}
    for index,figure in enumerate(report.figures):
        traces=[]
        for trace in figure.traces:
            if "x" not in trace or "y" not in trace: continue
            if trace.get("heatmap") and "z" in trace:
                traces.append({"type":"heatmap","x":encoded_array(trace["x"]),"y":encoded_array(trace["y"]),
                    "z":encoded_array(trace["z"]),"name":trace.get("label",""),"colorscale":"Viridis",
                    "colorbar":{"title":trace.get("colorbar_title",figure.metadata.get("colorbar_title","error"))},
                    "hovertemplate":"x=%{x:.5g}<br>y=%{y:.5g}<br>value=%{z:.5g}<extra></extra>"})
                continue
            if "z" in trace:
                values=np.asarray(trace.get("value",np.zeros_like(trace["x"])))
                colors=[f"hsl({int((np.angle(value)%(2*np.pi))*180/np.pi)},82%,{18+62*min(1.,abs(value)/(np.max(np.abs(values)) or 1.)):.0f}%)" for value in values] if trace.get("encoding")=="complex-domain-color" else values.tolist()
                traces.append({"type":"scatter3d","x":encoded_array(trace["x"]),"y":encoded_array(trace["y"]),"z":encoded_array(trace["z"]),"name":trace.get("label",""),"mode":"markers","marker":{"size":trace.get("marker_size",3),"color":colors,"colorscale":"RdBu","showscale":trace.get("encoding")!="complex-domain-color"}})
                continue
            if trace.get("colored_line") and "value" in trace:
                values=np.asarray(trace["value"]); scale=np.max(np.abs(values)) or 1.
                for point in range(len(values)-1):
                    value=(values[point]+values[point+1])/2
                    color=f"hsl({int((np.angle(value)%(2*np.pi))*180/np.pi)},82%,{18+62*min(1.,abs(value)/scale):.0f}%)" if trace.get("encoding")=="complex-domain-color" else float(value)
                    traces.append({"x":[float(trace["x"][point]),float(trace["x"][point+1])],"y":[float(trace["y"][point]),float(trace["y"][point+1])],"name":trace.get("label","") if point==0 else "","mode":"lines","line":{"color":color,"width":trace.get("width",1.5)},"showlegend":point==0})
                continue
            marker={"size":trace.get("marker_size",3)} if trace.get("marker") else None
            if marker is not None and "value" in trace:
                values=np.asarray(trace["value"])
                marker["color"]=[f"hsl({int((np.angle(value)%(2*np.pi))*180/np.pi)},82%,{18+62*min(1.,abs(value)/(np.max(np.abs(values)) or 1.)):.0f}%)" for value in values] if trace.get("encoding")=="complex-domain-color" else values.tolist()
                marker["colorscale"]="RdBu"; marker["showscale"]=trace.get("encoding")!="complex-domain-color"
            plotly_trace={"x":encoded_array(trace["x"]),"y":encoded_array(trace["y"]),"name":trace.get("label",""),"mode":"markers" if trace.get("marker") else "lines","line":{"color":trace.get("color","#3478c8"),"width":trace.get("width",1.5),"dash":trace.get("line_style","solid")},"opacity":trace.get("opacity",1.)}
            # Plotly treats an explicitly null marker as malformed trace data.
            # Line traces do not need marker configuration at all.
            if marker is not None: plotly_trace["marker"]=marker
            traces.append(plotly_trace)
        axes=figure.metadata.get("axes",{}); title_color=figure.metadata.get("title_color")
        layout_title=figure.title if title_color is None else {"text":figure.title,"font":{"color":str(title_color)}}
        layout={"title":layout_title,"xaxis":{"title":axes.get("x","")},"yaxis":{"title":axes.get("y","")},
                "legend":{"orientation":"h","x":0.,"xanchor":"left","y":1.02,"yanchor":"bottom","font":{"size":10}}}
        if axes.get("xlim") is not None: layout["xaxis"]["range"]=axes["xlim"]
        inset=figure.metadata.get("detail_inset")
        if inset:
            for trace_index in inset.get("trace_indices",[]):
                if trace_index>=len(figure.traces): continue
                source_trace=figure.traces[trace_index]
                if "x" not in source_trace or "y" not in source_trace: continue
                traces.append({"x":encoded_array(source_trace["x"]),"y":encoded_array(source_trace["y"]),
                               "name":source_trace.get("label",""),"mode":"lines","showlegend":False,
                               "hoverinfo":"skip","xaxis":"x2","yaxis":"y2",
                               "line":{"color":source_trace.get("color","#3478c8"),"width":source_trace.get("width",1.5),
                                       "dash":source_trace.get("line_style","solid")},"opacity":source_trace.get("opacity",1.)})
            layout["xaxis2"]={"domain":[.59,.97],"anchor":"y2","range":inset["xlim"],"showgrid":True,
                              "zeroline":False,"ticks":"inside","tickfont":{"size":8}}
            layout["yaxis2"]={"domain":[.56,.94],"anchor":"x2","range":inset["ylim"],"showgrid":True,
                              "zeroline":False,"ticks":"inside","tickfont":{"size":8}}
            layout.setdefault("annotations",[]).append({"text":"detail","xref":"paper","yref":"paper","x":.78,"y":.94,
                                                        "showarrow":False,"font":{"size":9},"bgcolor":"rgba(255,255,255,.88)"})
            layout.setdefault("shapes",[]).append({"type":"rect","xref":"paper","yref":"paper","x0":.59,"x1":.97,
                                                   "y0":.56,"y1":.94,"line":{"color":"#666","width":1},
                                                   "fillcolor":"rgba(255,255,255,.88)","layer":"below"})
        if figure.kind=="koopman_spectrum":
            # A complex spectrum is geometric data: unequal axis scales can
            # make the unit circle and modal locations visually misleading.
            layout["yaxis"].update({"scaleanchor":"x","scaleratio":1.})
            layout["legend"]={"orientation":"h","x":0.,"xanchor":"left","y":-0.2,"font":{"size":10}}
        if figure.kind in {"control_error","observable_heatmap","observability"} and any(trace.get("heatmap") for trace in figure.traces):
            # Scalar fields should fill their square display panel even when
            # the sampled state-coordinate spans differ substantially.  Equal
            # numerical units are reserved for genuinely geometric plots such
            # as the complex Koopman spectrum above.
            layout["xaxis"]["automargin"]=True
            layout["yaxis"]["automargin"]=True
            layout["legend"]={"orientation":"h","x":0.,"xanchor":"left","y":-0.2,"font":{"size":10}}
        if any(trace.get("type")=="scatter3d" for trace in traces): layout["scene"]={"xaxis_title":"x₁","yaxis_title":"x₂","zaxis_title":"x₃"}
        if any(trace.get("type")=="scatter3d" for trace in traces):
            # Standard centered, shallow-oblique view for legible 3-D data.
            system_settings=get_viewer_settings(figure.metadata.get("system"))
            layout["scene"]={"xaxis_title":"x1","yaxis_title":"x2","zaxis_title":"x3",
                             "aspectmode":"data",
                             "camera":(system_settings.camera if system_settings and system_settings.camera else {"center":{"x":0.,"y":0.,"z":0.},
                                       "eye":{"x":1.55,"y":-1.65,"z":0.82},"up":{"x":0.,"y":0.,"z":1.}})}
            layout["height"]=360
        figures.append({"id":f"viewer-figure-{index}","data":traces,"layout":layout,
                        "height":int(figure.metadata.get("height",360)),"grid_span":int(figure.metadata.get("grid_span",1)),
                        "square":bool(figure.metadata.get("square",False))})
    payload=json.dumps({"format":"dynsyslab.viewer-plotly/v2","shared_arrays":shared_arrays,"figures":figures},
                       ensure_ascii=False,default=lambda item: item.tolist() if isinstance(item,np.ndarray) else str(item))
    columns=report.tile_plan.columns
    body="".join(f"<section id='{item['id']}' style='{f'width:min(100%,{item["height"]}px);aspect-ratio:1/1;justify-self:center;' if item['square'] else f'height:{item["height"]}px;'}grid-column:span {item['grid_span']}'></section>" for item in figures)
    path.write_text(f"<!doctype html><html><head><meta charset='utf-8'><style>body{{margin:0;padding:1rem}}.grid{{width:100%;display:grid;grid-template-columns:repeat({columns},minmax(0,1fr));gap:1.25rem;align-items:start}}section{{width:100%;min-width:0;min-height:0;overflow:visible}}@media(max-width:700px){{.grid{{grid-template-columns:1fr}}section{{grid-column:span 1!important}}}}</style></head><body><div class='grid'>{body}</div><script src='https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js'></script><script>const bundle={payload};const figures=bundle.figures;const resolveArray=v=>(v&&v.__dynsyslab_array__!==undefined)?bundle.shared_arrays[v.__dynsyslab_array__]:v;figures.forEach(f=>{{f.data.forEach(t=>{{t.x=resolveArray(t.x);t.y=resolveArray(t.y);if(t.z!==undefined)t.z=resolveArray(t.z)}});const sizing=f.square?{{autosize:true}}:{{height:f.height,autosize:true}};Plotly.newPlot(f.id,f.data,{{...f.layout,...sizing,margin:{{l:72,r:28,t:88,b:64}}}},{{responsive:true,scrollZoom:true}})}});</script></body></html>",encoding="utf-8")
    return path

def _jet_rgb(position: float) -> tuple[float,float,float]:
    """Compact MATLAB-like jet interpolation for a normalized scalar."""
    p=float(np.clip(position,0.,1.))
    return (float(np.clip(1.5-abs(4*p-3),0.,1.)),
            float(np.clip(1.5-abs(4*p-2),0.,1.)),
            float(np.clip(1.5-abs(4*p-1),0.,1.)))

def _trace_color(value: complex | float, encoding: str, magnitude: float, scheme: str) -> tuple[float,float,float]:
    """Map scalar or complex eigenfunction values to the requested display colors."""
    if encoding=="complex-domain-color":
        phase=(np.angle(value)+np.pi)/(2*np.pi)
        if scheme=="jet": return _jet_rgb(phase)
        light=.18+.62*min(1.,abs(value)/magnitude)
        return colorsys.hsv_to_rgb(phase,.82,light/.8)
    normalized=float(value)/magnitude
    if scheme=="jet": return _jet_rgb((normalized+1.)/2.)
    return (1.,1.-.75*max(normalized,0),1.-.75*max(normalized,0)) if normalized>=0 else (1.-.75*(-normalized),1.-.75*(-normalized),1.)

def _svg_colorbar(trace: dict[str,Any]) -> str:
    """Return a compact phase or scalar colorbar for an SVG eigenfunction view."""
    if not trace.get("show_colorbar"): return ""
    encoding=trace.get("encoding",""); scheme=trace.get("color_scheme","domain")
    if encoding=="complex-domain-color":
        label="phase"; lower="−π"; upper="π"
        gradient=("#000080 0%;#0000ff 12%;#00bfff 30%;#00ff80 45%;#ffff00 60%;#ff8000 78%;#800000 100%" if scheme=="jet" else "#7a1f9b 0%;#2d5ed2 25%;#22a65a 50%;#e4c82e 75%;#7a1f9b 100%")
    else:
        label="value"; lower="−"; upper="+"
        gradient="#000080 0%;#0000ff 12%;#00bfff 30%;#00ff80 45%;#ffff00 60%;#ff8000 78%;#800000 100%" if scheme=="jet" else "#2855b1 0%;#ffffff 50%;#d54444 100%"
    stops="".join(f"<stop offset='{offset}' stop-color='{color}'/>" for color,offset in (item.rsplit(" ",1) for item in gradient.split(",")))
    return f"<defs><linearGradient id='viewer-colorbar' x1='0' x2='0' y1='1' y2='0'>{stops}</linearGradient></defs><rect x='486' y='20' width='10' height='212' fill='url(#viewer-colorbar)'/><text x='501' y='25' text-anchor='start'>{upper}</text><text x='501' y='232' text-anchor='start'>{lower}</text><text x='510' y='126' transform='rotate(-90 510 126)' text-anchor='middle'>{label}</text>"

def _svg_figure(figure: ViewerFigure) -> str:
    """Minimal dependency-free renderer for 2D line/scatter data; preserves plot specs for richer renderers."""
    if figure.kind in {"control_error","observable_heatmap","observability"} and figure.traces and figure.traces[0].get("heatmap"):
        trace=figure.traces[0]; x=np.asarray(trace["x"],dtype=float); y=np.asarray(trace["y"],dtype=float)
        values=np.asarray(trace["z"],dtype=float); stride=max(1,int(np.ceil(max(values.shape)/64)))
        x=x[::stride]; y=y[::stride]; values=values[::stride,::stride]
        minimum=float(np.nanmin(values)); maximum=float(np.nanmax(values)); span=max(maximum-minimum,np.finfo(float).eps)
        dx=432/max(1,x.size-1); dy=432/max(1,y.size-1); rects=[]
        for row in range(y.size):
            for column in range(x.size):
                value=values[row,column]
                if not np.isfinite(value): continue
                red,green,blue=_jet_rgb((float(value)-minimum)/span)
                rects.append(f"<rect x='{48+max(0,column-.5)*dx:.2f}' y='{20+max(0,y.size-1-row-.5)*dy:.2f}' width='{dx+0.5:.2f}' height='{dy+0.5:.2f}' fill='rgb({int(red*255)},{int(green*255)},{int(blue*255)})'/>")
        xticks=_axis_ticks(float(x.min()),float(x.max()),lambda value:48+(value-x.min())/(x.max()-x.min())*432,
                          axis="x",x_axis_y=452,x_label_y=469,tick_count=4)
        yticks=_axis_ticks(float(y.min()),float(y.max()),lambda value:452-(value-y.min())/(y.max()-y.min())*432,
                          axis="y",y_axis_x=48,y_label_x=39,tick_count=4)
        axes=figure.metadata.get("axes",{}); aria=escape(str(figure.metadata.get("colorbar_title","control approximation error")))
        return f"<section><h2>{escape(figure.title)}</h2><svg viewBox='0 0 512 510' role='img' aria-label='{aria} over phase space'><rect x='48' y='20' width='432' height='432' fill='none' stroke='currentColor'/>{''.join(rects)}{xticks}{yticks}<text x='264' y='500' text-anchor='middle'>{escape(str(axes.get('x','state x1')))}</text><text x='14' y='236' transform='rotate(-90 14 236)' text-anchor='middle'>{escape(str(axes.get('y','state x2')))}</text></svg></section>"
    if figure.kind=="eigenfunction":
        trace=figure.traces[0]; xx=np.asarray(trace["x"]); yy=np.asarray(trace["y"]); values=np.asarray(trace["value"])
        stride=max(1,int(np.ceil(max(xx.shape)/64))); xx=xx[::stride,::stride]; yy=yy[::stride,::stride]; values=values[::stride,::stride]
        encoding=trace.get("encoding","complex-domain-color")
        magnitude=np.abs(values); scale=np.nanmax(magnitude) or 1.; rects=[]
        scalar_limit=float(np.nanmax(np.abs(values))) or 1.
        dx=432/max(1,xx.shape[1]-1); dy=212/max(1,xx.shape[0]-1)
        for row in range(xx.shape[0]-1):
            for column in range(xx.shape[1]-1):
                value=values[row,column]
                if not np.isfinite(value.real) or not np.isfinite(value.imag): continue
                if encoding=="complex-domain-color":
                    hue=(np.angle(value)%(2*np.pi))/(2*np.pi); light=.18+.62*min(1.,abs(value)/scale); red,green,blue=colorsys.hsv_to_rgb(hue,.82,light/.8)
                else:
                    # A symmetric blue-white-red scale makes sign and relative
                    # strength of a scalar component visible at a glance.
                    normalized=float(value)/scalar_limit
                    if normalized>=0: red,green,blue=1.,1.-.75*normalized,1.-.75*normalized
                    else: red,green,blue=1.-.75*(-normalized),1.-.75*(-normalized),1.
                rects.append(f"<rect x='{48+column*dx:.2f}' y='{20+row*dy:.2f}' width='{dx+0.5:.2f}' height='{dy+0.5:.2f}' fill='rgb({int(red*255)},{int(green*255)},{int(blue*255)})'/>")
        xticks=_axis_ticks(float(xx.min()),float(xx.max()),lambda value:48+(value-xx.min())/(xx.max()-xx.min())*432,
                          axis="x",x_axis_y=232,x_label_y=249,tick_count=3)
        yticks=_axis_ticks(float(yy.min()),float(yy.max()),lambda value:232-(value-yy.min())/(yy.max()-yy.min())*212,
                          axis="y",y_axis_x=48,y_label_x=39,tick_count=3)
        aria="Complex-domain coloring" if encoding=="complex-domain-color" else f"{encoding} component"
        return f"<section><h2>{escape(figure.title)}</h2><svg viewBox='0 0 512 280' role='img' aria-label='{aria} of a Koopman eigenfunction'><rect x='48' y='20' width='432' height='212' fill='none' stroke='currentColor'/>{''.join(rects)}{xticks}{yticks}<text x='264' y='274' text-anchor='middle'>state x₁</text><text x='16' y='145' transform='rotate(-90 16 145)' text-anchor='middle'>state x₂</text></svg></section>"
    if figure.kind=="three_d":
        traces=[]
        for trace in figure.traces:
            x=np.asarray(trace["x"],float); y=np.asarray(trace["y"],float); z=np.asarray(trace["z"],float)
            def unit(values):
                lower,upper=values.min(),values.max(); return (values-lower)/(upper-lower) if upper>lower else np.full(values.size,.5)
            # Centered shallow-yaw projection: z stays vertical while a small
            # x/y rotation retains depth and reveals two-lobed attractors.
            system_settings=get_viewer_settings(figure.metadata.get("system"))
            yaw=system_settings.svg_yaw if system_settings else -.18
            depth_weight=system_settings.svg_depth_weight if system_settings else .13
            cosine=np.cos(yaw); sine=np.sin(yaw)
            horizontal=x*cosine-y*sine
            vertical=z+depth_weight*(x*sine+y*cosine)
            h0,h1=float(horizontal.min()),float(horizontal.max())
            v0,v1=float(vertical.min()),float(vertical.max())
            scale=min(412/max(h1-h0,np.finfo(float).eps),192/max(v1-v0,np.finfo(float).eps))
            px=264+(horizontal-(h0+h1)/2)*scale
            py=126-(vertical-(v0+v1)/2)*scale
            points=" ".join(f"{a:.2f},{b:.2f}" for a,b in zip(px,py))
            if "value" not in trace: traces.append(f"<polyline points='{points}' fill='none' stroke='{escape(str(trace.get('color','#3478c8')))}' stroke-width='{trace.get('width',1.5)}'/>")
            else:
                values=np.asarray(trace["value"]); scale=np.max(np.abs(values)) or 1.; circles=[]
                for a,b,value in zip(px,py,values):
                    if trace.get("encoding")=="complex-domain-color": hue=(np.angle(value)%(2*np.pi))/(2*np.pi); r,g,bl=colorsys.hsv_to_rgb(hue,.82,(.18+.62*min(1.,abs(value)/scale))/.8)
                    else:
                        q=float(value)/(np.max(np.abs(values)) or 1.); r,g,bl=(1.,1-.75*max(q,0),1-.75*max(q,0)) if q>=0 else (1-.75*(-q),1-.75*(-q),1.)
                    circles.append(f"<circle cx='{a:.2f}' cy='{b:.2f}' r='1.8' fill='rgb({int(r*255)},{int(g*255)},{int(bl*255)})'/>")
                traces.append("".join(circles))
        return f"<section><h2>{escape(figure.title)}</h2><svg viewBox='0 0 512 280' role='img' aria-label='Projected three-dimensional trajectory'><line x1='48' y1='232' x2='480' y2='232' stroke='currentColor'/><line x1='48' y1='232' x2='48' y2='20' stroke='currentColor'/><line x1='48' y1='232' x2='169' y2='147' stroke='currentColor'/><text x='44' y='248' text-anchor='end'>0</text><text x='480' y='248' text-anchor='middle'>1</text><text x='40' y='236' text-anchor='end'>0</text><text x='40' y='24' text-anchor='end'>1</text>{''.join(traces)}</svg></section>"
    traces=figure.traces; arrays=[(np.asarray(trace["x"],float),np.asarray(trace["y"],float)) for trace in traces if "x" in trace and "y" in trace]
    if not arrays: return f"<section><h2>{escape(figure.title)}</h2></section>"
    requested=figure.metadata.get("axes",{}).get("xlim"); xmin=float(requested[0]) if requested else min(float(x.min()) for x,_ in arrays); xmax=float(requested[1]) if requested else max(float(x.max()) for x,_ in arrays)
    clipped=[(x[(x>=xmin)&(x<=xmax)],y[(x>=xmin)&(x<=xmax)]) for x,y in arrays]
    ymin=min(float(y.min()) for _,y in clipped if y.size); ymax=max(float(y.max()) for _,y in clipped if y.size)
    if xmin==xmax: xmin-=1; xmax+=1
    if ymin==ymax: ymin-=1; ymax+=1
    square=figure.kind=="phase"; plot_bottom=452 if square else 252; plot_height=432 if square else 212; svg_height=480 if square else 280
    sx=lambda v:48+(v-xmin)/(xmax-xmin)*432; sy=lambda v:plot_bottom-(v-ymin)/(ymax-ymin)*plot_height
    paths=[]
    for trace,(x,y) in zip(traces,clipped):
        points=" ".join(f"{sx(a):.2f},{sy(b):.2f}" for a,b in zip(x,y)); color=escape(str(trace.get("color","#3478c8")))
        if trace.get("colored_line") and "value" in trace:
            values=np.asarray(trace["value"]); magnitude=np.max(np.abs(values)) or 1.; segments=[]
            for a,b,c,d,first,second in zip(x[:-1],y[:-1],x[1:],y[1:],values[:-1],values[1:]):
                value=(first+second)/2
                red,green,blue=_trace_color(value,trace.get("encoding",""),magnitude,trace.get("color_scheme","domain"))
                segments.append(f"<line x1='{sx(a):.2f}' y1='{sy(b):.2f}' x2='{sx(c):.2f}' y2='{sy(d):.2f}' stroke='rgb({int(red*255)},{int(green*255)},{int(blue*255)})' stroke-width='{trace.get('width',1.5)}'/>")
            paths.append("".join(segments))
        elif trace.get("marker"):
            # Spectra are point sets: connecting eigenvalues implies a false
            # ordering.  Render one marker for every value instead.
            if "value" in trace:
                values=np.asarray(trace["value"]); magnitude=np.max(np.abs(values)) or 1.; circles=[]
                for a,b,value in zip(x,y,values):
                    if trace.get("encoding")=="complex-domain-color":
                        hue=(np.angle(value)%(2*np.pi))/(2*np.pi); red,green,blue=colorsys.hsv_to_rgb(hue,.82,(.18+.62*min(1.,abs(value)/magnitude))/.8)
                    else:
                        scaled=float(value)/magnitude; red,green,blue=(1.,1-.75*max(scaled,0),1-.75*max(scaled,0)) if scaled>=0 else (1-.75*(-scaled),1-.75*(-scaled),1.)
                    circles.append(f"<circle cx='{sx(a):.2f}' cy='{sy(b):.2f}' r='{trace.get('marker_size',2.5)}' fill='rgb({int(red*255)},{int(green*255)},{int(blue*255)})'/>")
                paths.append("".join(circles))
            else: paths.append("".join(f"<circle cx='{sx(a):.2f}' cy='{sy(b):.2f}' r='{trace.get('marker_size',2.5)}' fill='{color}'/>" for a,b in zip(x,y)))
        else:
            # Deliberately do not decimate archive trajectories; Viewer is an
            # archive renderer and should always expose their stored resolution.
            dash=" stroke-dasharray='8 6'" if trace.get("line_style") in {"dash","dashed"} else ""
            paths.append(f"<polyline points='{points}' fill='none' stroke='{color}' stroke-width='{trace.get('width',1.5)}'{dash}/>")
    axes=figure.metadata.get("axes",{}); xlabel=escape(axes.get("x","")); ylabel=escape(axes.get("y",""))
    xticks=_axis_ticks(xmin,xmax,sx,axis="x",x_axis_y=plot_bottom,x_label_y=svg_height-11); yticks=_axis_ticks(ymin,ymax,sy,axis="y")
    colorbar=_svg_colorbar(traces[0]) if traces else ""
    inset_svg=""; inset=figure.metadata.get("detail_inset")
    if inset:
        ix,iy,iw,ih=286.,34.,178.,102.; x0,x1=map(float,inset["xlim"]); y0,y1=map(float,inset["ylim"])
        inset_id="viewer-detail-"+hashlib.blake2b(figure.title.encode("utf-8"),digest_size=5).hexdigest()
        isx=lambda value:ix+(value-x0)/(x1-x0)*iw; isy=lambda value:iy+ih-(value-y0)/(y1-y0)*ih
        inset_paths=[]
        for trace_index in inset.get("trace_indices",[]):
            if trace_index>=len(traces): continue
            trace=traces[trace_index]; tx=np.asarray(trace.get("x",[]),float); ty=np.asarray(trace.get("y",[]),float)
            selected=np.isfinite(tx)&np.isfinite(ty)&(tx>=x0)&(tx<=x1); tx=tx[selected]; ty=ty[selected]
            if not tx.size: continue
            points=" ".join(f"{isx(a):.2f},{isy(b):.2f}" for a,b in zip(tx,ty))
            color=escape(str(trace.get("color","#3478c8"))); dash=" stroke-dasharray='5 4'" if trace.get("line_style") in {"dash","dashed"} else ""
            inset_paths.append(f"<polyline points='{points}' fill='none' stroke='{color}' stroke-width='{trace.get('width',1.5)}'{dash}/>")
        inset_svg=(f"<g data-viewer-detail-inset='true'><rect x='{ix}' y='{iy}' width='{iw}' height='{ih}' fill='white' fill-opacity='.92' stroke='#666'/>"
                   f"<defs><clipPath id='{inset_id}'><rect x='{ix}' y='{iy}' width='{iw}' height='{ih}'/></clipPath></defs>"
                   f"<g clip-path='url(#{inset_id})'>{''.join(inset_paths)}</g><text x='{ix+iw-4}' y='{iy+11}' text-anchor='end' font-size='9'>detail</text>"
                   f"<text x='{ix+2}' y='{iy+ih-3}' font-size='8'>{x0:.3g}</text><text x='{ix+iw-2}' y='{iy+ih-3}' text-anchor='end' font-size='8'>{x1:.3g}</text>"
                   f"<text x='{ix+2}' y='{iy+10}' font-size='8'>{y1:.3g}</text><text x='{ix+2}' y='{iy+ih-13}' font-size='8'>{y0:.3g}</text></g>")
    vertical_center=(20+plot_bottom)/2
    return f"<section><h2>{escape(figure.title)}</h2><svg viewBox='0 0 530 {svg_height}' role='img'><line x1='48' y1='{plot_bottom}' x2='480' y2='{plot_bottom}' stroke='currentColor'/><line x1='48' y1='20' x2='48' y2='{plot_bottom}' stroke='currentColor'/>{''.join(paths)}{xticks}{yticks}{colorbar}{inset_svg}<text x='264' y='{svg_height-6}' text-anchor='middle'>{xlabel}</text><text x='15' y='{vertical_center}' transform='rotate(-90 15 {vertical_center})' text-anchor='middle'>{ylabel}</text></svg></section>"

def _axis_ticks(lower: float, upper: float, project, *, axis: str, tick_count: int=5,
                x_axis_y: float=252, x_label_y: float=269,
                y_axis_x: float=48, y_label_x: float=39) -> str:
    """Compact numeric ticks with geometry supplied by each SVG renderer."""
    if upper==lower: upper=lower+1.
    values=np.linspace(lower,upper,tick_count); labels=[]
    for value in values:
        number=f"{value:.3g}"; position=float(project(value))
        if axis=="x": labels.append(f"<line x1='{position:.2f}' y1='{x_axis_y:.2f}' x2='{position:.2f}' y2='{x_axis_y+5:.2f}' stroke='currentColor'/><text x='{position:.2f}' y='{x_label_y:.2f}' text-anchor='middle'>{number}</text>")
        else: labels.append(f"<line x1='{y_axis_x-5:.2f}' y1='{position:.2f}' x2='{y_axis_x:.2f}' y2='{position:.2f}' stroke='currentColor'/><text x='{y_label_x:.2f}' y='{position+4:.2f}' text-anchor='end'>{number}</text>")
    return "".join(labels)
