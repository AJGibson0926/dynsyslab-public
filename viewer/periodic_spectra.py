"""Presentation-only native Viewer pages for periodic spectral comparisons."""
from __future__ import annotations
import json
from pathlib import Path
import subprocess


def _display_folder(directory):
    from archiver.gateway import VIEWER_ARCHIVE
    folder=Path(directory).resolve()
    if folder.is_relative_to(Path(VIEWER_ARCHIVE.archive_root).resolve()):
        raise ValueError("Viewer graphics must stay outside the primary solution archive")
    return folder


def comparison_payload(comparison, *, title="Periodic observable spectra", x_axis="frequency"):
    """One separate page per observable, family and packet center; no re-fitting.

    Density uses signed symlog. Amplitudes are linear. Values are never clipped,
    renormalized, interpolated or smoothed here. Sorting is presentation only.
    """
    if comparison.get("method")!="periodic-spectral-comparison":
        raise ValueError("Expected a compiled periodic spectral comparison")
    if x_axis not in {"frequency","parameter"}:raise ValueError("Choose frequency or parameter")
    xs=comparison["validation_frequencies" if x_axis=="frequency" else "validation_parameters"]
    if any(x is None for x in xs):raise ValueError("Requested plot coordinate was not supplied")
    families=comparison["validation_families"];n=len(xs);pages=[];figures=[]
    palette=[("#b57b00","dash-dot"),("#0072b2","dash"),("#009e73","dot"),("#e69f00","solid")]
    def trace(name,x,y,color,style="solid",width=1.6):
        return dict(label=name,x=x,y=y,color=color,line_style=style,line_width_points=width)
    for ref in comparison["references"]:
        matches=[r for r in comparison["comparisons"] if r["channel"]==ref["channel"]]
        for family in dict.fromkeys(families):
            ids=sorted((j for j,f in enumerate(families) if f==family),key=lambda j:xs[j]);x=[xs[j] for j in ids]
            for k,center in enumerate(comparison["packet_centers"]):
                page=len(pages);name=f'{ref["name"]} | {family} | angle={center:g}'
                pages.append(name)
                density=[trace("Fourier reference",comparison["angles"],ref["density"],"#c52b24",width=2.1)]
                amplitude=[trace("Fourier reference",x,[ref["packet_rms"][k*n+j] for j in ids],"#c52b24",width=2.1)]
                for i,row in enumerate(matches):
                    color,style=palette[i%len(palette)];label=f'{row["model"]} (rank {row["rank"]})'
                    density.append(trace(label,comparison["angles"],row["density"],color,style))
                    amplitude.append(trace(label,x,[row["packet_rms"][k*n+j] for j in ids],color,style))
                figures.extend([
                    dict(kind="fourier",title=f'{ref["name"]}: spectral measure (training bank)',x_label="Spectral angle (rad/sample)",
                         y_label="Raw spectral density",value_scale="symlog",dashboard_page=page,traces=density),
                    dict(kind="phase",title=f'{ref["name"]}: packet at angle {center:g} | {family}',
                         x_label="Orbit angular frequency (rad/model time)" if x_axis=="frequency" else "Orbit parameter",
                         y_label="RMS packet amplitude",value_scale="linear",dashboard_page=page,traces=amplitude)])
    return dict(schema="dynsyslab.viewer.payload/v1",source="computed-periodic-spectral-comparison",view="phase",
                tile=dict(name=title,columns=2),page_labels=pages,figures=figures)


def write_comparison_dashboard(comparison, directory, *, title="Periodic observable spectra", x_axis="frequency"):
    """Write disposable Viewer payload/catalog in a new caller-owned directory."""
    folder=_display_folder(directory);folder.mkdir(parents=True,exist_ok=True)
    payload=comparison_payload(comparison,title=title,x_axis=x_axis)
    catalog=dict(schema="dynsyslab.viewer.catalog/v1",items=[dict(id="comparison",label=title,kind="trajectory",source="payload.json",views=["phase"])])
    for name in ["payload","catalog"]:
        if (folder/(name+".json")).exists():raise FileExistsError(folder/(name+".json"))
    for name,obj in [("payload",payload),("catalog",catalog)]:
        with (folder/(name+".json")).open("x",encoding="utf-8") as stream:json.dump(obj,stream,allow_nan=False)
    return folder/"payload.json"


def render_comparison_pages(payload, *, executable, pages=None, width=1500, height=650):
    """Render the existing native Viewer pages; wheel zoom/drag remain Viewer-owned."""
    path=Path(payload).resolve();data=json.loads(path.read_text(encoding="utf-8"));folder=_display_folder(path.parent)
    selected=list(range(len(data["page_labels"]))) if pages is None else list(pages)
    if any(not isinstance(p,int) or p<0 or p>=len(data["page_labels"]) for p in selected):raise ValueError("Invalid Viewer page")
    (folder/"archive").mkdir(exist_ok=True);outputs=[]
    for p in selected:
        output=folder/f"page-{p:02d}.png"
        subprocess.run([str(executable),"render","--source",str(path),"--catalog",str(folder/"catalog.json"),
            "--project",str(folder),"--archives",str(folder/"archive"),"--cache",str(folder/f"cache-{p}"),
            "--page",str(p),"--width",str(width),"--height",str(height),"--output",str(output)],check=True,capture_output=True,text=True)
        outputs.append(output)
    return outputs
