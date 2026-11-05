"""Build the portable dynsyslab architecture ('bacon') PDF."""
from __future__ import annotations

from pathlib import Path
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, PageBreak, KeepTogether, Preformatted, Image)
from reportlab.graphics.shapes import Drawing, Line, Polygon, Rect, String

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"docs"/"DYNSYSLAB_BACON_REPORT_2026-07-22.pdf"
MERMAID_IMAGE=ROOT/"docs"/"dynsyslab-architecture-mermaid.png"

def header_footer(canvas, document):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#235789")); canvas.line(0.65*inch,0.58*inch,7.85*inch,0.58*inch)
    canvas.setFont("Helvetica",8); canvas.setFillColor(colors.HexColor("#52616b"))
    canvas.drawString(0.65*inch,0.38*inch,"dynsyslab architecture report")
    canvas.drawRightString(7.85*inch,0.38*inch,f"Page {document.page}")
    canvas.restoreState()

def para(text, style): return Paragraph(text,style)

def flow_diagram():
    """Vector rendering of the Mermaid architecture flow for the PDF."""
    drawing=Drawing(520,260)
    navy=colors.HexColor("#173f5f"); blue=colors.HexColor("#235789"); pale=colors.HexColor("#eaf2f8"); green=colors.HexColor("#e8f5e9"); orange=colors.HexColor("#fff3e0")
    def node(x,y,w,h,label,fill=pale):
        drawing.add(Rect(x,y,w,h,rx=7,ry=7,fillColor=fill,strokeColor=blue,strokeWidth=1))
        lines=label.split("\n")
        for index,line in enumerate(lines):
            drawing.add(String(x+w/2,y+h/2+5-(index-(len(lines)-1)/2)*11,line,textAnchor="middle",fontName="Helvetica-Bold" if index==0 else "Helvetica",fontSize=7.6,fillColor=navy))
    def arrow(x1,y1,x2,y2):
        drawing.add(Line(x1,y1,x2,y2,strokeColor=blue,strokeWidth=1.2))
        dx=x2-x1; dy=y2-y1; length=max((dx*dx+dy*dy)**.5,1); ux=dx/length; uy=dy/length; px=-uy; py=ux
        drawing.add(Polygon([x2,y2,x2-7*ux+3*px,y2-7*uy+3*py,x2-7*ux-3*px,y2-7*uy-3*py],fillColor=blue,strokeColor=blue))
    node(8,205,72,30,"dynsyslib\nsystems")
    node(104,205,76,30,"Campaigner\nplans")
    node(204,205,80,30,"JobManager\nSQLite + jobs",orange)
    node(308,205,66,30,"Solver\none solve",green)
    node(398,205,76,30,"Archiver\ngateway")
    node(398,145,76,30,"Archives\narrays + metadata",green)
    node(200,85,84,30,"Postprocessing\nKoopman / POD")
    node(308,85,66,30,"ROM\nKROM/PROM")
    node(8,25,72,30,"Viewer\nread only")
    node(104,25,76,30,"Viewer compute\ntemporary")
    node(204,25,80,30,"Control +\nEstimation")
    node(398,25,76,30,"Tester\nregression")
    arrow(80,220,104,220); arrow(180,220,204,220); arrow(284,220,308,220); arrow(374,220,398,220); arrow(436,205,436,175)
    arrow(398,160,284,100); arrow(284,100,308,100); arrow(308,85,341,55); arrow(284,40,308,205)
    arrow(398,160,80,40); arrow(80,40,104,40); arrow(180,40,398,40)
    drawing.add(String(260,247,"dynsyslab architectural flow",textAnchor="middle",fontName="Helvetica-Bold",fontSize=10,fillColor=navy))
    drawing.add(String(260,4,"Solid arrows show authorized data/control flow; Viewer never writes to the primary archive.",textAnchor="middle",fontName="Helvetica",fontSize=7.2,fillColor=colors.HexColor("#52616b")))
    return drawing

def build():
    styles=getSampleStyleSheet()
    title=ParagraphStyle("Title",parent=styles["Title"],fontName="Helvetica-Bold",fontSize=25,leading=30,textColor=colors.HexColor("#173f5f"),alignment=TA_CENTER,spaceAfter=12)
    subtitle=ParagraphStyle("Subtitle",parent=styles["Normal"],fontSize=11,leading=15,textColor=colors.HexColor("#52616b"),alignment=TA_CENTER,spaceAfter=28)
    h1=ParagraphStyle("H1",parent=styles["Heading1"],fontName="Helvetica-Bold",fontSize=16,leading=20,textColor=colors.HexColor("#173f5f"),spaceBefore=12,spaceAfter=8)
    h2=ParagraphStyle("H2",parent=styles["Heading2"],fontName="Helvetica-Bold",fontSize=11,leading=14,textColor=colors.HexColor("#235789"),spaceBefore=8,spaceAfter=4)
    body=ParagraphStyle("Body",parent=styles["BodyText"],fontSize=9.4,leading=13,spaceAfter=7)
    small=ParagraphStyle("Small",parent=body,fontSize=8.1,leading=10.5)
    mono=ParagraphStyle("Mono",parent=small,fontName="Courier",backColor=colors.HexColor("#f2f5f7"),borderPadding=5,leading=11)
    story=[para("dynsyslab",title),para("Architecture, package boundaries, and computational flow",subtitle)]
    story += [para("Purpose",h1),para("dynsyslab is a local Python suite for defining dynamical systems, solving them, preserving solution data, postprocessing trajectories, scheduling large studies, and displaying existing results. Its design separates scientific equations, numerical execution, archive policy, orchestration, and visualization so that each can evolve without silently changing the others.",body)]
    story += [para("Suite structure",h1)]
    rows=[("Layer","Primary role","Boundary"),
          ("dynsyslib","Equations, presets, Jacobians, drivers","Never solves or stores outputs"),
          ("solver","One requested ODE/PDE/DDE/lattice solve","Writes through Archiver only"),
          ("archives + archiver","Arrays, metadata, indexes, audit, archive gateway","Archiver owns archive policy"),
          ("campaigner + job_manager","Plan campaigns; dispatch/recover jobs","Manifests only; no solution arrays"),
          ("postprocessing + rom","Koopman, DMD, POD, SINDy, reduced models","Derived records retain source provenance"),
          ("control + estimation","Optional live solve plug-ins","No effect unless explicitly configured"),
          ("viewer","Read-only display interface","Never solves or changes primary archive"),
          ("tester","Automated regression runner","Writes test reports only")]
    table=Table([[para(a,small),para(b,small),para(c,small)] for a,b,c in rows],colWidths=[1.35*inch,2.85*inch,2.9*inch],repeatRows=1)
    table.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),colors.HexColor("#173f5f")),("TEXTCOLOR",(0,0),(-1,0),colors.white),("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("GRID",(0,0),(-1,-1),0.35,colors.HexColor("#c9d6df")),("VALIGN",(0,0),(-1,-1),"TOP"),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f5f8fa")]),("LEFTPADDING",(0,0),(-1,-1),6),("RIGHTPADDING",(0,0),(-1,-1),6),("TOPPADDING",(0,0),(-1,-1),5),("BOTTOMPADDING",(0,0),(-1,-1),5)]))
    if not MERMAID_IMAGE.is_file():
        raise FileNotFoundError(f"Render {MERMAID_IMAGE.name} from DYNSYSLAB_ARCHITECTURE.mmd before building the PDF.")
    story += [table,PageBreak(),para("Architectural flow",h1),para("Exact Mermaid rendering of the dynsyslab architecture flow.",body),Image(str(MERMAID_IMAGE),width=5.62*inch,height=7.89*inch)]
    tree="""dynsyslab/
|- dynsyslib/                 equations, presets, drivers
|  |- lib/                    RPE, Duffing, pendula, Lorenz, Rossler, ...
|  |- ode.py                  plug-in interface
|  `- drivers.py
|- solver/                    one numerical solve at a time
|  |- ode.py, pde.py, dde.py, lattice.py
|  |- hierarchy.py, checkpoints.py, interpolation.py
|  `- pde_parallel/           future domain-decomposition pathway
|- archives/                  compressed solutions, hierarchy records, indexes
|- archiver/                  archive gateway, audit, diagnostics, storage
|- campaigner/                campaign and Job manifest designer
|- job_manager/               dispatch, claims, backends, pause/recovery
|  |- state_store.py          SQLite coordination metadata
|  |- archive_inspection.py, campaign_reporting.py
|  `- backends/               local, remote, Slurm/PBS/cloud pathways
|- postprocessing/            Koopman, POD/PCA, LCS, SSM and pipeline
|- rom/                       KROM and POD-Galerkin models
|- control/ and estimation/   optional live solve plug-ins
|- viewer/                    read-only rendering and private artifacts
|- gui/                       desktop interface
|- tester/                    test runner and tester/tests/ regression suite
`- time_series.py, hierarchy.py, app.py, project_paths.py
"""
    story += [PageBreak(),para("Architecture tree",h1),para("This is the source-level organization. Generated build output, caches, viewer artifacts, and test reports are deliberately excluded from the source tree.",body),Preformatted(tree,mono),PageBreak(),para("Data and execution lifecycle",h1),para("The canonical hierarchy is Batch -> System -> Cluster -> Trajectory -> TimeSeries. A trajectory holds an initial condition and, after a solve, a time series. Cluster and System pooling records relationships without duplicating raw child solution arrays. A completed solve stores its numerical data once; hierarchy and derived records preserve pointers and compact metadata.",body),para("Archive safety",h2),para("Archiver is a capability gateway: Viewer has read capability, Solver/Postprocessing have write capability, and JobManager has manage capability. JobManager coordination metadata is stored transactionally in SQLite. Manifests remain readable files and numerical arrays remain compressed file records. Generated build output is excluded from source discovery and wheel packaging.",body),para("Solver behavior",h2),para("The solver remains intentionally single-job. JobManager provides serial/parallel backend selection, dependency handling, claims, recovery planning, pause/resume, and timing estimates. Checkpoints and high-frequency history are optional; the default is lightweight completion data and telemetry.",body),para("Dynamical-system library",h2),para("The library includes RPE, Duffing, van der Pol, nonlinear/double/coupled pendula, Lorenz, Rossler, Chua, Henon-Heiles, Hindmarsh-Rose, Hodgkin-Huxley, FPUT, Kuramoto, delay systems, and statistical/lattice models. Individual ODEs live in dynsyslib/lib and expose a uniform contract including RHS, Jacobian when available, presets, initial state, period metadata, invalid-state checks, and control input interfaces where meaningful.",body)]
    story += [para("Analysis and visualization",h1),para("The postprocessing layer supports DMD, EDMD, HDMD, HAVOK, Koopman-with-SINDy, POD/PCA, delay-embedded POD, and extension layers for LCS and SSM. Koopman results are represented with spectra, eigenvalues, eigenfunctions, and KROM objects. Algorithms that require uniform time spacing must validate it; any resampling for analysis must not silently alter the source archive.",body),para("Viewer supports time, Fourier, phase, 3-D, bifurcation, spectrum, eigenfunction, and POD-mode displays. For a new viewing session it clears only its own temporary artifact archive. The primary archive remains read-only to Viewer at all times.",body)]
    story += [para("Current verification",h1),para("The coordination refactor has been verified with 118 automated tests. A packaged wheel was inspected and contains no build/lib entries. Archive-root imports outside the Archive/Archiver ownership layers were removed. New regression tests belong under tester/tests.",body),para("Recommended extension discipline",h2),para("Add equations in dynsyslib, solve logic in solver, durable I/O through Archiver, data analysis in postprocessing, and graphics in Viewer. Preserve source provenance rather than copying arrays. Add a focused regression test for any repaired failure, then run python -m tester --all before handoff.",body)]
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    SimpleDocTemplate(str(OUTPUT),pagesize=letter,rightMargin=0.65*inch,leftMargin=0.65*inch,topMargin=0.62*inch,bottomMargin=0.76*inch,title="dynsyslab architecture report",author="dynsyslab").build(story,onFirstPage=header_footer,onLaterPages=header_footer)
    print(OUTPUT)

if __name__=="__main__": build()
