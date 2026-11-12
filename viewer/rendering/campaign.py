"""CampaignRender methods extracted from Viewer."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import numpy as np

from archiver.gateway import VIEWER_ARCHIVE
from ..contracts import (_ViewerArtifact, BifurcationSource, CampaignSource,
                         PlotRequest, ViewerFigure)
from .helpers import _sequence

class CampaignRenderMixin:
    """Render existing campaign artifacts and passive campaign status."""

    @staticmethod
    def _campaign_view_state(metadata: dict[str,Any]) -> dict[str,Any]:
        """Read campaign/log status without modifying either primary record."""
        path=metadata.get("campaign_manifest")
        if not path or not Path(path).is_file(): return {"status":"complete_or_unknown"}
        campaign=VIEWER_ARCHIVE.read_json(path); rows=VIEWER_ARCHIVE.coordination_state().get("jobs",{})
        statuses=[rows.get(job_id,{}).get("status","pending") for job_id in campaign.get("jobs",[])]; active=any(item not in {"completed","failed"} for item in statuses); failed=any(item=="failed" for item in statuses)
        reporting=(campaign.get("metadata") or {}).get("progress_reporting",{}); log_path=Path(reporting["log_path"]) if reporting.get("log_path") else None; last=None
        if log_path and log_path.is_file():
            lines=log_path.read_text(encoding="utf-8").splitlines()
            if lines:
                try: last=json.loads(lines[-1])
                except json.JSONDecodeError: pass
        status="incomplete_with_failures" if active and failed else ("incomplete" if active else ("failed" if failed else "complete"))
        return {"status":status,"job_statuses":statuses,"last_campaign_log_event":last}

    def _campaign_trajectory_figure(self, request: PlotRequest) -> ViewerFigure:
        if not isinstance(request.source,_ViewerArtifact): raise RuntimeError("Internal campaign artifact was not materialized.")
        _,metadata=self._viewer_compute_record(request.source); paths=[Path(item) for item in metadata["source_archives"]]
        style=request.style or self.default_style; count=len(paths); colors=_sequence(style.colors,count); widths=_sequence(style.line_width,count); opacity=_sequence(style.opacity,count); traces=[]
        for index,path in enumerate(paths):
            series=VIEWER_ARCHIVE.read_time_series(path); time=np.asarray(series.time); values=np.asarray(series.values); label=f"trajectory-{index:03d}"
            if request.kind=="time":
                for component in request.state_components: traces.append({"label":f"{label} · x{component}","x":time,"y":values[:,component],"color":colors[index],"width":widths[index],"opacity":opacity[index]})
            elif request.kind=="phase":
                if len(request.state_components)!=2: raise ValueError("A phase portrait needs exactly two state components.")
                a,b=request.state_components; traces.append({"label":label,"x":values[:,a],"y":values[:,b],"color":colors[index],"width":widths[index],"opacity":opacity[index]})
            elif request.kind=="three_d":
                if len(request.state_components)!=3: raise ValueError("A 3D plot needs exactly three state components.")
                a,b,c=request.state_components; traces.append({"label":label,"x":values[:,a],"y":values[:,b],"z":values[:,c],"color":colors[index],"width":widths[index],"opacity":opacity[index]})
            else: raise ValueError("Campaign trajectory artifact supports time, phase, or 3D plots.")
        return ViewerFigure(request.kind,request.title or f"Campaign trajectories: {metadata.get('campaign_type','generic')}",traces,{"campaign":metadata,"axes":{"x":"t" if request.kind=="time" else "x0","y":"state" if request.kind=="time" else "x1"}})

    def _campaign_fourier_figure(self, request: PlotRequest) -> ViewerFigure:
        if request.kind!="fourier" or not isinstance(request.source,_ViewerArtifact): raise RuntimeError("Internal campaign Fourier artifact was not materialized.")
        arrays,metadata=self._viewer_compute_record(request.source); frequency=np.asarray(arrays["frequency"]); power=np.asarray(arrays["power"]); components=np.asarray(arrays["components"],dtype=int)
        trace_count=power.shape[1]*power.shape[2]; colors=_sequence(request.style.colors,trace_count); traces=[]; index=0
        for trajectory in range(power.shape[2]):
            for component_index,component in enumerate(components):
                traces.append({"label":f"trajectory-{trajectory:03d} · x{component}","x":frequency,"y":power[:,component_index,trajectory],"color":colors[index],"width":1.0}); index+=1
        return ViewerFigure("fourier",request.title or f"Campaign spectra: {metadata.get('campaign_type','generic')}",traces,{"campaign":metadata,"axes":{"x":"frequency","y":"spectral power"}})

    def _bifurcation_figure(self, request: PlotRequest) -> ViewerFigure:
        raise ValueError("Viewer is display-only and does not derive bifurcation points from archived trajectories.")
        if isinstance(request.source,CampaignSource): return self._campaign_bifurcation_figure(request)
        if not isinstance(request.source,BifurcationSource): raise ValueError("Bifurcation plots require BifurcationSource or CampaignSource.")
        source=request.source; manifest_path=Path(source.manifest_path)
        payload=VIEWER_ARCHIVE.read_json(manifest_path); runs=payload.get("runs",[])
        if not runs: raise ValueError("Bifurcation manifest has no runs.")
        x=[]; y=[]
        for run in runs:
            parameter=float(run[source.parameter_key]); series=VIEWER_ARCHIVE.read_time_series(run["solution_archive"]); start=int((1-source.tail_fraction)*series.sample_count); values=np.asarray(series.values[start:,0])
            maxima=np.flatnonzero((values[1:-1]>=values[:-2])&(values[1:-1]>=values[2:]))+1
            x.extend([parameter]*maxima.size); y.extend(values[maxima].tolist())
        return ViewerFigure("bifurcation",request.title or "Bifurcation diagram",[{"label":"post-transient maxima","x":np.asarray(x),"y":np.asarray(y),"color":request.style.colors,"marker":True}],{"manifest":str(manifest_path),"tail_fraction":source.tail_fraction})

    def _campaign_bifurcation_figure(self, request: PlotRequest) -> ViewerFigure:
        raise ValueError("Viewer is display-only and does not derive campaign bifurcation points from archived trajectories.")
        source=request.source; campaign_path=Path(source.campaign_path); campaign=VIEWER_ARCHIVE.read_json(campaign_path)
        if "bifurcation" not in str(campaign.get("goal","")).lower(): raise ValueError("CampaignSource bifurcation view requires a bifurcation campaign.")
        jobs=VIEWER_ARCHIVE.coordination_state().get("jobs",{})
        archives=[]
        for job_id in campaign.get("jobs",[]):
            row=jobs.get(job_id,{})
            if row.get("status")!="completed": continue
            archives.extend(row.get("solution_archives",[]))
        points=[]
        for path in sorted({Path(item) for item in archives if str(item).endswith(".npz") and "hierarchies" not in str(item)}):
            with VIEWER_ARCHIVE.open_read(path) as data: metadata=json.loads(str(data["metadata_json"]))
            driver=(metadata.get("system_parameters") or {}).get("driver") or {}
            amplitude=driver.get("amplitude")
            if amplitude is None: continue
            series=VIEWER_ARCHIVE.read_time_series(path); values=np.asarray(series.values[:,0],dtype=float); start=max(1,int((1-source.tail_fraction)*values.size))
            maxima=np.flatnonzero((values[start+1:-1]>=values[start:-2])&(values[start+1:-1]>=values[start+2:]))+start+1
            radius=float((metadata.get("system_parameters") or {}).get("equilibrium_radius_m",1.))
            if maxima.size:
                points.extend((float(amplitude)/1e3,float(values[index]/radius)) for index in maxima)
            else:
                # Short or strongly damped records can have no tail maximum.
                # Use tail-aligned Poincare samples rather than inventing a
                # maximum; longer campaigns naturally replace this fallback.
                frequency=float(driver.get("frequency_hz",0.))
                if frequency>0:
                    period=1./frequency; first=np.ceil(series.time[start]/period)*period
                    sample_times=np.arange(first,series.time[-1]+period*1e-10,period)
                    sampled=np.interp(sample_times,series.time,values)
                else: sampled=np.asarray((values[-1],))
                if not sampled.size: sampled=np.asarray((values[-1],))
                points.extend((float(amplitude)/1e3,float(value/radius)) for value in sampled)
        if not points: raise ValueError("This campaign has no completed forced solution archives yet.")
        array=np.asarray(points); return ViewerFigure("bifurcation",request.title or f"Bifurcation diagram: {campaign.get('label','campaign')}",[{"label":"post-transient radius maxima","x":array[:,0],"y":array[:,1],"color":request.style.colors,"marker":True}],{"campaign_manifest":str(campaign_path),"tail_fraction":source.tail_fraction,"axes":{"x":"acoustic pressure amplitude (kPa)","y":"radius maximum (R/R₀)"}})
