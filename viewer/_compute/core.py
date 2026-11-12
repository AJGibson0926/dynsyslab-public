"""Derive disposable Viewer artifacts without ever changing the primary archive.

``_ViewerCompute`` may read solution and postprocessing records from the
primary archive, but every output is written as a self-contained compressed
record below ``viewer_compute_archive``.  It is private implementation detail
of Viewer and is never part of Viewer's public API.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np

from archiver.gateway import VIEWER_ARCHIVE

VIEWER_COMPUTE_ARCHIVE_ROOT = VIEWER_ARCHIVE.archive_root.parent / "viewer_compute_archive"

def _stamp() -> str: return datetime.now(timezone.utc).isoformat()
def _safe(value: str) -> str: return "".join(item if item.isalnum() or item in "-_" else "-" for item in value)

class ViewerCompute:
    """Private, archive-safe producer of temporary Viewer plot records."""
    SESSION_FILE = ".active_viewing_session.json"

    def __init__(self, root: str | Path | None=None, *, session_id: str | None=None):
        self.root=Path(root or VIEWER_COMPUTE_ARCHIVE_ROOT).resolve()
        self.root.mkdir(parents=True,exist_ok=True)
        self.session_id=session_id or self.active_session(self.root)

    @classmethod
    def active_session(cls, root: str | Path | None=None) -> str | None:
        target=Path(root or VIEWER_COMPUTE_ARCHIVE_ROOT).resolve()/cls.SESSION_FILE
        if not target.is_file(): return None
        return json.loads(target.read_text(encoding="utf-8")).get("session_id")

    @classmethod
    def start_viewing_session(cls, root: str | Path | None=None) -> str:
        """Clear only temporary Viewer artifacts and return a fresh session id.

        This method never traverses or changes ``archives`` (the primary
        solution archive).  It is intentionally the sole cleanup route used
        at the beginning of a Viewer session.
        """
        target=Path(root or VIEWER_COMPUTE_ARCHIVE_ROOT).resolve(); target.mkdir(parents=True,exist_ok=True)
        for child in target.iterdir():
            if child.is_file(): child.unlink()
            elif child.is_dir():
                import shutil
                shutil.rmtree(child)
        session_id=uuid4().hex
        (target/cls.SESSION_FILE).write_text(json.dumps({"session_id":session_id,"created_utc":_stamp()},sort_keys=True),encoding="utf-8")
        return session_id

    def _save(self, kind: str, arrays: dict[str,np.ndarray], metadata: dict[str,Any], label: str) -> Path:
        target=self.root/f"{_safe(label)}__{uuid4().hex}.npz"
        payload={name:np.asarray(value) for name,value in arrays.items()}
        active=self.active_session(self.root)
        if self.session_id is None or self.session_id!=active:
            raise RuntimeError("ViewerCompute has no active viewing session. Initialize Viewer, then construct ViewerCompute(session_id=viewer.viewer_compute_session_id).")
        payload["metadata_json"]=np.asarray(json.dumps({"format":"dynsyslab.viewer-compute/v1","kind":kind,"created_utc":_stamp(),"viewing_session_id":self.session_id,**metadata},sort_keys=True))
        np.savez_compressed(target,**payload)
        return target

    def register_campaign_view(self, campaign_path: str | Path) -> Path:
        """Record a disposable Viewer subscription outside the primary archive."""
        if self.session_id is None or self.session_id!=self.active_session(self.root): raise RuntimeError("ViewerCompute has no active viewing session.")
        target=self.root/f"campaign-view-{self.session_id}.json"
        target.write_text(json.dumps({"viewing_session_id":self.session_id,"campaign_manifest":str(Path(campaign_path).resolve()),"registered_utc":_stamp()},sort_keys=True),encoding="utf-8")
        return target

    @staticmethod
    def _primary_path(path: str | Path) -> Path:
        target=Path(path).expanduser().resolve()
        if not target.is_file(): raise FileNotFoundError(target)
        # The condition documents and enforces the data flow: primary files
        # are opened only by read-only readers; no output shares their tree.
        try:
            target.relative_to(VIEWER_COMPUTE_ARCHIVE_ROOT.resolve())
        except ValueError:
            pass
        else:
            raise ValueError("ViewerCompute input must be a primary/archive source, not a viewer-compute artifact.")
        return target

    def fourier(self, solution_archive: str | Path, *, components: tuple[int,...] | None=None,
                label: str="fourier") -> Path:
        """Compute a power spectrum from a uniform archived trajectory."""
        source=self._primary_path(solution_archive)
        try: series=VIEWER_ARCHIVE.read_time_series(source)
        except (KeyError,ValueError) as error: raise ValueError(f"{source} is not a primary solution time-series archive.") from error
        time=np.asarray(series.time,dtype=float); values=np.asarray(series.values,dtype=float)
        spacing=np.diff(time)
        if time.size<2 or not np.allclose(spacing,spacing[0],rtol=1e-9,atol=max(abs(spacing[0])*1e-12,1e-15)):
            raise ValueError("Fourier computation requires a uniformly sampled archived time series; resample in postprocessing first.")
        selected=tuple(range(values.shape[1])) if components is None else tuple(int(item) for item in components)
        if not selected or min(selected)<0 or max(selected)>=values.shape[1]: raise ValueError("components are outside the archived state dimension.")
        selected_values=values[:,selected]; coefficients=np.fft.rfft(selected_values,axis=0)
        return self._save("fourier",{"frequency":np.fft.rfftfreq(time.size,d=float(spacing[0])),"power":np.abs(coefficients)**2/time.size,"components":np.asarray(selected,dtype=int)},
                          {"source_archive":str(source),"uniform_time_steps":True,"delta_t":float(spacing[0]),"component_count":len(selected)},label)

    def koopman_eigenfunction(self, spectrum: Any, eigenfunction_index: int, *, state_samples: np.ndarray,
                              view: str="phase_portrait", domain: tuple[tuple[float,float],tuple[float,float]] | None=None,
                              density: int=100, label: str="koopman-eigenfunction") -> Path:
        """Evaluate one eigenfunction and store plot-ready source or grid data."""
        function=spectrum.eigenfunctions[int(eigenfunction_index)]
        samples=np.asarray(state_samples,dtype=float)
        if samples.ndim!=2: raise ValueError("state_samples must be two-dimensional.")
        if view=="phase_portrait":
            states=samples; shape=None
        elif view=="domain":
            if samples.shape[1]!=2 or domain is None or density<2: raise ValueError("A 2D domain evaluation needs a domain and density >= 2.")
            (x0,x1),(y0,y1)=domain; xx,yy=np.meshgrid(np.linspace(x0,x1,density),np.linspace(y0,y1,density)); states=np.column_stack((xx.ravel(),yy.ravel())); shape=xx.shape
        else: raise ValueError("view must be 'phase_portrait' or 'domain'.")
        values=np.asarray([function.psi(state) for state in states],dtype=complex)
        arrays={"values":values}
        if shape is None: arrays["states"]=states
        system=getattr(spectrum,"dynamical_system",None)
        metadata={"eigenfunction_index":int(eigenfunction_index),"view":view,"state_dimension":int(states.shape[1]),"system_name":getattr(system,"name",None),"representation":getattr(function,"representation",None)}
        if shape is not None: metadata.update({"domain":[[float(x0),float(x1)],[float(y0),float(y1)]],"domain_density":int(density)})
        return self._save("koopman_eigenfunction",arrays,metadata,label)

    def pod_mode(self, prom: Any, mode_index: int, *, state_samples: np.ndarray,
                 view: str="phase_portrait", domain: tuple[tuple[float,float],tuple[float,float]] | None=None,
                 density: int=100, label: str="pod-mode") -> Path:
        """Evaluate one POD mode and store the resulting display values."""
        mode=prom.modes[int(mode_index)]; samples=np.asarray(state_samples,dtype=float)
        if samples.ndim!=2: raise ValueError("state_samples must be two-dimensional.")
        if view=="phase_portrait": states=samples; shape=None
        elif view=="domain":
            if samples.shape[1]!=2 or domain is None or density<2: raise ValueError("A 2D domain evaluation needs a domain and density >= 2.")
            (x0,x1),(y0,y1)=domain; xx,yy=np.meshgrid(np.linspace(x0,x1,density),np.linspace(y0,y1,density)); states=np.column_stack((xx.ravel(),yy.ravel())); shape=xx.shape
        else: raise ValueError("view must be 'phase_portrait' or 'domain'.")
        values=np.asarray([mode.phi(state) for state in states],dtype=float); arrays={"values":values}
        if shape is None: arrays["states"]=states
        metadata={"mode_index":int(mode_index),"view":view,"state_dimension":int(states.shape[1]),"singular_value":float(getattr(mode,"singular_value",np.nan))}
        if shape is not None: metadata.update({"domain":[[float(x0),float(x1)],[float(y0),float(y1)]],"domain_density":int(density)})
        return self._save("pod_mode",arrays,metadata,label)

    def bifurcation_from_campaign(self, campaign_path: str | Path, *, tail_fraction: float=.25,
                                  label: str="bifurcation") -> Path:
        """Extract post-transient radius maxima from completed campaign records."""
        path=Path(campaign_path).expanduser().resolve(); campaign,sources_paths=self._campaign_solution_paths(path); points=[]; sources=[]
        for source in sources_paths:
            with np.load(source,allow_pickle=False) as data: metadata=json.loads(str(data["metadata_json"]))
            driver=(metadata.get("system_parameters") or {}).get("driver") or {}; amplitude=driver.get("amplitude")
            if amplitude is None: continue
            series=VIEWER_ARCHIVE.read_time_series(source); values=np.asarray(series.values[:,0],dtype=float); start=max(1,int((1-float(tail_fraction))*values.size))
            maxima=np.flatnonzero((values[start+1:-1]>=values[start:-2])&(values[start+1:-1]>=values[start+2:]))+start+1
            radius=float((metadata.get("system_parameters") or {}).get("equilibrium_radius_m",1.))
            points.extend((float(amplitude)/1e3,float(values[index]/radius)) for index in maxima); sources.append(str(source))
        point_array=np.asarray(points,dtype=float).reshape((-1,2)) if points else np.empty((0,2),dtype=float)
        return self._save("bifurcation",{"pressure_kpa":point_array[:,0],"radius_ratio":point_array[:,1]},
                          {"campaign_manifest":str(path),"source_archives":sources,"tail_fraction":float(tail_fraction),"point_count":int(point_array.shape[0])},label)

    @staticmethod
    def _campaign_solution_paths(campaign_path: str | Path) -> tuple[dict[str,Any],list[Path]]:
        path=Path(campaign_path).expanduser().resolve(); campaign=VIEWER_ARCHIVE.read_json(path)
        indexed=VIEWER_ARCHIVE.campaign_paths(str(campaign["campaign_id"]))
        if indexed:
            return campaign,sorted({item for item in indexed if item.is_file() and "hierarchies" not in str(item)})
        jobs=VIEWER_ARCHIVE.coordination_state().get("jobs",{}); paths=[]
        for job_id in campaign.get("jobs",[]):
            row=jobs.get(job_id,{})
            if row.get("status")=="completed": paths.extend(Path(item) for item in row.get("solution_archives",[]) if str(item).endswith(".npz") and "hierarchies" not in str(item))
        return campaign,sorted({path for path in paths if path.is_file()})

    def campaign_trajectories(self, campaign_path: str | Path, *, label: str="campaign-trajectories") -> Path:
        """Create a pointer-only display index for completed campaign solutions."""
        campaign,paths=self._campaign_solution_paths(campaign_path)
        if not paths: raise ValueError("The campaign has no completed solution records.")
        # The primary archive remains the sole owner of trajectory arrays.
        return self._save("campaign_trajectories",{},
                          {"campaign_manifest":str(Path(campaign_path).resolve()),"campaign_type":self.campaign_type(campaign),"source_archives":[str(path) for path in paths],"trajectory_count":len(paths),"storage":"primary-archive-pointers-only"},label)

    def campaign_fourier(self, campaign_path: str | Path, *, components: tuple[int,...] | None=None,
                         label: str="campaign-fourier") -> Path:
        """Compute spectra for all uniformly sampled completed campaign trajectories."""
        campaign,paths=self._campaign_solution_paths(campaign_path)
        if not paths: raise ValueError("The campaign has no completed solution records.")
        series=[VIEWER_ARCHIVE.read_time_series(path) for path in paths]; first=series[0]; selected=tuple(range(first.values.shape[1])) if components is None else tuple(components)
        dt=float(np.diff(first.time)[0]); count=first.sample_count
        if any(item.sample_count!=count or not item.uniform_time_steps or not np.isclose(item.delta_t,dt) for item in series): raise ValueError("Campaign Fourier display requires equally sampled, uniform trajectories.")
        power=np.stack([np.abs(np.fft.rfft(item.values[:,selected],axis=0))**2/count for item in series],axis=2)
        return self._save("campaign_fourier",{"frequency":np.fft.rfftfreq(count,dt),"power":power,"components":np.asarray(selected,dtype=int)},
                          {"campaign_manifest":str(Path(campaign_path).resolve()),"campaign_type":self.campaign_type(campaign),"source_archives":[str(path) for path in paths],"trajectory_count":len(series),"delta_t":dt},label)

    @staticmethod
    def campaign_type(campaign: dict[str,Any]) -> str:
        text=" ".join(str(item).lower() for item in (campaign.get("goal",""),(campaign.get("metadata") or {}).get("campaign_type",""),(campaign.get("metadata") or {}).get("kind","")))
        if "bifurcation" in text: return "bifurcation"
        if "lyapunov" in text: return "lyapunov"
        if "ensemble" in text or "unforced" in text: return "trajectory_ensemble"
        if "sweep" in text: return "parameter_sweep"
        return "generic"

    def bifurcation_from_manifest(self, manifest_path: str | Path, *, parameter_key: str="parameter",
                                  tail_fraction: float=.25, label: str="bifurcation") -> Path:
        """Extract maxima from an existing simple run-to-solution manifest."""
        path=Path(manifest_path).expanduser().resolve(); payload=VIEWER_ARCHIVE.read_json(path); points=[]; sources=[]
        for run in payload.get("runs",[]):
            source=self._primary_path(run["solution_archive"]); series=VIEWER_ARCHIVE.read_time_series(source); values=np.asarray(series.values[:,0],dtype=float)
            start=max(1,int((1-float(tail_fraction))*values.size)); maxima=np.flatnonzero((values[start+1:-1]>=values[start:-2])&(values[start+1:-1]>=values[start+2:]))+start+1
            parameter=float(run[parameter_key]); points.extend((parameter,float(values[index])) for index in maxima); sources.append(str(source))
        point_array=np.asarray(points,dtype=float).reshape((-1,2)) if points else np.empty((0,2),dtype=float)
        return self._save("bifurcation",{"pressure_kpa":point_array[:,0],"radius_ratio":point_array[:,1]},
                          {"source_manifest":str(path),"source_archives":sources,"parameter_key":parameter_key,"tail_fraction":float(tail_fraction),"point_count":int(point_array.shape[0])},label)
