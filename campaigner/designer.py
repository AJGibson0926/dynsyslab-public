"""Generic, dynsyslib-backed campaign design with no system-specific logic."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Literal
from uuid import uuid4
import json, os
import numpy as np
from infrastructure.persistence import require_format, require_keys
from archiver.gateway import MANAGER_ARCHIVE
from job_manager import JobCollection, JobManifest

CampaignKind=Literal["bifurcation","lyapunov","trajectory_ensemble","parameter_sweep"]
BIFURCATION_PROFILES={"low":{"point_count":31,"systems_per_job":10},"medium":{"point_count":101,"systems_per_job":25},"high":{"point_count":301,"systems_per_job":50},"ultra_high":{"point_count":1001,"systems_per_job":100}}
# A bifurcation point can require a long transient before its retained samples
# are meaningful.  This policy is intentionally declarative: JobManager, not
# Solver, decides how it becomes solver execution settings.
CHECKPOINT_POLICY_PROFILES={
    # Fast is deliberately the default: ordinary solve_ivp, lightweight
    # telemetry, and only the completed solution in the primary archive.
    "fast": {},
    # Use this profile when avoiding restart after interruption is worth the
    # low-level stepper's additional bookkeeping cost.
    "recoverable": {
    "enabled": True,
    "if_estimated_seconds_over": 300,
    # New campaigns generally have no empirical runtime estimate yet.  Do not
    # silently lose recovery protection on their first pass.
    "when_estimate_unknown": "enable",
    "mode": "checkpoint",
    "checkpoint_interval": 600.0,
    "strategy": "low_level_stepper",
    "pause_at_checkpoint": True,
    },
}
EXECUTION_PROFILES={"fixed":{"adaptive_policy":"fixed","backend":"serial","retry_limit":0},"adaptive_local":{"adaptive_policy":"adaptive_local","backend":"local_process","retry_limit":1},"adaptive_portable":{"adaptive_policy":"adaptive_portable","backend":"local_process","retry_limit":1,"remote_allowed":True}}
def _stamp() -> str: return datetime.now(timezone.utc).isoformat()
def _safe(value: str) -> str: return "".join(c if c.isalnum() or c in "-_" else "-" for c in value)

def _checkpoint_policy(integration: dict[str, Any]) -> tuple[dict[str, Any], str]:
    """Return a named policy profile plus an optional mapping override."""
    raw=integration.get("checkpoint_policy",{})
    profile=str(integration.get("checkpoint_policy_profile", raw if isinstance(raw,str) else "fast"))
    if profile not in CHECKPOINT_POLICY_PROFILES:
        raise ValueError(f"Unknown checkpoint policy profile {profile!r}.")
    overrides={} if isinstance(raw,str) else dict(raw)
    return {**CHECKPOINT_POLICY_PROFILES[profile],**overrides}, profile

@dataclass(frozen=True)
class CampaignRequest:
    kind: CampaignKind; label: str; system_name: str; preset: str|None=None; quality: str="medium"
    parameters: dict[str,Any]=field(default_factory=dict); integration: dict[str,Any]=field(default_factory=dict)

@dataclass(frozen=True)
class CampaignManifest:
    label: str; goal: str; jobs: tuple[str,...]; metadata: dict[str,Any]
    campaign_id: str=field(default_factory=lambda:uuid4().hex); created_utc: str=field(default_factory=_stamp)
    format: str="dynsyslab.campaign/v1"
    def payload(self)->dict[str,Any]: return asdict(self)
    def write(self,path:str|Path)->Path:
        target=Path(path); target.parent.mkdir(parents=True,exist_ok=True); tmp=target.with_suffix(target.suffix+f".{uuid4().hex}.partial")
        with tmp.open("w",encoding="utf-8",newline="\n") as h: h.write(json.dumps(self.payload(),indent=2,sort_keys=True)); h.flush(); os.fsync(h.fileno())
        os.replace(tmp,target); return target
    @classmethod
    def read(cls,path:str|Path)->"CampaignManifest":
        data=json.loads(Path(path).read_text(encoding="utf-8")); require_format(data,"dynsyslab.campaign/"); require_keys(data,"campaign_id","label","jobs","metadata"); data["jobs"]=tuple(data["jobs"]); return cls(**data)

class Campaigner:
    """Design lightweight plans for registered systems; never solve or analyse."""
    def __init__(self,root:str|Path|None=None): self.root=Path(root or MANAGER_ARCHIVE.archive_root/"campaigns").resolve(); self.root.mkdir(parents=True,exist_ok=True)
    def save(self,campaign:CampaignManifest)->Path: return campaign.write(self.root/f"{_safe(campaign.label)}__{campaign.campaign_id}.json")
    @staticmethod
    def request_profiles(kind:str)->dict[str,dict[str,Any]]:
        if kind=="bifurcation": return {key:dict(value) for key,value in BIFURCATION_PROFILES.items()}
        if kind=="lyapunov": return {"screen":{"transient_cycles":2500,"measure_cycles":1024},"verification":{"transient_cycles":5000,"measure_cycles":4096}}
        return {"small":{"systems_per_job":10},"medium":{"systems_per_job":50},"large":{"systems_per_job":100}}
    @staticmethod
    def execution_profiles()->dict[str,dict[str,Any]]: return {key:dict(value) for key,value in EXECUTION_PROFILES.items()}
    @staticmethod
    def checkpoint_policy_profiles()->dict[str,dict[str,Any]]:
        """Named policy choices available to every campaign type."""
        return {key:dict(value) for key,value in CHECKPOINT_POLICY_PROFILES.items()}
    def _validate(self,system_name:str,preset:str|None):
        from dynsyslib import SYSTEM_LIBRARY,SYSTEM_PRESETS,load_system
        if system_name not in SYSTEM_LIBRARY: raise KeyError(f"Unknown dynsyslib system {system_name!r}.")
        if preset is not None:
            if preset not in SYSTEM_PRESETS.get(system_name,{}): raise KeyError(f"Unknown preset {preset!r} for {system_name!r}.")
            load_system(system_name,preset=preset)
    def _write_bundle(self,campaign:CampaignManifest,jobs:list[JobManifest],quality:str)->tuple[CampaignManifest,JobCollection]:
        path=self.save(campaign); plan_dir=self.root/"job_plans"/campaign.campaign_id; plan_dir.mkdir(parents=True,exist_ok=True); tagged=[]
        for job in jobs:
            item=replace(job,tags={**job.tags,"campaign_manifest":str(path)}); tagged.append(item); item.write(plan_dir/f"{_safe(item.label)}__{item.job_id}.json")
        updated=replace(campaign,metadata={**campaign.metadata,"job_plan_paths":[str(plan_dir/f"{_safe(j.label)}__{j.job_id}.json") for j in tagged]}); updated.write(path)
        if updated.metadata["campaign_type"]=="bifurcation": self.enable_bifurcation_updates(path)
        self.write_planning_report(path); return updated,JobCollection(label=f"{updated.label}-{quality}",jobs=tuple(tagged),tags={"campaign_id":updated.campaign_id,"campaign_manifest":str(path),"goal":updated.goal,"quality":quality})
    def enable_bifurcation_updates(self,path:str|Path,*,interval_percent:int=10)->Path:
        if not 1<=interval_percent<=100: raise ValueError("interval_percent must be between 1 and 100.")
        path=Path(path); campaign=CampaignManifest.read(path); log=self.root/"logs"/f"{campaign.campaign_id}.jsonl"; log.parent.mkdir(parents=True,exist_ok=True)
        updated=replace(campaign,metadata={**campaign.metadata,"progress_reporting":{"enabled":True,"interval_percent":interval_percent,"log_path":str(log)}}); updated.write(path)
        if not log.exists(): log.write_text(json.dumps({"event":"campaign_designed","time_utc":_stamp(),"campaign_id":campaign.campaign_id})+"\n",encoding="utf-8")
        return log
    def write_planning_report(self,path:str|Path)->tuple[Path,Path]:
        path=Path(path); campaign=CampaignManifest.read(path); plans=[Path(p) for p in campaign.metadata.get("job_plan_paths",[]) if Path(p).is_file()]
        count=sum(len(JobManifest.read(p).systems) for p in plans); first=JobManifest.read(plans[0]).systems[0] if plans else {}
        samples=campaign.metadata.get("estimated_uniform_output_samples_per_trajectory"); state_dim=len(first.get("initial_state",[])) or None
        report={"format":"dynsyslab.campaign-planning-report/v1","campaign_id":campaign.campaign_id,"scheduled":False,"executed":False,"campaign_type":campaign.metadata.get("campaign_type"),"system":campaign.metadata.get("system_name"),"preset":campaign.metadata.get("preset"),"planned_trajectory_count":count,"parameter_space":campaign.metadata.get("parameter_space"),"integration":campaign.metadata.get("integration"),"checkpoint_policy":campaign.metadata.get("checkpoint_policy"),"solver_options":campaign.metadata.get("solver_options"),"estimated_uniform_output_samples_per_trajectory":samples,"uncompressed_final_time_state_bytes":None if not(samples and state_dim) else count*samples*(state_dim+1)*8}
        directory=self.root/"reports"; directory.mkdir(exist_ok=True); base=directory/f"{_safe(campaign.label)}__{campaign.campaign_id}"; json_path=base.with_suffix(".json"); md_path=base.with_suffix(".md")
        json_path.write_text(json.dumps(report,indent=2,sort_keys=True),encoding="utf-8"); md_path.write_text("# Campaign planning report\n\n```json\n"+json.dumps(report,indent=2,sort_keys=True)+"\n```\n",encoding="utf-8"); return json_path,md_path
    def bifurcation_campaign(self,*,system_name:str,preset:str|None,driver:dict[str,Any],sweep:dict[str,Any],integration:dict[str,Any],label:str,quality:str="medium",solver_options:dict[str,Any]|None=None,metadata:dict[str,Any]|None=None)->tuple[CampaignManifest,JobCollection]:
        self._validate(system_name,preset)
        if quality not in BIFURCATION_PROFILES: raise ValueError(f"Unknown quality {quality!r}.")
        target=str(sweep.get("target",""));
        if not(target.startswith("driver.") or target.startswith("parameters.")): raise ValueError("sweep.target must address driver.<field> or parameters.<field>.")
        if "values" in sweep:
            values=list(sweep["values"])
        else:
            values=list(np.linspace(float(sweep["start"]),float(sweep["stop"]),int(sweep.get("points",BIFURCATION_PROFILES[quality]["point_count"]))))
        if len(values)<2: raise ValueError("A bifurcation campaign requires at least two values.")
        driver=dict(driver); frequency=float(driver.get("frequency_hz",0.)); integration=dict(integration)
        if "time_span" in integration: time_span=tuple(map(float,integration["time_span"])); delta_t=float(integration["delta_t"])
        else:
            if frequency<=0: raise ValueError("Cycle integration requires driver.frequency_hz.")
            total=int(integration["transient_cycles"])+int(integration["retained_cycles"]); delta_t=1/frequency/int(integration["samples_per_cycle"]); time_span=(0.,total/frequency)
        options={"method":"DOP853","rtol":1e-9,"atol":1e-12,**dict(solver_options or integration.get("solver_options",{}))}; chunk=int(integration.get("systems_per_job",BIFURCATION_PROFILES[quality]["systems_per_job"])); cid=uuid4().hex
        policy,policy_profile=_checkpoint_policy(integration); specs=[]
        for i,value in enumerate(values):
            d=dict(driver); params={}; (d if target.startswith("driver.") else params)[target.split(".",1)[1]]=value
            specs.append({"system_name":system_name,"preset":preset,"parameters":params,"driver":d,"time_span":list(time_span),"output":{"delta_t":delta_t},"solver":{"options":options,"history":{},"history_policy":policy,"interpolation_method":"pchip"},"hierarchy":{"clusters":[{"label":"forced-trajectory","initial_conditions":["__system_initial_state__"]}]},"label":f"{system_name}__{target.replace('.','-')}__{i:05d}"})
        jobs=[JobManifest(label=f"{label}-{i//chunk:03d}",systems=tuple(specs[i:i+chunk]),tags={"campaign_id":cid,"campaign_kind":"bifurcation-data","system_name":system_name,"sweep_target":target}) for i in range(0,len(specs),chunk)]
        meta={"campaign_id":cid,"campaign_type":"bifurcation","system_name":system_name,"preset":preset,"resolution_quality":quality,"parameter_space":{"swept_parameter":target,"values":values,"point_count":len(values),"driver":driver},"integration":{**integration,"time_span":list(time_span),"delta_t":delta_t},"checkpoint_policy":policy,"checkpoint_policy_profile":policy_profile,"solver_options":options,"estimated_uniform_output_samples_per_trajectory":int(round((time_span[1]-time_span[0])/delta_t))+1,"hierarchy":"Batch > System > Cluster > Trajectory","execution_profiles":EXECUTION_PROFILES,"default_execution_profile":"fixed",**dict(metadata or {})}
        return self._write_bundle(CampaignManifest(label=label,goal=f"Generate {quality} {system_name} bifurcation data.",jobs=tuple(j.job_id for j in jobs),metadata=meta,campaign_id=cid),jobs,quality)
    def trajectory_ensemble(self,*,system_name:str,preset:str|None,initial_conditions:Iterable[Iterable[float]],integration:dict[str,Any],label:str,parameters:dict[str,Any]|None=None,metadata:dict[str,Any]|None=None)->tuple[CampaignManifest,JobCollection]:
        self._validate(system_name,preset); values=[list(map(float,x)) for x in initial_conditions]; selected_parameters=dict(parameters or {}); cid=uuid4().hex; policy,policy_profile=_checkpoint_policy(integration); spec={"system_name":system_name,"preset":preset,"parameters":selected_parameters,"driver":{"kind":"zero"},"time_span":list(integration["time_span"]),"output":{"delta_t":float(integration["delta_t"])},"solver":{"options":dict(integration.get("solver_options",{})),"history":{},"history_policy":policy},"hierarchy":{"clusters":[{"label":"ensemble","initial_conditions":values}]},"label":label}; job=JobManifest(label=label,systems=(spec,),tags={"campaign_id":cid,"campaign_kind":"trajectory-ensemble","system_name":system_name}); meta={"campaign_id":cid,"campaign_type":"trajectory_ensemble","system_name":system_name,"preset":preset,"parameters":selected_parameters,"integration":dict(integration),"checkpoint_policy":policy,"checkpoint_policy_profile":policy_profile,"initial_condition_count":len(values),"hierarchy":"System > Cluster > Trajectory",**dict(metadata or {})}; return self._write_bundle(CampaignManifest(label=label,goal=f"Generate {system_name} trajectory ensemble.",jobs=(job.job_id,),metadata=meta,campaign_id=cid),[job],"ensemble")
    def design_request(self,request:CampaignRequest)->tuple[CampaignManifest,JobCollection]:
        if request.kind=="trajectory_ensemble": return self.trajectory_ensemble(system_name=request.system_name,preset=request.preset,initial_conditions=request.parameters["initial_conditions"],integration=request.integration,label=request.label,parameters=request.parameters.get("system_parameters"),metadata=request.parameters.get("metadata"))
        if request.kind in {"bifurcation","lyapunov","parameter_sweep"}:
            campaign,collection=self.bifurcation_campaign(system_name=request.system_name,preset=request.preset,driver=request.parameters["driver"],sweep=request.parameters["sweep"],integration=request.integration,label=request.label,quality=request.quality,metadata=request.parameters.get("metadata"));
            if request.kind!="bifurcation":
                path=self.root/f"{_safe(campaign.label)}__{campaign.campaign_id}.json"; campaign=replace(campaign,metadata={**campaign.metadata,"campaign_type":request.kind}); campaign.write(path)
            return campaign,collection
        raise ValueError(f"Unsupported campaign type {request.kind!r}.")
