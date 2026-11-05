"""Composable, archive-aware workflows spanning the full post-processing suite."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from time import perf_counter
import numpy as np
from .core import Trajectory
from .method_handlers import SUPPORTED_METHODS, handler_for
from .policies import TRAJECTORY_HDMD_SETTINGS, TRAJECTORY_SINDY_SETTINGS, workflow_policy
from .workflow import SourceResolver, WorkflowContext

@dataclass(frozen=True)
class ProcessingMethod:
    name: str
    settings: dict[str,Any]=field(default_factory=dict)
@dataclass
class ProcessingMethodResult:
    method: str; result: Any; archive: Path|None=None; timing: dict[str,float]=field(default_factory=dict)
@dataclass
class ProcessingRun:
    source: Trajectory; results: list[ProcessingMethodResult]; workflow_policy: str|None=None
    @property
    def spectrum(self): return next((getattr(item.result,"spectrum",None) for item in reversed(self.results) if getattr(item.result,"spectrum",None) is not None),None)
    @property
    def krom(self): return next((getattr(item.result,"krom",None) for item in reversed(self.results) if getattr(item.result,"krom",None) is not None),None)

class PostprocessingPipeline:
    """Run ordered Koopman, PCA/POD, SVD, and HAVOK processing workflows."""
    supported_methods=SUPPORTED_METHODS
    def __init__(self,sources:Any|None=None,*,archive:bool=True): self.sources=[] if sources is None else (list(sources) if isinstance(sources,(list,tuple)) else [sources]); self.methods=[]; self.archive=bool(archive); self.workflow_policy=None
    def add_source(self,source): self.sources.append(source); return self
    def add_method(self,name:str,**settings):
        name=name.lower()
        if name not in self.supported_methods: raise ValueError(f"Unsupported post-processing method '{name}'. Available: {', '.join(sorted(self.supported_methods))}.")
        self.methods.append(ProcessingMethod(name,dict(settings))); return self
    @classmethod
    def preset(cls,name="default",sources=None,*,archive:bool=True,sindy_settings:dict[str,Any]|None=None,**hdmd_settings):
        policy=workflow_policy(name)
        hdmd_defaults=policy.hdmd_settings()
        hdmd_defaults.update(hdmd_settings)
        configured_sindy=policy.sindy_settings(); configured_sindy.update(sindy_settings or {})
        pipeline=cls(sources,archive=archive); pipeline.workflow_policy=policy.name
        return pipeline.add_method("hdmd",**hdmd_defaults).add_method("sindy",**configured_sindy)
    def describe(self) -> dict[str,Any]:
        """Return a JSON-safe, nonexecuting description of this workflow."""
        def clean(value):
            if isinstance(value,Path): return str(value)
            if isinstance(value,np.generic): return value.item()
            if isinstance(value,(str,int,float,bool)) or value is None: return value
            if isinstance(value,(list,tuple)): return [clean(item) for item in value]
            if isinstance(value,dict): return {str(key):clean(item) for key,item in value.items()}
            return repr(value)
        return {"format":"dynsyslab.postprocessing-workflow/v1","workflow_policy":self.workflow_policy,"archive":self.archive,
                "sources":[clean(source) for source in self.sources],
                "methods":[{"method":item.name,"settings":clean(item.settings)} for item in self.methods]}
    def _sources(self):
        return SourceResolver().resolve(self.sources)
    def _archive(self,source,name,result):
        if not self.archive or source.source is None: return None
        return handler_for(name).publish(source,result)
    def _run_one(self,source):
        context=WorkflowContext(); outputs=[]
        for attachment in self.methods:
            method_started=perf_counter()
            name=attachment.name
            prepared=context.prepare(source,attachment.settings)
            result=handler_for(name).run(prepared.source,source,prepared.settings,context)
            algorithm_seconds=perf_counter()-method_started
            if hasattr(result,"metadata"):
                result.metadata["processor_attachment"]={"method":name,"workflow_policy":self.workflow_policy,
                    "settings":{k:repr(v) for k,v in attachment.settings.items()}}
                result.metadata["processor_timing"]={"algorithm_seconds":algorithm_seconds,"clock":"time.perf_counter"}
            archive_started=perf_counter(); archive=self._archive(source,name,result); archive_seconds=perf_counter()-archive_started
            timing={"algorithm_seconds":algorithm_seconds,"archive_seconds":archive_seconds,"total_seconds":perf_counter()-method_started}
            if hasattr(result,"metadata"): result.metadata["processor_timing"].update(timing)
            context.observe(result,archive)
            outputs.append(ProcessingMethodResult(name,result,archive,timing))
        return ProcessingRun(source,outputs,self.workflow_policy)
    def run(self):
        if not self.methods: raise ValueError("Attach at least one post-processing method before running the processor.")
        return [self._run_one(source) for source in self._sources()]


# Historical compatibility for callers that imported the original lowercase
# class. New code should use PostprocessingPipeline.
processor=PostprocessingPipeline


__all__=["PostprocessingPipeline","processor","ProcessingMethod","ProcessingMethodResult","ProcessingRun",
         "TRAJECTORY_HDMD_SETTINGS","TRAJECTORY_SINDY_SETTINGS"]
