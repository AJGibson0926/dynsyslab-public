"""Stable application-facing API for GUIs, notebooks, and automation."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from typing import Any
import numpy as np
from dynsyslib import SYSTEM_PRESETS, SinusoidalDriver, TabulatedDriver, ZeroDriver, load_system
from solver import HistoryOptions, InterpolationOptions, SolveRequest, SolverOptions, solve

@dataclass(frozen=True)
class RunConfiguration:
    system_name: str
    preset: str | None=None
    parameters: dict[str,Any]=field(default_factory=dict)
    driver: dict[str,Any]=field(default_factory=lambda:{"kind":"zero"})
    time_span: tuple[float,float]=(0.,1.)
    initial_state: tuple[float,...] | None=None
    solver: SolverOptions=field(default_factory=SolverOptions)
    history: HistoryOptions=field(default_factory=HistoryOptions)
    interpolation: InterpolationOptions=field(default_factory=InterpolationOptions)
    label: str="solve"

    def make_driver(self):
        spec=dict(self.driver); kind=spec.pop("kind","zero")
        if kind=="zero": return ZeroDriver()
        if kind=="sinusoid": return SinusoidalDriver(**spec)
        if kind=="tabulated": return TabulatedDriver(np.asarray(spec["time"]),np.asarray(spec["values"]))
        raise ValueError(f"Unknown driver kind {kind!r}.")

    def build_request(self, *, controller=None, estimator=None, postprocessors=(), on_event=None, cancel_token=None) -> SolveRequest:
        system=load_system(self.system_name,driver=self.make_driver(),preset=self.preset,**self.parameters)
        return SolveRequest(system=system,kind=getattr(system,"kind","ode"),time_span=self.time_span,
            initial_state=None if self.initial_state is None else np.asarray(self.initial_state,dtype=float),
            controller=controller,estimator=estimator if estimator is not None else SolveRequest.__dataclass_fields__["estimator"].default_factory(),
            options=self.solver,history=self.history,interpolation=self.interpolation,postprocessors=tuple(postprocessors),on_event=on_event,cancel_token=cancel_token,label=self.label)

class DynsysLabService:
    """Small public facade; front ends need not import package internals."""
    def available_systems(self): return {name:sorted(presets) for name,presets in SYSTEM_PRESETS.items()}
    def validate(self, config: RunConfiguration) -> list[str]:
        errors=[]
        try: request=config.build_request();
        except Exception as error: return [str(error)]
        if request.time_span is None or request.time_span[1]<=request.time_span[0]: errors.append("time_span must be increasing.")
        return errors
    def request(self, config: RunConfiguration, **kwargs) -> SolveRequest:
        errors=self.validate(config)
        if errors: raise ValueError("; ".join(errors))
        return config.build_request(**kwargs)
    def run(self, config: RunConfiguration, **kwargs): return solve(self.request(config,**kwargs))
