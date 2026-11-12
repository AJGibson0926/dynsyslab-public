"""Named, documented presets for individual ODE system modules."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any

@dataclass(frozen=True)
class SystemPreset:
    """A reproducible model setup plus optional post-processing recommendations.

    ``processing_metadata`` is keyed by processing-method name (for example,
    ``"hdmd"``, ``"havok"``, ``"dmd"``, or ``"sindy"``).  Each entry may
    retain method settings and provenance such as sampling, delay, and matrix
    geometry without passing those analysis settings into the system model.
    """
    key: str
    description: str
    parameters: dict[str, Any] = field(default_factory=dict)
    initial_state: tuple[float, ...] | None = None
    recommended_driver: dict[str, float] = field(default_factory=dict)
    reference: str = ""
    processing_metadata: dict[str, dict[str, Any]] = field(default_factory=dict)

    def construction_kwargs(self) -> dict[str, Any]:
        values=dict(self.parameters)
        if self.initial_state is not None: values["initial_state_override"]=self.initial_state
        return values

    def processing_settings(self, method: str) -> dict[str, Any]:
        """Return a copy of directly usable settings for one named method.

        Metadata-only fields such as paper matrix geometry remain available in
        :attr:`processing_metadata`; only the nested ``settings`` mapping is
        returned here, so it can be supplied safely to a processor call.
        """
        entry=self.processing_metadata.get(method.lower(),{})
        return dict(entry.get("settings",{}))
