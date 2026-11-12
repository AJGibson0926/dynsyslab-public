"""Compatibility gateway for the experiment-level Duffing search workflow."""
from __future__ import annotations

from importlib import import_module

__all__=["run","DEFAULT_CONTROL_SOURCE","DEFAULT_MANIFEST","DEFAULT_OUTPUT",
    "ACTIVE_WINDOW","FINAL_TIME","INITIAL_STATE","LEFT_EQUILIBRIUM",
    "TARGET_STATE","_jsonable"]


def __getattr__(name: str):
    """Lazily resolve legacy imports without coupling core ROM to experiments."""
    if name not in __all__:
        raise AttributeError(name)
    return getattr(import_module("operations.duffing_bilinear"),name)
