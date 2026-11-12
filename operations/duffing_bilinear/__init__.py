"""Experiment-level bounded Duffing bilinear-control search workflow."""
from .campaign import run
from .common import (DEFAULT_CONTROL_SOURCE, DEFAULT_MANIFEST, DEFAULT_OUTPUT,
                     ACTIVE_WINDOW, FINAL_TIME, INITIAL_STATE, LEFT_EQUILIBRIUM,
                     TARGET_STATE, _jsonable)

__all__ = ["run", "DEFAULT_CONTROL_SOURCE", "DEFAULT_MANIFEST", "DEFAULT_OUTPUT",
           "ACTIVE_WINDOW", "FINAL_TIME", "INITIAL_STATE", "LEFT_EQUILIBRIUM",
           "TARGET_STATE"]
