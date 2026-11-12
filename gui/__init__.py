"""Lazy compatibility desktop façade; the native workbench is preferred."""
from __future__ import annotations


def __getattr__(name: str):
    if name not in {"DynsysLabDesktop", "launch"}:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from . import desktop
    value = getattr(desktop, name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | {"DynsysLabDesktop", "launch"})


__all__ = ["DynsysLabDesktop", "launch"]
