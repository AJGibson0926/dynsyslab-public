"""Lazy convenience façade for the independently installable suite."""
from __future__ import annotations

from importlib import import_module


_MODULE_EXPORTS={
    "archives","app","control","campaigner","dynsyslib","estimation",
    "job_manager","postprocessing","solver","viewer",
}
_VALUE_EXPORTS={
    "TimeSeries":("time_series","TimeSeries"),
    "FourierTransform":("time_series","FourierTransform"),
    "Trajectory":("hierarchy","Trajectory"),
    "Cluster":("hierarchy","Cluster"),
    "System":("hierarchy","System"),
    "Batch":("hierarchy","Batch"),
}


def __getattr__(name: str):
    if name in _MODULE_EXPORTS:
        value=import_module(name)
    elif name in _VALUE_EXPORTS:
        module,attribute=_VALUE_EXPORTS[name]
        value=getattr(import_module(module),attribute)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name]=value
    return value


def __dir__():
    return sorted(set(globals())|_MODULE_EXPORTS|set(_VALUE_EXPORTS))


__all__=["dynsyslib","solver","archives","app","job_manager","campaigner","control","estimation","postprocessing","viewer","TimeSeries","FourierTransform","Trajectory","Cluster","System","Batch"]


if __name__ == "__main__":
    from operations.front_door import main
    raise SystemExit(main())
