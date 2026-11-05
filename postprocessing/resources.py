"""Conservative resource preflight helpers for dense postprocessing."""
from __future__ import annotations

import math

class ResourceLimitError(MemoryError):
    """Raised before allocating a dense matrix beyond an explicit budget."""

def estimate_dense_bytes(*shapes: tuple[int,...], itemsize: int=8, workspace_factor: float=3.0) -> int:
    """Estimate dense storage plus SVD/linear-algebra workspace conservatively."""
    elements=sum(math.prod(int(max(0,dimension)) for dimension in shape) for shape in shapes)
    return int(elements*itemsize*workspace_factor)

def require_memory_budget(operation: str, *shapes: tuple[int,...], max_memory_mb: float|None=1024.0,
                          workspace_factor: float=3.0, itemsize: int=8) -> int:
    """Fail early rather than allowing a workstation to thrash or terminate."""
    estimate=estimate_dense_bytes(*shapes,itemsize=itemsize,workspace_factor=workspace_factor)
    if max_memory_mb is not None:
        if max_memory_mb<=0: raise ValueError("max_memory_mb must be positive or None.")
        budget=int(float(max_memory_mb)*1024*1024)
        if estimate>budget:
            raise ResourceLimitError(f"{operation} would require about {estimate/1024**2:.1f} MiB of dense working memory, exceeding max_memory_mb={max_memory_mb}. Reduce rank/embedding/subsampling or explicitly raise the budget.")
    return estimate
