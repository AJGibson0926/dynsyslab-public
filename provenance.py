"""Shared compact serialization and provenance helpers for dynsyslab."""
from __future__ import annotations
from dataclasses import asdict, fields, is_dataclass
from typing import Any
import hashlib, json
import numpy as np

def serializable(value: Any) -> Any:
    if is_dataclass(value): return {field.name:serializable(getattr(value,field.name)) for field in fields(value) if field.init}
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,(list,tuple)): return [serializable(item) for item in value]
    if isinstance(value,dict): return {str(key):serializable(item) for key,item in value.items()}
    if isinstance(value,(str,int,float,bool)) or value is None: return value
    return {"class":type(value).__name__,"repr":repr(value)}

def serializable_system(system: Any) -> dict[str, Any]:
    """Serialize a system while retaining its reconstructable type identity."""
    archive_metadata=getattr(system,"archive_metadata",None)
    if callable(archive_metadata):
        payload=serializable(archive_metadata())
        payload["_dynsyslab_class"]=type(system).__name__
        payload["_dynsyslab_module"]=type(system).__module__
        return payload
    if is_dataclass(system):
        payload=serializable(system)
        payload["_dynsyslab_class"]=type(system).__name__
        payload["_dynsyslab_module"]=type(system).__module__
        return payload
    return {"repr":repr(system),"class":type(system).__name__}

def compact_json(value: Any) -> str: return json.dumps(serializable(value),separators=(",",":"),sort_keys=True)
def fingerprint(value: Any) -> str: return hashlib.sha256(compact_json(value).encode()).hexdigest()
