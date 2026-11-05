"""Optional external-package adapters; nothing external is vendored here."""
from __future__ import annotations

def optional_backend(name: str):
    try: return __import__(name)
    except ImportError as exc: raise ImportError(f"Optional Koopman backend '{name}' is not installed.") from exc
