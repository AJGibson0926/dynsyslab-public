"""Backend selection owned by reduced-order-model execution."""
from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from typing import Literal


ROMBackendPreference = Literal["auto", "python", "native"]


@dataclass(frozen=True)
class ROMBackend:
    name: str
    compiled: bool
    krom_rollout: bool
    auto_admitted: bool = False
    version: str | None = None
    reason: str | None = None

    def provenance(self, requested: ROMBackendPreference) -> dict:
        return {
            "requested": requested,
            "name": self.name,
            "compiled": self.compiled,
            "krom_rollout": self.krom_rollout,
            "auto_admitted": self.auto_admitted,
            "version": self.version,
        }


PYTHON_ROM_BACKEND = ROMBackend("numpy-reference", False, True, auto_admitted=True, version=None)


def native_rom_backend() -> ROMBackend:
    try:
        module = import_module("_dynsyslab_native")
    except ImportError as error:
        return ROMBackend("dynsyslab-native", True, False, reason=str(error))
    capabilities = module.capabilities()
    return ROMBackend(
        "dynsyslab-native",
        True,
        bool(capabilities.get("krom_rollout", False)),
        bool(capabilities.get("krom_rollout_auto", False)),
        str(capabilities.get("version", "unknown")),
    )


def select_krom_backend(preference: ROMBackendPreference) -> ROMBackend:
    if preference not in {"auto", "python", "native"}:
        raise ValueError("execution_backend must be 'auto', 'python', or 'native'.")
    if preference == "python":
        return PYTHON_ROM_BACKEND
    candidate = native_rom_backend()
    if preference == "native":
        if not candidate.krom_rollout:
            raise RuntimeError(candidate.reason or "The native backend does not support KROM rollout.")
        return candidate
    return candidate if candidate.krom_rollout and candidate.auto_admitted else PYTHON_ROM_BACKEND


__all__ = ["ROMBackend", "ROMBackendPreference", "select_krom_backend"]
