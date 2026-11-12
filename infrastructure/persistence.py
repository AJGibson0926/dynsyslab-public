"""Neutral validation contracts for versioned persisted records."""
from __future__ import annotations

from typing import Any


class ArchiveSchemaError(ValueError):
    """Raised when a persisted record is not compatible with its reader."""


def require_format(
    payload: dict[str, Any],
    expected_prefix: str,
    *,
    supported_suffixes: tuple[str, ...] = ("v1",),
) -> str:
    """Validate and return a versioned record format identifier."""
    value = str(payload.get("format", ""))
    if not value:
        raise ArchiveSchemaError("Persisted record has no format identifier.")
    if not value.startswith(expected_prefix):
        raise ArchiveSchemaError(
            f"Expected {expected_prefix!r} record, got {value!r}."
        )
    if not any(value.endswith(suffix) for suffix in supported_suffixes):
        raise ArchiveSchemaError(f"Unsupported format version {value!r}.")
    return value


def require_keys(payload: dict[str, Any], *keys: str) -> None:
    """Reject a persisted record that omits required keys."""
    missing = [key for key in keys if key not in payload]
    if missing:
        raise ArchiveSchemaError(
            f"Persisted record is missing required fields: {missing!r}."
        )


__all__ = ["ArchiveSchemaError", "require_format", "require_keys"]
