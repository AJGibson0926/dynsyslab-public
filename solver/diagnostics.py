"""Compatibility adapters for cross-package backend diagnostics."""
from __future__ import annotations

from importlib import import_module


def backend_diagnostics() -> dict:
    """Return the versioned suite-wide backend capability record."""
    return import_module("operations.backend_diagnostics").backend_diagnostics()


def main() -> int:
    """Run the suite-wide backend diagnostic command."""
    return import_module("operations.backend_diagnostics").main()


__all__=["backend_diagnostics","main"]


if __name__=="__main__":
    raise SystemExit(main())
