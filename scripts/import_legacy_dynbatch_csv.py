"""Command-line entry point for the reusable legacy MATLAB CSV decoder."""
from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dynsyslib.legacy_dynbatch_csv import _main


if __name__ == "__main__":
    raise SystemExit(_main())
