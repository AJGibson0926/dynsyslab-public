"""Link dynsyslab to the centrally managed Antenna repository.

This adapter owns no Antenna implementation. It resolves the shared checkout
through ``ANTENNA_HOME``, verifies its approved Git commit, then invokes the
upstream launcher's exact command line.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_POLICY = PROJECT_ROOT / "config" / "antenna_policy.json"
DEPENDENCY_MANIFEST = PROJECT_ROOT / "config" / "antenna_dependency.json"


class AntennaUnavailable(RuntimeError):
    """Raised when the centrally managed Antenna cannot be resolved safely."""


@dataclass(frozen=True)
class AntennaLink:
    """A resolved reference to the canonical, external Antenna launcher."""

    antenna_root: Path
    launcher: Path
    policy: Path
    commit: str

    @classmethod
    def resolve(cls, *, policy: Path = DEFAULT_POLICY) -> "AntennaLink":
        try:
            dependency = json.loads(DEPENDENCY_MANIFEST.read_text(encoding="utf-8"))
            approved_commit = dependency["commit"]
        except (OSError, json.JSONDecodeError, KeyError) as error:
            raise AntennaUnavailable(f"Antenna dependency manifest is invalid: {DEPENDENCY_MANIFEST}") from error
        configured_root = os.environ.get("ANTENNA_HOME")
        if not configured_root:
            raise AntennaUnavailable("ANTENNA_HOME is not set; point it to the shared the-antenna checkout.")
        antenna_root = Path(configured_root).expanduser().resolve()
        launcher = antenna_root / "tools" / "antenna.py"
        if not (antenna_root / "antenna" / "__init__.py").is_file() or not launcher.is_file():
            raise AntennaUnavailable(f"No valid shared Antenna found under {antenna_root}")
        current = subprocess.run(["git", "-C", str(antenna_root), "rev-parse", "HEAD"], text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if current.returncode or current.stdout.strip() != approved_commit:
            found = current.stdout.strip() or "unavailable"
            raise AntennaUnavailable(f"Antenna commit {found} does not match approved pin {approved_commit}")
        resolved_policy = policy.expanduser().resolve()
        if not resolved_policy.is_file():
            raise AntennaUnavailable(f"dynsyslab Antenna policy is missing: {resolved_policy}")
        return cls(antenna_root=antenna_root, launcher=launcher, policy=resolved_policy, commit=approved_commit)

    def command(self, arguments: Sequence[str]) -> list[str]:
        """Return a command that invokes the approved upstream launcher."""
        return [sys.executable, str(self.launcher), *arguments]

    def run(self, arguments: Sequence[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
        return subprocess.run(self.command(arguments), text=True, check=check)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the approved shared Antenna for dynsyslab.")
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="Arguments passed verbatim to tools/antenna.py.")
    args = parser.parse_args()
    if not args.arguments:
        parser.error("supply an Antenna command, for example: doctor")
    try:
        return AntennaLink.resolve().run(args.arguments).returncode
    except (AntennaUnavailable, OSError, subprocess.SubprocessError) as error:
        print(f"dynsyslab Antenna link failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
