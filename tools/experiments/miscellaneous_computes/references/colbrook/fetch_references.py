"""Plan or download pinned research sources without running or installing them."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess


def git(*args: str, cwd: Path | None = None) -> str:
    """Run Git with argument boundaries preserved and return its output."""
    return subprocess.check_output(
        ["git", *args], cwd=cwd, text=True, stderr=subprocess.STDOUT
    ).strip()


def main() -> None:
    """Print the acquisition plan by default, or fetch the selected pinned commits."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true",
                        help="fetch source checkouts, without installing or executing them")
    parser.add_argument("--repository", action="append", default=[],
                        help="owner/name to select, repeatable; default is all available")
    parser.add_argument("--destination", type=Path,
                        default=Path(__file__).resolve().parents[5]
                        / "build" / "reference_sources" / "colbrook")
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name("repositories.json").read_text(encoding="utf-8"))
    if manifest["schema"] != "dynsyslab.reference-sources/v1":
        raise ValueError("Unsupported reference manifest")
    entries = manifest["repositories"]
    known = {r["repository"] for r in entries}
    if set(args.repository) - known:
        parser.error("Unknown repository selection: " + ", ".join(sorted(set(args.repository) - known)))
    root = args.destination.resolve()
    for entry in entries:
        name = entry["repository"]
        if args.repository and name not in args.repository:
            continue
        if entry["status"] != "verified":
            print(f"UNAVAILABLE {name}: {entry['reason']}")
            continue
        sha = entry["commit"]
        if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", name):
            raise ValueError("Invalid repository name")
        if not re.fullmatch(r"[0-9a-f]{40}", sha):
            raise ValueError(f"Invalid pinned commit for {name}")
        url = f"https://github.com/{name}.git"
        if entry["clone_url"] != url:
            raise ValueError(f"Unexpected clone URL for {name}")
        target = root / (name.replace("/", "--") + "--" + sha[:12])
        if not target.resolve().is_relative_to(root):
            raise ValueError("Checkout destination escaped the requested root")
        print(f"{name} @ {sha}\n  {target}\n  license: {entry.get('license') or 'not declared'}")
        if not args.download:
            continue
        if target.exists():
            # Existing checkouts are evidence. Never reset, clean, or overwrite them.
            if not (target / ".git").is_dir():
                raise RuntimeError(f"Existing target is not a Git checkout: {target}")
            if git("rev-parse", "HEAD", cwd=target) != sha:
                raise RuntimeError(f"Existing checkout has a different revision: {target}")
            if git("remote", "get-url", "origin", cwd=target) != url:
                raise RuntimeError(f"Existing checkout has a different origin: {target}")
            if git("status", "--porcelain", cwd=target):
                raise RuntimeError(f"Existing checkout has local edits: {target}")
            print("  Already present and unchanged")
            continue
        target.mkdir(parents=True)
        git("init", "--quiet", str(target))
        git("remote", "add", "origin", url, cwd=target)
        git("-c", "fetch.recurseSubmodules=false", "fetch", "--depth=1", "origin", sha, cwd=target)
        empty_hooks = target / ".git" / "disabled-hooks"
        empty_hooks.mkdir()
        git("-c", "core.autocrlf=false", "-c", "core.hooksPath=" + str(empty_hooks),
            "checkout", "--quiet", "--detach", sha, cwd=target)
        if git("rev-parse", "HEAD", cwd=target) != sha:
            raise RuntimeError(f"Downloaded commit mismatch for {name}")
        print("  Downloaded. No example, dependency installer, or build was run.")
    if not args.download:
        print("Plan only. Add --download to fetch the selected reference source checkouts.")


if __name__ == "__main__":
    main()
