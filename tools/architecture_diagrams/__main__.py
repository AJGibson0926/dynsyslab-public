"""Architecture-diagram CLI.

    python tools/architecture_diagrams new module-sindy --title "Module: SINDy"   scaffold a spec
    python tools/architecture_diagrams check docs/diagrams/module-sindy.diagram.json
    python tools/architecture_diagrams render docs/diagrams/module-sindy.diagram.json --png
    python tools/architecture_diagrams all --png        re-render every spec in docs/diagrams
    python tools/architecture_diagrams icons

`render SPEC` may be shortened to `SPEC`. Output goes next to the spec unless -o is given.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

if __package__ in (None, ""):  # run as `python tools/architecture_diagrams`
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    __package__ = "architecture_diagrams"

from architecture_diagrams import ICONS, SpecError, load, render_file  # noqa: E402
from architecture_diagrams.templates import TEMPLATES, scaffold  # noqa: E402

DIAGRAMS = Path(__file__).resolve().parents[2] / "docs" / "diagrams"
EDGE = [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe"]


def export_png(svg: Path, png: Path, scale: float = 1.5) -> None:
    """Rasterize with a headless Chromium browser (Edge or Chrome); no Python image dependency."""
    browser = next((p for p in EDGE if Path(p).exists()), None) or shutil.which("msedge") or shutil.which("chrome")
    if not browser:
        raise SystemExit("PNG export needs Edge or Chrome; the SVG was written")
    text = svg.read_text(encoding="utf-8")
    width, height = (float(v) for v in re.search(r'width="([0-9.]+)" height="([0-9.]+)"', text).groups())
    if png.exists():
        png.unlink()
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as profile:
        subprocess.run([browser, "--headless", "--disable-gpu", "--no-first-run", "--hide-scrollbars",
                        f"--user-data-dir={profile}", f"--force-device-scale-factor={scale}",
                        f"--window-size={int(width)},{int(height)}", f"--screenshot={png.resolve()}",
                        svg.resolve().as_uri()], check=True, capture_output=True, timeout=120)
    deadline = time.monotonic() + 20  # headless Chromium can return before the file is flushed
    while time.monotonic() < deadline and not (png.exists() and png.stat().st_size > 0):
        time.sleep(0.25)
    if not png.exists():
        raise SystemExit("the browser did not write the PNG; the SVG was written")


def _svg_for(spec: Path, output: str | None) -> Path:
    return Path(output) if output else spec.with_name(spec.name.replace(".diagram.json", "") + ".svg")


def _render(spec: Path, output: str | None, png) -> int:
    svg = _svg_for(spec, output)
    try:
        render_file(spec, svg)
    except (SpecError, json.JSONDecodeError) as error:
        print(f"{spec}: invalid diagram spec: {error}", file=sys.stderr)
        return 2
    if png:
        export_png(svg, svg.with_suffix(".png") if png is True else Path(png))
    print(svg)
    return 0


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0].endswith(".json"):
        argv.insert(0, "render")
    parser = argparse.ArgumentParser(prog="architecture_diagrams", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    render = sub.add_parser("render", help="render one spec to SVG (and PNG)")
    render.add_argument("spec")
    render.add_argument("-o", "--output", help="SVG path (default: next to the spec)")
    render.add_argument("--png", nargs="?", const=True, help="also write a PNG (default: next to the SVG)")
    every = sub.add_parser("all", help="render every *.diagram.json in a folder")
    every.add_argument("--dir", default=str(DIAGRAMS))
    every.add_argument("--png", action="store_true")
    check = sub.add_parser("check", help="validate specs without writing anything")
    check.add_argument("specs", nargs="+")
    new = sub.add_parser("new", help="scaffold a spec from a template")
    new.add_argument("name", help="file stem, e.g. module-input-map")
    new.add_argument("--title", default="")
    new.add_argument("--template", choices=sorted(TEMPLATES), default="module")
    new.add_argument("--dir", default=str(DIAGRAMS))
    new.add_argument("--force", action="store_true", help="overwrite an existing spec")
    sub.add_parser("icons", help="list icon names")
    args = parser.parse_args(argv)

    if args.command == "icons":
        print("\n".join(sorted(ICONS)))
        return 0
    if args.command == "render":
        return _render(Path(args.spec), args.output, args.png)
    if args.command == "all":
        specs = sorted(Path(args.dir).glob("*.diagram.json"))
        return max([_render(spec, None, args.png) for spec in specs], default=0)
    if args.command == "check":
        status = 0
        for spec in args.specs:
            try:
                diagram = load(spec, set(ICONS))
                print(f"{spec}: ok ({len(diagram.blocks)} blocks, {len(diagram.edges)} edges)")
            except (SpecError, json.JSONDecodeError) as error:
                print(f"{spec}: {error}", file=sys.stderr)
                status = 2
        return status
    path = Path(args.dir) / f"{args.name}.diagram.json"
    if path.exists() and not args.force:
        print(f"{path} exists; pass --force to overwrite", file=sys.stderr)
        return 2
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(scaffold(args.template, args.title or args.name), indent=1, ensure_ascii=False) + "\n",
                    encoding="utf-8", newline="\n")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
