"""Standard-library architecture audit with JSON and Markdown output."""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import json
from pathlib import Path
import re
try:
    import tomllib
except ModuleNotFoundError:  # Historical reference lane supports Python 3.10.
    import tomli as tomllib
from typing import Iterable, Sequence

from project_paths import PROJECT_ROOT, source_files


PACKAGES = frozenset({
    "analysis_data", "archiver", "archives", "campaigner", "codex_harness", "compute_workbench",
    "control", "dynsyslib", "estimation", "gui", "infrastructure",
    "job_manager", "operations", "performance", "postprocessing", "quality",
    "rom", "scripts", "solver", "spectral", "tester", "viewer",
})
LINK = re.compile(r"\[[^]]+\]\(([^)]+)\)")


def _relative(path: Path, root: Path) -> str:
    return path.relative_to(root).as_posix()


def _python_metrics(paths: Iterable[Path], root: Path) -> dict:
    modules: list[dict] = []
    parse_errors: list[dict] = []
    public_total = public_documented = test_methods = 0
    package_public: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    edges: Counter[tuple[str, str]] = Counter()
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        relative = _relative(path, root)
        lines = len(text.splitlines())
        try:
            tree = ast.parse(text, filename=relative)
        except SyntaxError as error:
            parse_errors.append({"path": relative, "line": error.lineno,
                                 "message": error.msg})
            continue
        public = [node for node in tree.body
                  if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
                  and not node.name.startswith("_")]
        documented = sum(ast.get_docstring(node) is not None for node in public)
        source = relative.split("/", 1)[0]
        if not relative.startswith("tester/tests/"):
            public_total += len(public)
            public_documented += documented
            if source in PACKAGES:
                package_public[source][0] += len(public)
                package_public[source][1] += documented
        if relative.startswith("tester/tests/"):
            test_methods += sum(
                1 for node in ast.walk(tree)
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                and node.name.startswith("test_")
            )
        if source in PACKAGES:
            for node in ast.walk(tree):
                targets: list[str] = []
                if isinstance(node, ast.Import):
                    targets = [alias.name.split(".", 1)[0] for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    targets = [node.module.split(".", 1)[0]]
                for target in targets:
                    if target in PACKAGES and target != source:
                        edges[(source, target)] += 1
        modules.append({"path": relative, "lines": lines,
                        "module_docstring": ast.get_docstring(tree) is not None,
                        "public_definitions": len(public),
                        "documented_public_definitions": documented})
    package_edges = [
        {"source": source, "target": target, "imports": count}
        for (source, target), count in sorted(edges.items())
    ]
    return {
        "files": len(modules),
        "lines": sum(module["lines"] for module in modules),
        "over_1000_lines": sorted((m for m in modules if m["lines"] > 1000),
                                  key=lambda m: -m["lines"]),
        "over_500_lines": sorted((m for m in modules if m["lines"] > 500),
                                 key=lambda m: -m["lines"]),
        "largest_files": sorted(modules, key=lambda m: -m["lines"])[:10],
        "module_docstring_coverage": round(
            sum(m["module_docstring"] for m in modules) / max(1, len(modules)), 4),
        "public_docstring_coverage": round(public_documented / max(1, public_total), 4),
        "public_docstring_scope": "top-level definitions excluding tester/tests",
        "package_public_docstrings": {
            package: {"definitions": counts[0], "documented": counts[1],
                      "coverage": round(counts[1] / max(1, counts[0]), 4)}
            for package, counts in sorted(package_public.items())
        },
        "public_definitions": public_total,
        "test_methods": test_methods,
        "parse_errors": parse_errors,
        "package_edges": package_edges,
        "package_cycles": _package_cycles(package_edges),
    }


def _package_cycles(edges: Sequence[dict]) -> list[list[str]]:
    """Return deterministic strongly connected package components."""
    graph: dict[str, set[str]] = defaultdict(set)
    for edge in edges:
        graph[edge["source"]].add(edge["target"])
        graph.setdefault(edge["target"], set())
    index = 0
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    stack: list[str] = []
    active: set[str] = set()
    components: list[list[str]] = []

    def visit(node: str) -> None:
        nonlocal index
        indices[node] = lowlinks[node] = index
        index += 1
        stack.append(node)
        active.add(node)
        for target in sorted(graph[node]):
            if target not in indices:
                visit(target)
                lowlinks[node] = min(lowlinks[node], lowlinks[target])
            elif target in active:
                lowlinks[node] = min(lowlinks[node], indices[target])
        if lowlinks[node] == indices[node]:
            component: list[str] = []
            while True:
                member = stack.pop()
                active.remove(member)
                component.append(member)
                if member == node:
                    break
            if len(component) > 1:
                components.append(sorted(component))

    for package in sorted(graph):
        if package not in indices:
            visit(package)
    return sorted(components)


def _documentation_metrics(root: Path) -> dict:
    pointers = root / "docs" / "LOCAL_ASSET_POINTERS.json"
    local_assets: set[Path] = set()
    if pointers.is_file():
        manifest = json.loads(pointers.read_text(encoding="utf-8"))
        if manifest.get("schema") == "dynsyslab.local-asset-pointers/v1":
            for asset in manifest.get("assets", []):
                candidate = (root / asset["path"]).resolve()
                if (candidate.is_relative_to(root.resolve())
                        and re.fullmatch(r"[0-9a-f]{64}", asset.get("sha256", ""))
                        and isinstance(asset.get("bytes"), int) and asset["bytes"] >= 0):
                    local_assets.add(candidate)
    documents = sorted((root / "docs").rglob("*.md"))
    documents.extend(path for path in (root / "README.md", root / "ARCHITECTURE.md")
                     if path.is_file())
    broken: list[dict] = []
    local_links: list[dict] = []
    for document in documents:
        text = document.read_text(encoding="utf-8", errors="replace")
        for target in LINK.findall(text):
            clean = target.split("#", 1)[0]
            if not clean or re.match(r"^[a-zA-Z]+://", clean):
                continue
            path = (document.parent / clean).resolve()
            if path in local_assets:
                local_links.append({"document": _relative(document, root), "target": target})
            elif not path.exists():
                broken.append({"document": _relative(document, root), "target": target})
    index = root / "docs" / "README.md"
    index_text = index.read_text(encoding="utf-8", errors="replace")
    indexed = {(index.parent / target.split("#", 1)[0]).resolve()
               for target in LINK.findall(index_text)
               if target and not re.match(r"^[a-zA-Z]+://", target)}
    unindexed = [_relative(path, root) for path in documents
                 if path.is_relative_to(root / "docs") and path != index
                 and path.resolve() not in indexed]
    return {"documents": len(documents), "broken_links": broken,
            "local_asset_links": local_links,
            "unindexed_documents": sorted(unindexed)}


def _native_metrics(root: Path) -> dict:
    """Measure maintained C/C++ sources while excluding native build trees."""
    paths: set[Path] = set()
    for base in (root / "native", root / "native_ui", root / "viewer" / "native",
                 root / "compute_workbench" / "native"):
        if not base.is_dir():
            continue
        for pattern in ("*.cpp", "*.hpp", "*.h"):
            paths.update(path for path in base.rglob(pattern)
                         if "build" not in path.relative_to(base).parts)
        cmake = base / "CMakeLists.txt"
        if cmake.is_file():
            paths.add(cmake)
    files = [{"path": _relative(path, root),
              "lines": len(path.read_text(encoding="utf-8", errors="replace").splitlines())}
             for path in sorted(paths)]
    return {
        "files": len(files),
        "lines": sum(item["lines"] for item in files),
        "over_1000_lines": sorted((item for item in files if item["lines"] > 1000),
                                  key=lambda item: -item["lines"]),
        "over_500_lines": sorted((item for item in files if item["lines"] > 500),
                                 key=lambda item: -item["lines"]),
        "largest_files": sorted(files, key=lambda item: -item["lines"])[:10],
    }


def _boundary_violations(paths: Iterable[Path], root: Path) -> list[dict]:
    rules = (
        ("viewer dependency from a producer", {"solver", "archives", "archiver", "job_manager"},
         re.compile(r"(?:from|import)\s+viewer(?:\.|\s)")),
        ("low-level archive access from an application package",
         {"campaigner", "compute_workbench", "control", "estimation", "gui",
          "job_manager", "postprocessing", "rom", "solver", "viewer"},
         re.compile(r"(?:from|import)\s+archives(?:\.|\s)")),
        ("ROM dependency inside system definitions", {"dynsyslib"},
         re.compile(r"(?:from|import)\s+rom(?:\.|\s)")),
    )
    violations: list[dict] = []
    for path in paths:
        relative = _relative(path, root)
        package = relative.split("/", 1)[0]
        text = path.read_text(encoding="utf-8", errors="replace")
        for label, sources, pattern in rules:
            if package in sources and pattern.search(text):
                violations.append({"rule": label, "path": relative})
    return violations


def audit_project(root: str | Path = PROJECT_ROOT) -> dict:
    """Return a versioned, serializable audit of one DynSysLab checkout."""
    root = Path(root).expanduser().resolve()
    paths = source_files() if root == PROJECT_ROOT else [
        path for path in root.rglob("*.py") if not any(part in {
            ".git", ".test-deps", "__pycache__", "build", "dist", "tmp",
            "viewer_compute_archive", "viewer_outputs",
        } for part in path.relative_to(root).parts)
    ]
    python = _python_metrics(paths, root)
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    result = {
        "schema": "dynsyslab.architecture-audit/v1",
        "project_root": str(root),
        "policy": {"preferred_file_lines": 500, "maximum_file_lines": 1000,
                   "minimum_public_docstring_coverage": 0.70},
        "python": python,
        "native": _native_metrics(root),
        "documentation": _documentation_metrics(root),
        "boundary_violations": _boundary_violations(paths, root),
        "console_scripts": sorted(project["project"].get("scripts", {})),
    }
    result["passes_strict"] = not (
        python["parse_errors"] or python["over_1000_lines"] or python["package_cycles"] or
        result["native"]["over_1000_lines"] or
        result["documentation"]["broken_links"] or result["boundary_violations"] or
        python["public_docstring_coverage"] < result["policy"]["minimum_public_docstring_coverage"]
    )
    return result


def compare_baseline(result: dict, baseline: dict) -> dict:
    """Report new and resolved debt relative to an explicit checked baseline."""
    if baseline.get("schema") != "dynsyslab.architecture-baseline/v1":
        raise ValueError("unsupported architecture baseline schema")
    current_size = {item["path"] for item in result["python"]["over_1000_lines"]}
    allowed_size = set(baseline.get("allowed_hard_size_paths", []))
    current_boundary = {(item["path"], item["rule"])
                        for item in result["boundary_violations"]}
    allowed_boundary = {(item["path"], item["rule"])
                        for item in baseline.get("allowed_boundary_violations", [])}
    regressions = {
        "hard_size_paths": sorted(current_size - allowed_size),
        "boundary_violations": [
            {"path": path, "rule": rule}
            for path, rule in sorted(current_boundary - allowed_boundary)
        ],
    }
    resolved = {
        "hard_size_paths": sorted(allowed_size - current_size),
        "boundary_violations": [
            {"path": path, "rule": rule}
            for path, rule in sorted(allowed_boundary - current_boundary)
        ],
    }
    passes = not (
        regressions["hard_size_paths"] or regressions["boundary_violations"] or
        result["python"]["parse_errors"] or result["python"]["package_cycles"] or
        result["documentation"]["broken_links"] or
        result["native"]["over_1000_lines"] or
        result["python"]["public_docstring_coverage"] <
        result["policy"]["minimum_public_docstring_coverage"]
    )
    return {"schema": baseline["schema"], "passes": passes,
            "regressions": regressions, "resolved_debt": resolved}


def render_markdown(result: dict) -> str:
    """Render a concise human-readable audit summary."""
    python = result["python"]
    native = result["native"]
    docs = result["documentation"]
    lines = [
        "# DynSysLab architecture audit", "",
        f"Strict status: **{'PASS' if result['passes_strict'] else 'DEBT REMAINS'}**", "",
        f"- Python source: {python['files']} files, {python['lines']:,} lines",
        f"- Maintained C/C++ source: {native['files']} files, {native['lines']:,} lines",
        f"- Hard size violations (>1000): {len(python['over_1000_lines'])}",
        f"- Preferred-size exceptions (>500): {len(python['over_500_lines'])}",
        f"- Native files over 500/1000: {len(native['over_500_lines'])}/{len(native['over_1000_lines'])}",
        f"- Module docstring coverage: {python['module_docstring_coverage']:.1%}",
        f"- Public-definition docstring coverage: {python['public_docstring_coverage']:.1%} "
        f"({python['public_docstring_scope']})",
        f"- Discoverable test methods: {python['test_methods']}",
        f"- Documentation: {docs['documents']} files, {len(docs['broken_links'])} broken links, "
        f"{len(docs['unindexed_documents'])} unindexed files",
        f"- Boundary violations: {len(result['boundary_violations'])}", "",
        f"- Package dependency cycles: {len(python['package_cycles'])}", "",
        "## Hard size violations", "",
    ]
    lines.extend(f"- `{item['path']}` - {item['lines']} lines"
                 for item in python["over_1000_lines"])
    if not python["over_1000_lines"]:
        lines.append("- None")
    lines.extend(["", "## Boundary violations", ""])
    lines.extend(f"- `{item['path']}` - {item['rule']}"
                 for item in result["boundary_violations"])
    if not result["boundary_violations"]:
        lines.append("- None")
    lines.extend(["", "## Package dependency cycles", ""])
    lines.extend(f"- {' → '.join(cycle)}" for cycle in python["package_cycles"])
    if not python["package_cycles"]:
        lines.append("- None")
    lines.extend(["", "## Documentation gaps", ""])
    lines.extend(f"- Unindexed: `{item}`" for item in docs["unindexed_documents"])
    lines.extend(f"- Broken link: `{item['document']}` → `{item['target']}`"
                 for item in docs["broken_links"])
    if not docs["unindexed_documents"] and not docs["broken_links"]:
        lines.append("- None")
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    """Run the audit and optionally fail when strict policy is unmet."""
    parser = argparse.ArgumentParser(prog="dynsyslab-audit",
        description="Measure DynSysLab architecture and documentation policy conformance.")
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--baseline", type=Path,
                        help="allow only debt recorded in a versioned baseline")
    args = parser.parse_args(argv)
    result = audit_project(args.project_root)
    if args.baseline:
        baseline = json.loads(args.baseline.read_text(encoding="utf-8"))
        result["baseline"] = compare_baseline(result, baseline)
    rendered = (json.dumps(result, indent=2, sort_keys=True) + "\n"
                if args.format == "json" else render_markdown(result))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    passes = (result["baseline"]["passes"] if args.baseline
              else result["passes_strict"])
    return 1 if args.strict and not passes else 0


if __name__ == "__main__":
    raise SystemExit(main())
