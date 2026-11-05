"""A minimal, policy-bounded MCP server for local DynSysLab development.

The server intentionally exposes only read-only project-orientation tools.
It does not expose credentials, account bindings, arbitrary shell execution,
or archive mutation.  Future capability-specific tools should retain these
boundaries and route archive operations through Archiver.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Any

from mcp.server.mcpserver import MCPServer


PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SAFE_DOCUMENTS = {
    "README.md",
    "ARCHITECTURE.md",
    "SOLVER.md",
    "docs/AI_CONTEXT.md",
    "docs/AGENT_HARNESS.md",
    "docs/ANTENNA_SECURITY.md",
}

server = MCPServer(
    name="dynsyslab-harness",
    title="DynSysLab Local Harness",
    description="Read-only local orientation tools for the DynSysLab suite.",
    instructions=(
        "Use these tools for project orientation only. They never grant "
        "authority to modify archives, execute arbitrary commands, or reveal "
        "credentials/account data."
    ),
    version="0.1.0",
)


@server.tool(name="harness_status", description="Return safe local harness capability status.")
def harness_status() -> dict[str, Any]:
    """Report local, non-sensitive development capabilities."""
    docker = subprocess.run(
        ["docker", "version", "--format", "{{.Server.Version}}"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return {
        "project": "dynsyslab",
        "project_root": str(PROJECT_ROOT),
        "transport": "stdio",
        "archive_access": "not exposed by this harness",
        "docker_engine": "available" if docker.returncode == 0 else "unavailable",
        "mcp_version_line": "2.x",
        "sensitive_data_policy": "never exposed",
    }


@server.tool(name="project_structure", description="List the top-level DynSysLab source packages only.")
def project_structure() -> dict[str, list[str]]:
    """Return an intentionally small source-level project inventory."""
    packages = sorted(
        item.name
        for item in PROJECT_ROOT.iterdir()
        if item.is_dir() and not item.name.startswith((".", "__"))
        and item.name not in {"archives", "build", "dist", "logs", "test_reports", "viewer_compute_archive"}
    )
    return {"source_packages": packages}


@server.tool(name="read_project_document", description="Read a small, pre-approved public project document.")
def read_project_document(document: str) -> dict[str, str]:
    """Read one allow-listed document without accepting arbitrary file paths."""
    normalized = document.replace("\\", "/")
    if normalized not in _SAFE_DOCUMENTS:
        return {"error": "Document is not in the harness allow-list."}
    source = PROJECT_ROOT / normalized
    if not source.is_file():
        return {"error": "Requested document is not available locally."}
    return {"document": normalized, "content": source.read_text(encoding="utf-8")}


def main() -> None:
    """Run the local server using MCP stdio transport."""
    if os.environ.get("DYNSYSLAB_MCP_TRANSPORT", "stdio") != "stdio":
        raise RuntimeError("This initial harness supports only local stdio transport.")
    server.run("stdio")


if __name__ == "__main__":
    main()
