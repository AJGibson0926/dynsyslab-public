#pragma once
// `dynsyslab-cli mcp`: a Model Context Protocol server on stdio (newline-delimited JSON-RPC 2.0)
// for Claude Code, Codex and other agents. Its instructions ask the agent to open the DynSysLab
// Hub at the start of a DynSysLab conversation; its tools (hub, quickstart, list_systems, doctor,
// validate_job, plan_job, run_job, fourier, render_figure) run this CLI and the Viewer CLI as
// child processes, so protocol output on stdout stays clean. The `hub` prompt is the same
// greeting as a slash command.
namespace dynsyslab::cli {

int mcpServer();

}  // namespace dynsyslab::cli
