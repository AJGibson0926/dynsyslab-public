#pragma once
// `dynsyslab-cli agent-setup [--remove]`: connects DynSysLab to the user's coding agents so a plain
// chat session opens with the DynSysLab Hub. For each agent found it:
//   Claude Code  registers the `dynsyslab` MCP server (`claude mcp add --scope user`) and a
//                SessionStart hook (`dynsyslab-cli hub --hook`) in ~/.claude/settings.json;
//   Codex        registers the `dynsyslab` MCP server (`codex mcp add`, or [mcp_servers.dynsyslab]
//                in ~/.codex/config.toml when the codex command is absent);
// and writes AGENTS.md, CLAUDE.md and the .dynsyslab-workspace marker into the DynSysLab workspace
// folder. Every change is reported; --remove undoes them. Entries DynSysLab did not write are kept.
// --home and --no-agent-cli redirect the user folder and skip the agent commands (tests).
#include <filesystem>
#include <string>

namespace dynsyslab::cli {

struct AgentSetupOptions {
    bool remove{};
    bool useAgentCommands{true};
    std::filesystem::path home;       // default: %USERPROFILE%
    std::filesystem::path workspace;  // default: the installer's data folder (runtime.json project_root)
    std::string format{"markdown"};
};

int agentSetup(const AgentSetupOptions& options);

}  // namespace dynsyslab::cli
