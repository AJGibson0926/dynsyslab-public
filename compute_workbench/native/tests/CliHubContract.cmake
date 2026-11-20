# DynSysLab Hub: first login shows the overview and options and is recorded per version; the
# start-up hook stays quiet outside a workspace after the first login and greets inside one;
# agent-setup writes and removes the Claude hook, the Codex MCP entry and the workspace files while
# keeping the user's own entries; the MCP server answers initialize, tools/list and tools/call.
file(REMOVE_RECURSE "${ROOT}")
file(MAKE_DIRECTORY "${ROOT}/local" "${ROOT}/home/.claude" "${ROOT}/home/.codex" "${ROOT}/ws" "${ROOT}/elsewhere")
set(env ${CMAKE_COMMAND} -E env "LOCALAPPDATA=${ROOT}/local")

function(require text expected what)
  string(FIND "${text}" "${expected}" at)
  if(at EQUAL -1)
    message(FATAL_ERROR "${what}: missing '${expected}' in:\n${text}")
  endif()
endfunction()

execute_process(COMMAND ${env} "${CLI}" hub --archive-root "${ROOT}/archives" OUTPUT_VARIABLE first RESULT_VARIABLE code)
require("${first}" "# Welcome to DynSysLab" "first login")
require("${first}" "## What would you like to do?" "first login")
require("${first}" "1. **Run the quickstart**" "first login")
file(READ "${ROOT}/local/DynSysLab/hub.json" state)
string(JSON welcomed GET "${state}" welcomed_version)
if(welcomed STREQUAL "")
  message(FATAL_ERROR "hub did not record the first login: ${state}")
endif()

# Past the first-login window: a short greeting, and the hook is quiet outside a workspace.
string(JSON state SET "${state}" first_login_unix 0)
file(WRITE "${ROOT}/local/DynSysLab/hub.json" "${state}")
execute_process(COMMAND ${env} "${CLI}" hub --peek --archive-root "${ROOT}/archives" OUTPUT_VARIABLE later)
require("${later}" "# Welcome back to DynSysLab" "later login")
file(WRITE "${ROOT}/elsewhere.json" "{\"cwd\":\"${ROOT}/elsewhere\"}")
execute_process(COMMAND ${env} "${CLI}" hub --hook --archive-root "${ROOT}/archives"
  INPUT_FILE "${ROOT}/elsewhere.json" OUTPUT_VARIABLE quiet)
string(STRIP "${quiet}" quiet)
if(NOT quiet STREQUAL "")
  message(FATAL_ERROR "hook greeted outside a workspace after the first login: ${quiet}")
endif()

file(WRITE "${ROOT}/home/.codex/config.toml" "model = \"x\"\n\n[mcp_servers.other]\ncommand = \"y\"\n")
execute_process(COMMAND ${env} "${CLI}" agent-setup --home "${ROOT}/home" --workspace "${ROOT}/ws" --no-agent-cli
  RESULT_VARIABLE code OUTPUT_VARIABLE setup ERROR_VARIABLE error)
if(NOT code EQUAL 0)
  message(FATAL_ERROR "agent-setup failed: ${setup}${error}")
endif()
file(READ "${ROOT}/home/.claude/settings.json" settings)
require("${settings}" "hub --hook" "Claude start-up hook")
file(READ "${ROOT}/home/.codex/config.toml" codex)
require("${codex}" "[mcp_servers.dynsyslab]" "Codex MCP entry")
require("${codex}" "[mcp_servers.other]" "Codex keeps the user's entries")
foreach(file AGENTS.md CLAUDE.md .dynsyslab-workspace)
  if(NOT EXISTS "${ROOT}/ws/${file}")
    message(FATAL_ERROR "workspace file missing: ${file}")
  endif()
endforeach()
file(WRITE "${ROOT}/ws.json" "{\"cwd\":\"${ROOT}/ws\"}")
execute_process(COMMAND ${env} "${CLI}" hub --hook --archive-root "${ROOT}/archives"
  INPUT_FILE "${ROOT}/ws.json" OUTPUT_VARIABLE hook)
string(JSON context GET "${hook}" hookSpecificOutput additionalContext)
require("${context}" "What would you like to do?" "hook in the workspace")

file(WRITE "${ROOT}/mcp.jsonl"
  "{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"initialize\",\"params\":{\"protocolVersion\":\"2025-06-18\"}}\n"
  "{\"jsonrpc\":\"2.0\",\"method\":\"notifications/initialized\"}\n"
  "{\"jsonrpc\":\"2.0\",\"id\":2,\"method\":\"tools/list\"}\n"
  "{\"jsonrpc\":\"2.0\",\"id\":3,\"method\":\"tools/call\",\"params\":{\"name\":\"hub\",\"arguments\":{}}}\n")
execute_process(COMMAND ${env} "${CLI}" mcp INPUT_FILE "${ROOT}/mcp.jsonl" OUTPUT_VARIABLE mcp RESULT_VARIABLE code)
require("${mcp}" "\"instructions\"" "MCP initialize")
require("${mcp}" "\"render_figure\"" "MCP tools/list")
require("${mcp}" "What would you like to do?" "MCP hub tool")

execute_process(COMMAND ${env} "${CLI}" agent-setup --remove --home "${ROOT}/home" --workspace "${ROOT}/ws" --no-agent-cli
  RESULT_VARIABLE code OUTPUT_QUIET)
file(READ "${ROOT}/home/.claude/settings.json" settings)
file(READ "${ROOT}/home/.codex/config.toml" codex)
string(FIND "${settings}${codex}" "dynsyslab-cli" left)
if(NOT left EQUAL -1 OR EXISTS "${ROOT}/ws/AGENTS.md")
  message(FATAL_ERROR "agent-setup --remove left DynSysLab entries")
endif()
require("${codex}" "[mcp_servers.other]" "removal keeps the user's entries")
