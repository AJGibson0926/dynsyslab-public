#pragma once
// `dynsyslab-cli hub`: the DynSysLab Hub greeting. On the first login after an install (tracked
// in %LOCALAPPDATA%/DynSysLab/hub.json per installed version) it includes the welcome and the
// overview of the suite; later logins get a short greeting. It always ends with the options and
// "What would you like to do?". Markdown by default, dynsyslab.hub/v1 JSON with --format json.
//
// --hook is the Claude Code SessionStart hook form: it reads the hook's JSON from stdin and, in a
// DynSysLab workspace (a folder holding .dynsyslab-workspace) or on the first login after an
// install, prints hookSpecificOutput.additionalContext with the greeting; otherwise nothing.
// --peek shows the greeting without recording a login.
#include <filesystem>
#include <string>

namespace dynsyslab::cli {

struct HubOptions {
    std::filesystem::path archiveRoot;
    std::string format{"markdown"};
    bool hook{};
    bool peek{};
};

int hub(const HubOptions& options);

// The marker file agent-setup writes into the user's DynSysLab workspace folder.
inline constexpr const wchar_t* kWorkspaceMarker = L".dynsyslab-workspace";

}  // namespace dynsyslab::cli
