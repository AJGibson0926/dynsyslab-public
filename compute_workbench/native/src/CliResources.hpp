#pragma once
// Files the user-facing CLI commands (quickstart, hub, mcp) read from an installation or a
// source checkout, and the welcome document they present.
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::cli {

struct Resources {
    std::filesystem::path examples;   // quickstart examples folder
    std::filesystem::path welcome;    // docs/AGENT_WELCOME.md
    std::filesystem::path viewer;     // DynSysLabViewerCli.exe
};

// Installed layout (<install>/resources) first, then the source checkout of a development build.
[[nodiscard]] Resources locateResources();

struct Welcome {
    std::string overview;             // Markdown up to "## Get started"
    std::vector<std::string> next;    // bullets under "## Things you can ask"
};
[[nodiscard]] Welcome readWelcome(const std::filesystem::path& document);

// Typographic punctuation folded to ASCII, so text survives any console code page.
[[nodiscard]] std::string asciiText(std::string text);
// Native Windows form of a path for display.
[[nodiscard]] std::string shownPath(const std::filesystem::path& path);

}  // namespace dynsyslab::cli
