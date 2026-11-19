#pragma once
// `dynsyslab-cli quickstart`: the first-run experience as one command. It prints the DynSysLab
// welcome and overview (docs/AGENT_WELCOME.md), runs the bundled quickstart examples
// (resources/examples/quickstart), summarises each run with a Fourier analysis, renders
// figures with the Viewer CLI, and ends with things to try next. Output is Markdown, or
// dynsyslab.quickstart/v1 JSON with --format json. It never starts a desktop application.
#include <filesystem>
#include <string>

namespace dynsyslab::cli {

struct QuickstartOptions {
    std::filesystem::path archiveRoot;
    std::string archiveRootSource;
    std::string format{"markdown"};
};

// Prints the report to stdout; returns the process exit code.
int quickstart(const QuickstartOptions& options);

}  // namespace dynsyslab::cli
