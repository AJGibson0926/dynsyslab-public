#pragma once
// Headless figures for CLI reports: runs DynSysLabViewerCli (next to dynsyslab-cli) without a
// window, with its own output silenced, against the directory that holds the archive so the
// Viewer catalogs only that folder.
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::cli {

struct FigureRequest {
    std::filesystem::path archive;
    std::string view;                     // Viewer route: time, phase, fourier, ...
    std::vector<std::wstring> extra;      // further Viewer CLI options
    std::filesystem::path output;         // PNG
};

// Renders one figure; false (and a reason) when the Viewer CLI is missing or fails.
bool renderFigure(const std::filesystem::path& viewerCli, const FigureRequest& request,
                  const std::filesystem::path& cache, std::string& reason);

}  // namespace dynsyslab::cli
