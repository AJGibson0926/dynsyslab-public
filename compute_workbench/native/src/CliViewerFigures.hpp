#pragma once
// Headless figures for CLI reports: runs DynSysLabViewerCli (next to dynsyslab-cli) without a
// window and with its own output silenced, on the archive root whose index lists the archive
// (index new archives with core::indexArchiveFiles first).
#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::cli {

struct FigureRequest {
    std::filesystem::path archiveRoot;
    std::filesystem::path archive;
    std::string view;                     // Viewer route: time, phase, fourier, ...
    std::vector<std::wstring> extra;      // further Viewer CLI options
    std::filesystem::path output;         // PNG
};

// Renders one figure; false (and a reason) when the Viewer CLI is missing or fails.
bool renderFigure(const std::filesystem::path& viewerCli, const FigureRequest& request,
                  const std::filesystem::path& cache, std::string& reason);

}  // namespace dynsyslab::cli
