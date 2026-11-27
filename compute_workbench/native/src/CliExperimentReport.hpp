#pragma once
// dynsyslab-cli experiment report: an experiment manifest as a self-contained HTML or PDF report
// (stage graph, method tables, results, Viewer plots), rendered by DynSysLabViewerCli report.
#include "CliArguments.hpp"

#include <dynsyslab/core/Json.hpp>

#include <filesystem>

namespace dynsyslab::cli {

// The dynsyslab.report/v1 spec of an experiment manifest (dynsyslab.experiment/v1).
[[nodiscard]] core::json::Value experimentReportSpec(const core::json::Value& manifest,
                                                     const std::filesystem::path& manifestPath);
int experimentCommand(const Arguments& args);

}  // namespace dynsyslab::cli
