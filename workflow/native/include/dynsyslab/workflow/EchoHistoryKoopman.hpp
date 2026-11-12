#pragma once

#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow {

// Fits a causal echo-history coordinate map and a polynomial Koopman
// eigenfunction across independent, unforced training trajectories. The
// returned native archive is directly consumable by Viewer and by the
// approximate-control-map path.
[[nodiscard]] std::filesystem::path runEchoHistoryKoopman(
    const std::vector<std::filesystem::path>& trainingSources,
    const json::Value& settings,
    const std::filesystem::path& archiveRoot,
    const std::string& label);

}  // namespace dynsyslab::workflow
