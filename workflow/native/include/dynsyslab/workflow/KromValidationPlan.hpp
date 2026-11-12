#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <filesystem>

namespace dynsyslab::workflow {
// Read-only reuse-first lookup. Resolves physics, never fits or simulates.
// If nothing matches, returns an explicit offer requiring IC/horizon approval.
[[nodiscard]] dynsyslab::core::json::Value planAutonomousKoopmanValidation(
    const dynsyslab::core::KoopmanModel& model,
    const std::filesystem::path& modelPath,
    const std::filesystem::path& archiveRoot,
    std::size_t maximumDatasets = 5);
[[nodiscard]] dynsyslab::core::json::Value planAutonomousKoopmanValidation(
    const std::filesystem::path& modelPath,
    const std::filesystem::path& archiveRoot,
    std::size_t maximumDatasets = 5);
}
