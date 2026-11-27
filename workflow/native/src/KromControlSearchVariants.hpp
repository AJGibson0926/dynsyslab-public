#pragma once
// krom_control_search schedule, closure ordering, resume and map variants
// (KromControlSearchVariants.cpp).
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/KromScreening.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::search_detail {

struct MapVariant {
    std::string name;
    std::filesystem::path archive;
    double relativeResidual{};
    std::string failure;
};

[[nodiscard]] std::vector<double> controlledStates(const std::filesystem::path& source, std::size_t& samples);
[[nodiscard]] std::vector<dynsyslab::core::KromCandidateScore> scheduleCandidates(const json::Value& job,
    std::vector<dynsyslab::core::KromCandidateScore> ranked, const dynsyslab::core::KoopmanModel& parent,
    const std::vector<double>& states, std::size_t samples, json::Value& evidence);
[[nodiscard]] json::Value::Array resumedAttempts(const json::Value& job);
[[nodiscard]] bool attempted(const json::Value::Array& attempts, const std::vector<std::size_t>& indices);
[[nodiscard]] std::vector<MapVariant> fitMapVariants(const json::Value& job, const std::filesystem::path& extracted,
    const std::vector<double>& controlledStates, std::size_t controlledSamples, const std::string& priorExactArchive,
    const std::filesystem::path& directory, std::size_t attempt);

}  // namespace dynsyslab::workflow::search_detail
