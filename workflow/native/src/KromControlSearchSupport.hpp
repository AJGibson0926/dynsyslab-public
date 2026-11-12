#pragma once

#include <dynsyslab/core/KromScreening.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace dynsyslab::workflow::search_detail {

struct CaseMetrics {
    double releaseError{};
    double finalError{};
    double controlPeak{};
    double activeRmse{};
};

[[nodiscard]] std::vector<double> numbers(const json::Value& value);
[[nodiscard]] json::Value values(const std::vector<double>& source);
[[nodiscard]] json::Value indices(const std::vector<std::size_t>& source);
[[nodiscard]] json::Value candidateJson(
    const dynsyslab::core::KromCandidateScore& candidate);
[[nodiscard]] json::Value extractionJob(const json::Value& search,
    const std::filesystem::path& parent,const std::vector<std::size_t>& modes,
    std::size_t attempt);
[[nodiscard]] json::Value regressionJob(const json::Value& search,
    const std::filesystem::path& controlled,const std::filesystem::path& krom,
    std::size_t attempt);
[[nodiscard]] json::Value solveJob(const json::Value& search,
    const dynsyslab::core::KoopmanModel& parent,const std::filesystem::path& krom,
    const std::string& representation,std::size_t attempt);
[[nodiscard]] CaseMetrics metrics(const std::filesystem::path& solution,
    double releaseTime,const std::vector<double>& target,
    const std::vector<double>& finalState,double settlingStart);
[[nodiscard]] json::Value metricsJson(const CaseMetrics& source);
void writeJsonAtomic(const std::filesystem::path& path,const json::Value& value);

}  // namespace dynsyslab::workflow::search_detail
