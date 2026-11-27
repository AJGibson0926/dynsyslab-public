#pragma once
// Per-trajectory kernels (TrajectoryKernelsPostprocess.cpp): affine_dynamics, prehistory, clock_audit.
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runAffineDynamics(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                        const std::filesystem::path& root, const std::string& label);
std::filesystem::path runPrehistory(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                    const std::filesystem::path& root, const std::string& label);
std::filesystem::path runClockAudit(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                    const std::filesystem::path& root, const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
