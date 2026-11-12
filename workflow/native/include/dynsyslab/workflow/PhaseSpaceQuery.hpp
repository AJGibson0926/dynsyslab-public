#pragma once
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/PhaseSpace.hpp>
#include <filesystem>
namespace dynsyslab::workflow {
namespace json = dynsyslab::core::json;
json::Value phaseSpaceOptionsJson(const dynsyslab::core::PhaseSpaceOptions &options);
dynsyslab::core::PhaseSpaceOptions phaseSpaceOptionsFromJson(const json::Value &value,
                                                             std::size_t dimension);
// mode=lookup never computes or writes. mode=compute first returns a matching cache.
json::Value evaluatePhaseSpaceQuery(const json::Value &request,
                                    const std::filesystem::path &archiveRoot,
                                    const dynsyslab::core::ExecutionContext *execution = nullptr);
} // namespace dynsyslab::workflow
