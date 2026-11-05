#pragma once
#include "PostprocessInternal.hpp"
#include <optional>
namespace dynsyslab::workflow::postprocess_detail {
// Reads measured arrays directly. Analytic evaluation is an explicit simulator
// adapter used only when the requested measurement array is absent.
SourceData loadObservationSource(const std::filesystem::path& path,
    const std::string& arrayName,const std::string& analyticObservable = {});
std::optional<dynsyslab::core::SystemConfiguration> observationSystem(const SourceData& source);
void requireSameObservationSystem(const SourceData& expected,const SourceData& actual);
}
