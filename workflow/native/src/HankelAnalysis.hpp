#pragma once
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/Hankel.hpp>

namespace dynsyslab::workflow::postprocess_detail {
struct HankelAnalysis {
    dynsyslab::core::HankelMatrix matrix;
    std::vector<std::string> names;
    std::vector<double> scales;
    bool alignEnd{};
    json::Value provenance;
};
void archiveHankelLayout(dynsyslab::core::NpzRecord&,const HankelAnalysis&,bool saveMatrix);
HankelAnalysis prepareHankel(const SourceData&, const json::Value& settings, bool paired, bool windowHandled = false);
// Contiguous source rows with time in settings "time_window" [start, stop] (all rows when absent).
// "first" receives the index of the first kept row in the full source.
SourceData hankelTimeWindow(const SourceData&, const json::Value& settings, std::size_t& first);
HankelAnalysis readHankel(const SourceData&, const std::filesystem::path& archive, bool requirePair);
// The Hankel matrix as a data object (layout, matrices, label times, source hash), for publishing.
dynsyslab::core::NpzRecord hankelMatrixObject(const SourceData&,const HankelAnalysis&,const std::string& label);
std::filesystem::path runHankelMatrix(const SourceData&,const json::Value&,const std::filesystem::path& root,const std::string& label);
} // namespace dynsyslab::workflow::postprocess_detail
