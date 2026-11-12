#pragma once
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/KoopmanModel.hpp>

namespace dynsyslab::workflow::postprocess_detail {
std::string datasetDynamics(const SourceData& source);
void requireMatchingUnforcedDataset(const dynsyslab::core::KoopmanModel& model,
                                    const SourceData& source);
// Select existing samples only; never interpolates or edits a source archive.
SourceData validationSamples(const SourceData& source,
    const dynsyslab::core::KoopmanModel& model,const json::Value& settings);
}
