#pragma once
// Builds the H-DMD K-ROM (spectrum and eigenfunctions, core/KoopmanModel.hpp) from a fit.
#include "HankelAnalysis.hpp"
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/HankelRepresentations.hpp>
#include <dynsyslab/core/Spectral.hpp>

namespace dynsyslab::workflow::postprocess_detail {

struct HdmdFitContext {
    const dynsyslab::core::TimeSeriesRecord* record{};  // full input record, all channels
    std::size_t first{};          // first record row of the time window
    std::size_t labelOffset{};    // window label minus window start, in record rows
    const HankelAnalysis* hankel{};
    const json::Value* settings{};
};

// Window layout in record columns, from a prepared Hankel analysis (H-DMD and H-SVD share it).
[[nodiscard]] dynsyslab::core::HankelWindowLayout recordWindowLayout(const HankelAnalysis& hankel,
    const dynsyslab::core::TimeSeriesRecord& record, std::size_t labelOffset);

// Spectrum, samples, origin and representations, ordered by setting "representations":
// "interpolation" (the shared interpolation on setting "interpolation_domain", default the full
// state; SampledDomain.hpp) and "expansion" (the window map).
[[nodiscard]] dynsyslab::core::KoopmanModel hdmdEigenfunctions(
    const dynsyslab::core::HankelDMDResult& result, const HdmdFitContext& context);

}  // namespace dynsyslab::workflow::postprocess_detail
