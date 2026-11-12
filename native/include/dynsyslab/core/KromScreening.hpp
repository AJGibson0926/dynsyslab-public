#pragma once

#include <dynsyslab/core/KoopmanModel.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct KromCandidateScore {
    std::vector<std::size_t> modeIndices;
    double selectionScore{};
    double liftedTargetDisplacement{};
    double worstNormalizedControllability{};
    double minimumExactInputAuthority{};
    double frequencySpread{};
    std::size_t minimumControllabilityRank{};
    std::size_t splitDimension{};
};

struct KromScreeningOptions {
    std::size_t maximumModeCount{4};
    std::size_t maximumVirtualCandidates{5000};
    std::size_t retainedCandidates{100};
    std::size_t seedModeCount{12};
};

// Rank virtual sub-KROMs without fitting a control map or writing an archive.
// Expensive extraction/regression/control campaign accounting starts only
// after this read-only screen.
[[nodiscard]] std::vector<KromCandidateScore> rankKromControlCandidates(
    const KoopmanModel& parent, const std::vector<double>& initialState,
    const std::vector<double>& targetState,
    const KromScreeningOptions& options = {});

}  // namespace dynsyslab::core
