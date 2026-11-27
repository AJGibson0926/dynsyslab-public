#pragma once
// Pathwise diagnostics of an approximate K-ROM input map B_map(x) against the exact map
// B_exact(x) = grad psi(x) . g(x) along a state path (e.g. a controlled run): per point the
// relative error |B_map - B_exact| / |B_exact|, the direction cosine |<B_map, B_exact>| /
// (|B_map| |B_exact|), the norm ratio |B_map| / |B_exact|, and the truth-coupled frozen poles:
// the KLQR gain K designed on (Lambda, B_map) at that point applied through B_exact,
// eig(Lambda - B_exact K) in the real split coordinates. A gate passes when the worst relative
// error, the smallest cosine and the fraction of points with stable truth-coupled poles meet
// their thresholds (frozen-coefficient analysis of a scheduled design: Shamma & Athans 1990,
// "Analysis of gain scheduled control for nonlinear plants", IEEE TAC 35:898-907).
#include <dynsyslab/core/KoopmanModel.hpp>

#include <string>
#include <vector>

namespace dynsyslab::core {

struct ControlMapPathOptions {
    std::string controlMap{"constant"};
    double statePenalty{5000.0};
    double inputPenalty{0.01};
    double maximumRelativeError{0.25};
    double minimumCosine{0.9};
    double minimumStableFraction{1.0};
    std::size_t stride{1};
};

struct ControlMapPathResult {
    std::vector<double> relativeError, cosine, normRatio, maximumPoleReal;
    std::vector<std::size_t> sampleIndices;
    double stableFraction{};
    bool passed{};
    std::vector<std::string> failures;
};

[[nodiscard]] ControlMapPathResult analyzeControlMapPath(const KoopmanModel& model,
    const std::vector<double>& states, std::size_t samples, const ControlMapPathOptions& options);

}  // namespace dynsyslab::core
