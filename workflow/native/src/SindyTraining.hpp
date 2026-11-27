#pragma once
// The rows a Koopman-with-SINDy fit trains on (SindyEigenfunctions.cpp):
//   one trajectory   the step's source (default);
//   trajectories     the source and "additional_sources": [archives...], pooled with their
//                    boundaries kept (finite differences never cross one);
//   a domain         "domain": {"lower": [...], "upper": [...], "density": n, "time": t}, a grid of
//                    states, or "domain_points": [[...], ...]; rates from the source system's vector
//                    field at time t (no time series, so gamma_method must be system_rhs).
// Every trajectory takes the same "sample_stride" and "time_window" selection.
#include "PostprocessInternal.hpp"

#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

struct SindyTraining {
    std::vector<double> states, time, velocities;  // samples x dimension; velocities empty: finite differences
    std::size_t samples{}, dimension{};
    std::vector<std::size_t> segments;  // sample counts per trajectory, or ones for a domain; empty: one trajectory
    std::string input{"trajectory"};    // "trajectory", "trajectories", "domain" or "domain_points"
    std::vector<std::filesystem::path> sources;  // the trajectories read (empty for a domain)
    [[nodiscard]] bool oneTrajectory() const noexcept { return segments.empty(); }
    [[nodiscard]] json::Value describe() const;
};

// "sample_stride" and "time_window" applied to one record.
[[nodiscard]] SourceData sindyRows(const SourceData& full, const json::Value& settings);

// The training rows of `source` (already row-selected) and the settings. `exactRates` evaluates the
// system's vector field (gamma_method system_rhs) instead of differencing.
[[nodiscard]] SindyTraining sindyTraining(const SourceData& source, const json::Value& settings, bool exactRates,
                                          bool includeInputTerms);

}  // namespace dynsyslab::workflow::postprocess_detail
