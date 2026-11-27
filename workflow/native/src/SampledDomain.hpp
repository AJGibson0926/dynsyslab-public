#pragma once
// The domain and base points of an observable built from samples of a time series (any method:
// H-DMD, EDMD, regression, ...). The workflow half of core::interpolateSamples: it follows the
// series' origin back to the dynamical system's states and applies the requested coordinate map.
//
// Setting "interpolation_domain":
//   null / "state"   the full state x (default), recovered from the solver record the series points to
//   [names...]       catalog observables of the state, e.g. ["far_field_acoustic_echo", ...]
//   {expression}     an expression definition q(x) (Observable::fromExpressionDefinition)
//   "input"          the series' own channels: computed from x when every channel is a catalog
//                    observable, else the sampled values themselves (measured data)
// Without a producing system only "input" exists, and null falls back to it.
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/TimeSeriesRecord.hpp>

#include <optional>

namespace dynsyslab::workflow::postprocess_detail {

// The full state at each sample time, from the record's state source (unit-normalized, matched by time).
struct SampleStates {
    std::vector<double> values;  // samples x dimension
    std::size_t dimension{};
    dynsyslab::core::SystemConfiguration system;
};
[[nodiscard]] std::optional<SampleStates> statesAtSamples(const dynsyslab::core::TimeSeriesRecord& record,
                                                          const std::vector<double>& times);

// The full state at a K-ROM's training samples, through its origin: its input series, then the
// solver record they point to. Throws when the series has no producing system.
[[nodiscard]] SampleStates trainingStates(const dynsyslab::core::KoopmanModel& set);

struct SampledDomain {
    dynsyslab::core::FunctionDomain domain;  // bound when it has a state map
    std::vector<double> states;              // what interpolateSamples maps: x_k, or q_k for a supplied domain
};
// `rows` are the record rows of the samples; `columns` the record channels the method read (for "input").
[[nodiscard]] SampledDomain sampledDomain(const json::Value& request, const dynsyslab::core::TimeSeriesRecord& record,
                                          const std::vector<std::size_t>& rows, const std::vector<std::size_t>& columns);

}  // namespace dynsyslab::workflow::postprocess_detail
