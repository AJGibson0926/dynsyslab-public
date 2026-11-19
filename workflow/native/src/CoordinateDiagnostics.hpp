#pragma once
// Optional diagnostics recorded with a coordinate change (settings
// "diagnostics": true or {"crossings":{...},"closure":{...}}):
//
//  crossings  Self-intersections of the path in two plotted coordinates
//             ("axes", default [0,1]). Each crossing is lifted into all the
//             coordinates and, when available, the physical state. A crossing
//             whose lifted separation is below "tolerance" (fraction of the
//             RMS radius) is genuine in these coordinates; if the physical
//             states differ it witnesses a non-injective map, and if the
//             coordinate velocities also differ it obstructs any single-valued
//             law z' = f(z). Promoted from the miscellaneous-computes side-car.
//  closure    For an analytic map of the physical state, the native observable
//             closure search (CoordinateClosure) over the record's state box,
//             seeded with the trajectory's own witnesses.
//
// Summaries go to metadata "coordinate_diagnostics"; every crossing is kept in
// the arrays crossing_segments and crossing_metrics.
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/Observable.hpp>
#include <optional>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {
void appendCoordinateDiagnostics(const SourceData& source, const json::Value& settings,
    const dynsyslab::core::CoordinateSeries& series, const std::optional<std::vector<double>>& states,
    const dynsyslab::core::Observable* analyticMap, json::Value& meta, dynsyslab::core::NpzRecord& archive);
} // namespace dynsyslab::workflow::postprocess_detail
