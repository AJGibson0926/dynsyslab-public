#pragma once
// Coordinate change: evaluate named observables of a physical state record as
// scalar time-series channels. This is the only place post-processing turns
// physics (catalog observables, their time derivatives, KROM observables) into
// time series; the Hankel module only selects columns of what it is given.
#include "PostprocessInternal.hpp"
#include <string>
#include <vector>

namespace dynsyslab::workflow::postprocess_detail {

struct CoordinateChannel {
    std::string name;       // column name in the coordinate archive
    std::string base;       // base observable
    std::string derivative; // "none" or the derivative method
    std::string units;
    std::vector<double> values; // one per source sample
};

// A channel spec is a state index, an observable name, an object
// {"observable": name, "derivative_order": k} or a KROM observable reference
// {"artifact", "model_identity", "component", "coordinate_contract", ...}.
// Duplicate specs are evaluated once. The source must be a physical state record.
std::vector<CoordinateChannel> evaluateCoordinateChannels(const SourceData& source,
                                                          const json::Value::Array& specs);

// True when a Hankel channel selector needs physics: a name that is not a
// column of the source, a time derivative, or a KROM observable reference.
bool selectorNeedsCoordinateChange(const json::Value& item, const std::vector<std::string>& columns);

// Column names of a time-series archive: its recorded coordinate names, else
// the system's state names for a physical record, always aliased by state_k.
std::vector<std::string> timeSeriesColumnNames(const SourceData& source);

} // namespace dynsyslab::workflow::postprocess_detail
