#pragma once
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/CoordinateTransform.hpp>
#include <dynsyslab/core/HistoryObservation.hpp>
#include <filesystem>
#include <optional>

namespace dynsyslab::workflow::postprocess_detail {

// The scalar a history coordinate change embeds: a physics observable of a
// physical record, or one column of an earlier coordinate archive. With
// withRates, its exact time derivative too (analytic push-forward, or the
// archive's stored coordinate_derivatives column).
struct HistoryMeasurement { SourceData series; std::optional<std::vector<double>> rates; std::string name; };
HistoryMeasurement historyMeasurement(const std::filesystem::path& path, const json::Value& settings, bool withRates);
// transform_archive may be the basis itself or history coordinates built with it.
std::filesystem::path historyTransformArchive(const std::filesystem::path& path);

// kind "channels": evaluate a list of observables (or derivatives, KROM
// observables) of a physical record as one multi-column coordinate series.
dynsyslab::core::CoordinateSeries channelCoordinates(const SourceData& source, const json::Value& settings,
                                                     json::Value& meta);

// kind "history" without a transform_archive: fit the Hankel-SVD history basis
// on this record's scalar measurement (pooled with any "fit_sources") through
// the shared Hankel module, save it beside the coordinates, and project.
struct FittedHistory {
    dynsyslab::core::CoordinateSeries series;
    std::filesystem::path transformArchive;
    dynsyslab::core::FrozenHistoryTransform transform;
};
FittedHistory fitHistoryCoordinates(const SourceData& source, const json::Value& settings,
                                    const std::filesystem::path& root, const std::string& label);

// Physical states at each output time, when the source is a physical record
// and every output time is one of its samples. Lets downstream methods that
// need the physical state (H-DMD reconstruction, interpolated eigenfunctions)
// work from a coordinate archive.
std::optional<std::vector<double>> alignedSourceStates(const SourceData& source,
                                                       const dynsyslab::core::TimeSeries& coordinates);
} // namespace dynsyslab::workflow::postprocess_detail

namespace dynsyslab::workflow::postprocess_detail {
// A model fit handed an analytic coordinate-change archive (identity, a named
// observable map, or plain observable channels) fits on the physical record
// with that exact map, so the resulting model keeps an analytic coordinate
// Jacobian. Any other input is returned unchanged.
struct ResolvedCoordinateInput { SourceData source; json::Value settings; };
ResolvedCoordinateInput analyticCoordinateInput(const SourceData& given, const json::Value& settings);
} // namespace dynsyslab::workflow::postprocess_detail

namespace dynsyslab::workflow::postprocess_detail {
// kind "columns": keep a subset of an earlier coordinate archive's columns
// ("columns": indices or names), with their stored derivatives, labels and the
// provenance of the coordinates they came from (e.g. (xi_1, xi_2) of xi_1..xi_4).
dynsyslab::core::CoordinateSeries selectedColumns(const SourceData& source, const json::Value& settings,
                                                  const json::Value& original, json::Value& meta);
} // namespace dynsyslab::workflow::postprocess_detail

namespace dynsyslab::workflow::postprocess_detail {
// kind "resample": the source on a new clock ("times" or "time_step", optional "start"/"stop")
// with "method" pchip (default), linear, cubic or nearest (CoordinateResample.cpp).
dynsyslab::core::CoordinateSeries resampledCoordinates(const SourceData& source, const json::Value& settings,
                                                       json::Value& meta);
} // namespace dynsyslab::workflow::postprocess_detail

namespace dynsyslab::workflow::postprocess_detail {
// kind "delay_taps": delay coordinates at physical-time lags, linear or as quintic delay jets
// (CoordinateDelayTaps.cpp, core/DelaySampling.hpp).
dynsyslab::core::CoordinateSeries delayTapCoordinates(const SourceData& source, const json::Value& settings, json::Value& meta);
} // namespace dynsyslab::workflow::postprocess_detail
