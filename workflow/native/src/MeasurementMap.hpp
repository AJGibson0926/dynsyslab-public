#pragma once
// One loader for the frozen measurement maps an eigenfunction domain may name (its measurement map):
// history-basis archives and H-SVD archives alike. Every consumer (controllers, input maps, target
// lifts, certificates, audits) goes through here, so the identity check is the same everywhere.
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/HankelBasis.hpp>
#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

namespace dynsyslab::workflow::postprocess_detail {

// The causal projector for a Hankel basis: single channel, newest-sample labels, uncentered.
// Energy scaling is represented by unit singular values (q = U^T h), whitened by Sigma (xi).
[[nodiscard]] dynsyslab::core::FrozenHistoryTransform historyTransform(const dynsyslab::core::HankelBasis& basis,
    double sampleInterval, dynsyslab::core::HankelCoordinateScaling scaling);

// Reads a history-transform archive or an H-SVD archive.
[[nodiscard]] dynsyslab::core::FrozenHistoryTransform readMeasurementMap(const std::filesystem::path& path);

struct MeasurementMap {
    std::string path;
    dynsyslab::core::FrozenHistoryTransform transform;
};

// The map named by the model's domain (or settings "history_transform"), checked against the
// domain's identity. A legacy domain without a sha256 identity needs settings
// "unverified_transform_reason".
[[nodiscard]] MeasurementMap loadMeasurementMap(const dynsyslab::core::KoopmanModel& model,
    const std::filesystem::path& kromPath, const json::Value& settings);

}  // namespace dynsyslab::workflow::postprocess_detail
