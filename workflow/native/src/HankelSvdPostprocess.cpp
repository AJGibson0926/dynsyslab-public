// H-SVD (hankel_svd / hankel_pod): time-series records in, a Hankel basis plus a coordinate
// time series out. The output archive is itself a time-series record (time, values = coordinates),
// so it can feed H-DMD, SINDy or a controller directly. Broomhead & King, Physica D 20 (1986).
#include "HankelAnalysis.hpp"
#include "HdmdObjects.hpp"
#include "MeasurementMap.hpp"
#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/HankelBasis.hpp>
#include <dynsyslab/core/PCA.hpp>

#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

std::filesystem::path runHankelSvd(const std::string& method, const SourceData& full, const J& settings,
                                   const std::filesystem::path& root, const std::string& label) {
    const auto view = joinedTimeSeries(full, settings);
    const auto& record = view.record;
    std::size_t first = 0;
    const SourceData source = hankelTimeWindow(view, settings, first);
    auto input = prepareHankel(source, settings, false, true);
    if (!settings.at("time_window").isNull()) input.provenance.object()["time_window"] = settings.at("time_window");
    const bool center = settings.at("center").boolean(false);
    const double rankValue = settings.at("rank").number(0);
    if (!(rankValue >= 0) || rankValue != std::floor(rankValue)) throw std::invalid_argument("hankel_svd rank must be a count");
    const auto scaling = c::hankelCoordinateScaling(settings.at("coordinate_scaling").string("energy"));
    c::NpzRecord archive;
    archiveHankelLayout(archive, input, settings.at("save_matrix").boolean());
    const auto result = c::fitHankelPOD(input.matrix, std::size_t(rankValue), center);
    const auto& svd = result.decomposition;
    const std::size_t rank = svd.rank, rows = result.embeddedFeatures, windows = result.embeddedSamples;

    // U (window-space directions) is the right factor of the SVD of H^T.
    std::vector<double> modes(rows * rank);
    for (std::size_t k = 0; k < rank; ++k)
        for (std::size_t r = 0; r < rows; ++r) modes[r * rank + k] = svd.right[k * rows + r];
    const std::size_t labelOffset = input.alignEnd ? input.matrix.ends.at(0) - input.matrix.starts.at(0) : 0;
    c::DataOrigin origin{method, settings, record.origin.inputs};
    const c::HankelBasis basis(recordWindowLayout(input, record, labelOffset), modes, svd.singularValues,
                               center ? result.mean : std::vector<double>{}, origin);

    // Coordinates at every window, labelled by its record row.
    const auto& offsets = input.alignEnd ? input.matrix.ends : input.matrix.starts;
    std::vector<double> time, labels;
    std::vector<double> coordinates(windows * rank);
    for (std::size_t j = 0; j < windows; ++j) {
        time.push_back(source.time().at(offsets[j]));
        labels.push_back(double(first + offsets[j]));
        const double weight = scaling == c::HankelCoordinateScaling::energy ? 1.0 : 0.0;
        for (std::size_t k = 0; k < rank; ++k)  // left column j is V; q = sigma v, xi = v
            coordinates[j * rank + k] = svd.left[j * svd.rank + k] * (weight ? svd.singularValues[k] : 1.0);
    }
    archive.arrays["time"] = c::NpyArray::fromDoubles({windows}, time);
    archive.arrays["values"] = c::NpyArray::fromDoubles({windows, rank}, coordinates);
    archive.arrays["label_rows"] = c::NpyArray::fromDoubles({windows}, labels);
    archive.arrays["basis_modes"] = c::NpyArray::fromDoubles({rows, rank}, modes);
    archive.arrays["singular_values"] = c::NpyArray::fromDoubles({rank}, svd.singularValues);
    archive.arrays["mean"] = c::NpyArray::fromDoubles({rows}, result.mean);
    // Legacy names: the SVD is of H^T, so "left" vectors are temporal.
    archive.arrays["left_singular_vectors"] = c::NpyArray::fromDoubles({svd.rows, rank}, svd.left);
    archive.arrays["right_singular_vectors"] = c::NpyArray::fromDoubles({rank, svd.columns}, svd.right);
    archive.arrays["pod_modes"] = archive.arrays["basis_modes"];
    std::vector<double> scores(svd.left);
    for (std::size_t j = 0; j < windows; ++j)
        for (std::size_t k = 0; k < rank; ++k) scores[j * rank + k] *= svd.singularValues[k];
    archive.arrays["scores"] = c::NpyArray::fromDoubles({windows, rank}, scores);

    J::Array names;
    for (std::size_t k = 0; k < rank; ++k) names.emplace_back("hsvd[" + std::to_string(k) + "]");
    auto meta = json::parse(metadata(source, method, label,
        "\"hankel_window\":" + json::dump(input.provenance)
        + ",\"pod_modes_array\":\"basis_modes\",\"temporal_vectors_array\":\"left_singular_vectors\""
        + ",\"temporal_scores_array\":\"scores\",\"centered\":" + std::string(center ? "true" : "false")
        + ",\"rank\":" + std::to_string(rank)));
    auto block = basis.describe();
    block.object()["coordinate_scaling"] = c::hankelCoordinateScalingName(scaling);
    block.object()["arrays"] = J::Object{{"modes", "basis_modes"}, {"singular_values", "singular_values"},
                                         {"mean", "mean"}, {"coordinates", "values"}, {"label_rows", "label_rows"}};
    const double sampleInterval = uniformDeltaTime(source) * double(input.matrix.window.sampleStride);
    block.object()["sample_interval"] = sampleInterval;
    // A causal single-channel basis is a measurement map: fits on these coordinates (e.g. coordinate
    // SINDy) record its identity, and controllers project live windows through it.
    try {
        const auto map = historyTransform(basis, sampleInterval, scaling);
        meta.object()["coordinate_provider"] = "causal-history";
        meta.object()["coordinate_transform_identity"] = map.identity();
        block.object()["measurement_map_identity"] = map.identity();
    } catch (const std::invalid_argument&) {
        block.object()["measurement_map_identity"] = J();  // not causal / multichannel / centered
    }
    meta.object()["hankel_basis"] = std::move(block);
    // The output is a time series of history coordinates labelled by window.
    meta.object()["value_space"] = "history-coordinates";
    meta.object()["coordinate_transform_kind"] = "hankel_basis";
    meta.object()["coordinate_names"] = names;
    meta.object()["state_source_archive"] = record.stateSource ? J(record.stateSource->path) : J();
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(meta));
    return publishObject(root, source, dynsyslab::core::object_kind::timeSeries, method, std::move(archive));  // coordinates (with their basis)
}

}  // namespace dynsyslab::workflow::postprocess_detail
