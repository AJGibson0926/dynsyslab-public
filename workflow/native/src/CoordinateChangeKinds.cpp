#include "CoordinateChangeKinds.hpp"
#include "ObjectPublishing.hpp"
#include "CoordinateChannels.hpp"
#include "KromDatasetSupport.hpp"
#include "ObservationSource.hpp"
#include <dynsyslab/core/HistoryKoopman.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
namespace c = dynsyslab::core;
using J = json::Value;

std::size_t positive(const J& s, const char* key, std::size_t fallback) {
    const double n = s.at(key).number(double(fallback));
    if (!std::isfinite(n) || n < 1 || n != std::floor(n) || n > double(std::numeric_limits<std::size_t>::max() / 2))
        throw std::invalid_argument(std::string("Invalid positive count: ") + key);
    return static_cast<std::size_t>(n);
}

bool physical(const SourceData& s) {
    return (s.metadataJson.empty() ? J{} : json::parse(s.metadataJson)).at("coordinate_provider").isNull();
}

std::vector<double> decimate(const std::vector<double>& x, std::size_t stride) {
    std::vector<double> out;
    for (std::size_t k = 0; k < x.size(); k += stride) out.push_back(x[k]);
    return out;
}

std::string unit(const SourceData& source, const J& settings, const char* setting, const char* key, const char* fallback) {
    const auto meta = source.metadataJson.empty() ? J{} : json::parse(source.metadataJson);
    const auto declared = meta.at(key).string();
    const auto selected = settings.at(setting).string(declared.empty() ? fallback : declared);
    if (!declared.empty() && declared != selected)
        throw std::invalid_argument(std::string("Observation unit override conflicts with source: ") + setting);
    return selected;
}
}  // namespace

HistoryMeasurement historyMeasurement(const std::filesystem::path& path, const J& settings, bool withRates) {
    const auto& selector = settings.at("measurement");
    std::string name = selector.isNull() ? settings.at("observable").string() : selector.string();
    auto source = loadSource(path);
    HistoryMeasurement out;
    if (physical(source)) {
        if (name.empty()) throw std::invalid_argument("History coordinates need a measurement observable");
        out.series = loadObservationSource(path, settings.at("measurement_array").string(name), name);
        if (withRates) {  // analytic push-forward of the unforced vector field
            if (datasetDynamics(source) != "unforced")
                throw std::invalid_argument("Exact measurement rates require confirmed unforced data");
            const auto map = c::createSystem(systemConfiguration(source))->observableCatalog().at(name);
            const auto q = c::TimeSeries::view(source.time, source.values.data(), source.samples, source.dimension);
            const c::TimeSeries velocity(source.time, exactVelocities(source, false), source.dimension);
            const auto mapped = c::transformCoordinates(q, map, &velocity);
            out.rates = std::vector<double>(mapped.derivatives->valuesData(), mapped.derivatives->valuesData() + source.samples);
        }
        out.name = name;
        return out;
    }
    // A column of an earlier coordinate change, e.g. the echo of (e, e').
    const auto columns = timeSeriesColumnNames(source);
    std::size_t k = 0;
    if (std::holds_alternative<double>(selector.storage())) k = static_cast<std::size_t>(selector.number());
    else {
        const auto found = std::find(columns.begin(), columns.end(), name);
        if (found == columns.end()) throw std::invalid_argument("History measurement is not a column of the coordinate series: " + name);
        k = std::size_t(found - columns.begin());
    }
    if (k >= source.dimension) throw std::invalid_argument("History measurement column is out of range");
    out.name = columns[k];
    if (withRates) {
        const auto record = c::readNpz(path);
        if (!record.contains("coordinate_derivatives"))
            throw std::invalid_argument("Exact measurement rates need a coordinate archive with stored derivatives (system_derivatives)");
        const auto stored = record.at("coordinate_derivatives").doubles();
        std::vector<double> rates(source.samples);
        for (std::size_t i = 0; i < source.samples; ++i) rates[i] = stored[i * source.dimension + k];
        out.rates = std::move(rates);
    }
    std::vector<double> column(source.samples);
    for (std::size_t i = 0; i < source.samples; ++i) column[i] = source.values[i * source.dimension + k];
    source.values = std::move(column);
    source.dimension = 1;
    out.series = std::move(source);
    return out;
}

std::filesystem::path historyTransformArchive(const std::filesystem::path& path) {
    // A coordinate archive built with a fitted basis names it in its metadata.
    const auto record = c::readNpz(path);
    if (record.contains("left_singular_vectors")) return path;
    const auto meta = record.metadataJson().empty() ? J{} : json::parse(record.metadataJson());
    const auto named = meta.at("coordinate_transform_archive").string();
    if (named.empty()) throw std::invalid_argument("transform_archive is neither a history basis nor history coordinates");
    return std::filesystem::u8path(named);
}

// Observables are functions of the physical state. A coordinate-change archive (e.g. the identity
// coordinates a pipeline feeds to H-DMD) stores the physical states at its samples as source_states,
// so observable channels are evaluated on those; only a series without physical states is refused.
SourceData physicalStates(const SourceData& source) {
    if (physical(source)) return source;
    const auto record = c::readNpz(source.path);
    if (!record.contains("source_states"))
        throw std::invalid_argument("Coordinate channels need a physical state source "
                                    "(or a coordinate archive that stores source_states)");
    const auto& states = record.at("source_states");
    if (states.shape.size() != 2 || states.shape[0] != source.samples)
        throw std::invalid_argument("source_states are not aligned with the coordinate samples");
    SourceData out = source;
    out.values = states.doubles(); out.dimension = states.shape[1];
    auto meta = source.metadataJson.empty() ? J{} : json::parse(source.metadataJson);
    if (std::holds_alternative<J::Object>(meta.storage())) {
        meta.object().erase("coordinate_provider"); meta.object().erase("coordinate_names");
    }
    out.metadataJson = json::dump(meta);
    return out;
}

c::CoordinateSeries channelCoordinates(const SourceData& given, const J& settings, J& meta) {
    const auto source = physicalStates(given);
    if (settings.at("system_derivatives").boolean(false))
        throw std::invalid_argument("Channel coordinates carry values only; use kind observable for exact rates");
    const auto channels = evaluateCoordinateChannels(source, settings.at("channels").array());
    const std::size_t width = channels.size();
    std::vector<double> values(source.samples * width);
    J::Array names, provenance;
    for (std::size_t k = 0; k < width; ++k) {
        for (std::size_t i = 0; i < source.samples; ++i) values[i * width + k] = channels[k].values[i];
        names.emplace_back(channels[k].name);
        provenance.emplace_back(J::Object{{"name", channels[k].name}, {"base", channels[k].base},
            {"derivative_source", channels[k].derivative}, {"units", channels[k].units}});
    }
    meta.object()["coordinate_channels"] = std::move(provenance);
    if (settings.at("names").isNull()) meta.object()["coordinate_names"] = std::move(names);
    return {c::TimeSeries(source.time, std::move(values), width), std::nullopt};
}

FittedHistory fitHistoryCoordinates(const SourceData& source, const J& settings,
                                    const std::filesystem::path& root, const std::string& label) {
    const bool withRates = settings.at("system_derivatives").boolean(false);
    auto primary = historyMeasurement(source.path, settings, withRates);
    const auto name = primary.name;
    std::vector<SourceData> records{primary.series};
    for (const auto& extra : settings.at("fit_sources").array()) {
        auto other = historyMeasurement(std::filesystem::u8path(extra.string()), settings, false);
        if (other.name != name) throw std::invalid_argument("Pooled history records measure different observables");
        records.push_back(std::move(other.series));
    }
    const auto stride = positive(settings, "decimation", 1);
    const double dt = uniformDeltaTime(records.front()) * double(stride);
    std::vector<std::vector<double>> observations;
    std::vector<c::KoopmanDatasetReference> references;
    for (const auto& record : records) {
        if (record.dimension != 1) throw std::invalid_argument("History coordinates require a scalar measurement");
        if (datasetDynamics(record) != "unforced")
            throw std::invalid_argument("History basis training requires confirmed unforced observations");
        requireSameObservationSystem(records.front(), record);
        if (std::abs(uniformDeltaTime(record) * double(stride) - dt) > 1e-8 * dt)
            throw std::invalid_argument("Pooled history clocks differ");
        observations.push_back(decimate(record.values, stride));
        references.push_back({std::filesystem::absolute(record.path).generic_u8string(), "training", "unforced", {}});
    }
    c::PooledHankelOptions options;
    options.delays = positive(settings, "delays", 100);
    options.rank = positive(settings, "hankel_rank", 4);
    options.oversampling = positive(settings, "hankel_oversampling", 4);
    options.subspaceIterations = positive(settings, "hankel_subspace_iterations", 8);
    options.seed = positive(settings, "seed", 20261002);
    auto basis = c::fitPooledHankelBasis(observations, options);
    const auto rank = positive(settings, "coordinate_rank", options.rank);
    if (rank > basis.rank) throw std::invalid_argument("Coordinate rank exceeds the fitted basis");
    std::vector<std::size_t> selection;
    for (std::size_t j = 0; j < rank; ++j) selection.push_back(j);
    const auto& first = records.front();
    c::FrozenHistoryTransform transform(std::move(basis), dt, selection, name,
        unit(first, settings, "time_units", "measurement_time_units", "model-time"),
        unit(first, settings, "observation_units", "measurement_units", "unspecified"));
    auto meta = json::parse(metadata(first, "history-transform", label, ""));
    meta.object()["dataset_references"] = c::koopmanDatasetReferencesJson(references);
    meta.object()["history_decimation"] = double(stride);
    const auto directory = objectFolder(root, "coordinate-transform");
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(source.path.stem().u8string() + "__history-transform__" + timestamp()) + ".npz");
    c::saveFrozenHistoryTransform(path, transform, json::dump(meta));
    const auto measurement = c::TimeSeries(decimate(first.time, stride), decimate(first.values, stride), 1);
    if (!primary.rates) return {c::transformHistoryCoordinates(measurement, transform), path, transform};
    const c::TimeSeries rates(decimate(first.time, stride), decimate(*primary.rates, stride), 1);
    return {c::transformHistoryCoordinates(measurement, transform, &rates), path, transform};
}

std::optional<std::vector<double>> alignedSourceStates(const SourceData& source, const c::TimeSeries& coordinates) {
    // A physical record supplies its own states; a coordinate archive passes on
    // the physical states it carries, so chained coordinate changes keep them.
    std::vector<double> carried;
    std::size_t width = source.dimension;
    if (!physical(source)) {
        const auto archive = c::readNpz(source.path);
        if (!archive.contains("source_states")) return std::nullopt;
        carried = archive.at("source_states").doubles();
        width = archive.at("source_states").shape.at(1);
    }
    const auto& values = physical(source) ? source.values : carried;
    std::vector<double> states;
    states.reserve(coordinates.sampleCount() * width);
    std::size_t row = 0;
    for (std::size_t k = 0; k < coordinates.sampleCount(); ++k) {
        const double t = coordinates.timeData()[k];
        while (row < source.samples && source.time[row] < t) ++row;
        if (row == source.samples || source.time[row] != t) return std::nullopt;
        states.insert(states.end(), values.begin() + std::ptrdiff_t(row * width),
                      values.begin() + std::ptrdiff_t((row + 1) * width));
    }
    return states;
}

} // namespace dynsyslab::workflow::postprocess_detail

namespace dynsyslab::workflow::postprocess_detail {
ResolvedCoordinateInput analyticCoordinateInput(const SourceData& given, const json::Value& settings) {
    using J = json::Value;
    const auto meta = given.metadataJson.empty() ? J{} : json::parse(given.metadataJson);
    const auto kind = meta.at("coordinate_transform_kind").string();
    if (meta.at("method").string() != "coordinate-transform" || meta.at("value_space").string() != "observation-coordinates"
        || !meta.at("acoustic_receiver").isNull() || (kind != "identity" && kind != "observable" && kind != "channels"))
        return {given, settings};
    J::Array names;
    if (kind == "channels") {
        for (const auto& channel : meta.at("coordinate_channels").array()) {
            if (channel.at("derivative_source").string() != "none") return {given, settings};
            names.emplace_back(channel.at("base"));
        }
    }
    auto resolved = settings;
    const auto requested = settings.at("coordinate_observable").string();
    const auto fixed = kind == "identity" ? std::string("state") : meta.at("coordinate_observable").string();
    if ((!requested.empty() && requested != fixed) || !settings.at("coordinate_observables").array().empty()
        || !settings.at("coordinate_expression").isNull())
        throw std::invalid_argument("The coordinate change already fixes the coordinates; remove coordinate_observable from this method");
    if (kind == "channels") resolved.at("coordinate_observables") = names;
    else if (!fixed.empty()) resolved.at("coordinate_observable") = fixed;
    else resolved.at("coordinate_expression") = meta.at("coordinate_definition");
    auto physicalSource = loadSource(std::filesystem::u8path(meta.at("source_archive").string()));
    if (physicalSource.time != given.time)
        throw std::invalid_argument("Coordinate archive samples differ from its physical source");
    return {std::move(physicalSource), std::move(resolved)};
}
} // namespace dynsyslab::workflow::postprocess_detail

namespace dynsyslab::workflow::postprocess_detail {
dynsyslab::core::CoordinateSeries selectedColumns(const SourceData& source, const json::Value& settings,
                                                  const json::Value& original, json::Value& meta) {
    namespace c = dynsyslab::core;
    const auto names = timeSeriesColumnNames(source);
    std::vector<std::size_t> keep;
    for (const auto& item : settings.at("columns").array()) {
        std::size_t k = names.size();
        if (std::holds_alternative<double>(item.storage())) k = std::size_t(item.number());
        else k = std::size_t(std::find(names.begin(), names.end(), item.string()) - names.begin());
        if (k >= source.dimension) throw std::invalid_argument("columns names a coordinate the input does not have");
        keep.push_back(k);
    }
    if (keep.empty()) throw std::invalid_argument("columns needs a nonempty list");
    const std::size_t n = source.samples, d = source.dimension, w = keep.size();
    std::vector<double> values(n * w);
    for (std::size_t r = 0; r < n; ++r) for (std::size_t j = 0; j < w; ++j) values[r * w + j] = source.values[r * d + keep[j]];
    std::optional<c::TimeSeries> rates;
    const auto record = c::readNpz(source.path);
    if (record.contains("coordinate_derivatives")) {
        const auto stored = record.at("coordinate_derivatives").doubles();
        std::vector<double> picked(n * w);
        for (std::size_t r = 0; r < n; ++r) for (std::size_t j = 0; j < w; ++j) picked[r * w + j] = stored[r * d + keep[j]];
        rates = c::TimeSeries(source.time, std::move(picked), w);
    }
    for (const char* key : {"coordinate_provider", "value_space", "coordinate_transform_identity", "coordinate_transform_archive",
                            "history_span", "warmup_samples", "derivative_evidence", "controller_kind"})
        if (!original.at(key).isNull()) meta.object()[key] = original.at(key);
    json::Value::Array picked, labels;
    const auto& sourceLabels = original.at("coordinate_labels").array();
    for (const auto k : keep) {
        picked.emplace_back(names[k]);
        labels.emplace_back(sourceLabels.size() == d ? sourceLabels[k].string() : names[k]);
    }
    if (settings.at("names").isNull()) meta.object()["coordinate_names"] = std::move(picked);
    meta.object()["coordinate_labels"] = std::move(labels);
    return {c::TimeSeries(source.time, std::move(values), w), rates};
}
} // namespace dynsyslab::workflow::postprocess_detail
