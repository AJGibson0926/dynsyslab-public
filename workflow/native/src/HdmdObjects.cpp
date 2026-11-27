#include "HdmdObjects.hpp"
#include "SampledDomain.hpp"
#include "TimeSeriesSource.hpp"

#include <algorithm>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

std::vector<std::string> orderedNames(const J& requested, std::vector<std::string> available, const char* what) {
    if (requested.isNull()) return available;
    std::vector<std::string> out;
    for (const auto& item : requested.array()) {
        const auto name = item.string();
        if (std::find(available.begin(), available.end(), name) == available.end())
            throw std::invalid_argument(std::string("hdmd ") + what + " " + name + " is not available here");
        out.push_back(name);
    }
    if (out.empty()) throw std::invalid_argument(std::string("hdmd ") + what + " must not be empty");
    return out;
}

}  // namespace

c::HankelWindowLayout recordWindowLayout(const HankelAnalysis& hankel, const c::TimeSeriesRecord& record,
                                        std::size_t labelOffset) {
    c::HankelWindowLayout layout;
    layout.window = hankel.matrix.window;
    layout.alignEnd = hankel.alignEnd;
    layout.labelOffset = labelOffset;
    for (const auto& channel : hankel.provenance.at("channels").array()) {
        const auto column = std::size_t(channel.at("measurement_column").number());
        layout.channels.push_back({column, std::size_t(channel.at("lag_samples").number()), channel.at("scale").number(1.0)});
        layout.channelObservables.push_back(record.channels.at(column).observable);
    }
    return layout;
}

c::KoopmanModel hdmdEigenfunctions(const c::HankelDMDResult& result, const HdmdFitContext& context) {
    const auto& record = *context.record;
    const auto& settings = *context.settings;
    const auto& hankel = *context.hankel;
    const std::size_t modes = result.rank, n = result.samples;

    c::DataOrigin origin{"hdmd", settings, {}};
    for (const auto& input : record.origin.inputs) origin.inputs.push_back(input);

    auto samples = std::make_shared<c::EigenfunctionSamples>();
    samples->modes = modes;
    samples->values = result.coordinates;
    for (std::size_t k = 0; k < n; ++k) {
        const auto row = context.first + result.coordinateSampleOffsets.at(k);
        samples->sampleIndices.push_back(row);
        samples->time.push_back(record.series.time(row));
    }

    // Window layout in record columns (the Hankel builder indexed its own channel list).
    const auto layout = recordWindowLayout(hankel, record, context.labelOffset);
    std::vector<std::size_t> columns;
    for (const auto& channel : layout.channels)
        if (std::find(columns.begin(), columns.end(), channel.input) == columns.end()) columns.push_back(channel.input);

    // The interpolation is the shared one (SampledDomain.hpp, core::interpolateSamples): by default on
    // the full state the input series' origin leads to; H-DMD stores no copy of its input.
    J request = settings.at("interpolation_domain");
    if (request.isNull() && !settings.at("interpolation_domains").isNull())  // former setting: its first domain
        request = settings.at("interpolation_domains").array().at(0);
    const auto sampled = sampledDomain(request, record, samples->sampleIndices, columns);
    std::size_t neighbors = std::size_t(settings.at("interpolation_neighbors").number(8));
    const auto method = settings.at("interpolation_method").isNull() ? c::InterpolationMethod{}
        : c::interpolationMethodFromJson(settings.at("interpolation_method"), &neighbors);
    std::shared_ptr<const c::FunctionRepresentation> interpolation = c::interpolateSamples(sampled.states,
        n, sampled.domain, samples->values, method, neighbors, "hdmd eigenfunctions at their training samples");
    // The time-delay representation: the window map T = G^-1 Phi^H, an expansion on the window domain.
    std::shared_ptr<const c::FunctionRepresentation> window = c::windowMapFromModes(layout, result.modes, result.modeGram, modes);
    std::vector<std::shared_ptr<const c::FunctionRepresentation>> representations;
    J requested = settings.at("representations");
    if (!requested.isNull()) {  // former names of the window map
        J::Array names;
        for (const auto& item : requested.array())
            names.emplace_back(item.string() == "window_map" || item.string() == "time-delay" ? std::string("expansion") : item.string());
        requested = names;
    }
    for (const auto& name : orderedNames(requested, {"interpolation", "expansion"}, "representation"))
        representations.push_back(name == "interpolation" ? interpolation : window);

    auto spectrum = result.spectrum;
    spectrum.withOrigin(origin);
    return c::KoopmanModel(std::move(spectrum), *samples, std::move(origin), std::move(representations));
}

}  // namespace dynsyslab::workflow::postprocess_detail
