#include "SampledDomain.hpp"
#include "TimeSeriesSource.hpp"
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

std::optional<SampleStates> statesAtSamples(const c::TimeSeriesRecord& record, const std::vector<double>& times) {
    if (!record.stateSource) return std::nullopt;
    const auto source = loadSource(std::filesystem::u8path(record.stateSource->path));
    SampleStates out;
    out.dimension = source.dimension();
    out.system = systemConfiguration(source);
    out.values.reserve(times.size() * source.dimension());
    for (const double t : times) {
        const double tolerance = 1e-9 * std::max(1.0, std::abs(t));
        const auto found = std::lower_bound(source.time().begin(), source.time().end(), t - tolerance);
        const auto row = std::size_t(found - source.time().begin());
        if (row >= source.samples() || std::abs(source.time()[row] - t) > tolerance)
            throw std::invalid_argument("state source " + record.stateSource->path + " has no sample at a training time");
        out.values.insert(out.values.end(), source.values().begin() + std::ptrdiff_t(row * source.dimension()),
                          source.values().begin() + std::ptrdiff_t((row + 1) * source.dimension()));
    }
    return out;
}

SampleStates trainingStates(const c::KoopmanModel& set) {
    // An interpolation on the physical state stores the states at the samples themselves.
    for (const auto& rep : set.representations())
        if (const auto* fi = dynamic_cast<const c::FunctionInterpolation*>(rep.function.get());
            fi && fi->coordinates().isPhysicalState() && fi->data().samples == set.samples().samples() && set.samples().samples()) {
            SampleStates out;
            out.values = fi->points();
            out.dimension = fi->coordinateDimension();
            if (set.hasPhysicalSystem()) out.system = set.systemConfiguration();
            return out;
        }
    std::vector<std::filesystem::path> inputs;
    for (const auto& input : set.origin().inputs)
        if (input.role == "time_series") inputs.push_back(std::filesystem::u8path(input.path));
    if (inputs.empty()) throw std::invalid_argument("these eigenfunctions record no input time series");
    auto states = statesAtSamples(loadTimeSeriesRecords(inputs), set.samples().time);
    if (!states) throw std::invalid_argument("the eigenfunctions' input series has no producing system states");
    return std::move(*states);
}

namespace {

c::FunctionDomain catalogDomain(const std::vector<std::string>& names, const SampleStates& states) {
    const auto system = c::createSystem(states.system);
    const auto catalog = system->observableCatalog();
    std::size_t dimension = 0;
    for (const auto& name : names) dimension += catalog.at(name).outputDimension();
    return c::FunctionDomain::fromCatalog(names, dimension, &catalog);
}

SampledDomain inputDomain(const c::TimeSeriesRecord& record, const std::vector<std::size_t>& rows,
                          const std::vector<std::size_t>& columns, const std::optional<SampleStates>& states) {
    std::vector<std::string> names;
    bool fromState = states.has_value();
    for (const auto column : columns) {
        const auto& observable = record.channels.at(column).observable;
        names.push_back(observable);
        fromState = fromState && !observable.empty() && observable != "state" && observable.find("hsvd[") != 0;
    }
    if (fromState) {
        try {
            auto domain = catalogDomain(names, *states);
            if (domain.dimension == columns.size()) return {std::move(domain), states->values};
        } catch (const std::exception&) {}  // not catalog observables: fall through to the sampled values
    }
    // Measured channels: the sampled values themselves are the base points (stored by the interpolation).
    SampledDomain out{c::FunctionDomain::supplied(columns.size()), {}};
    out.domain.names = names;
    for (const auto row : rows)
        for (const auto column : columns) out.states.push_back(record.series.value(row, column));
    return out;
}

}  // namespace

SampledDomain sampledDomain(const J& request, const c::TimeSeriesRecord& record,
                            const std::vector<std::size_t>& rows, const std::vector<std::size_t>& columns) {
    std::vector<double> times;
    for (const auto row : rows) times.push_back(record.series.time(row));
    const auto states = statesAtSamples(record, times);
    const bool isString = std::holds_alternative<std::string>(request.storage());
    const auto name = isString ? request.string() : std::string();
    if ((request.isNull() && !states) || name == "input") return inputDomain(record, rows, columns, states);
    if (!states)
        throw std::invalid_argument("interpolation_domain needs the producing system's states; this series has none");
    if (request.isNull() || name == "state") {
        return {c::FunctionDomain::physicalState(states->dimension), states->values};
    }
    if (std::holds_alternative<J::Array>(request.storage())) {
        std::vector<std::string> names;
        for (const auto& item : request.array()) names.push_back(item.string());
        if (names.empty()) throw std::invalid_argument("interpolation_domain observables must not be empty");
        return {catalogDomain(names, *states), states->values};
    }
    if (std::holds_alternative<J::Object>(request.storage()))
        return {c::FunctionDomain::fromStateMap(c::Observable::fromExpressionDefinition(request)), states->values};
    throw std::invalid_argument("interpolation_domain must be null, \"state\", \"input\", observable names or an expression");
}

}  // namespace dynsyslab::workflow::postprocess_detail
