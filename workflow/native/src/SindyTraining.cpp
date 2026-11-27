// Training rows for Koopman-with-SINDy (SindyTraining.hpp): one trajectory, several pooled with their
// boundaries kept, or a state-space domain whose rates come from the system's vector field (the
// generator equation needs only the rate at each state, Kaiser, Kutz & Brunton, Mach. Learn.: Sci.
// Technol. 2:035023, 2021).
#include "SindyTraining.hpp"

#include <dynsyslab/core/System.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;

SourceData sindyRows(const SourceData& full, const json::Value& settings) {
    const double strideValue = settings.at("sample_stride").number(1);
    if (!std::isfinite(strideValue) || strideValue < 1 || std::floor(strideValue) != strideValue)
        throw std::invalid_argument("SINDy sample_stride must be a positive integer");
    const auto stride = static_cast<std::size_t>(strideValue);
    const auto& window = settings.at("time_window").array();
    if (!window.empty() && window.size() != 2) throw std::invalid_argument("SINDy time_window must be [start, stop]");
    if (stride == 1 && window.empty()) return full;
    const double start = window.empty() ? -INFINITY : window[0].number();
    const double stop = window.empty() ? INFINITY : window[1].number();
    std::vector<double> time, values, inputs;
    const std::size_t m = full.inputDimension(), inputRows = m ? full.inputs().size() / m : 0;
    std::size_t kept = 0;
    for (std::size_t sample = 0; sample < full.samples(); ++sample) {
        if (full.time()[sample] < start || full.time()[sample] > stop) continue;
        if (kept++ % stride) continue;
        time.push_back(full.time()[sample]);
        values.insert(values.end(), full.values().begin() + sample * full.dimension(),
                      full.values().begin() + (sample + 1) * full.dimension());
        // One input row per kept sample; an interval history holds its last interval at the end.
        const std::size_t row = std::min(sample, inputRows ? inputRows - 1 : 0);
        if (inputRows) inputs.insert(inputs.end(), full.inputs().begin() + row * m, full.inputs().begin() + (row + 1) * m);
    }
    auto rows = withSeries(full, time, std::move(values), full.dimension());
    rows.record.inputs.reset(); rows.record.inputsOnIntervals = false;
    if (inputRows) rows.record.inputs = inputSeries(time, std::move(inputs), m, false);
    if (rows.samples() < 3) throw std::invalid_argument("SINDy row selection leaves fewer than three samples");
    return rows;
}

namespace {

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.number());
    return out;
}

void append(SindyTraining& out, const SourceData& rows, bool exactRates, bool includeInputTerms) {
    if (out.dimension && rows.dimension() != out.dimension)
        throw std::invalid_argument("SINDy trajectories must share one state dimension");
    out.dimension = rows.dimension();
    out.states.insert(out.states.end(), rows.values().begin(), rows.values().end());
    out.time.insert(out.time.end(), rows.time().begin(), rows.time().end());
    if (exactRates) {
        const auto rates = exactVelocities(rows, includeInputTerms);
        out.velocities.insert(out.velocities.end(), rates.begin(), rates.end());
    }
    out.segments.push_back(rows.samples());
    out.samples += rows.samples();
    out.sources.push_back(rows.path);
}

// States on a domain with the system's rates at time t; each state its own one-sample segment.
SindyTraining domainTraining(const SourceData& source, const json::Value& settings, std::vector<double> states,
                             std::string input) {
    auto system = c::createSystem(systemConfiguration(source));
    const std::size_t d = system->descriptor().stateDimension;
    if (states.empty() || states.size() % d) throw std::invalid_argument("SINDy domain points need one value per state component");
    SindyTraining out;
    out.input = std::move(input);
    out.dimension = d;
    out.samples = states.size() / d;
    if (out.samples < 3) throw std::invalid_argument("a SINDy domain needs at least three states");
    const double t = settings.at("domain").at("time").number(settings.at("domain_time").number(0.0));
    out.velocities.resize(states.size());
    for (std::size_t k = 0; k < out.samples; ++k) {
        for (std::size_t i = 0; i < d; ++i)
            if (!std::isfinite(states[k*d+i])) throw std::invalid_argument("SINDy domain points must be finite");
        system->rhs(t, states.data() + k * d, out.velocities.data() + k * d);
        out.time.push_back(double(k));  // order only: no pair is differenced
    }
    out.states = std::move(states);
    out.segments.assign(out.samples, 1);
    return out;
}

}  // namespace

json::Value SindyTraining::describe() const {
    json::Value::Array paths, lengths;
    for (const auto& path : sources) paths.emplace_back(std::filesystem::absolute(path).generic_u8string());
    for (const auto length : segments) lengths.emplace_back(double(length));
    json::Value::Object out{{"input", input}, {"samples", double(samples)}, {"state_dimension", double(dimension)},
                            {"trajectory_boundaries_preserved", true}};
    if (input == "trajectories") { out["sources"] = paths; out["trajectory_samples"] = lengths; }
    return out;
}

SindyTraining sindyTraining(const SourceData& source, const json::Value& settings, bool exactRates, bool includeInputTerms) {
    const auto& domain = settings.at("domain");
    const auto& points = settings.at("domain_points");
    if (!domain.isNull() || !points.isNull()) {
        if (!domain.isNull() && !points.isNull()) throw std::invalid_argument("give SINDy domain or domain_points, not both");
        if (!exactRates) throw std::invalid_argument("a SINDy domain has no time series to difference; use gamma_method system_rhs");
        if (!points.isNull()) {
            std::vector<double> states;
            for (const auto& row : points.array()) for (const auto& value : row.array()) states.push_back(value.number());
            return domainTraining(source, settings, std::move(states), "domain_points");
        }
        // Korda & Mezic-style state-space sampling on a regular grid, first coordinate slowest.
        const auto lower = numbers(domain.at("lower")), upper = numbers(domain.at("upper"));
        const double densityValue = domain.at("density").number(20);
        if (lower.empty() || lower.size() != upper.size()) throw std::invalid_argument("SINDy domain needs lower and upper of one length");
        if (!(densityValue >= 3) || std::floor(densityValue) != densityValue)
            throw std::invalid_argument("SINDy domain density must be an integer of at least 3");
        const auto density = std::size_t(densityValue);
        double count = 1;
        for (std::size_t i = 0; i < lower.size(); ++i) {
            if (!(upper[i] > lower[i])) throw std::invalid_argument("SINDy domain upper must exceed lower");
            count *= double(density);
        }
        if (count > 250000) throw std::invalid_argument("SINDy domain grid exceeds 250,000 states; lower its density");
        std::vector<double> states;
        std::vector<std::size_t> index(lower.size(), 0);
        for (std::size_t k = 0; k < std::size_t(count); ++k) {
            for (std::size_t i = 0; i < lower.size(); ++i)
                states.push_back(lower[i] + (upper[i] - lower[i]) * double(index[i]) / double(density - 1));
            for (std::size_t i = lower.size(); i-- > 0;) { if (++index[i] < density) break; index[i] = 0; }
        }
        return domainTraining(source, settings, std::move(states), "domain");
    }
    SindyTraining out;
    append(out, source, exactRates, includeInputTerms);
    const auto& additional = settings.at("additional_sources");
    if (additional.isNull() || additional.array().empty()) {
        out.segments.clear();
        return out;
    }
    out.input = "trajectories";
    for (const auto& path : additional.array()) {
        auto loaded = loadSource(std::filesystem::u8path(path.string()));
        append(out, sindyRows(loaded, settings), exactRates, includeInputTerms);
    }
    return out;
}

}  // namespace dynsyslab::workflow::postprocess_detail
