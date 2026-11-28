// Orbit records against a system's known spectrum (OrbitAverages.hpp).
#include "OrbitAverages.hpp"

#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;

namespace {
std::vector<double> firstState(const SourceData& record) {
    if (!record.samples()) throw std::invalid_argument("an orbit record has no samples");
    return {record.values().begin(), record.values().begin() + std::ptrdiff_t(record.dimension())};
}
double invariantAt(const c::ODESystem& system, const std::string& id, const std::vector<double>& x) {
    const auto psi = c::libraryEigenfunctions(system, {id}).front();
    return psi.evaluate(x.data()).front().real();
}
}  // namespace

c::SpectralBand recordBand(const SourceData& record, const std::string& band) {
    const auto system = c::createSystem(systemConfiguration(record));
    const auto bands = c::libraryBands(*system, {band});
    if (bands.empty()) throw std::invalid_argument("the system's catalogue has no band " + band);
    return bands.front();
}

OrbitSeries orbitSeries(const SourceData& record, const c::SpectralBand& band, const std::vector<std::string>& observables,
                        int harmonics, const std::string& invariant) {
    if (!band.eigenfunctions.levelSet()) throw std::invalid_argument("orbit series need a band with level-set generalized eigenfunctions");
    if (harmonics < 0 || harmonics > 100000) throw std::invalid_argument("orbit series need 0..100000 harmonics");
    const auto system = c::createSystem(systemConfiguration(record));
    const auto catalog = system->observableCatalog();
    std::vector<const c::Observable*> g;
    for (const auto& name : observables) {
        const auto& o = catalog.at(name);
        if (o.outputDimension() != 1) throw std::invalid_argument("orbit series observables are scalar: " + name);
        g.push_back(&o);
    }
    const std::size_t K = record.samples(), d = record.dimension(), H = 2 * std::size_t(harmonics) + 1;
    if (d != system->descriptor().stateDimension) throw std::invalid_argument("orbit series need records of the full state");
    const auto x0 = firstState(record);
    OrbitSeries out;
    out.frequency = band.eigenfunctions.level->evaluate(x0.data()).front().real();
    if (!invariant.empty()) out.invariant = invariantAt(*system, invariant, x0);
    std::vector<std::vector<std::complex<double>>> sums(g.size(), std::vector<std::complex<double>>(H));
    double total = 0;
    for (std::size_t k = 1; k + 1 < K; ++k) {  // the window vanishes at both ends
        const double s = double(k) / double(K - 1), w = std::exp(-1 / (s * (1 - s)));
        if (w == 0) continue;
        const double* x = record.values().data() + k * d;
        auto z = band.eigenfunctions.functions->evaluate(x).front();
        if (!(std::abs(z) > 0)) continue;
        z = std::conj(z / std::abs(z));  // e^{-i phi}
        total += w;
        for (std::size_t j = 0; j < g.size(); ++j) {
            const auto value = g[j]->evaluate(x).front();
            // e^{-i n phi} for n = -N..N: start at n = -N (e^{+i N phi}) and step by e^{-i phi}.
            std::complex<double> power = std::pow(std::conj(z), harmonics);
            for (std::size_t n = 0; n < H; ++n, power *= z) sums[j][n] += w * value * power;
        }
    }
    if (!(total > 0)) throw std::invalid_argument("an orbit record is too short for a weighted average");
    for (std::size_t j = 0; j < g.size(); ++j) {
        math::FourierSeries series{out.frequency, {}};
        for (std::size_t n = 0; n < H; ++n) series.coefficients[int(n) - harmonics] = sums[j][n] / total;
        out.channels.push_back(std::move(series));
    }
    return out;
}

std::vector<double> recordWeights(const json::Value& rule, const std::vector<SourceData>& records) {
    const std::size_t n = records.size();
    std::vector<double> w;
    if (rule.isNull() || rule.string() == "uniform") w.assign(n, 1.0);
    else if (std::holds_alternative<json::Value::Array>(rule.storage())) {
        for (const auto& v : rule.array()) w.push_back(v.number());
        if (w.size() != n) throw std::invalid_argument("one weight per record");
    } else {
        const auto& area = rule.at("canonical_area");
        if (area.isNull()) throw std::invalid_argument("record weights: uniform, a list, or {\"canonical_area\": {\"invariant\", \"band\"}}");
        std::vector<double> energy(n), omega(n);
        for (std::size_t j = 0; j < n; ++j) {
            const auto system = c::createSystem(systemConfiguration(records[j]));
            const auto x0 = firstState(records[j]);
            energy[j] = invariantAt(*system, area.at("invariant").string(), x0);
            omega[j] = std::abs(recordBand(records[j], area.at("band").string()).eigenfunctions.level->evaluate(x0.data()).front().real());
        }
        if (n < 2) throw std::invalid_argument("canonical-area weights need two records or more");
        std::vector<std::size_t> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return energy[a] < energy[b]; });
        w.assign(n, 0.0);
        for (std::size_t k = 0; k < n; ++k) {
            const double lower = k ? energy[order[k]] - energy[order[k - 1]] : 0;
            const double upper = k + 1 < n ? energy[order[k + 1]] - energy[order[k]] : 0;
            w[order[k]] = (lower + upper) / 2 / omega[order[k]];
        }
    }
    double total = 0;
    for (const double v : w) {
        if (!(v >= 0) || !std::isfinite(v)) throw std::invalid_argument("record weights must be finite and nonnegative");
        total += v;
    }
    if (!(total > 0)) throw std::invalid_argument("record weights sum to zero");
    for (auto& v : w) v /= total;
    return w;
}

}  // namespace dynsyslab::workflow::postprocess_detail
