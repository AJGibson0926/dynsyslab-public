// Further mode-selection rules for krom_extract (settings.mode_selection.rule):
//
//  "state_observability": greedy grouped selection of the modes that best reconstruct the state
//    by least squares from [Re psi, Im psi] (conjugate pairs are one group), the forward greedy
//    regression of Efroymson 1960 ("Multiple regression analysis", in Mathematical Methods for
//    Digital Computers, Wiley). Keys: r2_threshold (0.99), max_modes (groups, 8),
//    minimum_relative_improvement (1e-4), control_authority_weight (0: scores weighted by
//    1 + w |B_psi| at the system's initial state, normalised), states_source (an archive with
//    "values"; default the K-ROM's training_states). Each step's R2 is reported.
//  "nearest_frequency": modes whose continuous rate lies in "band" [low, high] of |Im lambda|,
//    nearest first to "target_rate" [re, im] (or |Im| = "target_frequency"); "count" groups.
//  "validation_score": modes ranked by the anchored relative error of psi_j(x(t)) against
//    exp(lambda_j (t - t0)) psi_j(x(t0)) on "validation_sources" (Koopman eigenfunction
//    validation, Mezic 2005, Nonlinear Dyn. 41:309-325); optional "band"; "count" groups.
// Conjugate partners are always completed.
#include "KromModeRules.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {

namespace c = dynsyslab::core;
using Complex = std::complex<double>;

std::vector<Complex> rates(const c::KoopmanModel& model) {
    if (!model.spectrum().hasContinuousRates())
        throw std::invalid_argument("this mode_selection rule needs continuous eigenvalues");
    return model.spectrum().continuousRates();
}

// Groups of mode indices: a nonreal mode and its stored conjugate, or one real mode.
std::vector<std::vector<std::size_t>> groups(const std::vector<Complex>& spectrum, double tolerance = 1e-9) {
    std::vector<std::vector<std::size_t>> result;
    std::vector<bool> used(spectrum.size(), false);
    for (std::size_t j = 0; j < spectrum.size(); ++j) {
        if (used[j]) continue;
        used[j] = true;
        std::vector<std::size_t> group{j};
        if (std::abs(spectrum[j].imag()) > tolerance)
            for (std::size_t k = j + 1; k < spectrum.size(); ++k)
                if (!used[k] && std::abs(spectrum[k] - std::conj(spectrum[j])) <= 1e-9 + 1e-7 * std::abs(spectrum[j])) {
                    used[k] = true; group.push_back(k); break;
                }
        result.push_back(std::move(group));
    }
    return result;
}

std::pair<std::vector<double>, std::size_t> loadStates(const std::string& path, const char* array,
                                                       std::size_t dimension) {
    const auto archive = c::readNpz(c::wide(path));
    if (!archive.contains(array)) throw std::invalid_argument(std::string("mode_selection source lacks ") + array);
    const auto& values = archive.at(array);
    if (values.shape.size() != 2 || values.shape[1] != dimension)
        throw std::invalid_argument("mode_selection states do not match the K-ROM state dimension");
    return {values.doubles(), values.shape[0]};
}

bool inBand(const json::Value& rule, Complex rate) {
    const auto& band = rule.at("band").array();
    if (band.empty()) return true;
    if (band.size() != 2) throw std::invalid_argument("mode_selection band is [low, high]");
    return std::abs(rate.imag()) >= band[0].number() && std::abs(rate.imag()) <= band[1].number();
}

std::size_t count(const json::Value& rule, std::size_t fallback) {
    const double value = rule.at("count").number(double(fallback));
    if (!(value >= 1) || value != std::floor(value)) throw std::invalid_argument("mode_selection.count must be a positive integer");
    return static_cast<std::size_t>(value);
}

std::vector<std::size_t> flatten(const std::vector<std::vector<std::size_t>>& all, const std::vector<std::size_t>& chosen) {
    std::vector<std::size_t> result;
    for (auto g : chosen) result.insert(result.end(), all[g].begin(), all[g].end());
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<std::size_t> stateObservability(const json::Value& rule, const c::KoopmanModel& model,
                                            const c::NpzRecord& parentArchive, json::Value& report) {
    const std::size_t n = model.stateDimension();
    std::vector<double> states; std::size_t samples = 0;
    if (!rule.at("states_source").isNull()) std::tie(states, samples) = loadStates(rule.at("states_source").string(), "values", n);
    else if (parentArchive.contains("training_states")) {
        states = parentArchive.at("training_states").doubles(); samples = states.size() / n;
    } else throw std::invalid_argument("state_observability needs states_source or K-ROM training_states");
    const auto spectrum = model.spectrum().values();
    const auto all = groups(spectrum);
    // Features: [Re psi_j, Im psi_j] per mode at every state, centred; targets: centred states.
    const std::size_t modes = model.liftedDimension();
    std::vector<Complex> psi(samples * modes);
    for (std::size_t s = 0; s < samples; ++s) {
        const auto value = model.evaluate(states.data() + s * n);
        std::copy(value.begin(), value.end(), psi.begin() + s * modes);
    }
    std::vector<double> mean(n, 0.0), total(n, 0.0);
    for (std::size_t s = 0; s < samples; ++s) for (std::size_t i = 0; i < n; ++i) mean[i] += states[s * n + i] / double(samples);
    std::vector<double> target(samples * n);
    for (std::size_t s = 0; s < samples; ++s) for (std::size_t i = 0; i < n; ++i) {
        target[s * n + i] = states[s * n + i] - mean[i]; total[i] += target[s * n + i] * target[s * n + i]; }
    for (double t : total) if (!(t > 0.0)) throw std::invalid_argument("state_observability: every state component must vary");
    const auto columnsOf = [&](std::size_t g) {
        std::vector<std::vector<double>> cols;
        const auto j = all[g].front();
        std::vector<double> re(samples), im(samples);
        for (std::size_t s = 0; s < samples; ++s) { re[s] = psi[s * modes + j].real(); im[s] = psi[s * modes + j].imag(); }
        cols.push_back(re);
        if (all[g].size() > 1 || std::abs(spectrum[j].imag()) > 1e-9) cols.push_back(im);
        for (auto& col : cols) { double m = 0; for (double v : col) m += v / double(samples); for (double& v : col) v -= m; }
        return cols;
    };
    const auto r2For = [&](const std::vector<std::size_t>& chosen) {
        std::vector<std::vector<double>> cols;
        for (auto g : chosen) for (auto& col : columnsOf(g)) cols.push_back(std::move(col));
        const std::size_t p = cols.size();
        std::vector<double> design(samples * p);
        for (std::size_t s = 0; s < samples; ++s) for (std::size_t k = 0; k < p; ++k) design[s * p + k] = cols[k][s];
        numerics::LeastSquaresOptions options; options.singularValueTolerance = 1e-12;
        const auto W = numerics::solveLeastSquares(design, samples, p, target, n, options);
        std::vector<double> r2(n, 0.0), residual(n, 0.0);
        for (std::size_t s = 0; s < samples; ++s) for (std::size_t i = 0; i < n; ++i) {
            double fit = 0; for (std::size_t k = 0; k < p; ++k) fit += design[s * p + k] * W[k * n + i];
            residual[i] += (target[s * n + i] - fit) * (target[s * n + i] - fit);
        }
        for (std::size_t i = 0; i < n; ++i) r2[i] = 1.0 - residual[i] / total[i];
        return r2;
    };
    std::vector<double> authority(all.size(), 0.0);
    const double weight = rule.at("control_authority_weight").number(0.0);
    if (weight > 0.0) {
        const auto x0 = model.systemConfiguration().initialState;
        const auto map = model.controlMap(x0.data(), rule.at("control_representation").string("exact"));
        for (std::size_t g = 0; g < all.size(); ++g) authority[g] = std::abs(map.at(all[g].front()));
        const double scale = *std::max_element(authority.begin(), authority.end());
        if (scale > 0) for (double& a : authority) a /= scale;
    }
    const double threshold = rule.at("r2_threshold").number(0.99);
    const std::size_t maxModes = static_cast<std::size_t>(rule.at("max_modes").number(8));
    const double minimumImprovement = rule.at("minimum_relative_improvement").number(1e-4);
    std::vector<std::size_t> chosen;
    double previous = 0.0;
    json::Value::Array steps;
    std::string stop = "max_modes";
    while (chosen.size() < maxModes) {
        std::optional<std::size_t> best; double bestScore = -1e300, bestMean = 0; std::vector<double> bestR2;
        for (std::size_t g = 0; g < all.size(); ++g) {
            if (std::find(chosen.begin(), chosen.end(), g) != chosen.end()) continue;
            auto trial = chosen; trial.push_back(g);
            const auto r2 = r2For(trial);
            double mean = 0; for (double v : r2) mean += v / double(n);
            const double score = (mean - previous) * (1.0 + weight * authority[g]);
            if (score > bestScore) { bestScore = score; best = g; bestMean = mean; bestR2 = r2; }
        }
        if (!best || bestMean - previous < minimumImprovement) { stop = "minimum_relative_improvement"; break; }
        chosen.push_back(*best);
        previous = bestMean;
        json::Value::Array r2;
        for (double v : bestR2) r2.emplace_back(v);
        steps.emplace_back(json::Value::Object{{"step", double(chosen.size())}, {"mode", double(all[*best].front())},
            {"score", bestScore}, {"r2", json::Value(std::move(r2))}, {"authority", authority[*best]}});
        if (std::all_of(bestR2.begin(), bestR2.end(), [&](double v) { return v >= threshold; })) { stop = "r2_threshold"; break; }
    }
    if (chosen.empty()) throw std::invalid_argument("state_observability selected no mode");
    report = json::Value::Object{{"rule", "state_observability"}, {"steps", json::Value(std::move(steps))},
                                 {"stop_reason", stop}, {"samples", double(samples)}};
    return flatten(all, chosen);
}

std::vector<std::size_t> nearestFrequency(const json::Value& rule, const c::KoopmanModel& model, json::Value& report) {
    const auto spectrum = rates(model);
    const auto all = groups(spectrum);
    Complex target;
    const auto& rate = rule.at("target_rate").array();
    if (rate.size() == 2) target = {rate[0].number(), rate[1].number()};
    else if (!rule.at("target_frequency").isNull()) target = {0.0, rule.at("target_frequency").number()};
    else throw std::invalid_argument("nearest_frequency needs target_rate [re, im] or target_frequency");
    std::vector<std::pair<double, std::size_t>> ranked;
    for (std::size_t g = 0; g < all.size(); ++g) {
        const auto value = spectrum[all[g].front()];
        if (!inBand(rule, value)) continue;
        const Complex folded{value.real(), std::abs(value.imag())};
        ranked.emplace_back(std::abs(folded - Complex{target.real(), std::abs(target.imag())}), g);
    }
    std::stable_sort(ranked.begin(), ranked.end());
    const auto wanted = count(rule, 1);
    if (ranked.size() < wanted) throw std::invalid_argument("nearest_frequency found too few modes in the band");
    std::vector<std::size_t> chosen;
    json::Value::Array distances;
    for (std::size_t k = 0; k < wanted; ++k) { chosen.push_back(ranked[k].second); distances.emplace_back(ranked[k].first); }
    report = json::Value::Object{{"rule", "nearest_frequency"}, {"distances", json::Value(std::move(distances))}};
    return flatten(all, chosen);
}

std::vector<std::size_t> validationScore(const json::Value& rule, const c::KoopmanModel& model, json::Value& report) {
    const auto spectrum = rates(model);
    const auto all = groups(spectrum);
    const std::size_t n = model.stateDimension(), modes = model.liftedDimension();
    std::vector<double> squares(modes, 0.0), norms(modes, 0.0);
    for (const auto& source : rule.at("validation_sources").array()) {
        const auto archive = c::readNpz(c::wide(source.string()));
        const auto time = archive.at("time").doubles();
        const auto values = archive.at("values").doubles();
        if (archive.at("values").shape.size() != 2 || archive.at("values").shape[1] != n)
            throw std::invalid_argument("validation source does not match the K-ROM state dimension");
        const auto first = model.evaluate(values.data());
        for (std::size_t k = 0; k < time.size(); ++k) {
            const auto psi = model.evaluate(values.data() + k * n);
            for (std::size_t j = 0; j < modes; ++j) {
                const auto predicted = std::exp(spectrum[j] * (time[k] - time[0])) * first[j];
                squares[j] += std::norm(psi[j] - predicted);
                norms[j] += std::norm(psi[j]);
            }
        }
    }
    if (rule.at("validation_sources").array().empty()) throw std::invalid_argument("validation_score needs validation_sources");
    std::vector<std::pair<double, std::size_t>> ranked;
    for (std::size_t g = 0; g < all.size(); ++g) {
        const auto j = all[g].front();
        if (!inBand(rule, spectrum[j]) || !(norms[j] > 0)) continue;
        ranked.emplace_back(std::sqrt(squares[j] / norms[j]), g);
    }
    std::stable_sort(ranked.begin(), ranked.end());
    const auto wanted = count(rule, 1);
    if (ranked.size() < wanted) throw std::invalid_argument("validation_score found too few modes");
    std::vector<std::size_t> chosen;
    json::Value::Array errors;
    for (std::size_t k = 0; k < wanted; ++k) { chosen.push_back(ranked[k].second); errors.emplace_back(ranked[k].first); }
    report = json::Value::Object{{"rule", "validation_score"}, {"anchored_relative_errors", json::Value(std::move(errors))}};
    return flatten(all, chosen);
}

}  // namespace

bool isExtendedModeRule(const std::string& rule) {
    return rule == "state_observability" || rule == "nearest_frequency" || rule == "validation_score";
}

std::vector<std::size_t> extendedModeRule(const json::Value& rule, const c::KoopmanModel& model,
                                          const c::NpzRecord& parentArchive, json::Value& report) {
    const auto name = rule.at("rule").string();
    if (name == "state_observability") return stateObservability(rule, model, parentArchive, report);
    if (name == "nearest_frequency") return nearestFrequency(rule, model, report);
    if (name == "validation_score") return validationScore(rule, model, report);
    throw std::invalid_argument("unknown mode_selection rule: " + name);
}

}  // namespace dynsyslab::workflow::postprocess_detail
