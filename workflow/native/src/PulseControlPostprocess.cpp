// Pulse and input-map audits (docs/architecture/KROM_TOOLS.md). Ported from the echo pulse
// side-car (tools/experiments/echo_pulse_analyze.cpp) and
// scripts/generate_rpe_control_identification_data.py.
//
//   pulse_response     paired causal sensitivities of a K-ROM from pulse_branches solves: per anchor,
//                      G(t) = (psi(x_+(t)) - psi(x_-(t))) / (u_+ - u_-) of its first +/- pair, with
//                      the immediate and peak sensitivity, the semigroup error of
//                      exp(lambda (t - t_1)) G(t_1) overall and by time band, the free branch's
//                      drift, the even (nonlinear) over odd response, the change of G with the
//                      pulse amplitude, the post-pulse semigroup error of wider pulses and, with
//                      reference_sources (the same branches with a finer step), the integration-
//                      step change of G. The input-affine lifted model is that of Korda and Mezic
//                      2018, Automatica 93:149-160.
//     {"name": "pulse_response", "settings": {"krom": "model.krom.npz", "evaluate": "state",
//        "bands": [[0.01, 0.05], [0.05, 0.5]], "reference_sources": [...]}}
//
//   control_map_audit  representation error of stored input maps against the exact map along
//                      path states, and the rank/condition audit of controlled training designs
//                      (core/ControlMapAudit.hpp).
//     {"name": "control_map_audit", "settings": {"krom": "model.krom.npz", "maps": ["fitted"],
//        "training_sources": [...], "audit_stride": 10, "relative_floor": 0}}
//     The job's sources are the paths (for example a tracking target), in model state units.
#include "PulseControlPostprocess.hpp"
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/ControlMapAudit.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;
using Complex = std::complex<double>;

namespace {

struct Branch {
    double amplitude{}, width{};
    long long pair{-1};
    std::size_t index{};
    std::vector<double> time;
    std::vector<Complex> psi;  // samples x modes
};

using Anchor = std::map<std::size_t, Branch>;  // by branch index

std::map<std::string, Anchor> readBranches(const std::vector<std::filesystem::path>& sources, const c::KoopmanModel& model, bool coordinates) {
    std::map<std::string, Anchor> anchors;
    for (const auto& path : sources) {
        const auto source = loadSource(path);
        const auto meta = json::parse(source.metadataJson);
        const auto& design = meta.at("input_design");
        if (design.at("kind").string() != "pulse_branches")
            throw std::invalid_argument("pulse_response reads pulse_branches solves: " + path.generic_u8string());
        if (meta.at("control_evaluation").string() != "sampled")
            throw std::invalid_argument("pulse_response needs solves with solver.control_evaluation sampled (held pulses): " + path.generic_u8string());
        Branch b;
        b.amplitude = design.at("command").number(); b.width = design.at("width").number();
        b.index = std::size_t(design.at("branch").number());
        if (!design.at("pair_id").isNull()) b.pair = (long long)design.at("pair_id").number();
        b.time = source.time();
        const auto n = source.dimension();
        for (std::size_t k = 0; k < source.samples(); ++k) {
            const double* x = source.values().data() + k * n;
            const auto v = coordinates ? model.evaluateFromCoordinates(x) : model.evaluateFromState(x);
            b.psi.insert(b.psi.end(), v.begin(), v.end());
        }
        // Anchors are numbered within one solve job (its batch).
        const auto batch = meta.at("batch_id").string(path.parent_path().parent_path().generic_u8string());
        const auto key = batch + "#" + std::to_string(std::size_t(design.at("anchor").number()));
        if (!anchors[key].emplace(b.index, std::move(b)).second)
            throw std::invalid_argument("pulse_response: a branch appears twice: " + path.generic_u8string());
    }
    return anchors;
}

double ratio(double e, double w) { return std::sqrt(e / std::max(w, 1e-300)); }

// The pair partner of branch b (same pair_id), or null.
const Branch* partner(const Anchor& a, const Branch& b) {
    for (const auto& [i, other] : a) if (i != b.index && other.pair == b.pair) return &other;
    return nullptr;
}

// G(t) of a +/- pair, samples x modes.
std::vector<Complex> sensitivity(const Branch& plus, const Branch& minus) {
    std::vector<Complex> g(plus.psi.size());
    const double du = plus.amplitude - minus.amplitude;
    for (std::size_t i = 0; i < g.size(); ++i) g[i] = (plus.psi[i] - minus.psi[i]) / du;
    return g;
}

std::pair<const Branch*, const Branch*> orderedPair(const Branch& a, const Branch& b) {
    return a.amplitude > 0 ? std::pair{&a, &b} : std::pair{&b, &a};
}

}  // namespace

std::filesystem::path runPulseResponse(const std::vector<std::filesystem::path>& sources, const J& settings,
                                       const std::filesystem::path& root, const std::string& label) {
    const auto model = c::loadKoopmanModel(std::filesystem::u8path(settings.at("krom").string()));
    const auto rates = model.continuousEigenvalues();
    const bool coordinates = settings.at("evaluate").string("state") == "coordinates";
    const auto anchors = readBranches(sources, model, coordinates);
    std::vector<std::filesystem::path> referencePaths;
    for (const auto& r : settings.at("reference_sources").array()) referencePaths.push_back(std::filesystem::u8path(r.string()));
    const auto references = referencePaths.empty() ? std::map<std::string, Anchor>{} : readBranches(referencePaths, model, coordinates);
    if (!references.empty() && references.size() != anchors.size()) throw std::invalid_argument("pulse_response reference_sources must hold the same anchors");
    std::vector<std::pair<double, double>> bands;
    for (const auto& b : settings.at("bands").array()) bands.emplace_back(b.array().at(0).number(), b.array().at(1).number());
    const auto m = rates.size();
    J::Array rows;
    std::vector<Complex> table;
    std::vector<double> clock;
    auto ref = references.begin();
    for (const auto& [key, anchor] : anchors) {
        const Branch* first = nullptr;
        for (const auto& [i, b] : anchor) if (b.pair >= 0 && partner(anchor, b)) { first = &b; break; }
        if (!first) throw std::invalid_argument("pulse_response: anchor " + key + " has no +/- pair");
        const auto [plus, minus] = orderedPair(*first, *partner(anchor, *first));
        const auto S = plus->time.size();
        for (const auto& [i, b] : anchor)
            if (b.time.size() != S || b.psi.size() != S * m) throw std::invalid_argument("pulse_response: branches of an anchor share one clock");
        if (S < 3) throw std::invalid_argument("pulse_response needs at least three samples");
        const auto G = sensitivity(*plus, *minus);
        const auto& t = plus->time;
        const double dt = t[1] - t[0];
        if (clock.empty()) clock = t;
        else if (clock.size() != S) throw std::invalid_argument("pulse_response: anchors share one clock length");
        table.insert(table.end(), G.begin(), G.end());
        const Branch* zero = nullptr;
        for (const auto& [i, b] : anchor) if (b.amplitude == 0) zero = &b;
        J::Array modes;
        for (std::size_t j = 0; j < m; ++j) {
            const auto g = [&](std::size_t k) { return G[k * m + j]; };
            double peak = 0, peakTime = 0, energy = 0, semigroup = 0;
            for (std::size_t k = 1; k < S; ++k) {
                if (std::abs(g(k)) > peak) { peak = std::abs(g(k)); peakTime = t[k] - t[0]; }
                energy += std::norm(g(k));
                semigroup += std::norm(g(k) - std::exp(rates[j] * (t[k] - t[1])) * g(1));
            }
            J::Object row{{"mode", double(j)}, {"immediate_sensitivity_magnitude", std::abs(g(1))},
                {"peak_sensitivity_magnitude", peak}, {"time_to_peak", peakTime},
                {"peak_over_immediate", peak / std::max(std::abs(g(1)), 1e-300)},
                {"pulse_semigroup_relative_error", ratio(semigroup, energy)}};
            if (zero) {
                const auto z = [&](std::size_t k) { return zero->psi[k * m + j]; };
                const auto mu = std::exp(rates[j] * dt);
                row["initial_zero_drift_over_pulse_effect"] = std::abs(z(1) - mu * z(0)) / std::max(std::abs(g(1) * plus->amplitude), 1e-300);
                double even = 0;
                for (std::size_t k = 0; k < S; ++k)
                    even += std::norm((plus->psi[k * m + j] + minus->psi[k * m + j] - 2.0 * z(k)) / (plus->amplitude - minus->amplitude));
                row["even_nonlinear_over_odd_response"] = ratio(even, energy);
            }
            J::Array bandRows;
            for (const auto& [lo, hi] : bands) {
                double e = 0, w = 0, drift = 0, free = 0;
                for (std::size_t k = 1; k < S; ++k) {
                    const double age = t[k] - t[0];
                    if (age < lo - 1e-9 || age >= hi - 1e-9) continue;
                    e += std::norm(g(k) - std::exp(rates[j] * (t[k] - t[1])) * g(1)); w += std::norm(g(k));
                    if (zero) { drift += std::norm(zero->psi[k * m + j] - std::exp(rates[j] * dt) * zero->psi[(k - 1) * m + j]); free += std::norm(zero->psi[(k - 1) * m + j]); }
                }
                J::Object band{{"from", lo}, {"to", hi}, {"pulse_semigroup_relative_error", w > 0 ? J(ratio(e, w)) : J()}};
                if (zero) band["zero_branch_step_error"] = free > 0 ? J(ratio(drift, free)) : J();
                bandRows.push_back(J(std::move(band)));
            }
            if (!bandRows.empty()) row["bands"] = J(std::move(bandRows));
            // Other pairs: amplitude change at equal width, post-pulse semigroup at other widths.
            J::Array amplitudeRows, widthRows;
            for (const auto& [i, b] : anchor) {
                if (b.pair < 0 || &b == plus || &b == minus || b.amplitude < 0) continue;
                const auto* other = partner(anchor, b);
                if (!other) continue;
                const auto H = sensitivity(b, *other);
                if (b.width == plus->width) {
                    double change = 0;
                    for (std::size_t k = 0; k < S; ++k) change += std::norm(H[k * m + j] - g(k));
                    amplitudeRows.push_back(J::Object{{"amplitude", b.amplitude}, {"sensitivity_change", ratio(change, energy)}});
                } else {
                    const auto end = std::max<std::size_t>(1, std::size_t(std::llround(b.width / dt)));
                    if (end >= S) continue;
                    double e = 0, w = 0;
                    for (std::size_t k = end; k < S; ++k) {
                        e += std::norm(H[k * m + j] - std::exp(rates[j] * (t[k] - t[end])) * H[end * m + j]); w += std::norm(H[k * m + j]);
                    }
                    widthRows.push_back(J::Object{{"width", b.width}, {"post_pulse_semigroup_relative_error", ratio(e, w)}});
                }
            }
            if (!amplitudeRows.empty()) row["amplitude_checks"] = J(std::move(amplitudeRows));
            if (!widthRows.empty()) row["width_checks"] = J(std::move(widthRows));
            modes.push_back(J(std::move(row)));
        }
        J::Object anchorRow{{"anchor", key}, {"amplitude", plus->amplitude}, {"width", plus->width}, {"modes", J(std::move(modes))}};
        if (!references.empty()) {
            const auto& refAnchor = ref->second;
            const auto* rp = refAnchor.count(plus->index) ? &refAnchor.at(plus->index) : nullptr;
            const auto* rm = refAnchor.count(minus->index) ? &refAnchor.at(minus->index) : nullptr;
            if (!rp || !rm || rp->psi.size() != plus->psi.size()) throw std::invalid_argument("pulse_response: reference branches do not match");
            const auto R = sensitivity(*rp, *rm);
            J::Array change;
            for (std::size_t j = 0; j < m; ++j) {
                double e = 0, w = 0;
                for (std::size_t k = 0; k < S; ++k) { e += std::norm(R[k * m + j] - G[k * m + j]); w += std::norm(G[k * m + j]); }
                change.push_back(ratio(e, w));
            }
            anchorRow["reference_step_sensitivity_change"] = J(std::move(change));
            ++ref;
        }
        rows.push_back(J(std::move(anchorRow)));
    }
    const auto first = loadSource(sources.front());
    auto meta = json::parse(metadata(first, "pulse-response", label, "\"schema\":\"dynsyslab.pulse-response/v1\"")).object();
    meta["settings"] = settings;
    meta["anchors"] = J(std::move(rows));
    meta["references"] = J::Array{"Korda and Mezic 2018, Automatica 93:149-160"};
    c::NpzRecord archive;
    archive.arrays["time"] = c::NpyArray::fromDoubles({clock.size()}, clock);
    archive.arrays["causal_sensitivity"] = c::NpyArray::fromComplex({anchors.size(), clock.size(), m}, table);
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, first, c::object_kind::report, "pulse-response", std::move(archive));
}

std::filesystem::path runControlMapAudit(const std::vector<std::filesystem::path>& sources, const J& settings,
                                         const std::filesystem::path& root, const std::string& label) {
    if (sources.empty()) throw std::invalid_argument("control_map_audit needs path sources");
    const auto model = c::loadKoopmanModel(std::filesystem::u8path(settings.at("krom").string()));
    std::vector<std::string> maps;
    for (const auto& name : settings.at("maps").array()) maps.push_back(name.string());
    if (maps.empty()) for (const auto& [name, map] : model.controlMaps()) maps.push_back(name);
    if (maps.empty()) throw std::invalid_argument("control_map_audit: the K-ROM stores no input map");
    const double floor = settings.at("relative_floor").number(0);
    const auto n = model.stateDimension();
    J::Object paths;
    for (const auto& path : sources) {
        const auto s = loadSource(path);
        if (s.dimension() != n) throw std::invalid_argument("control_map_audit: path states differ from the model's: " + path.generic_u8string());
        J::Object byMap;
        for (const auto& map : maps) {
            const auto e = c::controlMapPathError(model, map, s.values(), s.samples(), floor);
            byMap[map] = J::Object{{"samples", double(e.samples)}, {"rmse", e.rmse}, {"relative_frobenius", e.relativeFrobenius},
                {"median_pointwise_relative", e.medianRelative}, {"p95_pointwise_relative", e.p95Relative},
                {"maximum_pointwise_relative", e.maximumRelative}};
        }
        paths[path.generic_u8string()] = J(std::move(byMap));
    }
    J::Object designs;
    const auto stride = std::size_t(settings.at("audit_stride").number(10));
    for (const auto& item : settings.at("training_sources").array()) {
        const auto s = loadSource(std::filesystem::u8path(item.string()));
        const auto m = s.inputDimension();
        if (!m || s.dimension() != n) throw std::invalid_argument("control_map_audit training sources are controlled records of the model's states");
        const auto rows = std::min(s.samples(), s.inputs().size() / m);
        const std::vector<double> states(s.values().begin(), s.values().begin() + std::ptrdiff_t(rows * n));
        const std::vector<double> inputs(s.inputs().begin(), s.inputs().begin() + std::ptrdiff_t(rows * m));
        const auto a = c::controlDesignAudit(model, states, inputs, rows, m, stride);
        const auto list = [](const std::vector<double>& v) { J::Array out; for (const auto x : v) out.emplace_back(x); return J(std::move(out)); };
        designs[item.string()] = J::Object{{"audit_stride", double(stride)}, {"design_rows", double(a.rows)},
            {"constant_design_rank", double(a.constantRank)}, {"bilinear_design_rank", double(a.bilinearRank)},
            {"bilinear_design_shape", J::Array{double(a.rows), double(a.bilinearColumns)}},
            {"bilinear_design_condition_number", std::isfinite(a.bilinearCondition) ? J(a.bilinearCondition) : J()},
            {"input_minimum", a.inputMinimum}, {"input_maximum", a.inputMaximum}, {"input_mean", a.inputMean}, {"input_rms", a.inputRms},
            {"state_minimum", list(a.stateMinimum)}, {"state_maximum", list(a.stateMaximum)}, {"state_p01", list(a.stateP01)}, {"state_p99", list(a.stateP99)}};
    }
    const auto first = loadSource(sources.front());
    auto meta = json::parse(metadata(first, "control-map-audit", label, "\"schema\":\"dynsyslab.control-map-audit/v1\"")).object();
    meta["settings"] = settings;
    meta["path_error"] = J(std::move(paths));
    if (!designs.empty()) meta["design"] = J(std::move(designs));
    meta["evaluation"] = "the exact map grad psi(x) B enters only this audit, never a fit";
    meta["references"] = J::Array{"Willems et al. 2005, Syst. Control Lett. 54:325-329", "Korda and Mezic 2018, Automatica 93:149-160"};
    c::NpzRecord archive;
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, first, c::object_kind::report, "control-map-audit", std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
