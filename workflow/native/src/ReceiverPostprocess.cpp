// Post-processing methods for receiver arrays (docs/architecture/RECEIVER_ARRAYS.md):
//
//   receiver_signals  {"receivers": [[0.005, 0, 0], [0, 0.006, 0]], "sound_speed": 1480, "stride": 40,
//                      "window": [a, b], "terms": ["monopole_pressure", ...], "range_corrected": true,
//                      "contrast": false}
//       The RPE receiver catalog at every receiver on the nearest receiver's emission clock
//       (core/ReceiverArray.hpp). values = the summed selected terms per receiver (times the range
//       when range_corrected), plus, with contrast, C = D_2 p_2(t) - D_1 p_1(t - Delta) (the near
//       receiver's pressure at the far receiver's emission). Keeps source_states, every term
//       (<term>, near_at_<term>), ranges_m, lag_tau and ring_bounds (radius maxima).
//   state_readout     {"training": [features.npz], "mode": "neighbors" | "affine" | "linear", "neighbors": 8, "ridge": 0,
//                      "scales": [..] | "shared", "equilibrium": [1, 0]}
//       State estimates from the source's features. neighbors: asinh-scaled k-NN with neighbour
//       rows and weights (core/NeighborReadout.hpp); "shared" scales every feature by the training
//       RMS of the first. linear: the delayed-output inverse about the equilibrium
//       (core/ReceiverArray.hpp linearDelayInverse) from the system's Jacobian and the receivers'
//       feature gradients. Training labels are used only to fit; the source's own states are
//       carried untouched for scoring.
//   ring_errors       {"center": [1, 0], "windows": {"early": [0, 5]}, "role": "validation"}
//       Per-oscillation (ring) errors of values against source_states over ring_bounds, the pooled
//       whole-record relative RMS and window scores (a report object).
#include "PostprocessInternal.hpp"
#include "KromDatasetSupport.hpp"
#include "ReceiverPostprocess.hpp"

#include <dynsyslab/core/NeighborReadout.hpp>
#include <dynsyslab/core/ReceiverArray.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;
using J = json::Value;

namespace {

std::vector<double> numbers(const J& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.number());
    return out;
}

std::vector<std::string> terms(const J& settings) {
    if (settings.at("terms").isNull()) return {"monopole_pressure", "spherical_near_field_pressure", "dipole_induction_pressure", "dipole_radiation_pressure"};
    std::vector<std::string> out;
    for (const auto& t : settings.at("terms").array()) out.push_back(t.string());
    return out;
}

std::vector<c::AcousticReceiver> receivers(const J& settings) {
    std::vector<c::AcousticReceiver> out;
    const double speed = settings.at("sound_speed").number(1480);
    for (const auto& p : settings.at("receivers").array()) {
        const auto v = numbers(p);
        if (v.size() != 3) throw std::invalid_argument("receiver positions need three components in metres");
        out.push_back({{v[0], v[1], v[2]}, speed});
    }
    if (out.empty()) throw std::invalid_argument("receiver_signals needs receivers");
    return out;
}

double rms(const std::vector<double>& v, std::size_t width, std::size_t column) {
    double s = 0;
    const auto n = v.size() / width;
    for (std::size_t k = 0; k < n; ++k) s += v[k * width + column] * v[k * width + column];
    return n ? std::sqrt(s / double(n)) : 0;
}

c::NpyArray copy(const c::NpzRecord& from, const char* key) {
    if (!from.contains(key)) throw std::invalid_argument(std::string("the features archive has no ") + key);
    return from.at(key);
}

}  // namespace

std::filesystem::path runReceiverSignals(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                         const std::string& label) {
    const auto source = loadSource(path);
    const auto system = c::createSystem(systemConfiguration(source));
    if (system->descriptor().identifier != "rpe") throw std::invalid_argument("receiver_signals needs an RPE trajectory");
    if (datasetDynamics(source) == "forced") throw std::invalid_argument("receiver_signals uses the free-bubble acoustic catalog; the source is forced");
    const auto array = receivers(settings);
    c::ReceiverArrayOptions options;
    options.stride = std::size_t(settings.at("stride").number(1));
    if (!settings.at("window").isNull()) { const auto w = numbers(settings.at("window")); options.begin = w.at(0); options.end = w.at(1); }
    const auto signals = c::rpeReceiverArray(*system, array, source.time(), source.values(), options);
    const auto chosen = terms(settings);
    const bool corrected = settings.at("range_corrected").boolean(true), contrast = settings.at("contrast").boolean(false);
    if (contrast && array.size() != 2) throw std::invalid_argument("the receiver contrast feature needs exactly two receivers");
    const auto n = signals.samples(), r = signals.receivers;
    const auto width = r + (contrast ? 1 : 0);
    std::vector<double> features(n * width, 0);
    c::NpzRecord archive;
    for (std::size_t k = 0; k < signals.terms.size(); ++k) {
        std::vector<double> p(signals.pressures.begin() + std::ptrdiff_t(k * n * r), signals.pressures.begin() + std::ptrdiff_t((k + 1) * n * r));
        std::vector<double> q(signals.nearAtEmission.begin() + std::ptrdiff_t(k * n * r), signals.nearAtEmission.begin() + std::ptrdiff_t((k + 1) * n * r));
        archive.arrays[signals.terms[k]] = c::NpyArray::fromDoubles({n, r}, p);
        archive.arrays["near_at_" + signals.terms[k]] = c::NpyArray::fromDoubles({n, r}, q);
        if (std::find(chosen.begin(), chosen.end(), signals.terms[k]) == chosen.end()) continue;
        for (std::size_t s = 0; s < n; ++s) {
            for (std::size_t i = 0; i < r; ++i) features[s * width + i] += p[s * r + i] * (corrected ? signals.ranges[i] : 1.0);
            if (contrast) {  // D_2 p_2(t) - D_1 p_1(t - Delta): the near pressure at the far emission
                const auto far = signals.nearest == 0 ? 1 : 0;
                features[s * width + r] += signals.ranges[far] * p[s * r + far] - signals.ranges[signals.nearest] * q[s * r + far];
            }
        }
    }
    for (const auto& t : chosen)
        if (std::find(signals.terms.begin(), signals.terms.end(), t) == signals.terms.end())
            throw std::invalid_argument("unknown receiver pressure term " + t);
    const auto bounds = c::oscillationMaxima(source.time(), source.values(), source.dimension(), 1);
    archive.arrays["time"] = c::NpyArray::fromDoubles({n}, signals.time);
    archive.arrays["values"] = c::NpyArray::fromDoubles({n, width}, features);
    archive.arrays["source_states"] = c::NpyArray::fromDoubles({n, signals.stateDimension}, signals.states);
    archive.arrays["ring_bounds"] = c::NpyArray::fromDoubles({bounds.size()}, bounds);
    archive.arrays["ranges_m"] = c::NpyArray::fromDoubles({r}, signals.ranges);
    archive.arrays["lag_tau"] = c::NpyArray::fromDoubles({r}, signals.lags);
    auto meta = json::parse(metadata(source, "receiver-signals", label, "\"schema\":\"dynsyslab.receiver-signals/v1\"")).object();
    J::Array names, selected;
    for (std::size_t i = 0; i < r; ++i) names.push_back("receiver_" + std::to_string(i + 1));
    if (contrast) names.push_back("contrast");
    for (const auto& t : chosen) selected.push_back(t);
    meta["coordinate_names"] = names;
    meta["value_space"] = "receiver-pressures";
    meta["units"] = corrected ? "Pa m" : "Pa";
    meta["clock"] = "nearest receiver emission time (model time)";
    meta["receiver_settings"] = settings;
    meta["selected_terms"] = selected;
    meta["nearest_receiver"] = double(signals.nearest);
    meta["references"] = J::Array{"Lighthill 1952, Proc. R. Soc. A 211:564-587", "Ffowcs Williams and Hawkings 1969, Phil. Trans. R. Soc. A 264:321-342"};
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::timeSeries, "receiver-signals", std::move(archive));
}

std::filesystem::path runStateReadout(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                      const std::string& label) {
    const auto query = c::readNpz(path);
    const auto source = loadSource(path);
    const auto width = source.dimension(), n = source.samples();
    std::vector<std::filesystem::path> training;
    const auto& given = settings.at("training");
    if (std::holds_alternative<std::string>(given.storage())) training.push_back(std::filesystem::u8path(given.string()));
    else for (const auto& t : given.array()) training.push_back(std::filesystem::u8path(t.string()));
    if (training.empty()) throw std::invalid_argument("state_readout needs training feature archives");
    std::vector<double> features, states;
    std::size_t stateWidth = 0;
    for (const auto& t : training) {
        if (std::filesystem::absolute(t) == std::filesystem::absolute(path)) throw std::invalid_argument("state_readout: the query is also a training record");
        const auto record = c::readNpz(t);
        const auto f = record.at("values"), s = copy(record, "source_states");
        if (f.shape.size() != 2 || f.shape[1] != width) throw std::invalid_argument("state_readout: training features differ in width");
        stateWidth = s.shape.at(1);
        const auto fv = f.doubles(), sv = s.doubles();
        features.insert(features.end(), fv.begin(), fv.end());
        states.insert(states.end(), sv.begin(), sv.end());
    }
    const auto mode = settings.at("mode").string("neighbors");
    std::vector<double> estimates(n * stateWidth, 0);
    c::NpzRecord archive;
    auto meta = json::parse(metadata(source, "state-readout", label, "")).object();
    if (mode == "neighbors") {
        std::vector<double> scales;
        const auto& s = settings.at("scales");
        if (s.isNull() || (std::holds_alternative<std::string>(s.storage()) && s.string() == "shared"))
            scales.assign(width, rms(features, width, 0));
        else scales = numbers(s);
        const auto k = std::size_t(settings.at("neighbors").number(8));
        const c::NeighborReadout readout(features, states, scales, stateWidth);
        std::vector<double> rows, weights, nearest;
        for (std::size_t q = 0; q < n; ++q) {
            const auto e = readout.predict(source.values().data() + q * width, k);
            std::copy(e.values.begin(), e.values.end(), estimates.begin() + std::ptrdiff_t(q * stateWidth));
            for (std::size_t j = 0; j < k; ++j) { rows.push_back(double(e.rows[j])); weights.push_back(e.weights[j]); }
            nearest.push_back(e.nearestDistance);
        }
        archive.arrays["neighbor_rows"] = c::NpyArray::fromDoubles({n, k}, rows);
        archive.arrays["neighbor_weights"] = c::NpyArray::fromDoubles({n, k}, weights);
        archive.arrays["nearest_distance"] = c::NpyArray::fromDoubles({n}, nearest);
        archive.arrays["scales"] = c::NpyArray::fromDoubles({width}, scales);
        meta["neighbors"] = double(k);
        meta["references"] = J::Array{"Cover and Hart 1967, IEEE Trans. Inf. Theory 13:21-27", "Shepard 1968, Proc. 23rd ACM National Conference",
                                      "Burbidge, Magee and Robb 1988, J. Am. Stat. Assoc. 83:123-127", "Bentley 1975, Commun. ACM 18:509-517"};
    } else if (mode == "affine") {
        // x = W (y - y_mean) + x_mean by least squares on the training rows, with an optional ridge
        // (Hoerl and Kennard 1970, Technometrics 12:55-67).
        const auto rows = features.size() / width;
        std::vector<double> ym(width, 0), xm(stateWidth, 0);
        for (std::size_t k = 0; k < rows; ++k) {
            for (std::size_t j = 0; j < width; ++j) ym[j] += features[k * width + j] / double(rows);
            for (std::size_t j = 0; j < stateWidth; ++j) xm[j] += states[k * stateWidth + j] / double(rows);
        }
        std::vector<double> Y(rows * width), X(rows * stateWidth);
        for (std::size_t k = 0; k < rows; ++k) {
            for (std::size_t j = 0; j < width; ++j) Y[k * width + j] = features[k * width + j] - ym[j];
            for (std::size_t j = 0; j < stateWidth; ++j) X[k * stateWidth + j] = states[k * stateWidth + j] - xm[j];
        }
        dynsyslab::numerics::LeastSquaresOptions ls;
        ls.method = dynsyslab::numerics::LeastSquaresMethod::divideAndConquerSvd;
        ls.ridge = settings.at("ridge").number(0);
        const auto W = dynsyslab::numerics::solveLeastSquares(Y, rows, width, X, stateWidth, ls);  // width x stateWidth
        for (std::size_t q = 0; q < n; ++q)
            for (std::size_t j = 0; j < stateWidth; ++j) {
                double e = xm[j];
                for (std::size_t i = 0; i < width; ++i) e += W[i * stateWidth + j] * (source.values()[q * width + i] - ym[i]);
                estimates[q * stateWidth + j] = e;
            }
        archive.arrays["readout_weights"] = c::NpyArray::fromDoubles({width, stateWidth}, W);
        std::vector<double> bias(stateWidth);
        for (std::size_t j = 0; j < stateWidth; ++j) { bias[j] = xm[j]; for (std::size_t i = 0; i < width; ++i) bias[j] -= W[i * stateWidth + j] * ym[i]; }
        archive.arrays["readout_bias"] = c::NpyArray::fromDoubles({stateWidth}, bias);
        meta["references"] = J::Array{"Hoerl and Kennard 1970, Technometrics 12:55-67"};
    } else if (mode == "linear") {
        // The receivers, terms and scaling the features were made with.
        const auto featureMeta = json::parse(query.metadataJson());
        const auto& receiverSettings = featureMeta.at("receiver_settings");
        if (receiverSettings.isNull() || receiverSettings.at("contrast").boolean(false))
            throw std::invalid_argument("linear readout needs receiver_signals features without the contrast column");
        const auto system = c::createSystem(systemConfiguration(source));
        const auto x0 = settings.at("equilibrium").isNull() ? system->controlledInitialState() : numbers(settings.at("equilibrium"));
        const auto m = x0.size();
        const auto A = system->jacobian(0, x0.data());
        const auto array = receivers(receiverSettings);
        const auto lags = copy(query, "lag_tau").doubles(), ranges = copy(query, "ranges_m").doubles();
        const bool corrected = receiverSettings.at("range_corrected").boolean(true);
        std::vector<double> H(array.size() * m, 0);
        for (std::size_t i = 0; i < array.size(); ++i) {
            const auto catalog = c::rpeReceiverPressures(*system, array[i]);
            for (const auto& term : terms(receiverSettings)) {
                const auto g = catalog.at(term).jacobian(x0.data());
                for (std::size_t j = 0; j < m; ++j) H[i * m + j] += g[j].real() * (corrected ? ranges[i] : 1.0);
            }
        }
        const auto inverse = c::linearDelayInverse(A, H, lags, m);  // m x receivers
        std::vector<double> y0(array.size(), 0);  // the features at the equilibrium
        for (std::size_t i = 0; i < array.size(); ++i) {
            const auto catalog = c::rpeReceiverPressures(*system, array[i]);
            for (const auto& term : terms(receiverSettings)) y0[i] += catalog.at(term).evaluateReal(x0.data()).front() * (corrected ? ranges[i] : 1.0);
        }
        stateWidth = m;
        estimates.assign(n * m, 0);
        for (std::size_t q = 0; q < n; ++q)
            for (std::size_t j = 0; j < m; ++j) {
                double e = x0[j];
                for (std::size_t i = 0; i < array.size(); ++i) e += inverse[j * array.size() + i] * (source.values()[q * width + i] - y0[i]);
                estimates[q * m + j] = e;
            }
        archive.arrays["linear_inverse"] = c::NpyArray::fromDoubles({m, array.size()}, inverse);
        meta["equilibrium"] = J::Array(x0.begin(), x0.end());
        meta["references"] = J::Array{"Kalman 1960, Bol. Soc. Mat. Mexicana 5:102-119"};
    } else {
        throw std::invalid_argument("state_readout mode is neighbors, affine or linear");
    }
    archive.arrays["time"] = c::NpyArray::fromDoubles({n}, source.time());
    archive.arrays["values"] = c::NpyArray::fromDoubles({n, stateWidth}, estimates);
    if (query.contains("source_states")) archive.arrays["source_states"] = query.at("source_states");
    if (query.contains("ring_bounds")) archive.arrays["ring_bounds"] = query.at("ring_bounds");
    J::Array sources;
    for (const auto& t : training) sources.push_back(t.generic_u8string());
    meta["mode"] = mode;
    meta["training_sources"] = sources;
    meta["value_space"] = "state-estimates";
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::timeSeries, "state-readout", std::move(archive));
}

std::filesystem::path runRingErrors(const std::filesystem::path& path, const J& settings, const std::filesystem::path& root,
                                    const std::string& label) {
    const auto record = c::readNpz(path);
    const auto source = loadSource(path);
    const auto reference = copy(record, "source_states");
    const auto d = reference.shape.at(1), n = source.samples();
    if (source.dimension() != d || reference.shape[0] != n) throw std::invalid_argument("ring_errors: estimates and reference states differ in shape");
    const auto x = reference.doubles();
    const auto& y = source.values();
    const auto center = settings.at("center").isNull() ? std::vector<double>(d, 0) : numbers(settings.at("center"));
    const auto bounds = copy(record, "ring_bounds").doubles();
    const auto rings = c::intervalErrors(source.time(), x, y, d, bounds, center);
    // Pooled scores: the whole record and named windows.
    const auto score = [&](double lo, double hi) {
        std::vector<double> se(d, 0), ss(d, 0);
        for (std::size_t k = 0; k < n; ++k) {
            if (source.time()[k] < lo || source.time()[k] > hi) continue;
            for (std::size_t j = 0; j < d; ++j) { se[j] += std::pow(y[k * d + j] - x[k * d + j], 2); ss[j] += std::pow(x[k * d + j] - center[j], 2); }
        }
        J::Array pct;
        for (std::size_t j = 0; j < d; ++j) pct.push_back(ss[j] > 0 ? J(100 * std::sqrt(se[j] / ss[j])) : J());
        return pct;
    };
    J::Object windows;
    for (const auto& [name, w] : settings.at("windows").object()) { const auto v = numbers(w); windows[name] = score(v.at(0), v.at(1)); }
    const auto role = settings.at("role").string("validation");
    auto meta = json::parse(metadata(source, "ring-errors", label, "")).object();
    meta[role] = J::Object{{"relative_rms_percent", score(-1e300, 1e300)}, {"windows", windows}, {"rings", double(bounds.size() - 1)}};
    meta["center"] = J::Array(center.begin(), center.end());
    c::NpzRecord archive;
    const auto m = bounds.size() - 1;
    archive.arrays["ring_bounds"] = c::NpyArray::fromDoubles({bounds.size()}, bounds);
    archive.arrays["count"] = c::NpyArray::fromDoubles({m}, rings.count);
    archive.arrays["partial"] = c::NpyArray::fromDoubles({m}, rings.partial);
    archive.arrays["error_sum_squares"] = c::NpyArray::fromDoubles({m, d}, rings.errorSquares);
    archive.arrays["reference_sum_squares"] = c::NpyArray::fromDoubles({m, d}, rings.referenceSquares);
    archive.arrays["relative_rms_percent"] = c::NpyArray::fromDoubles({m, d}, rings.relativeRmsPercent);
    archive.arrays["metadata_json"] = c::NpyArray::fromString(json::dump(J(meta)));
    return publishObject(root, source, c::object_kind::report, "ring-errors", std::move(archive));
}

}  // namespace dynsyslab::workflow::postprocess_detail
