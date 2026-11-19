// `empirical_observability`: local observability of the archived nonlinear system
// through a named observable (e.g. the far-field echo) at points along a recorded
// trajectory, with the run's own applied input. One immutable derived archive.
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/EmpiricalObservability.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
using namespace dynsyslab::core;

json::Value finiteOrNull(double v) { return std::isfinite(v) ? json::Value(v) : json::Value{}; }

// Piecewise-constant input recorded with the run: the controller-clock applied
// input when present, otherwise the control history on the output grid.
std::function<double(double)> recordedInput(const NpzRecord& archive) {
    std::vector<double> t, u;
    if (archive.contains("applied_input") && archive.contains("estimator_time")) {
        t = archive.at("estimator_time").doubles(); u = archive.at("applied_input").doubles();
    } else if (archive.contains("control_values")) {
        t = archive.at("time").doubles(); u = archive.at("control_values").doubles();
        if (u.size() != t.size()) throw std::invalid_argument("empirical_observability needs one control per output sample");
    } else return {};
    if (t.size() != u.size() || t.empty()) throw std::invalid_argument("recorded input and its clock differ in length");
    auto times = std::make_shared<std::vector<double>>(std::move(t));
    auto values = std::make_shared<std::vector<double>>(std::move(u));
    return [times, values](double time) {
        const auto it = std::upper_bound(times->begin(), times->end(), time);
        const std::size_t k = it == times->begin() ? 0 : std::size_t(it - times->begin()) - 1;
        return (*values)[k];
    };
}

std::vector<double> numbers(const json::Value& v) {
    std::vector<double> out;
    for (const auto& x : v.array()) out.push_back(x.number());
    return out;
}
}  // namespace

std::filesystem::path runEmpiricalObservability(const std::filesystem::path& sourcePath, const json::Value& s,
                                                const std::filesystem::path& root, const std::string& label) {
    auto source = loadSource(sourcePath);
    const auto archive = readNpz(sourcePath);
    const auto system = createSystem(systemConfiguration(source));
    const std::size_t n = system->descriptor().stateDimension;
    if (source.dimension != n) throw std::invalid_argument("empirical_observability: archive state dimension differs from the system");
    const auto name = s.at("observable").string("far_field_acoustic_echo");
    const auto catalog = system->observableCatalog();
    if (!catalog.contains(name)) throw std::invalid_argument("unknown observable for this system: " + name);
    const auto& observable = catalog.at(name);
    const std::size_t p = observable.outputDimension();

    const auto& t = source.time;
    const double start = s.at("start").number(t.front()), stop = s.at("stop").number(t.back());
    const std::size_t points = static_cast<std::size_t>(s.at("points").number(41));
    EmpiricalObservabilityOptions o;
    o.horizon = s.at("horizon").number(1.0);
    o.samples = static_cast<std::size_t>(s.at("samples").number(201));
    o.substeps = static_cast<std::size_t>(s.at("substeps").number(20));
    o.input = recordedInput(archive);
    const double V = s.at("noise_intensity").number(1e-4);
    if (points < 1 || !(stop > start) || !(o.horizon > 0) || !(V > 0))
        throw std::invalid_argument("empirical_observability needs points >= 1, stop > start, horizon > 0, noise_intensity > 0");
    if (stop + o.horizon > t.back() + 1e-12) throw std::invalid_argument("stop + horizon must lie inside the recorded run");

    // Per-state scale (orbit RMS deviation over the window unless given); perturbations relative to it.
    auto scale = numbers(s.at("state_scale"));
    if (scale.empty()) {
        scale.assign(n, 0.0);
        std::vector<double> mean(n, 0.0); std::size_t count = 0;
        for (std::size_t k = 0; k < t.size(); ++k) if (t[k] >= start && t[k] <= stop) {
            for (std::size_t i = 0; i < n; ++i) mean[i] += source.values[k * n + i]; ++count; }
        for (auto& m : mean) m /= double(std::max<std::size_t>(count, 1));
        for (std::size_t k = 0; k < t.size(); ++k) if (t[k] >= start && t[k] <= stop)
            for (std::size_t i = 0; i < n; ++i) { const double d = source.values[k * n + i] - mean[i]; scale[i] += d * d; }
        for (auto& v : scale) v = std::max(std::sqrt(v / double(std::max<std::size_t>(count, 1))), 1e-6);
    }
    if (scale.size() != n) throw std::invalid_argument("state_scale needs one value per state");
    const double relative = s.at("relative_perturbation").number(1e-5);
    o.perturbation.clear();
    for (double v : scale) o.perturbation.push_back(relative * v);

    std::vector<double> times, states, gram, eig, vec, cond, sEig, sVec, sCond, bound, rel, jac, valid;
    json::Value::Array weakDirection;
    for (std::size_t m = 0; m < points; ++m) {
        const double target = points == 1 ? start : start + (stop - start) * double(m) / double(points - 1);
        const std::size_t k = std::size_t(std::lower_bound(t.begin(), t.end(), target) - t.begin());
        const std::size_t idx = std::min(k, t.size() - 1);
        std::vector<double> x(source.values.begin() + idx * n, source.values.begin() + (idx + 1) * n);
        const auto r = empiricalObservability(*system, observable, x, t[idx], o);
        times.push_back(t[idx]); states.insert(states.end(), x.begin(), x.end());
        gram.insert(gram.end(), r.gramian.begin(), r.gramian.end());
        eig.insert(eig.end(), r.eigenvalues.begin(), r.eigenvalues.end());
        vec.insert(vec.end(), r.eigenvectors.begin(), r.eigenvectors.end());
        cond.push_back(r.conditionNumber);
        jac.insert(jac.end(), r.outputJacobian.begin(), r.outputJacobian.end());
        valid.push_back(r.valid ? 1.0 : 0.0);
        // Scaled gramian S W S (states in units of their orbit scale): comparable eigenvalues.
        const auto e2 = scaledObservability(r, scale);
        sEig.insert(sEig.end(), e2.eigenvalues.begin(), e2.eigenvalues.end());
        sVec.insert(sVec.end(), e2.eigenvectors.begin(), e2.eigenvectors.end());
        sCond.push_back(e2.conditionNumber);
        const auto b = estimationBoundStd(r, V);
        bound.insert(bound.end(), b.begin(), b.end());
        for (std::size_t i = 0; i < n; ++i) rel.push_back(b[i] / scale[i]);
    }
    const std::size_t N = times.size();
    NpzRecord out;
    out.arrays["sample_times"] = NpyArray::fromDoubles({N}, times);
    out.arrays["sample_states"] = NpyArray::fromDoubles({N, n}, states);
    out.arrays["gramians"] = NpyArray::fromDoubles({N, n, n}, gram);
    out.arrays["eigenvalues"] = NpyArray::fromDoubles({N, n}, eig);
    out.arrays["eigenvectors"] = NpyArray::fromDoubles({N, n, n}, vec);
    out.arrays["condition_number"] = NpyArray::fromDoubles({N}, cond);
    out.arrays["scaled_eigenvalues"] = NpyArray::fromDoubles({N, n}, sEig);
    out.arrays["scaled_eigenvectors"] = NpyArray::fromDoubles({N, n, n}, sVec);
    out.arrays["scaled_condition_number"] = NpyArray::fromDoubles({N}, sCond);
    out.arrays["estimation_bound_std"] = NpyArray::fromDoubles({N, n}, bound);
    out.arrays["relative_estimation_bound"] = NpyArray::fromDoubles({N, n}, rel);
    out.arrays["output_jacobian"] = NpyArray::fromDoubles({N, p, n}, jac);
    out.arrays["valid"] = NpyArray::fromDoubles({N}, valid);
    out.arrays["state_scale"] = NpyArray::fromDoubles({n}, scale);

    std::vector<double> sortedCond = sCond; std::sort(sortedCond.begin(), sortedCond.end());
    std::vector<double> medianRel(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::vector<double> col; for (std::size_t m = 0; m < N; ++m) col.push_back(rel[m * n + i]);
        std::sort(col.begin(), col.end()); medianRel[i] = col[col.size() / 2];
    }
    // Weakest scaled direction, averaged by absolute component (sign-free).
    std::vector<double> weak(n, 0.0);
    for (std::size_t m = 0; m < N; ++m) for (std::size_t i = 0; i < n; ++i)
        weak[i] += std::abs(sVec[m * n * n + i * n + (n - 1)]) / double(N);
    json::Value::Array names, medians, weakArray;
    for (std::size_t i = 0; i < n; ++i) {
        names.push_back(i < system->descriptor().stateNames.size() ? system->descriptor().stateNames[i] : "x" + std::to_string(i));
        medians.push_back(finiteOrNull(medianRel[i])); weakArray.push_back(weak[i]);
    }
    std::size_t invalid = 0; for (double v : valid) invalid += v == 0.0;
    json::Value::Object details{{"schema", "dynsyslab.empirical-observability/v1"}, {"observable", name},
        {"horizon", o.horizon}, {"samples", double(o.samples)}, {"substeps", double(o.substeps)},
        {"relative_perturbation", relative}, {"noise_intensity", V}, {"points", double(N)},
        {"window", json::Value::Array{start, stop}}, {"input", o.input ? "recorded applied input" : "unforced"},
        {"state_names", names}, {"median_relative_estimation_bound", medians},
        {"median_scaled_condition_number", finiteOrNull(sortedCond[N / 2])},
        {"max_scaled_condition_number", finiteOrNull(sortedCond.back())},
        {"weakest_direction_mean_abs_components", weakArray}, {"invalid_points", double(invalid)},
        {"definition", "W = int dy_i . dy_j dt over the horizon, dy_i central-difference output sensitivity; "
                       "bound = sqrt(V (W^-1)_ii); scaled = diag(state_scale) W diag(state_scale)"}};
    auto serialized = json::dump(json::Value(details));
    out.arrays["metadata_json"] = NpyArray::fromString(metadata(source, "empirical-observability", label,
        serialized.substr(1, serialized.size() - 2)));
    return publish(root, source, "empirical-observability", label, std::move(out));
}

}  // namespace dynsyslab::workflow::postprocess_detail
