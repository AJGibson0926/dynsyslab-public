// `linear_systems` post-processing: the classical LTI toolkit applied to a
// KROM (split realization) or to raw A/B/C/D arrays (DMDc, linear SINDy part,
// any external model). Read-only on the source; one new immutable archive.
#include "PostprocessInternal.hpp"
#include "LinearSystemsClassical.hpp"

#include <dynsyslab/core/BalancedTruncation.hpp>
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanStateSpace.hpp>
#include <dynsyslab/core/ReducedOrderControl.hpp>
#include <dynsyslab/core/ReducedOrderEstimation.hpp>
#include <dynsyslab/core/SystemNorms.hpp>
#include <dynsyslab/core/SystemStructure.hpp>

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
using namespace dynsyslab::core;

json::Value finite(double value) { return std::isfinite(value) ? json::Value(value) : json::Value{}; }

std::size_t count(const json::Value& value, const char* what) {
    const double x = value.number(0);
    if (!std::isfinite(x) || x < 0 || std::floor(x) != x || x > 1e7)
        throw std::invalid_argument(std::string(what) + " must be a bounded nonnegative integer");
    return static_cast<std::size_t>(x);
}

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> out;
    for (const auto& item : value.array()) out.push_back(item.number(std::nan("")));
    for (double x : out) if (!std::isfinite(x)) throw std::invalid_argument("numeric array settings must be finite");
    return out;
}

std::vector<double> scaledIdentity(std::size_t n, double weight) {
    std::vector<double> out(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) out[i * n + i] = weight;
    return out;
}

// M' X M for M (n x k) and X (n x n), or M X M' for M (k x n).
std::vector<double> congruence(const std::vector<double>& M, const std::vector<double>& X,
                               std::size_t n, std::size_t k, bool transposeLeft) {
    const auto m = [&](std::size_t a, std::size_t i) { return transposeLeft ? M[a * k + i] : M[i * n + a]; };
    std::vector<double> out(k * k, 0.0);
    for (std::size_t i = 0; i < k; ++i) for (std::size_t j = 0; j < k; ++j)
        for (std::size_t a = 0; a < n; ++a) for (std::size_t b = 0; b < n; ++b)
            out[i * k + j] += m(a, i) * X[a * n + b] * m(b, j);
    return out;
}
std::vector<double> projectedPenalty(const std::vector<double>& lift, const std::vector<double>& Q,
                                     std::size_t n, std::size_t k) { return congruence(lift, Q, n, k, true); }
std::vector<double> projectedCovariance(const std::vector<double>& projection, const std::vector<double>& W,
                                        std::size_t n, std::size_t k) { return congruence(projection, W, n, k, false); }

// Compensator (Ac, Bc, Cc) closed around the plant (A, B, C, D) with u = Cc xc:
// [x; xc]' = [[A, B Cc], [Bc C, Ac + Bc D Cc]] [x; xc].
StateSpaceModel closeCompensator(const StateSpaceModel& p, const StateSpaceModel& c) {
    const std::size_t n = p.states, k = c.states, N = n + k, m = p.inputs, q = p.outputs;
    std::vector<double> A(N * N, 0.0);
    for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) A[i * N + j] = p.A[i * n + j];
    for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < k; ++j)
        for (std::size_t u = 0; u < m; ++u) A[i * N + n + j] += p.B[i * m + u] * c.C[u * k + j];
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t y = 0; y < q; ++y) A[(n + i) * N + j] += c.B[i * q + y] * p.C[y * n + j];
        for (std::size_t j = 0; j < k; ++j) {
            double v = c.A[i * k + j];
            for (std::size_t y = 0; y < q; ++y) for (std::size_t u = 0; u < m; ++u)
                v += c.B[i * q + y] * p.D[y * m + u] * c.C[u * k + j];
            A[(n + i) * N + n + j] = v;
        }
    }
    return makeStateSpace(std::move(A), {}, {}, {}, N, 0, 0, p.domain, p.sampleTime);
}

void complexArray(NpzRecord& out, const std::string& name, const std::vector<std::complex<double>>& values) {
    std::vector<double> flat;
    for (const auto v : values) { flat.push_back(v.real()); flat.push_back(v.imag()); }
    out.arrays[name] = NpyArray::fromDoubles({values.size(), 2}, flat);
}

void realization(NpzRecord& out, const std::string& prefix, const StateSpaceModel& m) {
    out.arrays[prefix + "_A"] = NpyArray::fromDoubles({m.states, m.states}, m.A);
    out.arrays[prefix + "_B"] = NpyArray::fromDoubles({m.states, m.inputs}, m.B);
    out.arrays[prefix + "_C"] = NpyArray::fromDoubles({m.outputs, m.states}, m.C);
    out.arrays[prefix + "_D"] = NpyArray::fromDoubles({m.outputs, m.inputs}, m.D);
    if (m.states) complexArray(out, prefix + "_poles", systemPoles(m));
}

std::vector<double> symmetricEigenvalues(const std::vector<double>& W, std::size_t n) {
    // Principal gramian strengths (ellipse semi-axes are their square roots).
    StateSpaceModel probe = makeStateSpace(W, {}, {}, {}, n, 0, 0);
    auto values = systemPoles(probe);
    std::vector<double> out;
    for (auto v : values) out.push_back(std::max(v.real(), 0.0));
    std::sort(out.rbegin(), out.rend());
    return out;
}

StateSpaceModel sourceModel(const std::filesystem::path& path, const NpzRecord& archive,
                            const json::Value& s, json::Value::Object& details) {
    const bool discrete = s.at("time_domain").string("continuous") == "discrete";
    if (s.at("time_domain").string("continuous") != "continuous" && !discrete)
        throw std::invalid_argument("time_domain must be continuous or discrete");
    if (archive.contains("A")) {
        const auto& a = archive.at("A");
        if (a.shape.size() != 2 || a.shape[0] != a.shape[1]) throw std::invalid_argument("A must be square");
        const std::size_t n = a.shape[0];
        auto matrix = [&](const char* name, std::size_t rows, bool byRows) -> std::pair<std::vector<double>, std::size_t> {
            if (!archive.contains(name)) return {{}, 0};
            const auto& x = archive.at(name);
            if (x.shape.size() != 2 || (byRows ? x.shape[0] : x.shape[1]) != rows)
                throw std::invalid_argument(std::string(name) + " has inconsistent dimensions");
            return {x.doubles(), byRows ? x.shape[1] : x.shape[0]};
        };
        const auto [B, m] = matrix("B", n, true);
        const auto [C, p] = matrix("C", n, false);
        std::vector<double> D;
        if (archive.contains("D")) D = archive.at("D").doubles();
        details["source_kind"] = "matrices";
        return makeStateSpace(a.doubles(), B, C, D, n, m, p,
            discrete ? TimeDomain::discrete : TimeDomain::continuous, s.at("sample_time").number());
    }
    std::vector<std::size_t> selected;
    for (const auto& v : s.at("selected_modes").array()) selected.push_back(count(v, "selected mode"));
    const auto model = loadKoopmanModel(path, selected);
    KoopmanStateSpaceOptions o;
    for (const auto& v : s.at("measurements").array()) o.measurements.push_back(v.string());
    if (!s.at("measurement").string().empty()) o.measurements.push_back(s.at("measurement").string());
    o.controlMap = s.at("control_map").string();
    o.operatingState = numbers(s.at("operating_state"));
    o.operatingLifted = numbers(s.at("operating_lifted"));
    o.domain = discrete ? TimeDomain::discrete : TimeDomain::continuous;
    o.sampleTime = s.at("sample_time").number();
    details["source_kind"] = "krom";
    details["source_krom"] = path.generic_u8string();
    json::Value::Array names; for (const auto& n : o.measurements) names.push_back(n);
    details["measurements"] = names; details["control_map"] = o.controlMap;
    auto system = koopmanStateSpace(model, o);
    const auto balancing = s.at("balancing_input").string(o.controlMap.empty() ? "none" : "control_map");
    if (balancing == "process_noise") {  // unit disturbance on every lifted coordinate (estimation view)
        if (!o.controlMap.empty()) throw std::invalid_argument("balancing_input process_noise excludes control_map");
        system.inputs = system.states; system.B.assign(system.states * system.states, 0.0);
        for (std::size_t i = 0; i < system.states; ++i) system.B[i * system.states + i] = 1.0;
        system.D.assign(system.outputs * system.inputs, 0.0);
    } else if (balancing != "none" && balancing != "control_map")
        throw std::invalid_argument("balancing_input must be control_map, process_noise or none");
    details["balancing_input"] = balancing;
    return system;
}

json::Value structureSummary(const StateSpaceModel& m, double tolerance) {
    const auto r = analyzeStructure(m, tolerance);
    return json::Value::Object{{"controllable_dimension", double(r.controllableDimension)},
        {"observable_dimension", double(r.observableDimension)}, {"minimal_order", double(r.minimalOrder)},
        {"controllable", r.controllable}, {"observable", r.observable}, {"stabilizable", r.stabilizable},
        {"detectable", r.detectable}, {"pbh_agrees", r.pbhAgrees},
        {"pbh", json::Value::Object{{"controllable", r.pbhControllable}, {"observable", r.pbhObservable},
            {"stabilizable", r.pbhStabilizable}, {"detectable", r.pbhDetectable},
            {"min_controllability_margin", r.pbhControllabilityMargin},
            {"min_observability_margin", r.pbhObservabilityMargin}}}};
}

void frequencyArrays(NpzRecord& out, const std::string& prefix, const StateSpaceModel& m,
                     const std::vector<double>& grid) {
    const auto f = frequencyResponse(m, grid);
    const std::size_t k = std::min(m.inputs, m.outputs);
    out.arrays[prefix + "_singular_values"] = NpyArray::fromDoubles({grid.size(), k}, f.singularValues);
    if (m.inputs == 1 && m.outputs == 1) {
        std::vector<double> phase;
        for (auto g : f.response) phase.push_back(std::arg(g));
        out.arrays[prefix + "_phase"] = NpyArray::fromDoubles({grid.size()}, phase);
    }
}

}  // namespace

std::filesystem::path runLinearSystems(const std::filesystem::path& sourcePath, const json::Value& s,
                                       const std::filesystem::path& root, const std::string& label) {
    const auto archive = readNpz(sourcePath);
    SourceData source; source.path = sourcePath; source.metadataJson = archive.metadataJson();
    json::Value::Object details;
    details["schema"] = "dynsyslab.linear-systems/v1";
    const auto model = sourceModel(sourcePath, archive, s, details);
    const std::size_t n = model.states;
    details["time_domain"] = timeDomainName(model.domain);
    details["sample_time"] = model.sampleTime;
    details["states"] = double(n); details["inputs"] = double(model.inputs); details["outputs"] = double(model.outputs);
    NpzRecord out;
    realization(out, "system", model);
    std::vector<double> damping;
    for (const auto& d : modalDamping(model))
        damping.insert(damping.end(), {d.pole.real(), d.pole.imag(), d.naturalFrequency, d.dampingRatio});
    out.arrays["modal_damping"] = NpyArray::fromDoubles({n, 4}, damping);
    if (model.inputs && model.outputs) details["structure"] = structureSummary(model, s.at("structure_tolerance").number(1e-10));

    BalancedReductionOptions reduction;
    reduction.order = count(s.at("reduction_order"), "reduction_order");
    reduction.errorTolerance = s.at("reduction_error_tolerance").number();
    const auto method = s.at("reduction_method").string("truncation");
    if (method != "truncation" && method != "residualization") throw std::invalid_argument("reduction_method must be truncation or residualization");
    reduction.method = method == "residualization" ? BalancedReductionMethod::residualization : BalancedReductionMethod::truncation;
    reduction.stabilityMargin = s.at("stability_margin").number(1e-8);
    const bool io = model.inputs && model.outputs;
    std::optional<BalancedReduction> reduced;
    if (io) {
        reduced = balancedReduction(model, reduction);
        const auto& r = *reduced;
        details["reduction"] = json::Value::Object{{"method", balancedReductionMethodName(r.method)},
            {"order", double(r.order)}, {"stable_order", double(r.stableOrder)},
            {"unstable_order", double(r.unstableOrder)}, {"error_bound", finite(r.errorBound)},
            {"lower_bound", finite(r.lowerBound)}, {"stability_margin", reduction.stabilityMargin}};
        const auto& hsv = r.hankelSingularValues;
        out.arrays["hankel_singular_values"] = NpyArray::fromDoubles({hsv.size()}, hsv);
        std::vector<double> bound(hsv.size() + 1, 0.0);  // bound[k] = 2 sum_{i>=k} sigma_i
        for (std::size_t k = hsv.size(); k-- > 0;) bound[k] = bound[k + 1] + 2.0 * hsv[k];
        out.arrays["hankel_error_bound"] = NpyArray::fromDoubles({bound.size()}, bound);
        if (r.order) realization(out, "reduced", r.reduced);
        out.arrays["reduction_projection"] = NpyArray::fromDoubles({r.order, n}, r.projection);
        out.arrays["reduction_lift"] = NpyArray::fromDoubles({n, r.order}, r.lift);
        if (r.stablePart.states) {  // gramians of the stable part (the whole model when stable)
            const auto Wc = controllabilityGramian(r.stablePart), Wo = observabilityGramian(r.stablePart);
            const std::size_t ns = r.stablePart.states;
            out.arrays["controllability_gramian"] = NpyArray::fromDoubles({ns, ns}, Wc);
            out.arrays["observability_gramian"] = NpyArray::fromDoubles({ns, ns}, Wo);
            const auto ec = symmetricEigenvalues(Wc, ns), eo = symmetricEigenvalues(Wo, ns);
            out.arrays["controllability_gramian_eigenvalues"] = NpyArray::fromDoubles({ns}, ec);
            out.arrays["observability_gramian_eigenvalues"] = NpyArray::fromDoubles({ns}, eo);
            details["gramian_coordinates"] = r.unstableOrder ? "stable part of the stable/unstable split" : "model coordinates";
        }
        if (r.stablePart.states) {
            StateSpaceModel empty; empty.inputs = model.inputs; empty.outputs = model.outputs;
            empty.D = r.reduced.D; empty.domain = model.domain; empty.sampleTime = model.sampleTime;
            const auto error = parallelConnection(r.stablePart, r.reducedStablePart.states ? r.reducedStablePart : empty, -1.0);
            details["reduction_error_hinf"] = finite(hInfinityNorm(error).value);
        }
        const auto hinf = hInfinityNorm(model);
        json::Value::Object norms{{"h2", finite(h2Norm(model))}, {"hinf", finite(hinf.value)},
            {"hinf_peak_frequency", finite(hinf.peakFrequency)}};
        if (r.unstableOrder && r.stablePart.states) {  // full norms are infinite; report the stable part
            const auto stableHinf = hInfinityNorm(r.stablePart);
            norms["stable_part_h2"] = finite(h2Norm(r.stablePart));
            norms["stable_part_hinf"] = finite(stableHinf.value);
            norms["stable_part_hinf_peak_frequency"] = finite(stableHinf.peakFrequency);
        }
        details["norms"] = norms;
        const auto grid = logarithmicFrequencyGrid(model, count(s.at("frequency_points"), "frequency_points") ?
            count(s.at("frequency_points"), "frequency_points") : 400);
        out.arrays["frequencies"] = NpyArray::fromDoubles({grid.size()}, grid);
        frequencyArrays(out, "full", model, grid);
        if (r.order) frequencyArrays(out, "reduced", r.reduced, grid);
        if (model.inputs == model.outputs) complexArray(out, "transmission_zeros", transmissionZeros(model));
        std::vector<std::string> names;  // channel labels: K-ROM measurement names, indexed when multi-row
        for (const auto& v : details["measurements"].array()) names.push_back(v.string());
        if (names.size() == 1 && model.outputs > 1) {
            const auto base = names.front(); names.clear();
            for (std::size_t k = 0; k < model.outputs; ++k) names.push_back(base + "[" + std::to_string(k) + "]");
        }
        appendTransferFunctions(out, details, model, names, grid);
    }
    const double qWeight = s.at("lqr_state_weight").number(), rWeight = s.at("lqr_input_weight").number();
    if (qWeight > 0 && rWeight > 0 && model.inputs) {
        const auto Q = scaledIdentity(n, qWeight), R = scaledIdentity(model.inputs, rWeight);
        const auto full = designLQR(model, Q, R);
        out.arrays["lqr_gain"] = NpyArray::fromDoubles({model.inputs, n}, full.gain);
        complexArray(out, "lqr_closed_loop_poles", full.closedLoopPoles);
        if (reduced && reduced->order) {
            const auto red = designReducedLQR(model, reduction, Q, R);
            out.arrays["reduced_lqr_lifted_gain"] = NpyArray::fromDoubles({model.inputs, n}, red.liftedGain);
            complexArray(out, "reduced_lqr_lifted_poles", red.liftedClosedLoopPoles);
            details["reduced_lqr_lifted_stable"] = red.liftedClosedLoopStable;
        }
    }
    const double wWeight = s.at("lqe_process_weight").number(), vWeight = s.at("lqe_measurement_weight").number();
    if (wWeight > 0 && vWeight > 0 && model.outputs) {
        const auto W = scaledIdentity(n, wWeight), V = scaledIdentity(model.outputs, vWeight);
        const auto full = designLQE(model, W, V);
        out.arrays["lqe_gain"] = NpyArray::fromDoubles({n, model.outputs}, full.gain);
        complexArray(out, "lqe_observer_poles", full.observerPoles);
        if (reduced && reduced->order) {
            const auto red = designReducedLQE(model, reduction, W, V);
            out.arrays["reduced_lqe_lifted_gain"] = NpyArray::fromDoubles({n, model.outputs}, red.liftedGain);
            complexArray(out, "reduced_lqe_lifted_poles", red.liftedObserverPoles);
            details["reduced_lqe_lifted_stable"] = red.liftedObserverStable;
        }
    }
    // Coordinated estimation + control: the LQG output-feedback compensator (KLQE
    // feeding KLQR), at full order and, with a reduction, at reduced order closed
    // around the FULL model.
    if (qWeight > 0 && rWeight > 0 && wWeight > 0 && vWeight > 0 && model.inputs && model.outputs) {
        const auto lqg = designLQG(model, scaledIdentity(n, qWeight), scaledIdentity(model.inputs, rWeight),
                                   scaledIdentity(n, wWeight), scaledIdentity(model.outputs, vWeight));
        complexArray(out, "lqg_closed_loop_poles", lqg.closedLoopPoles);
        out.arrays["lqg_compensator_A"] = NpyArray::fromDoubles({n, n}, lqg.compensator.A);
        out.arrays["lqg_compensator_B"] = NpyArray::fromDoubles({n, model.outputs}, lqg.compensator.B);
        out.arrays["lqg_compensator_C"] = NpyArray::fromDoubles({model.inputs, n}, lqg.compensator.C);
        json::Value::Object feedback{{"compensator_order", double(n)},
            {"closed_loop_stable", isAsymptoticallyStable(lqg.closedLoop)}};
        if (reduced && reduced->order) {
            const auto& red = *reduced;
            const std::size_t k = red.order;
            const auto Q = projectedPenalty(red.lift, scaledIdentity(n, qWeight), n, k);
            const auto W = projectedCovariance(red.projection, scaledIdentity(n, wWeight), n, k);
            const auto small = designLQG(red.reduced, Q, scaledIdentity(model.inputs, rWeight), W,
                                         scaledIdentity(model.outputs, vWeight));
            const auto loop = closeCompensator(model, small.compensator);
            const auto poles = systemPoles(loop);
            complexArray(out, "reduced_lqg_full_closed_loop_poles", poles);
            out.arrays["reduced_lqg_compensator_A"] = NpyArray::fromDoubles({k, k}, small.compensator.A);
            out.arrays["reduced_lqg_compensator_B"] = NpyArray::fromDoubles({k, model.outputs}, small.compensator.B);
            out.arrays["reduced_lqg_compensator_C"] = NpyArray::fromDoubles({model.inputs, k}, small.compensator.C);
            feedback["reduced_compensator_order"] = double(k);
            feedback["reduced_closed_loop_stable"] = isAsymptoticallyStable(loop);
        }
        details["output_feedback"] = std::move(feedback);
    }
    appendLoopRobustness(out, details, model, s);
    if (!s.at("readouts").array().empty()) {  // K-ROM readouts (e.g. state) seen through the measurement
        if (archive.contains("A")) throw std::invalid_argument("readouts need a K-ROM source");
        std::vector<std::size_t> selected;
        for (const auto& v : s.at("selected_modes").array()) selected.push_back(count(v, "selected mode"));
        KoopmanStateSpaceOptions ro;
        std::vector<std::string> readoutNames;
        for (const auto& v : s.at("readouts").array()) ro.measurements.push_back(v.string());
        const auto krom = loadKoopmanModel(sourcePath, selected);
        const auto readoutModel = koopmanStateSpace(krom, ro);
        for (const auto& name : ro.measurements) {
            const auto rows = krom.reconstruction(name).outputs;
            for (std::size_t k = 0; k < rows; ++k) readoutNames.push_back(rows > 1 ? name + "[" + std::to_string(k) + "]" : name);
        }
        appendReadoutObservability(out, details, model, readoutModel, readoutNames, s);
    }  // classical margins of LQR, Kalman, LQG and LTR loops
    auto serialized = json::dump(json::Value(details));
    out.arrays["metadata_json"] = NpyArray::fromString(metadata(source, "linear-systems", label,
        serialized.substr(1, serialized.size() - 2)));
    return publish(root, source, "linear-systems", label, std::move(out));
}

}  // namespace dynsyslab::workflow::postprocess_detail
