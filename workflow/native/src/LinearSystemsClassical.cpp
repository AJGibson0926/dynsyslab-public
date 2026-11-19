#include "LinearSystemsClassical.hpp"

#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/LoopMargins.hpp>
#include <dynsyslab/core/TransferFunction.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <cmath>
#include <complex>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace {
using namespace dynsyslab::core;
namespace json = dynsyslab::core::json;

json::Value finiteOrNull(double v) { return std::isfinite(v) ? json::Value(v) : json::Value{}; }

void storeComplex(NpzRecord& out, const std::string& name, const std::vector<std::complex<double>>& v) {
    std::vector<double> flat;
    for (const auto& x : v) { flat.push_back(x.real()); flat.push_back(x.imag()); }
    out.arrays[name] = NpyArray::fromDoubles({v.size(), 2}, flat);
}

std::vector<double> identity(std::size_t n, double w) {
    std::vector<double> m(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) m[i * n + i] = w;
    return m;
}

json::Value::Object summary(const LoopRobustness& r) {
    return json::Value::Object{
        {"closed_loop_stable", r.closedLoopStable},
        {"gain_margin_upper", finiteOrNull(r.gainMarginUpper)},
        {"gain_margin_lower", r.gainMarginLower},
        {"gain_margin_upper_db", finiteOrNull(20.0 * std::log10(r.gainMarginUpper))},
        {"phase_margin_deg", finiteOrNull(r.phaseMargin)},
        {"phase_crossover", finiteOrNull(r.phaseCrossover)},
        {"gain_crossover", finiteOrNull(r.gainCrossover)},
        {"minimum_return_difference", r.minimumReturnDifference},
        {"return_difference_frequency", finiteOrNull(r.returnDifferenceFrequency)},
        {"disk_gain_lower", r.diskGainLower},
        {"disk_gain_upper", finiteOrNull(r.diskGainUpper)},
        {"disk_phase_margin_deg", r.diskPhaseMargin},
        {"sensitivity_peak", finiteOrNull(r.sensitivityPeak)},
        {"complementary_sensitivity_peak", finiteOrNull(r.complementarySensitivityPeak)}};
}

void storeLoop(NpzRecord& out, json::Value::Object& loops, json::Value::Array& order, const std::string& name,
               const StateSpaceModel& loop, std::size_t points, json::Value::Object extra = {}) {
    const auto r = loopRobustness(loop, points);
    const std::string prefix = "loop__" + name + "__";
    out.arrays[prefix + "frequencies"] = NpyArray::fromDoubles({r.frequencies.size()}, r.frequencies);
    storeComplex(out, prefix + "response", r.response);
    std::vector<double> d;
    for (const auto& l : r.response) d.push_back(std::abs(1.0 + l));
    out.arrays[prefix + "return_difference"] = NpyArray::fromDoubles({d.size()}, d);
    auto s = summary(r);
    for (auto& [k, v] : extra) s[k] = v;
    loops[name] = std::move(s);
    order.push_back(name);
}
}  // namespace

void appendTransferFunctions(NpzRecord& out, json::Value::Object& details, const StateSpaceModel& model,
                             const std::vector<std::string>& outputNames, const std::vector<double>& w) {
    if (!model.inputs || !model.outputs || !model.states) return;
    const std::size_t channels = model.inputs * model.outputs, F = w.size();
    std::vector<double> magnitude(F * channels), phase(F * channels);
    json::Value::Array list, names;
    for (std::size_t i = 0; i < model.outputs; ++i) {
        for (std::size_t j = 0; j < model.inputs; ++j) {
            const auto tf = transferFunction(model, i, j);
            const std::string y = i < outputNames.size() ? outputNames[i] : "y" + std::to_string(i);
            const std::string tag = "y" + std::to_string(i) + "_u" + std::to_string(j);
            storeComplex(out, "transfer_function_zeros__" + tag, tf.zeros);
            out.arrays["transfer_function_numerator__" + tag] = NpyArray::fromDoubles({tf.numerator.size()}, tf.numerator);
            json::Value::Array cancelled;
            for (const auto& [p, z] : tf.cancellations) cancelled.push_back(json::Value::Array{double(p), double(z)});
            list.push_back(json::Value::Object{{"channel", tag}, {"output", y}, {"input", "u" + std::to_string(j)},
                {"gain", tf.gain}, {"direct_gain", tf.directGain}, {"relative_degree", double(tf.relativeDegree)},
                {"zeros", double(tf.zeros.size())}, {"poles", double(tf.poles.size())},
                {"cancellations", cancelled}});
            names.push_back(y + " / u" + std::to_string(j));
            const auto b = bode(model, i, j, w);
            const std::size_t c = i * model.inputs + j;
            for (std::size_t k = 0; k < F; ++k) { magnitude[k * channels + c] = b.magnitudeDb[k]; phase[k * channels + c] = b.phaseDeg[k]; }
            if (i == 0 && j == 0) {
                out.arrays["transfer_function_denominator"] = NpyArray::fromDoubles({tf.denominator.size()}, tf.denominator);
                storeComplex(out, "transfer_function_poles", tf.poles);
            }
        }
    }
    out.arrays["bode_frequencies"] = NpyArray::fromDoubles({F}, w);
    out.arrays["bode_magnitude_db"] = NpyArray::fromDoubles({F, channels}, magnitude);
    out.arrays["bode_phase_deg"] = NpyArray::fromDoubles({F, channels}, phase);
    details["transfer_functions"] = list;
    details["bode_channels"] = names;
}

void appendReadoutObservability(NpzRecord& out, json::Value::Object& details, const StateSpaceModel& model,
                                const StateSpaceModel& readouts, const std::vector<std::string>& names,
                                const json::Value& s) {
    const std::size_t n = model.states, q = readouts.outputs;
    if (!n || !model.outputs || !q || readouts.states != n) return;
    if (!isAsymptoticallyStable(model)) throw std::invalid_argument("readout observability needs a stable model (infinite gramian)");
    const auto Wo = observabilityGramian(model);
    const auto eig = numerics::symmetricEigensystem(Wo, n);  // ascending
    std::vector<double> sv(n), dirs(n * n), signature(q * n);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t src = n - 1 - k;
        sv[k] = std::max(eig.eigenvalues[src], 0.0);
        for (std::size_t i = 0; i < n; ++i) dirs[i * n + k] = eig.eigenvectors[i * n + src];
    }
    for (std::size_t r = 0; r < q; ++r)  // readout signature of each direction: C_r v_k
        for (std::size_t k = 0; k < n; ++k) {
            double v = 0;
            for (std::size_t i = 0; i < n; ++i) v += readouts.C[r * n + i] * dirs[i * n + k];
            signature[r * n + k] = v;
        }
    out.arrays["observability_gramian_full"] = NpyArray::fromDoubles({n, n}, Wo);
    out.arrays["observability_singular_values"] = NpyArray::fromDoubles({n}, sv);
    out.arrays["observability_directions"] = NpyArray::fromDoubles({n, n}, dirs);
    out.arrays["readout_direction_signature"] = NpyArray::fromDoubles({q, n}, signature);
    const double top = sv.front();
    std::size_t rank = 0;
    for (double v : sv) rank += v > 1e-10 * top;
    // fraction of each readout row lying in the numerically unobservable subspace
    json::Value::Array readoutNames, unobservedFraction;
    for (std::size_t r = 0; r < q; ++r) {
        double total = 0, hidden = 0;
        for (std::size_t k = 0; k < n; ++k) {
            const double c = signature[r * n + k];
            total += c * c; if (k >= rank) hidden += c * c;
        }
        readoutNames.push_back(r < names.size() ? names[r] : "readout " + std::to_string(r));
        unobservedFraction.push_back(total > 0 ? std::sqrt(hidden / total) : 0.0);
    }
    json::Value::Object summary{{"readouts", readoutNames}, {"numerical_rank", double(rank)},
        {"condition_observed_subspace", rank ? json::Value(sv.front() / sv[rank - 1]) : json::Value{}},
        {"readout_norm_fraction_unobservable", unobservedFraction}};
    const double w = s.at("lqe_process_weight").number(), v = s.at("lqe_measurement_weight").number();
    if (w > 0 && v > 0) {  // the model's own claim: steady-state Kalman readout uncertainty
        const auto lqe = designLQE(model, identity(n, w), identity(model.outputs, v));
        std::vector<double> stdv(q);
        json::Value::Array list;
        for (std::size_t r = 0; r < q; ++r) {
            double acc = 0;
            for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j)
                acc += readouts.C[r * n + i] * lqe.covariance[i * n + j] * readouts.C[r * n + j];
            stdv[r] = std::sqrt(std::max(acc, 0.0)); list.push_back(stdv[r]);
        }
        out.arrays["lqe_readout_std"] = NpyArray::fromDoubles({q}, stdv);
        out.arrays["lqe_covariance"] = NpyArray::fromDoubles({n, n}, lqe.covariance);
        summary["lqe_readout_std"] = list;
        // Kalman covariance inside the observable subspace (V_o^T P V_o): its spread shows how
        // unevenly the tuned filter pins the observable directions (large W trusts each sample
        // and leaves the quadrature direction loose).
        if (rank) {
            std::vector<double> Pobs(rank * rank, 0.0);
            for (std::size_t a = 0; a < rank; ++a) for (std::size_t b = 0; b < rank; ++b)
                for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j)
                    Pobs[a * rank + b] += dirs[i * n + a] * lqe.covariance[i * n + j] * dirs[j * n + b];
            const auto pe = numerics::symmetricEigensystem(Pobs, rank).eigenvalues;  // ascending
            out.arrays["lqe_observed_covariance_eigenvalues"] = NpyArray::fromDoubles({rank}, pe);
            summary["lqe_observed_covariance_anisotropy"] = pe.front() > 0 ? json::Value(pe.back() / pe.front()) : json::Value{};
        }
    }
    details["readout_observability"] = summary;
}

void appendLoopRobustness(NpzRecord& out, json::Value::Object& details, const StateSpaceModel& model,
                          const json::Value& s) {
    if (model.inputs != 1 || !model.outputs) return;  // classical loops need a single control input
    const std::size_t n = model.states;
    const double q = s.at("lqr_state_weight").number(), r = s.at("lqr_input_weight").number();
    const double w = s.at("lqe_process_weight").number(), v = s.at("lqe_measurement_weight").number();
    const std::size_t points = static_cast<std::size_t>(s.at("loop_frequency_points").number(2000));
    if (points < 16 || points > 1000000) throw std::invalid_argument("loop_frequency_points must be 16..1e6");
    const bool control = q > 0 && r > 0, estimation = w > 0 && v > 0, siso = model.outputs == 1;
    json::Value::Object loops; json::Value::Array order;
    if (control) {
        const auto lqr = designLQR(model, identity(n, q), {r});
        storeLoop(out, loops, order, "lqr_input", stateFeedbackInputLoop(model, lqr.gain), points);
    }
    if (estimation && siso) {
        const auto lqe = designLQE(model, identity(n, w), {v});
        storeLoop(out, loops, order, "kalman_output", estimatorOutputLoop(model, lqe.gain), points);
    }
    if (control && estimation) {
        const auto V = identity(model.outputs, v);
        const auto lqg = designLQG(model, identity(n, q), {r}, identity(n, w), V);
        storeLoop(out, loops, order, "lqg_input", outputFeedbackInputLoop(model, lqg.compensator), points);
        if (siso) storeLoop(out, loops, order, "lqg_output", outputFeedbackOutputLoop(model, lqg.compensator), points);
        std::size_t k = 0;
        for (const auto& item : s.at("ltr_weights").array()) {  // Doyle-Stein recovery at the input
            const double q2 = item.number(-1);
            if (!(q2 >= 0) || !std::isfinite(q2)) throw std::invalid_argument("ltr_weights must be nonnegative");
            const auto W = loopTransferRecoveryCovariance(model, identity(n, w), q2);
            const auto design = designLQG(model, identity(n, q), {r}, W, V);
            storeLoop(out, loops, order, "ltr_" + std::to_string(k++) + "_input",
                      outputFeedbackInputLoop(model, design.compensator), points,
                      json::Value::Object{{"ltr_weight", q2}});
        }
    }
    if (!order.empty()) { details["loop_robustness"] = loops; details["loop_names"] = order; }
}

}  // namespace dynsyslab::workflow::postprocess_detail
