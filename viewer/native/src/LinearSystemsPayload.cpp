#include "LinearSystemsPayload.hpp"
#include "LinearSystemsLoopPayload.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace dsl {
namespace {

using dynsyslab::core::NpyArray;
using dynsyslab::core::NpzRecord;
namespace json = dynsyslab::core::json;

constexpr const char* kFull = "#2563eb";     // blue
constexpr const char* kReduced = "#ef4444";  // red
constexpr const char* kGold = "#d4a017";     // gold
constexpr double kPi = 3.14159265358979323846;

std::string number(double value) {
    std::ostringstream text; text.precision(3); text << value; return text.str();
}

std::string label(const json::Value& meta, const std::filesystem::path& source) {
    return meta.at("label").string(source.stem().u8string());
}

Trace points(std::string name, std::string color, bool marker = true) {
    Trace t; t.label = std::move(name); t.color = std::move(color); t.marker = marker; return t;
}

// Poles and zeros are a scatter: markers only, never joined in storage order.
Trace complexPoints(const NpzRecord& record, const char* name, std::string text, std::string color) {
    Trace t = points(std::move(text), std::move(color));
    t.lineStyle = "none";
    if (!record.contains(name)) return t;
    const auto v = record.at(name).doubles();
    for (std::size_t k = 0; k + 1 < v.size(); k += 2) { t.x.push_back(v[k]); t.y.push_back(v[k + 1]); }
    return t;
}

Payload hankel(const std::filesystem::path& source, const std::string& view, const NpzRecord& r,
               const json::Value& meta) {
    const auto hsv = r.at("hankel_singular_values").doubles();
    const auto bound = r.at("hankel_error_bound").doubles();
    const auto& red = meta.at("reduction");
    const double order = red.at("stable_order").number();
    Figure f{"hankel_spectrum", "Hankel singular values and balanced-truncation bound",
             "index k (stable part)", "value", false};
    f.valueScale = "log";
    // Values below 1e-12 of the largest are round-off; plotting them would stretch
    // the log axis over dozens of empty decades.
    double floor = 0;
    for (double v : hsv) floor = std::max(floor, v);
    floor *= 1e-12;
    Trace s = points("Hankel singular value sigma_k", kFull);
    for (std::size_t k = 0; k < hsv.size(); ++k)
        if (hsv[k] > floor) { s.x.push_back(double(k + 1)); s.y.push_back(hsv[k]); }
    Trace b = points("error bound 2 sum_{i>r} sigma_i at order r", kGold, false);
    for (std::size_t k = 0; k + 1 < bound.size(); ++k)
        if (bound[k] > floor) { b.x.push_back(double(k)); b.y.push_back(bound[k]); }
    f.traces.push_back(std::move(s)); f.traces.push_back(std::move(b));
    const auto chosenIndex = static_cast<std::size_t>(order);
    if (chosenIndex < bound.size() && bound[chosenIndex] > 0) {
        // Vertical line at the chosen order, from the smallest plotted value to the largest bound.
        double low = bound[chosenIndex], high = bound[chosenIndex];
        for (double v : hsv) if (v > floor) low = std::min(low, v);
        for (double v : bound) if (v > 0) high = std::max(high, v);
        Trace chosen = points("reduced order r = " + number(order) + " (" + red.at("method").string() + ")", kReduced, false);
        chosen.lineStyle = "dash";
        chosen.x = {order, order}; chosen.y = {low, high};
        f.traces.push_back(std::move(chosen));
    }
    const double unstable = red.at("unstable_order").number();
    if (unstable > 0) f.title += " (" + number(unstable) + " unstable/marginal modes kept exactly)";
    Payload p{source.u8string(), view, label(meta, source), 1};
    p.figures.push_back(std::move(f));
    return p;
}

// Controllability: states reachable with unit input energy, semi-axes sqrt(lambda).
// Observability: initial states producing unit output energy, semi-axes 1/sqrt(lambda).
Trace ellipse(const std::vector<double>& eig, bool inverse, std::string text, std::string color) {
    Trace t = points(std::move(text), std::move(color), false);
    auto axis = [&](std::size_t i) {
        const double v = i < eig.size() ? eig[i] : 0.0;
        return inverse ? (v > 0 ? 1.0 / std::sqrt(v) : 0.0) : std::sqrt(std::max(v, 0.0));
    };
    const double a = axis(0), b = axis(1);
    if (!(a > 0) || !std::isfinite(a) || !std::isfinite(b)) return t;
    for (int k = 0; k <= 360; ++k) {
        const double th = 2.0 * kPi * k / 360.0;
        t.x.push_back(a * std::cos(th)); t.y.push_back(b * std::sin(th));
    }
    return t;
}

Trace indexed(const NpzRecord& r, const char* name, std::string text, std::string color) {
    Trace t = points(std::move(text), std::move(color));
    if (!r.contains(name)) return t;
    const auto v = r.at(name).doubles();
    double floor = 0;  // drop round-off values below 1e-12 of this series' largest
    for (double x : v) floor = std::max(floor, x);
    floor *= 1e-12;
    for (std::size_t k = 0; k < v.size(); ++k) if (v[k] > floor) { t.x.push_back(double(k + 1)); t.y.push_back(v[k]); }
    return t;
}

Payload gramians(const std::filesystem::path& source, const std::string& view, const NpzRecord& r,
                 const json::Value& meta) {
    Payload p{source.u8string(), view, label(meta, source), 2};
    Figure f{"gramian_spectrum", "Gramian spectra (" + meta.at("gramian_coordinates").string("model coordinates") + ")",
             "principal direction", "eigenvalue", false};
    f.valueScale = "log";
    f.traces.push_back(indexed(r, "controllability_gramian_eigenvalues", "controllability gramian Wc", kFull));
    f.traces.push_back(indexed(r, "observability_gramian_eigenvalues", "observability gramian Wo", kReduced));
    f.traces.push_back(indexed(r, "hankel_singular_values", "Hankel singular values sqrt(eig(Wc Wo))", kGold));
    p.figures.push_back(std::move(f));
    if (r.contains("controllability_gramian_eigenvalues") && r.contains("observability_gramian_eigenvalues")) {
        Figure e{"gramian_ellipses", "Gramian ellipses in each gramian's two strongest principal directions",
                 "principal direction 1", "principal direction 2", true};
        e.traces.push_back(ellipse(r.at("controllability_gramian_eigenvalues").doubles(), false,
                                   "reachable, unit input energy (Wc)", kFull));
        e.traces.push_back(ellipse(r.at("observability_gramian_eigenvalues").doubles(), true,
                                   "unit output energy (Wo)", kReduced));
        p.figures.push_back(std::move(e));
    }
    return p;
}

Payload frequency(const std::filesystem::path& source, const std::string& view, const NpzRecord& r,
                  const json::Value& meta) {
    const auto w = r.at("frequencies").doubles();
    std::vector<double> x;
    for (double v : w) x.push_back(std::log10(v));
    Payload p{source.u8string(), view, label(meta, source), 1};
    const double bound = meta.at("reduction").at("error_bound").number();
    const double actual = meta.at("reduction_error_hinf").number();
    Figure f{"frequency_response", "Singular values: full vs reduced (H-infinity error " + number(actual) +
             ", balanced bound " + number(bound) + ")", "log10 omega (rad per model time unit)", "singular value", false};
    f.valueScale = "log";
    auto add = [&](const char* name, const std::string& text, const char* color, const char* style) {
        if (!r.contains(name)) return;
        const auto& a = r.at(name);
        const std::size_t k = a.shape.size() > 1 ? a.shape[1] : 1;
        const auto v = a.doubles();
        for (std::size_t j = 0; j < k; ++j) {
            Trace t = points(k > 1 ? text + " sigma_" + std::to_string(j + 1) : text, color, false);
            t.lineStyle = style; t.x = x;
            for (std::size_t i = 0; i < w.size(); ++i) t.y.push_back(std::max(v[i * k + j], 1e-300));
            f.traces.push_back(std::move(t));
        }
    };
    add("full_singular_values", "full model", kFull, "solid");
    add("reduced_singular_values", "reduced model", kReduced, "dash");
    p.figures.push_back(std::move(f));
    if (r.contains("full_phase")) {
        Figure ph{"phase_response", "Phase (single input, single output)", "log10 omega (rad per model time unit)",
                  "phase (deg)", false};
        auto phase = [&](const char* name, const char* text, const char* color, const char* style) {
            if (!r.contains(name)) return;
            Trace t = points(text, color, false); t.lineStyle = style; t.x = x;
            for (double v : r.at(name).doubles()) t.y.push_back(v * 180.0 / kPi);
            ph.traces.push_back(std::move(t));
        };
        phase("full_phase", "full model", kFull, "solid");
        phase("reduced_phase", "reduced model", kReduced, "dash");
        p.figures.push_back(std::move(ph));
        p.columns = 2;
    }
    return p;
}

Payload poleZero(const std::filesystem::path& source, const std::string& view, const NpzRecord& r,
                 const json::Value& meta) {
    const std::string axis = meta.at("time_domain").string() == "discrete" ? "mu" : "lambda";
    Payload p{source.u8string(), view, label(meta, source), 1};
    Figure f{"pole_zero", "Poles and transmission zeros", "Re(" + axis + ")", "Im(" + axis + ")", false};
    f.traces.push_back(complexPoints(r, "system_poles", "full-model poles", kFull));
    f.traces.push_back(complexPoints(r, "reduced_poles", "reduced-model poles", kReduced));
    f.traces.push_back(complexPoints(r, "transmission_zeros", "transmission zeros", kGold));
    if (!r.contains("transmission_zeros"))  // non-square: zeros of each input/output channel
        for (const auto& tf : meta.at("transfer_functions").array()) {
            const auto name = "transfer_function_zeros__" + tf.at("channel").string();
            auto z = complexPoints(r, name.c_str(), [&] { auto t = "zeros of " + tf.at("output").string() + " / " + tf.at("input").string();
                std::replace(t.begin(), t.end(), '_', ' '); return t; }(), kGold);
            z.opacity = 0.7; f.traces.push_back(std::move(z));
        }
    p.figures.push_back(std::move(f));
    if (r.contains("lqr_closed_loop_poles") || r.contains("lqe_observer_poles")) {
        Figure g{"feedback_poles", "Feedback poles on the full model", "Re(" + axis + ")", "Im(" + axis + ")", false};
        g.traces.push_back(complexPoints(r, "lqr_closed_loop_poles", "full LQR closed loop", kFull));
        g.traces.push_back(complexPoints(r, "reduced_lqr_lifted_poles", "reduced LQR gain lifted to the full model", kReduced));
        g.traces.push_back(complexPoints(r, "lqe_observer_poles", "full LQE observer", kGold));
        auto lifted = complexPoints(r, "reduced_lqe_lifted_poles", "reduced LQE gain lifted to the full model", kGold);
        lifted.opacity = 0.5;
        g.traces.push_back(std::move(lifted));
        p.figures.push_back(std::move(g));
        p.columns = 2;
    }
    if (r.contains("lqg_closed_loop_poles")) {
        Figure h{"output_feedback_poles", "Output feedback (LQG: LQE feeding LQR) closed loop on the full model",
                 "Re(" + axis + ")", "Im(" + axis + ")", false};
        h.traces.push_back(complexPoints(r, "lqg_closed_loop_poles", "full-order compensator", kFull));
        h.traces.push_back(complexPoints(r, "reduced_lqg_full_closed_loop_poles",
                                         "reduced-order compensator on the full model", kReduced));
        p.figures.push_back(std::move(h));
        p.columns = 2;
    }
    for (auto& figure : p.figures)
        figure.traces.erase(std::remove_if(figure.traces.begin(), figure.traces.end(),
            [](const Trace& t) { return t.x.empty(); }), figure.traces.end());
    return p;
}

}  // namespace

Payload buildLinearSystemsPayload(const std::filesystem::path& source, const std::string& view,
                                  const NpzRecord& record, const json::Value& metadata) {
    if (metadata.at("schema").string() != "dynsyslab.linear-systems/v1")
        throw std::runtime_error("archive is not a linear-systems analysis");
    if (view == "hankel_spectrum") return hankel(source, view, record, metadata);
    if (view == "gramian_spectrum") return gramians(source, view, record, metadata);
    if (view == "frequency_response") return frequency(source, view, record, metadata);
    if (view == "pole_zero") return poleZero(source, view, record, metadata);
    if (isClassicalLinearSystemsView(view)) return buildClassicalLinearSystemsPayload(source, view, record, metadata);
    throw std::runtime_error("unknown linear-systems view: " + view);
}

}  // namespace dsl
