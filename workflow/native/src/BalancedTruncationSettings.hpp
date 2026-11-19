#pragma once
// Shared parser for the KLQE/KLQR stage option "balanced_truncation".
//
//   "balanced_truncation": null | false          -> full-order design (default)
//   "balanced_truncation": true                  -> numerically minimal order
//   "balanced_truncation": 4                     -> reduced order 4
//   "balanced_truncation": {"order": 4, "error_tolerance": 0.01,
//                           "method": "truncation" | "residualization",
//                           "stability_margin": 1e-8}
//
// The order counts retained unstable/marginal modes. With "error_tolerance" and
// no order, the smallest order whose bound 2 * sum(discarded HSV) meets it is used.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/BalancedTruncation.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>

namespace dynsyslab::workflow {

inline std::optional<core::BalancedReductionOptions> parseBalancedTruncation(const json::Value& v) {
    if (v.isNull()) return std::nullopt;
    if (std::holds_alternative<bool>(v.storage())) {
        if (!v.boolean()) return std::nullopt;
        return core::BalancedReductionOptions{};
    }
    core::BalancedReductionOptions o;
    const auto order = [](double value) {
        if (!(value >= 1) || value != std::floor(value) || !std::isfinite(value))
            throw std::invalid_argument("balanced_truncation order must be a positive integer");
        return static_cast<std::size_t>(value);
    };
    if (!std::holds_alternative<json::Value::Object>(v.storage())) {
        o.order = order(v.number(-1));
        return o;
    }
    if (!v.at("order").isNull()) o.order = order(v.at("order").number(-1));
    o.errorTolerance = v.at("error_tolerance").number(0.0);
    if (!(o.errorTolerance >= 0) || !std::isfinite(o.errorTolerance))
        throw std::invalid_argument("balanced_truncation error_tolerance must be finite and nonnegative");
    const auto method = v.at("method").string("truncation");
    if (method == "residualization") o.method = core::BalancedReductionMethod::residualization;
    else if (method != "truncation")
        throw std::invalid_argument("balanced_truncation method must be truncation or residualization");
    o.stabilityMargin = v.at("stability_margin").number(o.stabilityMargin);
    if (!(o.stabilityMargin >= 0) || !std::isfinite(o.stabilityMargin))
        throw std::invalid_argument("balanced_truncation stability_margin must be finite and nonnegative");
    return o;
}

inline json::Value balancedTruncationEvidence(const core::BalancedReduction& r) {
    json::Value::Array hsv;
    for (double s : r.hankelSingularValues) hsv.push_back(json::Value(s));
    return json::Value(json::Value::Object{
        {"method", json::Value(std::string(core::balancedReductionMethodName(r.method)))},
        {"order", json::Value(double(r.order))}, {"full_order", json::Value(double(r.fullOrder))},
        {"stable_order", json::Value(double(r.stableOrder))}, {"unstable_order", json::Value(double(r.unstableOrder))},
        {"error_bound", json::Value(r.errorBound)}, {"lower_bound", json::Value(r.lowerBound)},
        {"hankel_singular_values", json::Value(std::move(hsv))}});
}

// The requested options, for evidence recorded before any gain is solved (KLQR
// solves its gain at run time, when the input map is known).
inline json::Value balancedTruncationRequest(const core::BalancedReductionOptions& o) {
    return json::Value(json::Value::Object{
        {"method", json::Value(std::string(core::balancedReductionMethodName(o.method)))},
        {"order", o.order ? json::Value(double(o.order)) : json::Value("numerically minimal or error_tolerance")},
        {"error_tolerance", json::Value(o.errorTolerance)}, {"stability_margin", json::Value(o.stabilityMargin)}});
}

// Archive arrays written by a KLQR controller with balanced_truncation.
inline void storeKlqrReduction(const core::KoopmanLQRController& controller, core::NpzRecord& archive) {
    if (!controller.balancedReduction() || !controller.lastReducedOrder()) return;
    archive.arrays["klqr_reduced_order"] = core::NpyArray::fromDoubles({1}, {double(controller.lastReducedOrder())});
    archive.arrays["klqr_reduction_error_bound"] =
        core::NpyArray::fromDoubles({1}, {controller.lastReductionErrorBound()});
    const auto& last = controller.lastSolution();
    archive.arrays["klqr_gain"] = core::NpyArray::fromDoubles(
        {last.inputDimension, last.stateDimension}, last.gain);
    archive.arrays["klqr_closed_loop_poles"] = core::NpyArray::fromComplex(
        {last.closedLoopPoles.size()}, last.closedLoopPoles);
}

}  // namespace dynsyslab::workflow
