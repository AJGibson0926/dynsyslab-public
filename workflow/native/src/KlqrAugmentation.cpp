#include "KlqrAugmentation.hpp"

#include <dynsyslab/core/Observability.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace dynsyslab::workflow {
namespace {

std::vector<double> numbers(const json::Value& v) {
    std::vector<double> result;
    if (const auto* list = std::get_if<json::Value::Array>(&v.storage()))
        for (const auto& item : *list) result.push_back(item.number());
    else if (!v.isNull()) result.push_back(v.number());
    return result;
}

std::vector<std::size_t> indices(const json::Value& v, const char* name) {
    std::vector<std::size_t> result;
    for (double x : numbers(v)) {
        if (!(x >= 0) || x != std::floor(x)) throw std::invalid_argument(std::string(name) + " must be nonnegative integers");
        result.push_back(static_cast<std::size_t>(x));
    }
    return result;
}

// Scalar, diagonal or full rows x columns matrix.
std::vector<double> matrix(const json::Value& v, std::size_t rows, std::size_t columns, double fallback, const char* name) {
    auto values = numbers(v);
    if (values.size() == rows * columns && !(rows == 1 && columns == 1 && values.empty())) return values;
    std::vector<double> result(rows * columns, 0.0);
    if (values.empty()) values.assign(1, fallback);
    if (values.size() == 1 && rows == columns) { for (std::size_t k = 0; k < rows; ++k) result[k * columns + k] = values[0]; return result; }
    if (values.size() == 1) { std::fill(result.begin(), result.end(), values[0]); return result; }
    if (rows == columns && values.size() == rows) { for (std::size_t k = 0; k < rows; ++k) result[k * columns + k] = values[k]; return result; }
    throw std::invalid_argument(std::string(name) + " must be scalar, diagonal or " + std::to_string(rows) + " x " + std::to_string(columns));
}

core::LiftedOutputChannel channel(const core::KoopmanModel& model, const json::Value& spec, const char* name) {
    const auto reconstruction = spec.at("reconstruction").string("state");
    if (!model.hasReconstruction(reconstruction))
        throw std::invalid_argument(std::string(name) + ": the KROM has no reconstruction " + reconstruction);
    const auto& readout = model.reconstruction(reconstruction);
    auto components = indices(spec.at("components"), "components");
    if (components.empty()) components.push_back(0);
    auto targets = indices(spec.at("target_components"), "target_components");
    if (targets.empty()) targets = components;
    if (targets.size() != components.size())
        throw std::invalid_argument(std::string(name) + ": target_components must match components");
    const auto full = core::reconstructionSplitMatrix(model, readout);
    const auto n = model.splitDimension();
    core::LiftedOutputChannel result;
    for (auto k : components) {
        if (k >= readout.outputs) throw std::invalid_argument(std::string(name) + ": component is out of range");
        result.outputMatrix.insert(result.outputMatrix.end(), full.begin() + std::ptrdiff_t(k * n), full.begin() + std::ptrdiff_t((k + 1) * n));
        result.outputOffset.push_back(readout.offset.empty() ? 0.0 : readout.offset[k]);
    }
    result.targetComponents = targets;
    return result;
}

std::pair<double, double> limits(const json::Value& v) {
    const auto values = numbers(v);
    const double inf = std::numeric_limits<double>::infinity();
    if (values.empty()) return {-inf, inf};
    if (values.size() == 1) return {-std::abs(values[0]), std::abs(values[0])};
    if (values.size() == 2) return {values[0], values[1]};
    throw std::invalid_argument("output_limits must be L or [low, high]");
}

}  // namespace

bool klqrIsStateful(const json::Value& settings) {
    return !settings.at("integral_action").isNull() || !settings.at("pid_correction").isNull();
}

void applyKlqrAugmentations(core::KoopmanLQRController& controller, const core::KoopmanModel& model,
                            const json::Value& settings) {
    const auto m = model.inputDimension();
    if (const auto& spec = settings.at("integral_action"); !spec.isNull()) {
        core::KoopmanIntegralAction action;
        action.channel = channel(model, spec, "integral_action");
        const auto p = action.channel.targetComponents.size();
        action.penalty = matrix(spec.at("penalty"), p, p, 1.0, "integral_action penalty");
        action.integratorLimit = numbers(spec.at("integrator_limit"));
        controller.setIntegralAction(std::move(action));
    }
    if (const auto& spec = settings.at("pid_correction"); !spec.isNull()) {
        core::KoopmanPIDCorrection pid;
        pid.channel = channel(model, spec, "pid_correction");
        const auto p = pid.channel.targetComponents.size();
        pid.proportional = matrix(spec.at("Kp"), m, p, 0.0, "pid_correction Kp");
        pid.integral = matrix(spec.at("Ki"), m, p, 0.0, "pid_correction Ki");
        pid.derivative = matrix(spec.at("Kd"), m, p, 0.0, "pid_correction Kd");
        pid.derivativeFilterTime = spec.at("derivative_filter_time").number(0.0);
        pid.integratorLimit = numbers(spec.at("integrator_limit"));
        std::tie(pid.outputMinimum, pid.outputMaximum) = limits(spec.at("output_limits"));
        pid.conditionalIntegration = spec.at("conditional_integration").boolean(true);
        controller.setPIDCorrection(std::move(pid));
    }
}

json::Value klqrAugmentationEvidence(const json::Value& settings) {
    json::Value::Object evidence;
    if (!settings.at("integral_action").isNull()) evidence["integral_action"] = settings.at("integral_action");
    if (!settings.at("pid_correction").isNull()) evidence["pid_correction"] = settings.at("pid_correction");
    return json::Value(std::move(evidence));
}

}  // namespace dynsyslab::workflow
