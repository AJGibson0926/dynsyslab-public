#include "KromDomainPolicy.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/TrainingDomain.hpp>

#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

struct Tally {
    std::size_t evaluations{}, stateOutside{}, targetOutside{};
    double stateExcursion{}, targetExcursion{};
};

}  // namespace

std::pair<std::vector<double>, std::vector<double>> kromTrainingDomain(
    const std::string& path, const json::Value& settings, std::size_t n) {
    const auto& given = settings.at("training_domain");
    if (!given.isNull()) {
        auto low = numbers(given.at("minimum")), high = numbers(given.at("maximum"));
        if (low.size() != n || high.size() != n)
            throw std::invalid_argument("training_domain needs one minimum and maximum per state");
        return {low, high};
    }
    if (path.empty()) return {};
    const auto archive = dynsyslab::core::readNpz(dynsyslab::core::wide(path));
    if (archive.contains("training_minimum") && archive.contains("training_maximum")) {
        auto low = archive.at("training_minimum").doubles(), high = archive.at("training_maximum").doubles();
        if (low.size() == n && high.size() == n) return {low, high};
    }
    if (!archive.contains("training_states")) return {};
    const auto& states = archive.at("training_states");
    if (states.shape.size() != 2 || states.shape[1] != n || !states.shape[0]) return {};
    const auto values = states.doubles();
    std::vector<double> low(n, std::numeric_limits<double>::infinity()), high(n, -low[0]);
    for (std::size_t row = 0; row < states.shape[0]; ++row)
        for (std::size_t axis = 0; axis < n; ++axis) {
            low[axis] = std::min(low[axis], values[row * n + axis]);
            high[axis] = std::max(high[axis], values[row * n + axis]);
        }
    return {low, high};
}

void applyTargetDomainPolicy(ControllerSetup& setup, const json::Value& settings, const std::string& path,
                             std::size_t n, dynsyslab::core::StateTargetFunction target, std::size_t count) {
    const auto policy = settings.at("target_domain_policy").string("warn");
    if (policy != "warn" && policy != "error" && policy != "ignore")
        throw std::invalid_argument("target_domain_policy must be error, warn or ignore");
    auto& evidence = setup.modelBindingEvidence;
    if (!std::holds_alternative<json::Value::Object>(evidence.storage())) evidence = json::Value(json::Value::Object{});
    evidence.object()["target_domain_policy"] = policy;
    if (policy == "ignore" || !(setup.policy.value || setup.policy.vector)) return;
    auto [low, high] = kromTrainingDomain(path, settings, n);
    evidence.object()["training_domain_available"] = !low.empty();
    if (low.empty()) return;
    dynsyslab::core::AxisAlignedStateDomain domain{low, high, "model"};
    auto tallies = std::make_shared<std::vector<Tally>>(count);
    const bool fail = policy == "error";
    const auto guard = [domain, target, tallies, fail, n](std::size_t trajectory, double time, const double* state) {
        auto& tally = tallies->at(trajectory);
        ++tally.evaluations;
        const auto check = [&](const double* x, const char* what, std::size_t& outside, double& worst) {
            const auto report = domain.report(x);
            if (report.inside) return;
            ++outside;
            worst = std::max(worst, report.maximumRelativeExcursion);
            if (fail) {
                std::string axes;
                for (auto axis : report.violatingAxes) axes += (axes.empty() ? "" : ",") + std::to_string(axis);
                throw std::runtime_error(std::string("Koopman control ") + what + " leaves the K-ROM training domain at t="
                    + std::to_string(time) + " (axes " + axes + ", relative excursion "
                    + std::to_string(report.maximumRelativeExcursion)
                    + "); set target_domain_policy warn or ignore only for an explicit extrapolation study");
            }
        };
        check(state, "state", tally.stateOutside, tally.stateExcursion);
        if (target) {
            const auto reference = target(time);
            if (reference.size() == n) check(reference.data(), "target", tally.targetOutside, tally.targetExcursion);
        }
    };
    if (setup.policy.value)
        setup.policy.value = [law = std::move(setup.policy.value), guard](std::size_t k, double time, const double* state) {
            guard(k, time, state);
            return law(k, time, state);
        };
    else
        setup.policy.vector = [law = std::move(setup.policy.vector), guard](std::size_t k, double time,
                                                                           const double* state, double* input) {
            guard(k, time, state);
            law(k, time, state, input);
        };
    auto previous = std::move(setup.appendDiagnostics);
    setup.appendDiagnostics = [previous, tallies, low, high](std::size_t trajectory, dynsyslab::core::NpzRecord& archive) {
        if (previous) previous(trajectory, archive);
        const auto& t = tallies->at(trajectory);
        archive.arrays["domain_report_counts"] = dynsyslab::core::NpyArray::fromDoubles(
            {3}, {double(t.evaluations), double(t.stateOutside), double(t.targetOutside)});
        archive.arrays["domain_report_maximum_excursion"] = dynsyslab::core::NpyArray::fromDoubles(
            {2}, {t.stateExcursion, t.targetExcursion});
        archive.arrays["domain_report_bounds"] = dynsyslab::core::NpyArray::fromDoubles(
            {2, low.size()}, [&] { auto v = low; v.insert(v.end(), high.begin(), high.end()); return v; }());
    };
}

}  // namespace dynsyslab::workflow::detail
