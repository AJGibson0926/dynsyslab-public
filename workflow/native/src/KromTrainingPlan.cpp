#include <dynsyslab/workflow/KromTrainingPlan.hpp>

#include <dynsyslab/core/FileDigest.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/Targets.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/core/TrainingDomain.hpp>
#include <dynsyslab/workflow/ControlTargets.hpp>

#include "KromDomainPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::workflow {
namespace {

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

json::Value vector(const std::vector<double>& values) {
    json::Value::Array result;
    for (double v : values) result.emplace_back(v);
    return json::Value(std::move(result));
}

json::Value domainJson(const core::AxisAlignedStateDomain& d) {
    return json::Value::Object{{"minimum", vector(d.minimum)}, {"maximum", vector(d.maximum)},
                               {"coordinate_system", d.coordinateSystem}};
}

void include(core::AxisAlignedStateDomain& domain, const std::vector<double>& low, const std::vector<double>& high) {
    if (domain.minimum.empty()) { domain.minimum = low; domain.maximum = high; return; }
    if (low.size() != domain.minimum.size()) throw std::invalid_argument("training requirement dimension mismatch");
    for (std::size_t i = 0; i < low.size(); ++i) {
        domain.minimum[i] = std::min(domain.minimum[i], low[i]);
        domain.maximum[i] = std::max(domain.maximum[i], high[i]);
    }
}

}  // namespace

json::Value designKromTrainingPlan(const json::Value& request) {
    core::SystemConfiguration system;
    const auto kromPath = request.at("krom").string();
    if (!request.at("system").isNull()) {
        system.identifier = request.at("system").at("name").string();
        system.preset = request.at("system").at("preset").string();
        for (const auto& [name, value] : request.at("system").at("parameters").object()) system.parameters[name] = value.number();
    } else if (!kromPath.empty()) {
        system = core::loadKoopmanModel(core::wide(kromPath)).systemConfiguration();
    } else throw std::invalid_argument("a training plan needs a krom or a system");
    auto ode = core::createSystem(system);
    const std::size_t n = ode->descriptor().stateDimension;
    const std::string coordinates = ode->descriptor().preferredCoordinates;
    core::AxisAlignedStateDomain current{{}, {}, coordinates};
    if (!request.at("current_domain").isNull()) {
        current.minimum = numbers(request.at("current_domain").at("minimum"));
        current.maximum = numbers(request.at("current_domain").at("maximum"));
    } else {
        auto [low, high] = detail::kromTrainingDomain(kromPath, json::Value{}, n);
        if (low.empty()) throw std::invalid_argument("the K-ROM records no training domain; give current_domain");
        current.minimum = low; current.maximum = high;
    }
    // The requirements: explicit domains, target envelopes with probes, the RPE preset.
    json::Value::Array requirements = request.at("requirements").array();
    if (request.at("preset").string() == "rpe_targets") {
        std::vector<std::string> names;
        for (const auto& v : request.at("targets").array()) names.push_back(v.string());
        if (names.empty()) for (const auto& t : core::rpeControlTargets()) names.push_back(t.name);
        for (const auto& name : names)
            requirements.emplace_back(json::Value::Object{{"name", name},
                {"target", json::Value::Object{{"kind", "rpe_canonical"}, {"name", name}}},
                {"probes", request.at("probes_per_target").number(33)}});
    }
    core::AxisAlignedStateDomain required{{}, {}, coordinates};
    std::vector<double> anchors;
    std::size_t anchorCount = 0;
    for (const auto& row : request.at("anchors").array()) {
        const auto a = numbers(row);
        if (a.size() != n) throw std::invalid_argument("anchor dimension mismatch");
        anchors.insert(anchors.end(), a.begin(), a.end()); ++anchorCount;
    }
    json::Value::Array probes, summaries;
    constexpr double twoPi = 6.283185307179586476925286766559;
    for (const auto& item : requirements) {
        const auto name = item.at("name").string();
        if (!item.at("domain").isNull()) {
            include(required, numbers(item.at("domain").at("minimum")), numbers(item.at("domain").at("maximum")));
            summaries.emplace_back(json::Value::Object{{"name", name}, {"kind", "domain"}});
            continue;
        }
        const auto& target = item.at("target");
        if (target.isNull()) throw std::invalid_argument("a training requirement needs a domain or a target");
        const double start = item.at("start").number(0.0), stop = item.at("stop").number(twoPi);
        const auto envelope = controlTargetEnvelope(target, n, start, stop);
        include(required, envelope.lower, envelope.upper);
        const auto evaluate = resolveControlTarget(target, n);
        const auto count = static_cast<std::size_t>(item.at("probes").number(33));
        if (count < 5) throw std::invalid_argument("a target requirement needs at least five probes");
        for (std::size_t k = 0; k < count; ++k) {
            const double t = start + (stop - start) * double(k) / double(count);
            const auto state = evaluate(t);
            for (double v : state) if (!std::isfinite(v)) throw std::invalid_argument("target probe is not finite");
            anchors.insert(anchors.end(), state.begin(), state.end()); ++anchorCount;
            probes.emplace_back(json::Value::Object{{"requirement", name}, {"parameter", t}, {"state", vector(state)}});
        }
        summaries.emplace_back(json::Value::Object{{"name", name}, {"kind", "target"}, {"envelope", envelope.method},
                                                   {"minimum", vector(envelope.lower)}, {"maximum", vector(envelope.upper)}});
    }
    if (required.minimum.empty()) throw std::invalid_argument("a training plan needs at least one requirement");
    core::StateDomainExpansionOptions options;
    for (const auto& v : request.at("points_per_axis").array()) options.pointsPerAxis.push_back(static_cast<std::size_t>(v.number()));
    options.marginFraction = request.at("margin_fraction").number(0.1);
    options.anchors = anchors; options.anchorCount = anchorCount;
    const auto plan = core::planStateDomainExpansion(current, required, options);
    std::vector<double> positive = numbers(request.at("positive_axes"));
    if (positive.empty() && request.at("preset").string() == "rpe_targets") positive.push_back(0.0);
    for (double axis : positive)
        if (plan.planned.minimum.at(static_cast<std::size_t>(axis)) <= 0.0)
            throw std::invalid_argument("planned domain crosses the singular boundary on axis "
                                        + std::to_string(static_cast<std::size_t>(axis)));
    // The solve job that runs and screens the plan.
    json::Value::Array initial;
    for (std::size_t k = 0; k < plan.initialConditionCount; ++k)
        initial.push_back(vector(std::vector<double>(plan.initialConditions.begin() + k * n,
                                                     plan.initialConditions.begin() + (k + 1) * n)));
    json::Value::Array anchorRows;
    for (std::size_t k = 0; k < anchorCount; ++k)
        anchorRows.push_back(vector(std::vector<double>(anchors.begin() + k * n, anchors.begin() + (k + 1) * n)));
    json::Value::Array points;
    for (auto v : plan.pointsPerAxis) points.emplace_back(double(v));
    const auto& time = request.at("time");
    json::Value solve = json::Value::Object{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "solve"},
        {"label", request.at("label").string("krom-training-expansion")},
        {"system", json::Value::Object{{"name", system.identifier}, {"preset", system.preset},
                                       {"parameters", request.at("system").at("parameters")}}},
        {"time", json::Value::Object{{"start", 0.0}, {"stop", time.at("stop").number(25.0)},
                                     {"output_step", time.at("output_step").number(0.001)}, {"coordinate_system", "automatic"}}},
        {"initial_conditions", json::Value::Object{{"mode", "explicit"}, {"count", double(plan.initialConditionCount)},
            {"clusters", 1.0}, {"seed", 0.0}, {"values", json::Value(std::move(initial))}, {"ranges", json::Value::Array{}}}},
        {"controller", json::Value::Object{{"kind", "none"}, {"settings", json::Value::Object{}}}},
        {"solver", request.at("solver").isNull() ? json::Value(json::Value::Object{{"method", "DOP853"}, {"rtol", 1e-9}, {"atol", 1e-12}})
                                                  : request.at("solver")},
        {"preflight", json::Value::Object{{"enabled", true}, {"kind", "krom_training_domain"}, {"source_krom", kromPath},
            {"current_domain", domainJson(current)}, {"required_domain", domainJson(required)},
            {"points_per_axis", json::Value(std::move(points))}, {"margin_fraction", options.marginFraction},
            {"anchors", json::Value(std::move(anchorRows))},
            {"screening", request.at("screening").isNull() ? json::Value(json::Value::Object{}) : request.at("screening")}}},
        {"archive", json::Value::Object{{"enabled", true}, {"data_role", "free_training"}}}};
    json::Value::Object result{{"schema", "dynsyslab.krom-training-plan/v1"}, {"source_krom", kromPath},
        {"system", system.identifier}, {"current", domainJson(current)}, {"required", domainJson(required)},
        {"planned", domainJson(plan.planned)}, {"requirements", json::Value(std::move(summaries))},
        {"probes", json::Value(std::move(probes))}, {"anchor_count", double(anchorCount)},
        {"initial_condition_count", double(plan.initialConditionCount)}, {"scheduled", false}, {"executed", false}};
    result["plan_fingerprint_sha256"] = core::sha256Text(json::dump(json::Value(result)));
    result["solve_job"] = std::move(solve);
    return json::Value(std::move(result));
}

}  // namespace dynsyslab::workflow
