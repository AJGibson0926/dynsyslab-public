#include "TrainingPreflight.hpp"

#include <dynsyslab/core/TrainingDomain.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {
namespace json = dynsyslab::core::json;

std::vector<double> numbers(const json::Value& value) {
    std::vector<double> result;
    for (const auto& item : value.array()) result.push_back(item.number());
    return result;
}

std::vector<std::size_t> sizes(const json::Value& value) {
    std::vector<std::size_t> result;
    for (const auto& item : value.array()) {
        const double number = item.number();
        if (number < 0.0 || number != std::floor(number))
            throw std::invalid_argument("points_per_axis values must be non-negative integers");
        result.push_back(static_cast<std::size_t>(number));
    }
    return result;
}

dynsyslab::core::AxisAlignedStateDomain domain(
    const json::Value& source, const std::string& coordinates) {
    return {numbers(source.at("minimum")), numbers(source.at("maximum")),
            source.at("coordinate_system").string(coordinates)};
}

std::vector<double> flattenedRows(
    const json::Value& value, std::size_t dimension, std::size_t& rows) {
    std::vector<double> result;
    rows = 0;
    for (const auto& row : value.array()) {
        const auto values = numbers(row);
        if (values.size() != dimension)
            throw std::invalid_argument("training-domain anchor dimension mismatch");
        result.insert(result.end(), values.begin(), values.end());
        ++rows;
    }
    return result;
}

dynsyslab::core::IntegrationMethod method(const std::string& value) {
    return dynsyslab::core::integrationMethodFromName(value);
}

json::Value vectorValue(const std::vector<double>& values) {
    json::Value::Array result;
    for (const double value : values) result.emplace_back(value);
    return json::Value(std::move(result));
}

json::Value entryValue(const dynsyslab::core::StateDomainScreeningEntry& entry,
                       std::size_t index) {
    json::Value::Object result;
    result["index"] = static_cast<double>(index);
    result["accepted"] = entry.accepted;
    result["reason"] = entry.reason;
    result["samples"] = static_cast<double>(entry.samples);
    result["minimum"] = vectorValue(entry.minimum);
    result["maximum"] = vectorValue(entry.maximum);
    result["maximum_absolute"] = vectorValue(entry.maximumAbsolute);
    result["final_state"] = vectorValue(entry.finalState);
    result["accepted_steps"] = static_cast<double>(entry.statistics.acceptedSteps);
    result["rejected_steps"] = static_cast<double>(entry.statistics.rejectedSteps);
    result["rhs_calls"] = static_cast<double>(entry.statistics.rightHandSideCalls);
    return json::Value(std::move(result));
}

struct Prepared {
    dynsyslab::core::StateDomainExpansionPlan plan;
    dynsyslab::core::StateDomainScreeningPolicy screening;
};

Prepared prepare(const json::Value& job, std::size_t dimension,
                 const std::string& coordinates) {
    const auto& source = job.at("preflight");
    if (source.at("kind").string() != "krom_training_domain")
        throw std::invalid_argument("unknown native solve preflight kind");
    if (job.at("initial_conditions").at("mode").string() != "explicit")
        throw std::invalid_argument("KROM training-domain preflight requires explicit initial conditions");
    const std::string requested = source.at("coordinate_system").string(coordinates);
    if (!requested.empty() && requested != coordinates)
        throw std::invalid_argument("training preflight coordinates differ from system coordinates");
    dynsyslab::core::StateDomainExpansionOptions options;
    options.pointsPerAxis = sizes(source.at("points_per_axis"));
    options.marginFraction = source.at("margin_fraction").number(0.1);
    options.outsideCurrentOnly = source.at("outside_current_only").boolean(true);
    options.anchors = flattenedRows(source.at("anchors"), dimension, options.anchorCount);
    const double maximum = source.at("maximum_initial_conditions").number(100000.0);
    if (maximum < 1.0 || maximum != std::floor(maximum))
        throw std::invalid_argument("maximum_initial_conditions must be a positive integer");
    options.maximumInitialConditions = static_cast<std::size_t>(maximum);
    Prepared result;
    result.plan = dynsyslab::core::planStateDomainExpansion(
        domain(source.at("current_domain"), requested),
        domain(source.at("required_domain"), requested), options);
    const auto& screen = source.at("screening");
    result.screening.integration.initialTime = 0.0;
    result.screening.integration.finalTime = screen.at("time_end").number(0.25);
    result.screening.integration.outputStep = screen.at("output_step").number(
        screen.at("maximum_step").number(0.002));
    result.screening.integration.maximumStep = screen.at("maximum_step").number(0.002);
    result.screening.integration.relativeTolerance = screen.at("rtol").number(1.0e-9);
    result.screening.integration.absoluteTolerance = screen.at("atol").number(1.0e-12);
    result.screening.integration.method = method(screen.at("method").string("RK45"));
    result.screening.stateMinimum = numbers(screen.at("state_minimum"));
    result.screening.stateMaximum = numbers(screen.at("state_maximum"));
    result.screening.maximumAbsoluteState = numbers(screen.at("maximum_absolute_state"));
    result.screening.requireAll = screen.at("require_all").boolean(true);
    if (!result.screening.requireAll)
        throw std::invalid_argument("native publication preflight currently requires all initial conditions");
    if (!(result.screening.integration.finalTime > 0.0)
        || !(result.screening.integration.outputStep > 0.0)
        || result.screening.integration.outputStep > result.screening.integration.finalTime)
        throw std::invalid_argument("training screening requires a valid positive time horizon and output step");
    return result;
}

void validateInitials(const dynsyslab::core::StateDomainExpansionPlan& plan,
                      const std::vector<double>& initialStates,
                      std::size_t count) {
    if (plan.initialConditionCount != count
        || plan.initialConditions.size() != initialStates.size())
        throw std::invalid_argument("solve initial conditions no longer match the training-domain plan");
    for (std::size_t index = 0; index < initialStates.size(); ++index) {
        const double scale = std::max({1.0, std::abs(initialStates[index]),
                                      std::abs(plan.initialConditions[index])});
        if (std::abs(initialStates[index] - plan.initialConditions[index]) > 1.0e-12 * scale)
            throw std::invalid_argument("solve initial conditions no longer match the training-domain plan");
    }
}

}  // namespace

void validateTrainingPreflightPlan(const json::Value& job,
                                   std::size_t stateDimension,
                                   const std::string& coordinateSystem) {
    const auto& source = job.at("preflight");
    if (!source.at("enabled").boolean(false)) return;
    const auto prepared = prepare(job, stateDimension, coordinateSystem);
    const double count = job.at("initial_conditions").at("count").number();
    if (count != static_cast<double>(prepared.plan.initialConditionCount))
        throw std::invalid_argument("planned training-domain count differs from solve count");
}

TrainingPreflightResult executeTrainingPreflight(
    const json::Value& job, dynsyslab::core::ODESystem& system,
    const std::vector<double>& initialStates, std::size_t initialConditionCount,
    TrainingPreflightProgress progress,
    const dynsyslab::core::ExecutionContext* execution) {
    TrainingPreflightResult result;
    const auto& source = job.at("preflight");
    result.enabled = source.at("enabled").boolean(false);
    if (!result.enabled) { result.report = json::Value(json::Value::Object{}); return result; }
    const auto prepared = prepare(job, system.descriptor().stateDimension,
                                  system.descriptor().preferredCoordinates);
    validateInitials(prepared.plan, initialStates, initialConditionCount);
    const auto report = dynsyslab::core::screenInitialConditions(
        system, initialStates, initialConditionCount, prepared.screening,
        [&](std::size_t completed, std::size_t planned) {
            if (progress) progress(static_cast<double>(completed) / planned,
                completed, planned, "Screened initial condition "
                    + std::to_string(completed) + " of " + std::to_string(planned));
        }, execution);
    json::Value::Object encoded;
    encoded["schema"] = "dynsyslab.training-screen/v1";
    encoded["kind"] = "krom_training_domain";
    encoded["source_krom"] = source.at("source_krom").string();
    encoded["accepted_count"] = static_cast<double>(report.acceptedIndices.size());
    encoded["planned_count"] = static_cast<double>(initialConditionCount);
    encoded["all_accepted"] = report.allAccepted;
    json::Value::Array entries;
    for (std::size_t index = 0; index < report.entries.size(); ++index)
        entries.push_back(entryValue(report.entries[index], index));
    encoded["entries"] = json::Value(std::move(entries));
    result.report = json::Value(std::move(encoded));
    if (!report.allAccepted) {
        const auto rejected = std::find_if(report.entries.begin(), report.entries.end(),
            [](const auto& entry) { return !entry.accepted; });
        const std::size_t index = static_cast<std::size_t>(
            std::distance(report.entries.begin(), rejected));
        throw std::runtime_error("training-domain preflight rejected initial condition "
            + std::to_string(index) + ": " + rejected->reason);
    }
    return result;
}

}  // namespace dynsyslab::workflow::detail
