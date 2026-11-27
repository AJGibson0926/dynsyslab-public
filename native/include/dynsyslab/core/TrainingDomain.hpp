#pragma once

#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/Integrator.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

// Membership of one state in an axis-aligned training domain: the axes it leaves (beyond a
// 64-ulp tolerance) and the excursion per axis relative to the domain width.
struct DomainReport {
    bool inside{true};
    std::vector<std::size_t> violatingAxes;
    std::vector<double> relativeExcursion;
    double maximumRelativeExcursion{};
};

struct AxisAlignedStateDomain {
    std::vector<double> minimum;
    std::vector<double> maximum;
    std::string coordinateSystem;

    [[nodiscard]] std::size_t dimension() const noexcept { return minimum.size(); }
    [[nodiscard]] bool contains(const double* state) const;
    [[nodiscard]] DomainReport report(const double* state) const;
};

struct StateDomainExpansionOptions {
    std::vector<std::size_t> pointsPerAxis;
    double marginFraction{0.1};
    bool outsideCurrentOnly{true};
    std::vector<double> anchors;
    std::size_t anchorCount{};
    std::size_t maximumInitialConditions{100000};
};

struct StateDomainExpansionPlan {
    AxisAlignedStateDomain current;
    AxisAlignedStateDomain required;
    AxisAlignedStateDomain planned;
    std::vector<std::size_t> pointsPerAxis;
    std::vector<double> initialConditions;
    std::size_t initialConditionCount{};
    bool outsideCurrentOnly{};
};

[[nodiscard]] StateDomainExpansionPlan planStateDomainExpansion(
    const AxisAlignedStateDomain& current,
    const AxisAlignedStateDomain& required,
    const StateDomainExpansionOptions& options);

struct StateDomainScreeningPolicy {
    IntegrationOptions integration;
    std::vector<double> stateMinimum;
    std::vector<double> stateMaximum;
    std::vector<double> maximumAbsoluteState;
    bool requireAll{true};
};

struct StateDomainScreeningEntry {
    bool accepted{};
    std::string reason;
    std::size_t samples{};
    std::vector<double> minimum;
    std::vector<double> maximum;
    std::vector<double> maximumAbsolute;
    std::vector<double> finalState;
    IntegrationStatistics statistics;
};

struct StateDomainScreeningReport {
    std::vector<StateDomainScreeningEntry> entries;
    std::vector<std::size_t> acceptedIndices;
    bool allAccepted{};
};

using StateDomainScreeningProgress =
    std::function<void(std::size_t completed, std::size_t planned)>;

[[nodiscard]] StateDomainScreeningReport screenInitialConditions(
    ODESystem& system, const std::vector<double>& initialStates,
    std::size_t initialConditionCount,
    const StateDomainScreeningPolicy& policy,
    StateDomainScreeningProgress progress = {},
    const ExecutionContext* execution = nullptr);

}  // namespace dynsyslab::core
