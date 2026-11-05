#include <dynsyslab/core/TrainingDomain.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

void validateDomain(const AxisAlignedStateDomain& domain, const char* name) {
    if (domain.minimum.empty() || domain.minimum.size() != domain.maximum.size())
        throw std::invalid_argument(std::string(name) + " domain has inconsistent dimensions");
    for (std::size_t axis = 0; axis < domain.minimum.size(); ++axis) {
        if (!std::isfinite(domain.minimum[axis]) || !std::isfinite(domain.maximum[axis])
            || domain.minimum[axis] > domain.maximum[axis])
            throw std::invalid_argument(std::string(name) + " domain has invalid bounds");
    }
}

bool samePoint(const double* left, const double* right, std::size_t dimension) {
    for (std::size_t axis = 0; axis < dimension; ++axis)
        if (left[axis] != right[axis]) return false;
    return true;
}

void validateOptionalBounds(const std::vector<double>& values,
                            std::size_t dimension, const char* name) {
    if (!values.empty() && values.size() != dimension)
        throw std::invalid_argument(std::string(name) + " must cover every state component");
    for (const double value : values)
        if (!std::isfinite(value))
            throw std::invalid_argument(std::string(name) + " must contain finite values");
}

std::string boundReason(const StateDomainScreeningPolicy& policy,
                        const double* state, std::size_t dimension) {
    for (std::size_t axis = 0; axis < dimension; ++axis) {
        if (!std::isfinite(state[axis])) return "non-finite state";
        if (!policy.stateMinimum.empty() && state[axis] < policy.stateMinimum[axis])
            return "state below screening minimum on component " + std::to_string(axis);
        if (!policy.stateMaximum.empty() && state[axis] > policy.stateMaximum[axis])
            return "state above screening maximum on component " + std::to_string(axis);
        if (!policy.maximumAbsoluteState.empty()
            && std::abs(state[axis]) > policy.maximumAbsoluteState[axis])
            return "absolute state exceeded screening limit on component "
                + std::to_string(axis);
    }
    return {};
}

}  // namespace

bool AxisAlignedStateDomain::contains(const double* state) const {
    for (std::size_t axis = 0; axis < dimension(); ++axis)
        if (state[axis] < minimum[axis] || state[axis] > maximum[axis]) return false;
    return true;
}

StateDomainExpansionPlan planStateDomainExpansion(
    const AxisAlignedStateDomain& current,
    const AxisAlignedStateDomain& required,
    const StateDomainExpansionOptions& options) {
    validateDomain(current, "current");
    validateDomain(required, "required");
    if (current.dimension() != required.dimension())
        throw std::invalid_argument("current and required domains have different dimensions");
    if (!current.coordinateSystem.empty() && !required.coordinateSystem.empty()
        && current.coordinateSystem != required.coordinateSystem)
        throw std::invalid_argument("current and required domains use different coordinates");
    const std::size_t dimension = current.dimension();
    if (options.pointsPerAxis.size() != dimension)
        throw std::invalid_argument("points_per_axis must cover every state component");
    if (!std::isfinite(options.marginFraction) || options.marginFraction < 0.0)
        throw std::invalid_argument("margin_fraction must be finite and non-negative");
    if (options.anchorCount * dimension != options.anchors.size())
        throw std::invalid_argument("anchor dimensions are inconsistent");
    if (!options.maximumInitialConditions)
        throw std::invalid_argument("maximum_initial_conditions must be positive");

    StateDomainExpansionPlan result;
    result.current = current;
    result.required = required;
    result.planned.coordinateSystem = current.coordinateSystem.empty()
        ? required.coordinateSystem : current.coordinateSystem;
    result.pointsPerAxis = options.pointsPerAxis;
    result.outsideCurrentOnly = options.outsideCurrentOnly;
    result.planned.minimum.resize(dimension);
    result.planned.maximum.resize(dimension);
    std::size_t gridCount = 1;
    for (std::size_t axis = 0; axis < dimension; ++axis) {
        if (options.pointsPerAxis[axis] < 2)
            throw std::invalid_argument("each points_per_axis value must be at least two");
        if (gridCount > options.maximumInitialConditions / options.pointsPerAxis[axis])
            throw std::invalid_argument("training-domain grid exceeds its initial-condition limit");
        gridCount *= options.pointsPerAxis[axis];
        double low = std::min(current.minimum[axis], required.minimum[axis]);
        double high = std::max(current.maximum[axis], required.maximum[axis]);
        const double margin = options.marginFraction
            * std::max(high - low, std::numeric_limits<double>::epsilon());
        if (required.minimum[axis] < current.minimum[axis]) low -= margin;
        if (required.maximum[axis] > current.maximum[axis]) high += margin;
        result.planned.minimum[axis] = low;
        result.planned.maximum[axis] = high;
    }
    result.initialConditions.reserve((gridCount + options.anchorCount) * dimension);
    for (std::size_t point = 0; point < gridCount; ++point) {
        std::vector<double> state(dimension);
        std::size_t stride = gridCount;
        for (std::size_t axis = 0; axis < dimension; ++axis) {
            stride /= options.pointsPerAxis[axis];
            const std::size_t index = (point / stride) % options.pointsPerAxis[axis];
            const double fraction = static_cast<double>(index)
                / static_cast<double>(options.pointsPerAxis[axis] - 1);
            state[axis] = result.planned.minimum[axis]
                + fraction * (result.planned.maximum[axis] - result.planned.minimum[axis]);
        }
        if (options.outsideCurrentOnly && current.contains(state.data())) continue;
        result.initialConditions.insert(
            result.initialConditions.end(), state.begin(), state.end());
    }
    for (std::size_t anchor = 0; anchor < options.anchorCount; ++anchor) {
        const double* state = options.anchors.data() + anchor * dimension;
        if (!result.planned.contains(state))
            throw std::invalid_argument("training-domain anchor lies outside the planned domain");
        bool duplicate = false;
        for (std::size_t row = 0; row < result.initialConditions.size() / dimension; ++row)
            if (samePoint(result.initialConditions.data() + row * dimension,
                          state, dimension)) { duplicate = true; break; }
        if (!duplicate)
            result.initialConditions.insert(result.initialConditions.end(), state, state + dimension);
    }
    result.initialConditionCount = result.initialConditions.size() / dimension;
    if (result.initialConditionCount > options.maximumInitialConditions)
        throw std::invalid_argument("training-domain plan exceeds its initial-condition limit");
    return result;
}

StateDomainScreeningReport screenInitialConditions(
    ODESystem& system, const std::vector<double>& initialStates,
    std::size_t initialConditionCount, const StateDomainScreeningPolicy& policy,
    StateDomainScreeningProgress progress, const ExecutionContext* execution) {
    const std::size_t dimension = system.descriptor().stateDimension;
    if (!initialConditionCount || initialStates.size() != initialConditionCount * dimension)
        throw std::invalid_argument("screening initial-condition dimensions are inconsistent");
    validateOptionalBounds(policy.stateMinimum, dimension, "state_minimum");
    validateOptionalBounds(policy.stateMaximum, dimension, "state_maximum");
    validateOptionalBounds(policy.maximumAbsoluteState, dimension, "maximum_absolute_state");
    for (std::size_t axis = 0; axis < dimension; ++axis) {
        if (!policy.stateMinimum.empty() && !policy.stateMaximum.empty()
            && policy.stateMinimum[axis] > policy.stateMaximum[axis])
            throw std::invalid_argument("screening state minimum exceeds maximum");
        if (!policy.maximumAbsoluteState.empty() && policy.maximumAbsoluteState[axis] < 0.0)
            throw std::invalid_argument("maximum_absolute_state must be non-negative");
    }
    StateDomainScreeningReport report;
    report.entries.reserve(initialConditionCount);
    ExecutionContext child;
    if (execution) { child.limits = execution->limits; child.cancellation = execution->cancellation; }
    for (std::size_t index = 0; index < initialConditionCount; ++index) {
        if (execution) execution->checkpoint();
        StateDomainScreeningEntry entry;
        entry.minimum.assign(dimension, std::numeric_limits<double>::infinity());
        entry.maximum.assign(dimension, -std::numeric_limits<double>::infinity());
        entry.maximumAbsolute.assign(dimension, 0.0);
        const auto begin = initialStates.begin() + static_cast<std::ptrdiff_t>(index * dimension);
        const std::vector<double> initial(begin, begin + static_cast<std::ptrdiff_t>(dimension));
        entry.reason = boundReason(policy, initial.data(), dimension);
        try {
            if (entry.reason.empty()) {
                const auto solution = integrate(system, initial, 1, policy.integration,
                    {}, {}, execution ? &child : nullptr);
                entry.samples = solution.sampleCount;
                entry.statistics = solution.statistics;
                for (std::size_t row = 0; row < solution.sampleCount; ++row) {
                    const double* state = solution.state(0, row);
                    if (entry.reason.empty()) entry.reason = boundReason(policy, state, dimension);
                    for (std::size_t axis = 0; axis < dimension; ++axis) {
                        entry.minimum[axis] = std::min(entry.minimum[axis], state[axis]);
                        entry.maximum[axis] = std::max(entry.maximum[axis], state[axis]);
                        entry.maximumAbsolute[axis] = std::max(
                            entry.maximumAbsolute[axis], std::abs(state[axis]));
                    }
                }
                entry.finalState.assign(solution.state(0, solution.sampleCount - 1),
                                        solution.state(0, solution.sampleCount - 1) + dimension);
            }
        } catch (const ExecutionCancelled&) {
            throw;
        } catch (const std::exception& error) {
            entry.reason = error.what();
        }
        entry.accepted = entry.reason.empty();
        if (entry.accepted) report.acceptedIndices.push_back(index);
        report.entries.push_back(std::move(entry));
        if (progress) progress(index + 1, initialConditionCount);
    }
    report.allAccepted = report.acceptedIndices.size() == initialConditionCount;
    return report;
}

}  // namespace dynsyslab::core
