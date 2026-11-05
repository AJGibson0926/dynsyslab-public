#include <dynsyslab/core/System.hpp>
#include <dynsyslab/core/TrainingDomain.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        using namespace dynsyslab::core;
        const AxisAlignedStateDomain current{{-1.0, -0.5}, {1.0, 0.5}, "dimensionless"};
        const AxisAlignedStateDomain required{{-2.0, -0.5}, {1.0, 1.5}, "dimensionless"};
        StateDomainExpansionOptions options;
        options.pointsPerAxis = {4, 3};
        options.anchors = {0.0, 0.0, -2.2, 1.7};
        options.anchorCount = 2;
        const auto plan = planStateDomainExpansion(current, required, options);
        require(plan.planned.minimum[0] < -2.0 && plan.planned.maximum[1] > 1.5,
                "expansion margin");
        require(plan.initialConditionCount > 0 && plan.initialConditionCount <= 14,
                "outside-current grid and anchors");
        bool foundOrigin = false;
        for (std::size_t row = 0; row < plan.initialConditionCount; ++row)
            foundOrigin = foundOrigin
                || (plan.initialConditions[row * 2] == 0.0
                    && plan.initialConditions[row * 2 + 1] == 0.0);
        require(foundOrigin, "anchor retained even inside current domain");
        StateDomainExpansionOptions orderingOptions;
        orderingOptions.pointsPerAxis = {2,3};
        orderingOptions.outsideCurrentOnly = false;
        const auto ordered = planStateDomainExpansion(
            {{0,0},{1,2},"dimensionless"},{{0,0},{1,2},"dimensionless"},
            orderingOptions);
        require(ordered.initialConditions == std::vector<double>({
                    0,0, 0,1, 0,2, 1,0, 1,1, 1,2}),
                "Cartesian ordering matches ij-indexed row-major grids");
        orderingOptions.outsideCurrentOnly = true;
        require(planStateDomainExpansion(
                    {{0,0},{1,2},"dimensionless"},{{0,0},{1,2},"dimensionless"},
                    orderingOptions).initialConditionCount == 0,
                "a fully covered requested domain needs no expansion points");

        const AxisAlignedStateDomain lorenzCurrent{{-1,-1,-1},{1,1,1},"dimensionless"};
        const AxisAlignedStateDomain lorenzRequired{{-2,-2,-2},{2,2,2},"dimensionless"};
        StateDomainExpansionOptions lorenzOptions;
        lorenzOptions.pointsPerAxis = {3,3,3};
        lorenzOptions.outsideCurrentOnly = false;
        require(planStateDomainExpansion(lorenzCurrent, lorenzRequired, lorenzOptions)
                    .initialConditionCount == 27,
                "system-agnostic three-dimensional planner");

        SystemConfiguration configuration;
        configuration.identifier = "duffing";
        configuration.preset = "double_well_reference";
        auto system = createSystem(configuration);
        StateDomainScreeningPolicy screening;
        screening.integration.initialTime = 0.0;
        screening.integration.finalTime = 0.02;
        screening.integration.outputStep = 0.01;
        screening.integration.maximumStep = 0.002;
        screening.maximumAbsoluteState = {5.0, 5.0};
        std::size_t updates = 0;
        const auto accepted = screenInitialConditions(
            *system, {-1.0, 0.0, 1.0, 0.0}, 2, screening,
            [&](std::size_t, std::size_t) { ++updates; });
        require(accepted.allAccepted && accepted.acceptedIndices.size() == 2,
                "stable Duffing initial conditions accepted");
        require(updates == 2, "screening progress is trajectory-granular");
        screening.maximumAbsoluteState = {0.5, 5.0};
        const auto rejected = screenInitialConditions(*system, {-1.0, 0.0}, 1, screening);
        require(!rejected.allAccepted && !rejected.entries[0].reason.empty(),
                "component limit rejects unsafe initial condition");
        std::cout << "dynsyslab_training_domain_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_training_domain_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
