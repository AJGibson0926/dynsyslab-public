#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/PDEParallel.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void decayFieldParity() {
    dynsyslab::core::StructuredMesh mesh{{2, 3}, {1.0, 1.0}, {0.0, 0.0}};
    dynsyslab::core::PDEIntegrationOptions options;
    options.finalTime = 1.0;
    options.outputStep = 0.01;
    options.internalStep = 0.01;
    options.method = dynsyslab::core::PDEIntegrationMethod::rk4;
    const std::vector<double> initial{1, 2, 3, 4, 5, 6};
    const auto result = dynsyslab::core::integratePDE(
        mesh, {2, 3}, initial, options,
        [](double, const auto&, const double* field, std::size_t count,
           double* derivative) {
            for (std::size_t index = 0; index < count; ++index)
                derivative[index] = -field[index];
        });
    require(result.sampleCount == 101 && result.fieldShape == std::vector<std::size_t>({2, 3}),
            "PDE decay-field shape parity");
    for (std::size_t index = 0; index < initial.size(); ++index)
        require(std::abs(result.values[(result.sampleCount - 1) * initial.size() + index]
            - initial[index] * std::exp(-1.0)) < 2.0e-9,
            "PDE RK4 decay-field accuracy");
}

void heatEquationMode() {
    dynsyslab::core::SystemConfiguration request;
    request.identifier = "heat_equation_1d";
    request.preset = "periodic_decay";
    request.parameters = {{"points", 32}, {"diffusivity", 0.01},
                          {"decay_rate", 0.1}, {"wave_number", 2}};
    const auto configuration = dynsyslab::core::resolveHeatEquation(request);
    const auto initial = dynsyslab::core::heatEquationInitialField(configuration);
    dynsyslab::core::PDEIntegrationOptions options;
    options.finalTime = 0.1;
    options.outputStep = 0.01;
    options.internalStep = 0.001;
    const auto mesh = dynsyslab::core::heatEquationMesh(configuration);
    const auto result = dynsyslab::core::integratePDE(
        mesh, {initial.size()}, initial, options,
        dynsyslab::core::heatEquationRightHandSide(configuration));
    const double discreteRate = 0.1 + 4.0 * 0.01
        * std::pow(std::sin(2.0 * 3.14159265358979323846 / 32.0), 2)
        / std::pow(mesh.spacing[0], 2);
    const double scale = std::exp(-discreteRate * options.finalTime);
    for (std::size_t index = 0; index < initial.size(); ++index)
        require(std::abs(result.values[(result.sampleCount - 1) * initial.size() + index]
            - initial[index] * scale) < 2.0e-10,
            "heat-equation discrete Fourier mode evolution");
    require(dynsyslab::core::heatEquationStableStep(configuration) > options.internalStep,
            "heat-equation stability guidance");
}

void partitionedHeatParity() {
    dynsyslab::core::SystemConfiguration request;
    request.identifier = "heat_equation_1d";
    request.preset = "periodic_decay";
    request.parameters = {{"points", 31}, {"diffusivity", 0.01},
                          {"decay_rate", 0.05}, {"wave_number", 3}};
    const auto configuration = dynsyslab::core::resolveHeatEquation(request);
    const auto mesh = dynsyslab::core::heatEquationMesh(configuration);
    const auto initial = dynsyslab::core::heatEquationInitialField(configuration);
    dynsyslab::core::PDEIntegrationOptions options;
    options.finalTime = 0.04;
    options.outputStep = 0.01;
    options.internalStep = 0.001;
    const auto serial = dynsyslab::core::integratePDE(
        mesh, {initial.size()}, initial, options,
        dynsyslab::core::heatEquationRightHandSide(configuration));
    const auto plan = dynsyslab::core::buildParallelMeshPlan(mesh, {4}, 1, 3);
    require(plan.partitions.size() == 4
        && plan.partitions.front().interior.front()
            == std::pair<std::size_t, std::size_t>(0, 8)
        && plan.partitions.back().interior.front()
            == std::pair<std::size_t, std::size_t>(24, 31),
        "uneven structured partition cuts");
    require(plan.partitions[1].neighbors.at("axis0-") == 0
        && plan.partitions[1].neighbors.at("axis0+") == 2,
        "structured partition neighbors");
    auto overlapping = plan;
    overlapping.partitions[1].interior[0].first = 3;
    bool overlapRejected = false;
    try {
        (void)dynsyslab::core::partitionedPDERightHandSide(
            overlapping, dynsyslab::core::heatEquationSubdomainRightHandSide(configuration));
    } catch (const std::invalid_argument&) {
        overlapRejected = true;
    }
    require(overlapRejected, "overlapping structured partitions are rejected");
    const auto parallel = dynsyslab::core::integratePDE(
        mesh, {initial.size()}, initial, options,
        dynsyslab::core::partitionedPDERightHandSide(
            plan, dynsyslab::core::heatEquationSubdomainRightHandSide(configuration)));
    require(serial.time == parallel.time && serial.values.size() == parallel.values.size(),
            "partitioned heat output grid parity");
    for (std::size_t index = 0; index < serial.values.size(); ++index)
        require(std::abs(serial.values[index] - parallel.values[index]) < 1.0e-14,
                "partitioned heat numerical parity");
}

void multidimensionalComponentParity() {
    dynsyslab::core::StructuredMesh mesh{
        {4, 5}, {1.0, 1.0}, {0.0, 0.0},
        {dynsyslab::core::PDEBoundary::edge, dynsyslab::core::PDEBoundary::edge}};
    std::vector<double> initial(4 * 5 * 2);
    for (std::size_t index = 0; index < initial.size(); ++index)
        initial[index] = 0.25 + static_cast<double>(index);
    dynsyslab::core::PDEIntegrationOptions options;
    options.finalTime = 0.02;
    options.outputStep = 0.01;
    options.internalStep = 0.002;
    const auto serialRhs = [](double, const auto&, const double* field,
        std::size_t count, double* derivative) {
        for (std::size_t index = 0; index < count; ++index)
            derivative[index] = -field[index];
    };
    const auto serial = dynsyslab::core::integratePDE(
        mesh, {4, 5, 2}, initial, options, serialRhs);
    const auto plan = dynsyslab::core::buildParallelMeshPlan(mesh, {2, 3}, 2, 4);
    require(plan.localWorkers == 4, "parallel worker count records effective workers");
    const auto subdomainRhs = [](double, const auto&,
        const dynsyslab::core::ParallelMeshPlan& activePlan,
        const dynsyslab::core::MeshPartition& partition, const double* ghosted,
        const std::vector<std::size_t>& ghostShape, std::size_t components,
        double* derivative) {
        const auto localShape = partition.shape();
        for (std::size_t row = 0; row < localShape[0]; ++row)
            for (std::size_t column = 0; column < localShape[1]; ++column) {
                const std::size_t center = ((row + activePlan.ghostWidth)
                    * ghostShape[1] + column + activePlan.ghostWidth) * components;
                const std::size_t target = (row * localShape[1] + column) * components;
                for (std::size_t component = 0; component < components; ++component)
                    derivative[target + component] = -ghosted[center + component];
            }
    };
    const auto parallel = dynsyslab::core::integratePDE(
        mesh, {4, 5, 2}, initial, options,
        dynsyslab::core::partitionedPDERightHandSide(plan, subdomainRhs));
    require(serial.values == parallel.values,
            "multidimensional component partition parity");
}

void resourceAndCatalogContracts() {
    dynsyslab::core::StructuredMesh mesh{{4}, {1.0}, {0.0}};
    dynsyslab::core::PDEIntegrationOptions options;
    dynsyslab::core::ExecutionContext execution;
    execution.limits.maximumMemoryBytes = 1;
    bool rejected = false;
    try {
        (void)dynsyslab::core::integratePDE(mesh, {4}, {1, 2, 3, 4}, options,
            [](double, const auto&, const double*, std::size_t count, double* derivative) {
                std::fill(derivative, derivative + count, 0.0);
            }, &execution);
    } catch (const dynsyslab::core::ExecutionLimitExceeded&) {
        rejected = true;
    }
    require(rejected, "PDE memory budget enforcement");
    const auto systems = dynsyslab::core::availableSystems();
    const auto found = std::find_if(systems.begin(), systems.end(), [](const auto& item) {
        return item.identifier == "heat_equation_1d" && item.family == "pde";
    });
    require(found != systems.end(), "PDE system catalog registration");
}

}  // namespace

int main() {
    try {
        decayFieldParity();
        heatEquationMode();
        partitionedHeatParity();
        multidimensionalComponentParity();
        resourceAndCatalogContracts();
        std::cout << "dynsyslab_pde_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_pde_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
