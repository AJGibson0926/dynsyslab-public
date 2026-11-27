// Halo exchange, MPI availability and adaptive rebalance of partitioned PDE solves.
#include <dynsyslab/core/PDE.hpp>
#include <dynsyslab/core/PDEAdaptivePartition.hpp>
#include <dynsyslab/core/PDEHaloExchange.hpp>
#include <dynsyslab/core/PDEParallel.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace core = dynsyslab::core;

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

// Sum over the full (2w+1)^2 window, corners included, times a per-component weight.
void boxStencil(double, const core::StructuredMesh&, const core::ParallelMeshPlan& plan,
                const core::MeshPartition& partition, const double* ghosted,
                const std::vector<std::size_t>& ghostShape, std::size_t components, double* derivative) {
    const auto shape = partition.shape();
    const auto w = static_cast<std::ptrdiff_t>(plan.ghostWidth);
    for (std::size_t row = 0; row < shape[0]; ++row)
        for (std::size_t column = 0; column < shape[1]; ++column)
            for (std::size_t c = 0; c < components; ++c) {
                double sum = 0.0;
                for (std::ptrdiff_t dr = -w; dr <= w; ++dr)
                    for (std::ptrdiff_t dc = -w; dc <= w; ++dc) {
                        const auto r = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(row) + w + dr);
                        const auto k = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(column) + w + dc);
                        sum += (1.0 + 0.1 * static_cast<double>(dr * 3 + dc)) * ghosted[(r * ghostShape[1] + k) * components + c];
                    }
                derivative[(row * shape[1] + column) * components + c] = (c + 1.0) * sum;
            }
}

void haloExchangeMatchesSharedMemory() {
    using B = core::PDEBoundary;
    for (const auto& boundary : std::vector<std::vector<B>>{{B::edge, B::periodic}, {B::periodic, B::zeroValue},
                                                            {B::zeroValue, B::edge}, {B::periodic, B::periodic}})
        for (const auto& [parts, width] : std::vector<std::pair<std::vector<std::size_t>, std::size_t>>{
                 {{2, 3}, 1}, {{3, 2}, 2}, {{1, 3}, 2}, {{2, 1}, 1}}) {
            core::StructuredMesh mesh{{6, 7}, {1.0, 1.0}, {0.0, 0.0}, boundary};
            std::vector<double> field(6 * 7 * 2);
            for (std::size_t index = 0; index < field.size(); ++index)
                field[index] = std::sin(0.37 * static_cast<double>(index)) + 0.01 * static_cast<double>(index);
            auto shared = core::buildParallelMeshPlan(mesh, parts, width, 3);
            auto halo = shared;
            halo.communication = core::PDECommunication::haloExchange;
            std::vector<double> a(field.size()), b(field.size());
            core::partitionedPDERightHandSide(shared, boxStencil)(0.0, mesh, field.data(), field.size(), a.data());
            core::partitionedPDERightHandSide(halo, boxStencil)(0.0, mesh, field.data(), field.size(), b.data());
            for (std::size_t index = 0; index < a.size(); ++index)
                require(std::abs(a[index] - b[index]) < 1.0e-13, "halo exchange equals the shared-memory ghost read");
        }
}

void mpiAvailability() {
    core::StructuredMesh mesh{{8}, {1.0}, {0.0}};
    auto plan = core::buildParallelMeshPlan(mesh, {2});
    plan.communication = core::PDECommunication::mpi;
    bool rejected = false;
    try { (void)core::partitionedPDERightHandSide(plan, boxStencil); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "communication mpi is refused in process");
    if (!core::mpiHaloExchangeAvailable()) {
        bool thrown = false;
        const std::vector<double> interior(4, 1.0);
        try { (void)core::mpiExchangeRankField(plan, interior.data(), 1); }
        catch (const std::runtime_error&) { thrown = true; }
        require(thrown, "mpiExchangeRankField reports a build without MPI");
    }
    require(core::pdeCommunicationFromName("halo_exchange") == core::PDECommunication::haloExchange
            && std::string(core::pdeCommunicationName(core::PDECommunication::mpi)) == "mpi",
            "communication names round-trip");
}

void rebalanceMovesCutsTowardTheHeavyPartition() {
    core::StructuredMesh mesh{{40}, {1.0}, {0.0}};
    auto plan = core::buildParallelMeshPlan(mesh, {4});
    plan.adaptive.enabled = true;
    const std::vector<double> heavy{4.0, 1.0, 1.0, 1.0};
    require(std::abs(core::loadImbalance(plan, heavy) - 4.0 / 1.75) < 1.0e-12, "imbalance is max over mean");
    require(core::needsRebalance(plan, heavy), "an imbalanced load needs a rebalance");
    const auto [next, record] = core::rebalanceParallelMeshPlan(plan, heavy);
    require(record.repartitioned && record.predictedImbalance < record.imbalance, "the rebalance improves the balance");
    require(next.partitions[0].shape()[0] < 10 && next.partitions[3].shape()[0] > 10, "the heavy partition shrinks");
    require(record.checkpointRequired && record.json().find("dynsyslab.pde-repartition/v1") != std::string::npos,
            "the record carries the policy and schema");
    require(!core::needsRebalance(plan, {1.0, 1.0, 1.0, 1.0}), "a balanced load needs no rebalance");
    const auto [same, unchanged] = core::rebalanceParallelMeshPlan(plan, {1.0, 1.0, 1.0, 1.0});
    require(!unchanged.repartitioned && unchanged.cuts == unchanged.priorCuts
            && same.partitions[2].interior == plan.partitions[2].interior, "a balanced plan keeps its cuts");
    bool rejected = false;
    try { (void)core::loadImbalance(plan, {1.0, 0.0, 1.0, 1.0}); } catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "weights must be positive");
}

void adaptiveSolveKeepsTheNumbers() {
    core::SystemConfiguration request;
    request.identifier = "heat_equation_1d";
    request.preset = "periodic_decay";
    request.parameters = {{"points", 33}, {"diffusivity", 0.01}, {"decay_rate", 0.05}, {"wave_number", 2}};
    const auto configuration = core::resolveHeatEquation(request);
    const auto mesh = core::heatEquationMesh(configuration);
    const auto initial = core::heatEquationInitialField(configuration);
    core::PDEIntegrationOptions options;
    options.finalTime = 0.05;
    options.outputStep = 0.01;
    options.internalStep = 0.001;
    const auto serial = core::integratePDE(mesh, {initial.size()}, initial, options,
                                           core::heatEquationRightHandSide(configuration));
    auto plan = core::buildParallelMeshPlan(mesh, {4}, 1, 2);
    plan.adaptive.enabled = true;
    plan.adaptive.rebalanceImbalance = 1.0;
    plan.adaptive.rebalanceIntervalCalls = 7;
    plan.communication = core::PDECommunication::haloExchange;
    auto log = std::make_shared<core::RepartitionLog>();
    const auto adaptive = core::integratePDE(mesh, {initial.size()}, initial, options,
        core::adaptivePartitionedPDERightHandSide(plan, core::heatEquationSubdomainRightHandSide(configuration), log));
    for (std::size_t index = 0; index < serial.values.size(); ++index)
        require(std::abs(serial.values[index] - adaptive.values[index]) < 1.0e-14, "repartitioning keeps the solution");
    for (const auto& record : log->snapshot())
        require(record.cuts.size() == 1 && record.cuts[0].size() == 4 && record.cuts[0].back().second == 33,
                "every record holds complete cuts");
}

}  // namespace

int main() {
    try {
        haloExchangeMatchesSharedMemory();
        mpiAvailability();
        rebalanceMovesCutsTowardTheHeavyPartition();
        adaptiveSolveKeepsTheNumbers();
        std::cout << "dynsyslab_pde_parallel_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_pde_parallel_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
