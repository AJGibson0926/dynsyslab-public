#pragma once

#include <dynsyslab/core/PDE.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::core {

struct MeshPartition {
    std::size_t rank{};
    std::vector<std::pair<std::size_t, std::size_t>> interior;
    std::map<std::string, std::size_t> neighbors;
    double weight{1.0};   // relative load estimate (adaptive rebalance, PDEAdaptivePartition.hpp)

    [[nodiscard]] std::vector<std::size_t> shape() const;
};

// How subdomains obtain their ghost (halo) cells.
//   sharedMemory: every partition reads its ghosted window from the global field.
//   haloExchange: every partition owns only its interior; face halos move rank to rank as
//     packed buffers (PDEHaloExchange.hpp), the in-process equivalent of an MPI run.
//   mpi: halo buffers move through MPI (optional DYNSYSLAB_WITH_MPI build; one rank per
//     partition, launched by mpiexec).
enum class PDECommunication { sharedMemory, haloExchange, mpi };

// Load-balance policy of a partitioned solve (port of the Python AdaptiveMeshPolicy).
// Rebalance when max(weight) / mean(weight) exceeds rebalanceImbalance; the refinement and
// coarsening thresholds are recorded for a future refining mesh adapter.
struct AdaptiveMeshPolicy {
    bool enabled{false};
    double refinementThreshold{0.8};
    double coarseningThreshold{0.1};
    double rebalanceImbalance{1.25};
    bool checkpointBeforeRepartition{true};
    std::size_t rebalanceIntervalCalls{50};   // right-hand-side calls between load checks
};

struct ParallelMeshPlan {
    StructuredMesh mesh;
    std::vector<std::size_t> partitionShape;
    std::size_t ghostWidth{1};
    std::size_t localWorkers{1};
    std::vector<MeshPartition> partitions;
    PDECommunication communication{PDECommunication::sharedMemory};
    AdaptiveMeshPolicy adaptive;

    [[nodiscard]] std::size_t workerCount() const;
    void validate() const;
};

[[nodiscard]] ParallelMeshPlan buildParallelMeshPlan(
    StructuredMesh mesh, std::vector<std::size_t> partitionShape,
    std::size_t ghostWidth = 1, std::size_t localWorkers = 1);

// A plan from explicit per-axis cuts (half-open, contiguous, covering each axis).
[[nodiscard]] ParallelMeshPlan buildParallelMeshPlanFromCuts(
    StructuredMesh mesh, std::vector<std::vector<std::pair<std::size_t, std::size_t>>> axisCuts,
    std::size_t ghostWidth = 1, std::size_t localWorkers = 1);

[[nodiscard]] const char* pdeCommunicationName(PDECommunication communication);
[[nodiscard]] PDECommunication pdeCommunicationFromName(const std::string& name);

using PDESubdomainRightHandSide = std::function<void(
    double time, const StructuredMesh& mesh, const ParallelMeshPlan& plan,
    const MeshPartition& partition, const double* ghostedField,
    const std::vector<std::size_t>& ghostedShape, std::size_t componentCount,
    double* interiorDerivative)>;

// localWorkers > 1 may invoke rhs concurrently for distinct partitions. The
// callback must therefore be reentrant; every invocation owns its input/output
// buffers and no callback writes another partition's state.
[[nodiscard]] PDERightHandSide partitionedPDERightHandSide(
    ParallelMeshPlan plan, PDESubdomainRightHandSide rhs);

[[nodiscard]] PDESubdomainRightHandSide heatEquationSubdomainRightHandSide(
    const SystemConfiguration& configuration);

}  // namespace dynsyslab::core
