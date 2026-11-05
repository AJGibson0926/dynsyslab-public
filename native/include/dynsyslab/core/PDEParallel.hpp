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

    [[nodiscard]] std::vector<std::size_t> shape() const;
};

struct ParallelMeshPlan {
    StructuredMesh mesh;
    std::vector<std::size_t> partitionShape;
    std::size_t ghostWidth{1};
    std::size_t localWorkers{1};
    std::vector<MeshPartition> partitions;

    [[nodiscard]] std::size_t workerCount() const;
    void validate() const;
};

[[nodiscard]] ParallelMeshPlan buildParallelMeshPlan(
    StructuredMesh mesh, std::vector<std::size_t> partitionShape,
    std::size_t ghostWidth = 1, std::size_t localWorkers = 1);

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
