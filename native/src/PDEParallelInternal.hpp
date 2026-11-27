#pragma once

// Internal helpers shared by the partitioned PDE right-hand sides (PDEParallel.cpp) and the
// adaptive load balancer (PDEAdaptivePartition.cpp).

#include <dynsyslab/core/PDEParallel.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core::detail {

inline constexpr const char* mpiInProcessMessage =
    "communication \"mpi\" runs one partition per MPI rank: launch the ranks with mpiexec "
    "and exchange halos with mpiExchangeRankField (DYNSYSLAB_WITH_MPI build)";

// One evaluation of the partitioned right-hand side; `seconds` (optional) receives each
// partition's wall time, the measured load used by the adaptive rebalance.
void evaluatePartitioned(const ParallelMeshPlan& plan, const PDESubdomainRightHandSide& rhs,
                         double time, const StructuredMesh& mesh, const double* field,
                         std::size_t count, double* derivative, std::vector<double>* seconds);

// The partition's interior cells of a global row-major field (cells x components).
[[nodiscard]] std::vector<double> gatherInterior(const ParallelMeshPlan& plan,
                                                 const MeshPartition& partition,
                                                 const double* field, std::size_t components);

}  // namespace dynsyslab::core::detail
