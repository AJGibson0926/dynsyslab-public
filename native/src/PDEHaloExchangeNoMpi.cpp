// Build without DYNSYSLAB_WITH_MPI: the in-process halo exchange (PDEHaloExchange.cpp) is the
// parallel path; the MPI transport reports that it is not built.
#include <dynsyslab/core/PDEHaloExchange.hpp>

#include <stdexcept>

namespace dynsyslab::core {

bool mpiHaloExchangeAvailable() { return false; }

std::vector<double> mpiExchangeRankField(const ParallelMeshPlan&, const double*, std::size_t) {
    throw std::runtime_error("this DynSysLab build has no MPI: configure with -DDYNSYSLAB_WITH_MPI=ON, "
                             "or use communication \"halo_exchange\" for the in-process equivalent");
}

}  // namespace dynsyslab::core
