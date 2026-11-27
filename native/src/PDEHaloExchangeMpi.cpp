// MPI transport of the rank-local halo exchange (built only with DYNSYSLAB_WITH_MPI).
// The same packed regions as the in-process exchange (PDEHaloExchange.cpp) move through
// nonblocking point-to-point messages, one axis at a time (MPI Forum, MPI: A Message-Passing
// Interface Standard 4.0, 2021, sec. 3.7; Gropp, Lusk & Skjellum, Using MPI, 3rd ed., ch. 4).
// Port of the Python MPIExchange.exchange_rank_field.
#include <dynsyslab/core/PDEHaloExchange.hpp>

#include <mpi.h>

#include <stdexcept>

namespace dynsyslab::core {

bool mpiHaloExchangeAvailable() { return true; }

std::vector<double> mpiExchangeRankField(const ParallelMeshPlan& plan, const double* interior,
                                         std::size_t components) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) throw std::runtime_error("mpiExchangeRankField needs MPI_Init first");
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (static_cast<std::size_t>(size) != plan.workerCount())
        throw std::invalid_argument("MPI world size must equal the parallel mesh plan's worker count");
    const auto& partition = plan.partitions.at(static_cast<std::size_t>(rank));
    const auto shape = paddedShape(plan, partition);
    auto padded = detail::padInterior(plan, partition, interior, components);
    const auto tag = [](std::size_t axis, bool destinationHigh) {
        return static_cast<int>(2 * axis + (destinationHigh ? 1 : 0));
    };
    for (std::size_t axis = 0; axis < plan.mesh.shape.size(); ++axis) {
        std::vector<std::vector<double>> sends, receives;
        std::vector<HaloRegion> targets;
        std::vector<MPI_Request> requests;
        sends.reserve(2); receives.reserve(2);
        for (const bool high : {false, true}) {
            const auto source = haloSource(plan, static_cast<std::size_t>(rank), axis, high);
            if (!source) {
                detail::fillPhysicalBoundary(plan, partition, padded, components, axis, high);
                continue;
            }
            // Receive this rank's `high` halo from `source`; send `source`'s opposite halo.
            const auto target = detail::haloDestination(plan, partition, axis, high);
            std::size_t cells = 1;
            for (const auto& [begin, end] : target) cells *= end - begin;
            receives.emplace_back(cells * components);
            targets.push_back(target);
            requests.emplace_back();
            MPI_Irecv(receives.back().data(), static_cast<int>(receives.back().size()), MPI_DOUBLE,
                      static_cast<int>(*source), tag(axis, high), MPI_COMM_WORLD, &requests.back());
            sends.push_back(detail::packRegion(padded, shape,
                detail::haloSourceRegion(plan, partition, axis, !high), components));
            requests.emplace_back();
            MPI_Isend(sends.back().data(), static_cast<int>(sends.back().size()), MPI_DOUBLE,
                      static_cast<int>(*source), tag(axis, !high), MPI_COMM_WORLD, &requests.back());
        }
        if (!requests.empty())
            MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        for (std::size_t index = 0; index < targets.size(); ++index)
            detail::unpackRegion(padded, shape, targets[index], components, receives[index].data());
    }
    return padded;
}

}  // namespace dynsyslab::core
