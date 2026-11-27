#pragma once

// Rank-local ghost (halo) exchange for partitioned PDE solves.
//
// Each partition owns only its interior; its padded field (interior + ghostWidth layers per
// side) is completed axis by axis: along axis a it receives the neighbour's ghostWidth
// boundary layers, spanning the already-completed ghost layers of the axes before a, so
// edge and corner ghosts equal those of a global read (the dimension-by-dimension halo
// exchange of distributed stencil codes; Gropp, Lusk & Skjellum, Using MPI, 3rd ed., MIT
// Press 2014, ch. 4). Physical boundaries follow StructuredMesh::boundary: edge replicates the
// boundary layer, zeroValue fills zeros, periodic exchanges with the partition at the far end.
//
// exchangeHalos moves the packed buffers in process (the thread-parallel equivalent);
// mpiExchangeRankField moves the same buffers with MPI_Isend/MPI_Irecv when the library is
// built with DYNSYSLAB_WITH_MPI (MPI Forum, MPI: A Message-Passing Interface Standard 4.0,
// 2021). Every partition needs at least ghostWidth interior cells along every axis.

#include <dynsyslab/core/PDEParallel.hpp>

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace dynsyslab::core {

using HaloRegion = std::vector<std::pair<std::size_t, std::size_t>>;   // half-open per axis

// Interior extent + 2 ghostWidth along every axis.
[[nodiscard]] std::vector<std::size_t> paddedShape(const ParallelMeshPlan& plan,
                                                   const MeshPartition& partition);

// The rank that supplies `rank`'s halo on side `high` of `axis` (the periodic wrap included);
// nullopt at a non-periodic physical boundary.
[[nodiscard]] std::optional<std::size_t> haloSource(const ParallelMeshPlan& plan, std::size_t rank,
                                                    std::size_t axis, bool high);

// In process: interiors[rank] (row-major interior cells x components) -> padded fields.
[[nodiscard]] std::vector<std::vector<double>> exchangeHalos(
    const ParallelMeshPlan& plan, const std::vector<std::vector<double>>& interiors,
    std::size_t components);

// True when the library was built with DYNSYSLAB_WITH_MPI.
[[nodiscard]] bool mpiHaloExchangeAvailable();

// One MPI rank (MPI_COMM_WORLD rank r owns partition r; the world size must equal the plan's
// worker count; MPI must be initialised by the caller): the padded field of this rank.
// Throws std::runtime_error in a build without MPI.
[[nodiscard]] std::vector<double> mpiExchangeRankField(const ParallelMeshPlan& plan,
                                                       const double* interior,
                                                       std::size_t components);

namespace detail {
// Building blocks shared by the in-process and MPI transports.
[[nodiscard]] std::vector<double> padInterior(const ParallelMeshPlan& plan, const MeshPartition& partition,
                                              const double* interior, std::size_t components);
// Ghost layers of `partition` on side `high` of `axis` in its padded field.
[[nodiscard]] HaloRegion haloDestination(const ParallelMeshPlan& plan, const MeshPartition& partition,
                                         std::size_t axis, bool high);
// Layers of `source`'s padded field that fill a neighbour's halo on side `destinationHigh`.
[[nodiscard]] HaloRegion haloSourceRegion(const ParallelMeshPlan& plan, const MeshPartition& source,
                                          std::size_t axis, bool destinationHigh);
[[nodiscard]] std::vector<double> packRegion(const std::vector<double>& padded,
                                             const std::vector<std::size_t>& shape,
                                             const HaloRegion& region, std::size_t components);
void unpackRegion(std::vector<double>& padded, const std::vector<std::size_t>& shape,
                  const HaloRegion& region, std::size_t components, const double* data);
// Edge or zero ghosts at a non-periodic physical boundary.
void fillPhysicalBoundary(const ParallelMeshPlan& plan, const MeshPartition& partition,
                          std::vector<double>& padded, std::size_t components, std::size_t axis,
                          bool high);
}  // namespace detail

}  // namespace dynsyslab::core
