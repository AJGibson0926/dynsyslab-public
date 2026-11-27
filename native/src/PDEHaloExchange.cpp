// Rank-local halo exchange for partitioned PDE solves (see PDEHaloExchange.hpp for the
// method and its references). The MPI transport is PDEHaloExchangeMpi.cpp.
#include <dynsyslab/core/PDEHaloExchange.hpp>

#include <stdexcept>
#include <string>

namespace dynsyslab::core {
namespace {

std::size_t cellsIn(const HaloRegion& region) {
    std::size_t count = 1;
    for (const auto& [begin, end] : region) count *= end - begin;
    return count;
}

// Visit every cell of `region` (row-major) as its linear index in an array of `shape`.
template <typename Visit>
void forEachCell(const std::vector<std::size_t>& shape, const HaloRegion& region, Visit visit) {
    const std::size_t count = cellsIn(region), dimension = shape.size();
    std::vector<std::size_t> coordinate(dimension);
    for (std::size_t cell = 0; cell < count; ++cell) {
        std::size_t rest = cell;
        for (std::size_t axis = dimension; axis-- > 0;) {
            const std::size_t extent = region[axis].second - region[axis].first;
            coordinate[axis] = region[axis].first + rest % extent;
            rest /= extent;
        }
        std::size_t linear = 0;
        for (std::size_t axis = 0; axis < dimension; ++axis) linear = linear * shape[axis] + coordinate[axis];
        visit(cell, linear, coordinate);
    }
}

PDEBoundary boundaryOf(const ParallelMeshPlan& plan, std::size_t axis) {
    return plan.mesh.boundary.empty() ? PDEBoundary::edge : plan.mesh.boundary[axis];
}

HaloRegion layers(const ParallelMeshPlan& plan, const MeshPartition& partition, std::size_t axis,
                  std::pair<std::size_t, std::size_t> alongAxis) {
    const auto interior = partition.shape();
    const std::size_t w = plan.ghostWidth;
    HaloRegion region(interior.size());
    for (std::size_t b = 0; b < interior.size(); ++b)
        region[b] = b == axis ? alongAxis
                  : b < axis ? std::pair<std::size_t, std::size_t>{0, interior[b] + 2 * w}
                             : std::pair<std::size_t, std::size_t>{w, interior[b] + w};
    return region;
}

}  // namespace

std::vector<std::size_t> paddedShape(const ParallelMeshPlan& plan, const MeshPartition& partition) {
    auto shape = partition.shape();
    for (auto& extent : shape) extent += 2 * plan.ghostWidth;
    return shape;
}

std::optional<std::size_t> haloSource(const ParallelMeshPlan& plan, std::size_t rank,
                                      std::size_t axis, bool high) {
    const auto& partition = plan.partitions.at(rank);
    const auto found = partition.neighbors.find("axis" + std::to_string(axis) + (high ? "+" : "-"));
    if (found != partition.neighbors.end()) return found->second;
    if (boundaryOf(plan, axis) != PDEBoundary::periodic) return std::nullopt;
    // Periodic wrap: the partition at the far end of this axis (rank order is row-major
    // over partitionShape, as buildParallelMeshPlanFromCuts assigns it).
    const auto& shape = plan.partitionShape;
    std::vector<std::size_t> coordinate(shape.size());
    std::size_t rest = rank;
    for (std::size_t b = shape.size(); b-- > 0;) { coordinate[b] = rest % shape[b]; rest /= shape[b]; }
    coordinate[axis] = high ? 0 : shape[axis] - 1;
    std::size_t result = 0;
    for (std::size_t b = 0; b < shape.size(); ++b) result = result * shape[b] + coordinate[b];
    return result;
}

namespace detail {

std::vector<double> padInterior(const ParallelMeshPlan& plan, const MeshPartition& partition,
                                const double* interior, std::size_t components) {
    const auto extents = partition.shape();
    for (const auto extent : extents)
        if (extent < plan.ghostWidth)
            throw std::invalid_argument("halo exchange needs at least ghost_width interior cells per axis");
    const auto shape = paddedShape(plan, partition);
    std::size_t cells = 1;
    for (const auto extent : shape) cells *= extent;
    std::vector<double> padded(cells * components, 0.0);
    HaloRegion inner;
    for (const auto extent : extents) inner.emplace_back(plan.ghostWidth, plan.ghostWidth + extent);
    unpackRegion(padded, shape, inner, components, interior);
    return padded;
}

HaloRegion haloDestination(const ParallelMeshPlan& plan, const MeshPartition& partition,
                           std::size_t axis, bool high) {
    const std::size_t w = plan.ghostWidth, n = partition.shape().at(axis);
    return layers(plan, partition, axis, high ? std::pair<std::size_t, std::size_t>{n + w, n + 2 * w}
                                              : std::pair<std::size_t, std::size_t>{0, w});
}

HaloRegion haloSourceRegion(const ParallelMeshPlan& plan, const MeshPartition& source,
                            std::size_t axis, bool destinationHigh) {
    const std::size_t w = plan.ghostWidth, n = source.shape().at(axis);
    // A neighbour's high halo is this partition's first layers; its low halo our last ones.
    return layers(plan, source, axis, destinationHigh ? std::pair<std::size_t, std::size_t>{w, 2 * w}
                                                      : std::pair<std::size_t, std::size_t>{n, n + w});
}

std::vector<double> packRegion(const std::vector<double>& padded, const std::vector<std::size_t>& shape,
                               const HaloRegion& region, std::size_t components) {
    std::vector<double> data(cellsIn(region) * components);
    forEachCell(shape, region, [&](std::size_t cell, std::size_t linear, const auto&) {
        for (std::size_t c = 0; c < components; ++c) data[cell * components + c] = padded[linear * components + c];
    });
    return data;
}

void unpackRegion(std::vector<double>& padded, const std::vector<std::size_t>& shape,
                  const HaloRegion& region, std::size_t components, const double* data) {
    forEachCell(shape, region, [&](std::size_t cell, std::size_t linear, const auto&) {
        for (std::size_t c = 0; c < components; ++c) padded[linear * components + c] = data[cell * components + c];
    });
}

void fillPhysicalBoundary(const ParallelMeshPlan& plan, const MeshPartition& partition,
                          std::vector<double>& padded, std::size_t components, std::size_t axis,
                          bool high) {
    const auto shape = paddedShape(plan, partition);
    const auto region = haloDestination(plan, partition, axis, high);
    const bool zero = boundaryOf(plan, axis) == PDEBoundary::zeroValue;
    const std::size_t edge = high ? shape[axis] - plan.ghostWidth - 1 : plan.ghostWidth;
    forEachCell(shape, region, [&](std::size_t, std::size_t linear, const auto& coordinate) {
        std::size_t source = 0;
        for (std::size_t b = 0; b < shape.size(); ++b) source = source * shape[b] + (b == axis ? edge : coordinate[b]);
        for (std::size_t c = 0; c < components; ++c)
            padded[linear * components + c] = zero ? 0.0 : padded[source * components + c];
    });
}

}  // namespace detail

std::vector<std::vector<double>> exchangeHalos(const ParallelMeshPlan& plan,
                                               const std::vector<std::vector<double>>& interiors,
                                               std::size_t components) {
    if (interiors.size() != plan.partitions.size())
        throw std::invalid_argument("halo exchange needs one interior field per partition");
    std::vector<std::vector<double>> padded;
    std::vector<std::vector<std::size_t>> shapes;
    for (const auto& partition : plan.partitions) {
        std::size_t cells = 1;
        for (const auto extent : partition.shape()) cells *= extent;
        if (interiors[partition.rank].size() != cells * components)
            throw std::invalid_argument("halo exchange interior field has the wrong size");
        padded.push_back(detail::padInterior(plan, partition, interiors[partition.rank].data(), components));
        shapes.push_back(paddedShape(plan, partition));
    }
    // Axis by axis: a pass writes only ghost layers of its axis and reads interior layers of
    // it, so passes over different ranks never conflict.
    for (std::size_t axis = 0; axis < plan.mesh.shape.size(); ++axis)
        for (const auto& partition : plan.partitions)
            for (const bool high : {false, true}) {
                const auto source = haloSource(plan, partition.rank, axis, high);
                if (!source) {
                    detail::fillPhysicalBoundary(plan, partition, padded[partition.rank], components, axis, high);
                    continue;
                }
                const auto data = detail::packRegion(padded[*source], shapes[*source],
                    detail::haloSourceRegion(plan, plan.partitions[*source], axis, high), components);
                const auto target = detail::haloDestination(plan, partition, axis, high);
                if (data.size() != cellsIn(target) * components)
                    throw std::invalid_argument("halo exchange neighbours do not share their cut");
                detail::unpackRegion(padded[partition.rank], shapes[partition.rank], target, components, data.data());
            }
    return padded;
}

}  // namespace dynsyslab::core
