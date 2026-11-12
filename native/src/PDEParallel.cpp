#include <dynsyslab/core/PDEParallel.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace dynsyslab::core {
namespace {

std::size_t product(const std::vector<std::size_t>& shape) {
    std::size_t result = 1;
    for (const auto value : shape) {
        if (!value || result > (std::numeric_limits<std::size_t>::max)() / value)
            throw std::invalid_argument("parallel PDE extent is empty or too large");
        result *= value;
    }
    return result;
}

std::vector<std::size_t> decode(std::size_t linear,
                                const std::vector<std::size_t>& shape) {
    std::vector<std::size_t> coordinate(shape.size());
    for (std::size_t axis = shape.size(); axis-- > 0;) {
        coordinate[axis] = linear % shape[axis];
        linear /= shape[axis];
    }
    return coordinate;
}

std::size_t encode(const std::vector<std::size_t>& coordinate,
                   const std::vector<std::size_t>& shape) {
    std::size_t result = 0;
    for (std::size_t axis = 0; axis < shape.size(); ++axis)
        result = result * shape[axis] + coordinate[axis];
    return result;
}

std::vector<std::pair<std::size_t, std::size_t>> cuts(
    std::size_t cells, std::size_t parts) {
    if (!parts || parts > cells)
        throw std::invalid_argument("parallel PDE has more partitions than cells");
    const std::size_t base = cells / parts, remainder = cells % parts;
    std::vector<std::pair<std::size_t, std::size_t>> result;
    result.reserve(parts);
    std::size_t start = 0;
    for (std::size_t index = 0; index < parts; ++index) {
        const std::size_t end = start + base + (index < remainder ? 1 : 0);
        result.emplace_back(start, end);
        start = end;
    }
    return result;
}

PDEBoundary boundary(const StructuredMesh& mesh, std::size_t axis) {
    return mesh.boundary.empty() ? PDEBoundary::edge : mesh.boundary[axis];
}

std::vector<double> ghostedField(const ParallelMeshPlan& plan,
                                 const MeshPartition& partition,
                                 const double* field, std::size_t components,
                                 std::vector<std::size_t>& shape) {
    const auto interiorShape = partition.shape();
    shape.resize(interiorShape.size());
    for (std::size_t axis = 0; axis < shape.size(); ++axis) {
        if (plan.ghostWidth > (std::numeric_limits<std::size_t>::max)() / 2)
            throw std::invalid_argument("parallel PDE ghost extent is too large");
        if (interiorShape[axis] > (std::numeric_limits<std::size_t>::max)()
            - 2 * plan.ghostWidth)
            throw std::invalid_argument("parallel PDE ghost extent is too large");
        shape[axis] = interiorShape[axis] + 2 * plan.ghostWidth;
    }
    const std::size_t cells = product(shape);
    if (components > (std::numeric_limits<std::size_t>::max)() / cells)
        throw std::invalid_argument("parallel PDE ghost field is too large");
    std::vector<double> result(cells * components, 0.0);
    for (std::size_t local = 0; local < cells; ++local) {
        const auto coordinate = decode(local, shape);
        std::vector<std::size_t> global(shape.size());
        bool zero = false;
        for (std::size_t axis = 0; axis < shape.size(); ++axis) {
            const auto relative = static_cast<std::ptrdiff_t>(coordinate[axis])
                - static_cast<std::ptrdiff_t>(plan.ghostWidth);
            auto index = static_cast<std::ptrdiff_t>(partition.interior[axis].first)
                + relative;
            const auto axisCells = static_cast<std::ptrdiff_t>(plan.mesh.shape[axis]);
            if (index < 0 || index >= axisCells) {
                const auto mode = boundary(plan.mesh, axis);
                if (mode == PDEBoundary::zeroValue) zero = true;
                else if (mode == PDEBoundary::periodic)
                    index = (index % axisCells + axisCells) % axisCells;
                else index = std::clamp<std::ptrdiff_t>(index, 0, axisCells - 1);
            }
            global[axis] = static_cast<std::size_t>(index);
        }
        if (zero) continue;
        const std::size_t source = encode(global, plan.mesh.shape) * components;
        std::copy(field + source, field + source + components,
                  result.begin() + static_cast<std::ptrdiff_t>(local * components));
    }
    return result;
}

void scatter(const ParallelMeshPlan& plan, const MeshPartition& partition,
             const std::vector<double>& local, std::size_t components,
             double* derivative) {
    const auto shape = partition.shape();
    for (std::size_t cell = 0; cell < product(shape); ++cell) {
        auto global = decode(cell, shape);
        for (std::size_t axis = 0; axis < shape.size(); ++axis)
            global[axis] += partition.interior[axis].first;
        const std::size_t target = encode(global, plan.mesh.shape) * components;
        std::copy(local.begin() + static_cast<std::ptrdiff_t>(cell * components),
                  local.begin() + static_cast<std::ptrdiff_t>((cell + 1) * components),
                  derivative + target);
    }
}

}  // namespace

std::vector<std::size_t> MeshPartition::shape() const {
    std::vector<std::size_t> result;
    result.reserve(interior.size());
    for (const auto& range : interior) {
        if (range.second <= range.first)
            throw std::invalid_argument("parallel PDE partition range is empty");
        result.push_back(range.second - range.first);
    }
    return result;
}

std::size_t ParallelMeshPlan::workerCount() const { return product(partitionShape); }

void ParallelMeshPlan::validate() const {
    mesh.validate();
    if (partitionShape.size() != mesh.shape.size() || !ghostWidth || !localWorkers)
        throw std::invalid_argument("parallel PDE plan dimensions or worker counts are invalid");
    for (std::size_t axis = 0; axis < partitionShape.size(); ++axis)
        if (!partitionShape[axis] || partitionShape[axis] > mesh.shape[axis])
            throw std::invalid_argument("parallel PDE partition shape is invalid");
    if (partitions.size() != workerCount())
        throw std::invalid_argument("parallel PDE plan has incomplete partitions");
    std::vector<unsigned char> coverage(mesh.cellCount(), 0);
    for (std::size_t rank = 0; rank < partitions.size(); ++rank) {
        const auto& partition = partitions[rank];
        if (partition.rank != rank || partition.interior.size() != mesh.shape.size())
            throw std::invalid_argument("parallel PDE partition identity is invalid");
        for (std::size_t axis = 0; axis < mesh.shape.size(); ++axis)
            if (partition.interior[axis].second > mesh.shape[axis])
                throw std::invalid_argument("parallel PDE partition exceeds its mesh");
        const auto localShape = partition.shape();
        for (std::size_t cell = 0; cell < product(localShape); ++cell) {
            auto global = decode(cell, localShape);
            for (std::size_t axis = 0; axis < global.size(); ++axis)
                global[axis] += partition.interior[axis].first;
            auto& count = coverage[encode(global, mesh.shape)];
            if (count++) throw std::invalid_argument("parallel PDE partitions overlap");
        }
    }
    if (std::any_of(coverage.begin(), coverage.end(),
                    [](unsigned char count) { return count != 1; }))
        throw std::invalid_argument("parallel PDE partitions do not cover the mesh");
}

ParallelMeshPlan buildParallelMeshPlan(
    StructuredMesh mesh, std::vector<std::size_t> partitionShape,
    std::size_t ghostWidth, std::size_t localWorkers) {
    mesh.validate();
    if (partitionShape.size() != mesh.shape.size())
        throw std::invalid_argument("parallel PDE partition shape must match the mesh");
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> axisCuts;
    for (std::size_t axis = 0; axis < mesh.shape.size(); ++axis)
        axisCuts.push_back(cuts(mesh.shape[axis], partitionShape[axis]));
    ParallelMeshPlan result{std::move(mesh), std::move(partitionShape),
                            ghostWidth, localWorkers, {}};
    const std::size_t count = result.workerCount();
    result.localWorkers = std::min(result.localWorkers, count);
    result.partitions.reserve(count);
    for (std::size_t rank = 0; rank < count; ++rank) {
        const auto coordinate = decode(rank, result.partitionShape);
        MeshPartition partition; partition.rank = rank;
        for (std::size_t axis = 0; axis < coordinate.size(); ++axis) {
            partition.interior.push_back(axisCuts[axis][coordinate[axis]]);
            for (const int offset : {-1, 1}) {
                const auto neighbor = static_cast<std::ptrdiff_t>(coordinate[axis]) + offset;
                if (neighbor < 0 || neighbor >= static_cast<std::ptrdiff_t>(result.partitionShape[axis]))
                    continue;
                auto nearby = coordinate; nearby[axis] = static_cast<std::size_t>(neighbor);
                partition.neighbors["axis" + std::to_string(axis)
                    + (offset < 0 ? "-" : "+")] = encode(nearby, result.partitionShape);
            }
        }
        result.partitions.push_back(std::move(partition));
    }
    result.validate();
    return result;
}

PDERightHandSide partitionedPDERightHandSide(
    ParallelMeshPlan plan, PDESubdomainRightHandSide rhs) {
    plan.validate();
    if (!rhs) throw std::invalid_argument("parallel PDE right-hand side is empty");
    return [plan = std::move(plan), rhs = std::move(rhs)](
        double time, const StructuredMesh& mesh, const double* field,
        std::size_t count, double* derivative) {
        if (mesh.shape != plan.mesh.shape || mesh.spacing != plan.mesh.spacing
            || mesh.origin != plan.mesh.origin || mesh.boundary != plan.mesh.boundary
            || count % mesh.cellCount() != 0)
            throw std::invalid_argument("parallel PDE field does not match its plan");
        const std::size_t components = count / mesh.cellCount();
        std::atomic<std::size_t> next{0};
        std::atomic<bool> failed{false};
        std::exception_ptr error;
        std::mutex errorMutex;
        const auto worker = [&] {
            while (!failed.load(std::memory_order_relaxed)) {
                const std::size_t index = next.fetch_add(1);
                if (index >= plan.partitions.size()) return;
                try {
                    const auto& partition = plan.partitions[index];
                    std::vector<std::size_t> ghostShape;
                    const auto ghosted = ghostedField(
                        plan, partition, field, components, ghostShape);
                    std::vector<double> local(product(partition.shape()) * components);
                    rhs(time, mesh, plan, partition, ghosted.data(), ghostShape,
                        components, local.data());
                    if (std::any_of(local.begin(), local.end(),
                        [](double value) { return !std::isfinite(value); }))
                        throw std::runtime_error("parallel PDE right-hand side is non-finite");
                    scatter(plan, partition, local, components, derivative);
                } catch (...) {
                    failed.store(true, std::memory_order_relaxed);
                    std::lock_guard<std::mutex> lock(errorMutex);
                    if (!error) error = std::current_exception();
                }
            }
        };
        const std::size_t workers = std::min(plan.localWorkers, plan.partitions.size());
        std::vector<std::thread> threads;
        threads.reserve(workers > 0 ? workers - 1 : 0);
        for (std::size_t index = 1; index < workers; ++index)
            threads.emplace_back(worker);
        worker();
        for (auto& thread : threads) thread.join();
        if (error) std::rethrow_exception(error);
    };
}

PDESubdomainRightHandSide heatEquationSubdomainRightHandSide(
    const SystemConfiguration& configuration) {
    const auto resolved = resolveHeatEquation(configuration);
    const double diffusivity = resolved.parameters.at("diffusivity");
    const double decay = resolved.parameters.at("decay_rate");
    const int boundaryCode = static_cast<int>(resolved.parameters.at("boundary_code"));
    return [diffusivity, decay, boundaryCode](double, const StructuredMesh& mesh,
        const ParallelMeshPlan& plan, const MeshPartition& partition,
        const double* field, const std::vector<std::size_t>& ghostShape,
        std::size_t components, double* derivative) {
        if (mesh.shape.size() != 1 || partition.interior.size() != 1
            || ghostShape.size() != 1 || components != 1)
            throw std::invalid_argument("partitioned heat equation requires a scalar 1-D field");
        const std::size_t cells = partition.shape()[0];
        const double scale = diffusivity / (mesh.spacing[0] * mesh.spacing[0]);
        for (std::size_t index = 0; index < cells; ++index) {
            const std::size_t global = partition.interior[0].first + index;
            if (boundaryCode == 2 && (global == 0 || global + 1 == mesh.shape[0])) {
                derivative[index] = 0.0;
                continue;
            }
            const std::size_t center = plan.ghostWidth + index;
            derivative[index] = scale * (field[center - 1] - 2.0 * field[center]
                + field[center + 1]) - decay * field[center];
        }
    };
}

}  // namespace dynsyslab::core
