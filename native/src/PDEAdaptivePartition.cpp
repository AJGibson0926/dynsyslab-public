// Adaptive load balancing of partitioned PDE solves; method and references in
// PDEAdaptivePartition.hpp.
#include <dynsyslab/core/PDEAdaptivePartition.hpp>

#include "PDEParallelInternal.hpp"
#include <dynsyslab/core/Json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

std::size_t cellsOf(const MeshPartition& partition) {
    std::size_t count = 1;
    for (const auto extent : partition.shape()) count *= extent;
    return count;
}

std::string utcNow() {
    const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm parts{};
#ifdef _WIN32
    gmtime_s(&parts, &now);
#else
    gmtime_r(&now, &parts);
#endif
    char text[32];
    std::strftime(text, sizeof(text), "%Y-%m-%dT%H:%M:%SZ", &parts);
    return text;
}

// Per-cell cost: each partition's weight spread uniformly over its cells.
std::vector<double> cellCosts(const ParallelMeshPlan& plan, const std::vector<double>& weights) {
    std::vector<double> cost(plan.mesh.cellCount(), 0.0);
    const std::size_t dimension = plan.mesh.shape.size();
    for (const auto& partition : plan.partitions) {
        const double density = weights[partition.rank] / static_cast<double>(cellsOf(partition));
        const auto shape = partition.shape();
        std::vector<std::size_t> local(dimension);
        for (std::size_t cell = 0; cell < cellsOf(partition); ++cell) {
            std::size_t rest = cell, linear = 0;
            for (std::size_t axis = dimension; axis-- > 0;) { local[axis] = rest % shape[axis]; rest /= shape[axis]; }
            for (std::size_t axis = 0; axis < dimension; ++axis)
                linear = linear * plan.mesh.shape[axis] + partition.interior[axis].first + local[axis];
            cost[linear] = density;
        }
    }
    return cost;
}

std::vector<double> partitionLoads(const ParallelMeshPlan& plan, const std::vector<double>& cost) {
    std::vector<double> loads(plan.partitions.size(), 0.0);
    const std::size_t dimension = plan.mesh.shape.size();
    std::vector<std::size_t> coordinate(dimension);
    for (std::size_t linear = 0; linear < cost.size(); ++linear) {
        std::size_t rest = linear;
        for (std::size_t axis = dimension; axis-- > 0;) { coordinate[axis] = rest % plan.mesh.shape[axis]; rest /= plan.mesh.shape[axis]; }
        for (const auto& partition : plan.partitions) {
            bool inside = true;
            for (std::size_t axis = 0; axis < dimension && inside; ++axis)
                inside = coordinate[axis] >= partition.interior[axis].first
                      && coordinate[axis] < partition.interior[axis].second;
            if (inside) { loads[partition.rank] += cost[linear]; break; }
        }
    }
    return loads;
}

// Equal-cost cuts of one axis from its marginal cost profile; every part keeps >= 1 cell.
std::vector<std::pair<std::size_t, std::size_t>> recut(const std::vector<double>& marginal, std::size_t parts) {
    const std::size_t cells = marginal.size();
    std::vector<double> cumulative(cells + 1, 0.0);
    for (std::size_t index = 0; index < cells; ++index) cumulative[index + 1] = cumulative[index] + marginal[index];
    std::vector<std::pair<std::size_t, std::size_t>> result;
    std::size_t start = 0;
    for (std::size_t part = 1; part < parts; ++part) {
        const double target = cumulative[cells] * static_cast<double>(part) / static_cast<double>(parts);
        std::size_t best = start + 1;
        for (std::size_t cut = start + 1; cut <= cells - (parts - part); ++cut)
            if (std::abs(cumulative[cut] - target) < std::abs(cumulative[best] - target)) best = cut;
        result.emplace_back(start, best);
        start = best;
    }
    result.emplace_back(start, cells);
    return result;
}

json::Value cutsJson(const AxisCuts& cuts) {
    json::Value::Array axes;
    for (const auto& axis : cuts) {
        json::Value::Array ranges;
        for (const auto& [begin, end] : axis)
            ranges.emplace_back(json::Value::Array{static_cast<double>(begin), static_cast<double>(end)});
        axes.emplace_back(std::move(ranges));
    }
    return json::Value(std::move(axes));
}

}  // namespace

std::string RepartitionRecord::json() const {
    json::Value::Array weights;
    for (const auto weight : priorWeights) weights.emplace_back(weight);
    return json::dump(json::Value(json::Value::Object{
        {"schema", "dynsyslab.pde-repartition/v1"}, {"created_utc", createdUtc}, {"reason", reason},
        {"prior_weights", json::Value(std::move(weights))}, {"imbalance", imbalance},
        {"predicted_imbalance", predictedImbalance}, {"checkpoint_required", checkpointRequired},
        {"repartitioned", repartitioned}, {"prior_cuts", cutsJson(priorCuts)}, {"cuts", cutsJson(cuts)}}));
}

AxisCuts axisCuts(const ParallelMeshPlan& plan) {
    AxisCuts result(plan.mesh.shape.size());
    for (const auto& partition : plan.partitions)
        for (std::size_t axis = 0; axis < result.size(); ++axis)
            if (std::find(result[axis].begin(), result[axis].end(), partition.interior[axis]) == result[axis].end())
                result[axis].push_back(partition.interior[axis]);
    for (auto& axis : result) std::sort(axis.begin(), axis.end());
    return result;
}

double loadImbalance(const ParallelMeshPlan& plan, const std::vector<double>& weights) {
    if (weights.size() != plan.partitions.size()
        || std::any_of(weights.begin(), weights.end(), [](double w) { return !(w > 0.0) || !std::isfinite(w); }))
        throw std::invalid_argument("one positive load weight per mesh partition is required");
    const double mean = std::accumulate(weights.begin(), weights.end(), 0.0) / static_cast<double>(weights.size());
    return *std::max_element(weights.begin(), weights.end()) / mean;
}

bool needsRebalance(const ParallelMeshPlan& plan, const std::vector<double>& weights) {
    return loadImbalance(plan, weights) > plan.adaptive.rebalanceImbalance;
}

std::pair<ParallelMeshPlan, RepartitionRecord> rebalanceParallelMeshPlan(
    const ParallelMeshPlan& plan, const std::vector<double>& weights, const std::string& reason) {
    plan.validate();
    RepartitionRecord record;
    record.createdUtc = utcNow();
    record.reason = reason;
    record.priorWeights = weights;
    record.imbalance = loadImbalance(plan, weights);
    record.checkpointRequired = plan.adaptive.checkpointBeforeRepartition;
    record.priorCuts = axisCuts(plan);
    const auto cost = cellCosts(plan, weights);
    const std::size_t dimension = plan.mesh.shape.size();
    AxisCuts cuts(dimension);
    for (std::size_t axis = 0; axis < dimension; ++axis) {
        std::vector<double> marginal(plan.mesh.shape[axis], 0.0);
        std::vector<std::size_t> coordinate(dimension);
        for (std::size_t linear = 0; linear < cost.size(); ++linear) {
            std::size_t rest = linear;
            for (std::size_t b = dimension; b-- > 0;) { coordinate[b] = rest % plan.mesh.shape[b]; rest /= plan.mesh.shape[b]; }
            marginal[coordinate[axis]] += cost[linear];
        }
        cuts[axis] = recut(marginal, plan.partitionShape[axis]);
    }
    auto next = buildParallelMeshPlanFromCuts(plan.mesh, cuts, plan.ghostWidth, plan.localWorkers);
    next.communication = plan.communication;
    next.adaptive = plan.adaptive;
    const auto loads = partitionLoads(next, cost);
    record.predictedImbalance = loadImbalance(next, loads);
    record.repartitioned = record.predictedImbalance < record.imbalance * (1.0 - 1.0e-12);
    if (!record.repartitioned) {
        record.cuts = record.priorCuts;
        record.predictedImbalance = record.imbalance;
        auto same = plan;
        for (auto& partition : same.partitions) partition.weight = weights[partition.rank];
        return {std::move(same), std::move(record)};
    }
    record.cuts = cuts;
    for (auto& partition : next.partitions) partition.weight = loads[partition.rank];
    return {std::move(next), std::move(record)};
}

std::vector<RepartitionRecord> RepartitionLog::snapshot() const {
    std::lock_guard<std::mutex> lock(mutex);
    return records;
}

PDERightHandSide adaptivePartitionedPDERightHandSide(
    ParallelMeshPlan plan, PDESubdomainRightHandSide rhs, std::shared_ptr<RepartitionLog> log) {
    plan.validate();
    if (!rhs) throw std::invalid_argument("parallel PDE right-hand side is empty");
    if (plan.communication == PDECommunication::mpi)
        throw std::invalid_argument(detail::mpiInProcessMessage);
    if (!plan.adaptive.rebalanceIntervalCalls)
        throw std::invalid_argument("adaptive rebalance interval must be positive");
    struct State {
        ParallelMeshPlan plan;
        std::vector<double> accumulated;
        std::size_t calls{};
    };
    auto state = std::make_shared<State>(State{std::move(plan), {}, 0});
    state->accumulated.assign(state->plan.partitions.size(), 0.0);
    return [state, rhs = std::move(rhs), log = std::move(log)](
        double time, const StructuredMesh& mesh, const double* field, std::size_t count, double* derivative) {
        std::vector<double> seconds;
        detail::evaluatePartitioned(state->plan, rhs, time, mesh, field, count, derivative, &seconds);
        // A cell-proportional floor keeps weights positive below the clock resolution.
        for (const auto& partition : state->plan.partitions)
            state->accumulated[partition.rank] += seconds[partition.rank] + 1.0e-9 * static_cast<double>(cellsOf(partition));
        if (!state->plan.adaptive.enabled || ++state->calls % state->plan.adaptive.rebalanceIntervalCalls) return;
        if (needsRebalance(state->plan, state->accumulated)) {
            auto [next, record] = rebalanceParallelMeshPlan(state->plan, state->accumulated);
            state->plan = std::move(next);
            if (log) { std::lock_guard<std::mutex> lock(log->mutex); log->records.push_back(std::move(record)); }
        }
        state->accumulated.assign(state->plan.partitions.size(), 0.0);
    };
}

}  // namespace dynsyslab::core
