#pragma once

// Adaptive load balancing of a partitioned PDE solve (port of the Python
// AdaptiveMeshCoordinator and RepartitionRecord, with a working rebalance).
//
// The load of partition r is a positive weight w_r (measured right-hand-side wall time, or a
// caller's estimate). The imbalance is max(w) / mean(w) (Hendrickson & Devine, "Dynamic load
// balancing in computational mechanics", Comput. Methods Appl. Mech. Engrg. 184 (2000)
// 485-500). A rebalance spreads each partition's load uniformly over its cells, sums the cell
// costs along every axis and recuts each axis so its parts carry equal cumulative cost, the
// tensor-product form of recursive coordinate bisection (Berger & Bokhari, IEEE Trans.
// Comput. C-36 (1987) 570-580). The new cuts are kept only when the predicted imbalance
// improves; either way a RepartitionRecord documents the decision.

#include <dynsyslab/core/PDEParallel.hpp>

#include <cstddef>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace dynsyslab::core {

using AxisCuts = std::vector<std::vector<std::pair<std::size_t, std::size_t>>>;

struct RepartitionRecord {
    std::string createdUtc;
    std::string reason;
    std::vector<double> priorWeights;
    double imbalance{};            // max / mean of priorWeights
    double predictedImbalance{};   // of the new cuts under the same cell costs
    bool checkpointRequired{};     // policy.checkpointBeforeRepartition
    bool repartitioned{};          // false when the new cuts would not improve the balance
    AxisCuts priorCuts;
    AxisCuts cuts;

    [[nodiscard]] std::string json() const;   // a dynsyslab.pde-repartition/v1 object
};

[[nodiscard]] AxisCuts axisCuts(const ParallelMeshPlan& plan);

// max(weights) / mean(weights); one positive finite weight per partition is required.
[[nodiscard]] double loadImbalance(const ParallelMeshPlan& plan, const std::vector<double>& weights);

[[nodiscard]] bool needsRebalance(const ParallelMeshPlan& plan, const std::vector<double>& weights);

// The rebalanced plan (same partition shape, ghost width, workers, communication and policy;
// partition weights set to the predicted loads) and its record.
[[nodiscard]] std::pair<ParallelMeshPlan, RepartitionRecord> rebalanceParallelMeshPlan(
    const ParallelMeshPlan& plan, const std::vector<double>& weights,
    const std::string& reason = "adaptive load imbalance");

struct RepartitionLog {
    mutable std::mutex mutex;
    std::vector<RepartitionRecord> records;
    [[nodiscard]] std::vector<RepartitionRecord> snapshot() const;
};

// A partitioned right-hand side that measures each partition's wall time and, every
// plan.adaptive.rebalanceIntervalCalls calls when plan.adaptive.enabled, rebalances when the
// accumulated imbalance exceeds plan.adaptive.rebalanceImbalance. Repartitioning never
// changes the computed derivative. Records are appended to `log` (may be null).
[[nodiscard]] PDERightHandSide adaptivePartitionedPDERightHandSide(
    ParallelMeshPlan plan, PDESubdomainRightHandSide rhs, std::shared_ptr<RepartitionLog> log);

}  // namespace dynsyslab::core
