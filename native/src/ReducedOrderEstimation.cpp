#include <dynsyslab/core/ReducedOrderEstimation.hpp>

#include <stdexcept>
#include <string>

namespace dynsyslab::core {
namespace {
// Row-major product of (r x k) and (k x c).
std::vector<double> multiply(const std::vector<double>& a, const std::vector<double>& b,
                             std::size_t r, std::size_t k, std::size_t c) {
    std::vector<double> out(r * c, 0.0);
    for (std::size_t i = 0; i < r; ++i)
        for (std::size_t l = 0; l < k; ++l)
            for (std::size_t j = 0; j < c; ++j) out[i * c + j] += a[i * k + l] * b[l * c + j];
    return out;
}
std::vector<double> transpose(const std::vector<double>& a, std::size_t r, std::size_t c) {
    std::vector<double> out(a.size());
    for (std::size_t i = 0; i < r; ++i) for (std::size_t j = 0; j < c; ++j) out[j * r + i] = a[i * c + j];
    return out;
}
}  // namespace

ReducedEstimatorDesign designReducedLQE(const StateSpaceModel& model, const BalancedReductionOptions& options,
    const std::vector<double>& W, const std::vector<double>& V) {
    validateStateSpace(model);
    const std::size_t n = model.states, p = model.outputs;
    if (!p) throw std::invalid_argument("reduced LQE needs at least one measurement");
    if (W.size() != n * n) throw std::invalid_argument("process covariance must be full states x states");
    ReducedEstimatorDesign out;
    out.reduction = balancedReduction(model, options);
    const std::size_t r = out.reduction.order;
    if (!r) throw std::invalid_argument("balanced reduction removed every state");
    const auto& P = out.reduction.projection;  // r x n
    const auto Wr = multiply(multiply(P, W, r, n, n), transpose(P, r, n), r, n, r);
    out.estimator = designLQE(out.reduction.reduced, Wr, V);
    out.liftedGain = multiply(out.reduction.lift, out.estimator.gain, n, r, p);
    auto error = model;
    const auto LC = multiply(out.liftedGain, model.C, n, p, n);
    for (std::size_t k = 0; k < n * n; ++k) error.A[k] -= LC[k];
    out.liftedObserverPoles = systemPoles(error);
    out.liftedObserverStable = isAsymptoticallyStable(error);
    return out;
}

ReducedEstimatorDesign designReducedKoopmanLQE(const KoopmanModel& model,
    const KoopmanStateSpaceOptions& realization, const BalancedReductionOptions& options,
    const std::vector<double>& W, const std::vector<double>& V) {
    auto system = koopmanStateSpace(model, realization);
    if (realization.controlMap.empty()) {  // process noise on every lifted coordinate
        system.inputs = system.states; system.B.assign(system.states * system.states, 0.0);
        for (std::size_t i = 0; i < system.states; ++i) system.B[i * system.states + i] = 1.0;
        system.D.assign(system.outputs * system.inputs, 0.0);
    }
    return designReducedLQE(system, options, W, V);
}

ReducedEstimatorDesign reduceKoopmanLQEDesign(KoopmanLQEDesign& design,
    const BalancedReductionOptions& options, const std::vector<double>& W, const std::vector<double>& V) {
    const std::size_t n = design.solution.stateDimension, p = design.solution.measurementDimension;
    std::vector<double> identity(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) identity[i * n + i] = 1.0;
    const auto plant = makeStateSpace(design.generator, identity, design.measurementMatrix, {}, n, n, p);
    auto reduced = designReducedLQE(plant, options, W, V);
    if (!reduced.liftedObserverStable)
        throw std::invalid_argument("balanced_truncation: the order-" + std::to_string(reduced.reduction.order)
            + " KLQE gain does not give a stable observer on the full KROM; raise balanced_truncation.order");
    const std::size_t r = reduced.reduction.order;
    const auto& L = reduced.reduction.lift;  // n x r
    design.solution.gain = reduced.liftedGain;
    design.solution.observerPoles = reduced.liftedObserverPoles;
    design.solution.covariance = multiply(multiply(L, reduced.estimator.covariance, n, r, r), transpose(L, n, r), n, r, n);
    design.solution.relativeResidual = reduced.estimator.relativeResidual;
    return reduced;
}

}  // namespace dynsyslab::core
