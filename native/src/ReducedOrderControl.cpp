#include <dynsyslab/core/ReducedOrderControl.hpp>

#include <stdexcept>

namespace dynsyslab::core {
namespace {
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

ReducedRegulatorDesign designReducedLQR(const StateSpaceModel& model, const BalancedReductionOptions& options,
    const std::vector<double>& Q, const std::vector<double>& R) {
    validateStateSpace(model);
    const std::size_t n = model.states, m = model.inputs;
    if (!m) throw std::invalid_argument("reduced LQR needs at least one input");
    if (Q.size() != n * n) throw std::invalid_argument("state penalty must be full states x states");
    ReducedRegulatorDesign out;
    out.reduction = balancedReduction(model, options);
    const std::size_t r = out.reduction.order;
    if (!r) throw std::invalid_argument("balanced reduction removed every state");
    const auto& L = out.reduction.lift;  // n x r
    const auto Qr = multiply(multiply(transpose(L, n, r), Q, r, n, n), L, r, n, r);
    out.regulator = designLQR(out.reduction.reduced, Qr, R);
    out.liftedGain = multiply(out.regulator.gain, out.reduction.projection, m, r, n);
    const auto closed = stateFeedbackClosedLoop(model, out.liftedGain);
    out.liftedClosedLoopPoles = systemPoles(closed);
    out.liftedClosedLoopStable = isAsymptoticallyStable(closed);
    return out;
}

ReducedRegulatorDesign designReducedKoopmanLQR(const KoopmanModel& model,
    const KoopmanStateSpaceOptions& realization, const BalancedReductionOptions& options,
    const std::vector<double>& Q, const std::vector<double>& R) {
    if (realization.controlMap.empty()) throw std::invalid_argument("reduced KLQR needs a control map");
    auto system = koopmanStateSpace(model, realization);
    if (realization.measurements.empty()) {  // performance output: every lifted coordinate
        system.outputs = system.states; system.C.assign(system.states * system.states, 0.0);
        for (std::size_t i = 0; i < system.states; ++i) system.C[i * system.states + i] = 1.0;
        system.D.assign(system.outputs * system.inputs, 0.0);
    }
    return designReducedLQR(system, options, Q, R);
}

}  // namespace dynsyslab::core
