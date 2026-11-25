// Discrete LQ tracking on the measured delay state of an ARX model. The steady gain comes from
// the suite's DARE (designLQR); the feedforward is the backward adjoint recursion for a known
// reference (B. D. O. Anderson and J. B. Moore, Optimal Control: Linear Quadratic Methods,
// Prentice Hall, 1990, ch. 4). The companion structure makes each backward step O(n).
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/IoDelayModel.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

IoDelayTracker::IoDelayTracker(const IoDelayModel& model, double stateWeight, double inputWeight,
                               const std::vector<double>& reference)
    : model_(model) {
    if (!(stateWeight > 0) || !(inputWeight > 0) || !std::isfinite(stateWeight) || !std::isfinite(inputWeight))
        throw std::invalid_argument("io_delay_lqr needs positive finite Q and R");
    if (reference.size() < 2) throw std::invalid_argument("io_delay_lqr needs a reference over its horizon");
    const auto realization = ioDelayRealization(model_);
    const auto& system = realization.system;
    const std::size_t n = system.states, ny = model_.ny, nu = model_.nu;
    std::vector<double> Q(n * n, 0.0);
    Q[0] = stateWeight;
    const auto design = designLQR(system, Q, {inputWeight});
    gain_ = design.gain; poles_ = design.closedLoopPoles; residual_ = design.relativeResidual;

    double BPB = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            BPB += system.B[i] * design.riccati[i * n + j] * system.B[j];
    const double S = inputWeight + BPB;
    std::vector<double> Pd(n);
    for (std::size_t i = 0; i < n; ++i) Pd[i] = design.riccati[i * n] * model_.c;

    const std::size_t horizon = reference.size() - 1;
    feedforward_.resize(horizon);
    std::vector<double> v(n, 0.0), w(n), next(n);
    v[0] = stateWeight * reference[horizon];
    const double* firstRow = system.A.data();
    for (std::size_t k = horizon; k-- > 0;) {
        for (std::size_t i = 0; i < n; ++i) w[i] = v[i] - Pd[i];
        const double Bw = model_.b[0] * w[0] + (nu > 1 ? w[ny] : 0.0);
        feedforward_[k] = Bw / S;
        // (A - B K)^T w = A^T w - K^T (B^T w), with A = dense first row + two shift chains.
        for (std::size_t j = 0; j < n; ++j) {
            const bool shifted = j + 1 < n && j + 1 != ny;
            next[j] = firstRow[j] * w[0] + (shifted ? w[j + 1] : 0.0) - gain_[j] * Bw;
        }
        next[0] += stateWeight * reference[k];
        v.swap(next);
    }
}

IoDelayTrackerStream::IoDelayTrackerStream(std::shared_ptr<const IoDelayTracker> tracker, double warmupInput,
                                           std::optional<std::pair<double, double>> limits)
    : tracker_(std::move(tracker)), limits_(limits), warmup_(warmupInput) {
    if (!tracker_ || !std::isfinite(warmupInput))
        throw std::invalid_argument("io_delay_lqr stream needs a tracker and a finite warmup input");
    const auto& m = tracker_->model();
    z_.assign(m.stateDimension(), 0.0);
    for (std::size_t j = m.ny; j < z_.size(); ++j) z_[j] = warmup_;
    required_ = std::max(m.ny, m.nu);
}

double IoDelayTrackerStream::step(std::size_t tick, double measurement) {
    if (started_ && tick == lastTick_) return last_;
    if (started_ && tick != lastTick_ + 1)
        throw std::runtime_error("io_delay_lqr needs consecutive controller ticks");
    if (!std::isfinite(measurement)) throw std::runtime_error("io_delay_lqr measurement is not finite");
    const auto& m = tracker_->model();
    std::copy_backward(z_.begin(), z_.begin() + static_cast<std::ptrdiff_t>(m.ny - 1),
                       z_.begin() + static_cast<std::ptrdiff_t>(m.ny));
    z_[0] = measurement;
    ++samples_;
    double command = warmup_;
    if (samples_ >= required_) {
        const auto& feedforward = tracker_->feedforward();
        if (tick >= feedforward.size())
            throw std::runtime_error("io_delay_lqr tick is past the reference horizon");
        command = feedforward[tick];
        const auto& gain = tracker_->gain();
        for (std::size_t i = 0; i < z_.size(); ++i) command -= gain[i] * z_[i];
        if (!std::isfinite(command)) throw std::runtime_error("io_delay_lqr command is not finite");
    }
    if (limits_) command = std::clamp(command, limits_->first, limits_->second);
    if (m.nu > 1) {  // the history stores the command actually applied, after the cap
        std::copy_backward(z_.begin() + static_cast<std::ptrdiff_t>(m.ny),
                           z_.end() - 1, z_.end());
        z_[m.ny] = command;
    }
    started_ = true; lastTick_ = tick; last_ = command;
    return command;
}

}  // namespace dynsyslab::core
