#pragma once
// Several acoustic receivers of one source, synchronized on the nearest receiver's emission clock
// (docs/architecture/RECEIVER_ARRAYS.md).
//
// A receiver at range D records at time t the pressure emitted at t - D/c (retarded time;
// Lighthill 1952, Proc. R. Soc. A 211:564-587; Ffowcs Williams and Hawkings 1969, Phil. Trans.
// R. Soc. A 264:321-342). On the nearest receiver's emission clock tau, receiver i therefore holds
// the pressure emitted at tau - (D_i - D_near)/(c T), T seconds per model time unit; each term of
// the receiver catalog (core/AcousticReceiver.hpp) is evaluated at the two bracketing source
// samples and interpolated linearly. The record starts once every path has history.
//
// The linearized delay inverse estimates the state deviation from one sample of n receivers:
// near an equilibrium x0 with A = Df(x0) and feature gradients H_i, receiver i observes
// H_i exp(-A delta_i) (x - x0), so the stacked rows O give x - x0 = O^+ y (an observability
// matrix of delayed outputs; Kalman 1960, Bol. Soc. Mat. Mexicana 5:102-119).
#include <dynsyslab/core/AcousticReceiver.hpp>
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct ReceiverArraySignals {
    std::vector<double> time;            // samples: the nearest receiver's emission clock (model time)
    std::vector<double> states;          // samples x state: the source state at that clock
    std::vector<std::string> terms;      // pressure terms of the receiver catalog
    std::vector<double> pressures;       // terms x samples x receivers, Pa
    std::vector<double> nearAtEmission;  // terms x samples x receivers: the near receiver at each column's emission
    std::vector<double> ranges;          // receivers, m
    std::vector<double> lags;            // receivers: (D_i - D_near) / (c T), model time
    std::size_t receivers{}, stateDimension{}, nearest{};
    [[nodiscard]] std::size_t samples() const noexcept { return time.size(); }
    [[nodiscard]] double pressure(std::size_t term, std::size_t sample, std::size_t receiver) const {
        return pressures[(term * samples() + sample) * receivers + receiver];
    }
};

struct ReceiverArrayOptions {
    std::size_t stride{1};               // keep every stride-th source sample of the clock
    double begin{-1e300}, end{1e300};    // emission-time window (model time)
};

// The RPE receiver catalog at every receiver along a recorded trajectory (time, states row-major).
[[nodiscard]] ReceiverArraySignals rpeReceiverArray(const ODESystem& system, const std::vector<AcousticReceiver>& receivers,
                                                    const std::vector<double>& time, const std::vector<double>& states,
                                                    const ReceiverArrayOptions& options = {});

// O^+ for the delayed-output rows H_i exp(-A delta_i): A is n x n, H receivers x n (row-major),
// delays per receiver; the result is n x receivers. Throws when O has rank below n.
[[nodiscard]] std::vector<double> linearDelayInverse(const std::vector<double>& A, const std::vector<double>& H,
                                                     const std::vector<double>& delays, std::size_t n);

}  // namespace dynsyslab::core
