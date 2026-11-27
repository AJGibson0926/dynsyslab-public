#pragma once
// Checks and repairs of recorded trajectories (docs/architecture/TRAJECTORY_KERNELS.md).
//
//   backwardHistory  states before a record's first sample, by integrating the unforced field
//                    backward in time (classical RK4, steps <= maximumStep), then forward again to
//                    verify the return to the initial state; causal history coordinates of a
//                    simulated record need such a prehistory to start at t0. Simulation only:
//                    measured data cannot be extended this way.
//   compareClocks    a sampled record against a dense reference at the reference's samples in
//                    their common window, by linear interpolation of the sampled record (no
//                    extrapolation): max |error|, RMS and the worst time per component.
#include <dynsyslab/core/System.hpp>

#include <cstddef>
#include <vector>

namespace dynsyslab::core {

struct Prehistory {
    std::vector<double> states;  // times x state, oldest first
    double roundtripError{};     // max |x(t0) after backward and forward integration - x0|
};

// `times` increasing and strictly before t0.
[[nodiscard]] Prehistory backwardHistory(const ODESystem& system, const std::vector<double>& x0, double t0,
                                         const std::vector<double>& times, double maximumStep);

struct ClockComparison {
    double first{}, last{};               // compared window
    std::size_t referenceCount{}, sampleCount{};
    std::vector<double> maximumAbsolute, rms, worstTime;  // per component
};

[[nodiscard]] ClockComparison compareClocks(const std::vector<double>& referenceTime, const std::vector<double>& reference,
                                            const std::vector<double>& sampleTime, const std::vector<double>& samples,
                                            std::size_t width, double windowBegin = -1e300, double windowEnd = 1e300);

}  // namespace dynsyslab::core
