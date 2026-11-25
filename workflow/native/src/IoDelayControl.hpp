#pragma once
#include <dynsyslab/core/Integrator.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/System.hpp>

#include <optional>
#include <string>
#include <utility>

namespace dynsyslab::workflow::detail {

struct IoDelayControl {
    dynsyslab::core::ScalarControlLaw law;
    dynsyslab::core::StateTargetFunction stateTarget;
    dynsyslab::core::StateTargetFunction measuredTarget;
    std::string observable;
    double sampleTime{};
    dynsyslab::core::json::Value evidence;
};

// Controller kind "io_delay_lqr". The target is lifted to y_ref once, offline, through the
// from_state observable path. Online the law reads one scalar measurement and its own past
// commands; the state argument only simulates the sensor. See docs/architecture/IO_DELAY_CONTROL.md.
[[nodiscard]] IoDelayControl buildIoDelayControl(const dynsyslab::core::json::Value& settings,
    const dynsyslab::core::json::Value& target, const dynsyslab::core::ODESystem& system,
    std::size_t streams, double controlSampleTime,
    std::optional<std::pair<double, double>> limits);

// The model clock, read from the archive, for resolving the sampled controller clock.
[[nodiscard]] double ioDelayModelSampleTime(const dynsyslab::core::json::Value& settings);

}  // namespace dynsyslab::workflow::detail
