#pragma once
// When a solve writes a checkpoint (partial archive) and when it pauses: model-time interval
// (history.checkpoint_interval), wall-clock interval (history.checkpoint_wall_seconds), a
// wall-clock deadline (history.deadline_unix_seconds, set by the queue's deadline policy) and
// the pause file (history.pause_file). Without the wall-clock keys the integrator gates calls
// by model time and every call writes, as before.

#include "SolveSettings.hpp"

#include <chrono>
#include <cstddef>

namespace dynsyslab::workflow::detail {

struct CheckpointDecision {
    bool write{};      // write the partial archive for this sample
    bool pause{};      // stop after this sample (NativeSolution::paused)
    bool deadline{};   // the pause is the wall-clock deadline
};

class CheckpointClock {
public:
    explicit CheckpointClock(const CheckpointSettings& settings);
    // Called after each output sample the integrator hands to the checkpoint callback.
    [[nodiscard]] CheckpointDecision at(std::size_t trajectory, double time);
    [[nodiscard]] bool deadlineReached() const noexcept { return deadlineReached_; }

private:
    [[nodiscard]] bool pauseFileExists() const;

    CheckpointSettings settings_;
    bool gated_{};
    std::size_t trajectory_{static_cast<std::size_t>(-1)};
    double lastModelTime_{};
    std::chrono::steady_clock::time_point lastWrite_, lastPauseCheck_;
    bool deadlineReached_{};
};

}  // namespace dynsyslab::workflow::detail
