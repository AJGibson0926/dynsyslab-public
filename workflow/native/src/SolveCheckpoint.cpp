#include "SolveCheckpoint.hpp"

#include <dynsyslab/core/LongPath.hpp>

#include <filesystem>

namespace dynsyslab::workflow::detail {
namespace {
// A pause file is polled at most this often between due checkpoints (a stat per sample
// would dominate small integrations).
constexpr double pausePollSeconds = 0.2;

double seconds(std::chrono::steady_clock::duration value) {
    return std::chrono::duration<double>(value).count();
}
}  // namespace

CheckpointClock::CheckpointClock(const CheckpointSettings& settings)
    : settings_(settings), gated_(settings.wallInterval > 0.0 || settings.deadlineUnix > 0.0),
      lastWrite_(std::chrono::steady_clock::now()), lastPauseCheck_(lastWrite_) {}

bool CheckpointClock::pauseFileExists() const {
    return !settings_.pauseFile.empty()
        && std::filesystem::exists(core::extendedLengthPath(settings_.pauseFile));
}

CheckpointDecision CheckpointClock::at(std::size_t trajectory, double time) {
    CheckpointDecision decision;
    if (!gated_) {  // the integrator already applied the model-time interval
        decision.write = true;
        decision.pause = pauseFileExists();
        return decision;
    }
    const auto now = std::chrono::steady_clock::now();
    if (trajectory != trajectory_) { trajectory_ = trajectory; lastModelTime_ = time; }
    const double interval = settings_.interval;
    if (interval > 0.0 && time - lastModelTime_ >= interval * (1.0 - 1.0e-12)) decision.write = true;
    if (settings_.wallInterval > 0.0 && seconds(now - lastWrite_) >= settings_.wallInterval)
        decision.write = true;
    if (settings_.deadlineUnix > 0.0) {
        const double nowUnix = std::chrono::duration<double>(
            std::chrono::system_clock::now().time_since_epoch()).count();
        if (nowUnix >= settings_.deadlineUnix) {
            decision.write = decision.pause = decision.deadline = true;
            deadlineReached_ = true;
        }
    }
    if (!decision.pause && (decision.write || seconds(now - lastPauseCheck_) >= pausePollSeconds)) {
        lastPauseCheck_ = now;
        if (pauseFileExists()) decision.write = decision.pause = true;
    }
    if (decision.write) { lastWrite_ = now; lastModelTime_ = time; }
    return decision;
}

}  // namespace dynsyslab::workflow::detail
