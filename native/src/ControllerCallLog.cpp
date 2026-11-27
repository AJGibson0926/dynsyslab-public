#include <dynsyslab/core/ControllerCallLog.hpp>

namespace dynsyslab::core {

void ControllerCallLog::record(std::size_t trajectoryIndex, double at, const double* state,
                               const double* input, Kind callKind, Status callStatus) {
    if (maximumEntries && size() >= maximumEntries) { ++droppedEntries; return; }
    trajectory.push_back(trajectoryIndex);
    time.push_back(at);
    states.insert(states.end(), state, state + stateDimension);
    inputs.insert(inputs.end(), input, input + inputDimension);
    kind.push_back(callKind);
    status.push_back(callStatus);
}

void ControllerCallLog::resolvePending(Status outcome) {
    for (std::size_t index = pendingBegin_; index < status.size(); ++index)
        if (status[index] == pending) status[index] = outcome;
    pendingBegin_ = status.size();
}

std::vector<std::size_t> ControllerCallLog::entriesOf(std::size_t trajectoryIndex) const {
    std::vector<std::size_t> result;
    for (std::size_t index = 0; index < trajectory.size(); ++index)
        if (trajectory[index] == trajectoryIndex) result.push_back(index);
    return result;
}

}  // namespace dynsyslab::core
