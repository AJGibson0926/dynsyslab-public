#pragma once

// Opt-in log of every controller evaluation of an integration (IntegrationOptions::
// controllerCallLog). A continuous-RHS controller is called at every stage of every trial
// step, including steps the error control later rejects; each such call is tagged with the
// fate of its step once the stepper decides it. Port of the Python solver's controller call
// log (solver/ode.py); the accepted/rejected tagging follows the step-acceptance test of
// embedded Runge-Kutta pairs (Hairer, Norsett & Wanner, Solving ODEs I, 2nd ed., Springer
// 1993, sec. II.4).

#include <cstddef>
#include <cstdint>
#include <vector>

namespace dynsyslab::core {

struct ControllerCallLog {
    enum Kind : std::int8_t { rightHandSide = 0, sampleTick = 1, outputSample = 2 };
    // 1 accepted, 0 rejected, -1 unknown (the stepper does not report trial outcomes, as
    // for CVODE, or the trajectory ended before its step was decided).
    enum Status : std::int8_t { rejected = 0, accepted = 1, unknown = -1, pending = -2 };

    std::size_t stateDimension{};
    std::size_t inputDimension{};
    std::size_t maximumEntries{0};   // 0: unlimited; past it calls are counted, not stored
    std::size_t droppedEntries{};

    std::vector<std::size_t> trajectory;
    std::vector<double> time;
    std::vector<double> states;      // entries x stateDimension
    std::vector<double> inputs;      // entries x inputDimension
    std::vector<std::int8_t> kind;
    std::vector<std::int8_t> status;

    [[nodiscard]] std::size_t size() const noexcept { return time.size(); }
    void record(std::size_t trajectoryIndex, double at, const double* state, const double* input,
                Kind callKind, Status callStatus);
    // Tags every pending entry (calls since the last decision) with the step's outcome.
    void resolvePending(Status outcome);
    // Indices of one trajectory's entries, in call order.
    [[nodiscard]] std::vector<std::size_t> entriesOf(std::size_t trajectoryIndex) const;

private:
    std::size_t pendingBegin_{};
};

}  // namespace dynsyslab::core
