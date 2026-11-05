#include <dynsyslab/core/ClassicalControl.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void fixedLqr() {
    dynsyslab::core::StateSpaceLQRController controller(
        {-1.0}, {1.0}, {1.0}, {1.0}, 1, 1);
    const double state[] = {1.0};
    const auto value = controller.value(0.0, state);
    require(!controller.stateDependent(), "fixed LQR identity");
    require(std::abs(controller.lastSolution().gain.front()
        - (std::sqrt(2.0) - 1.0)) < 1.0e-12, "fixed LQR gain");
    require(std::abs(value.front() + std::sqrt(2.0) - 1.0) < 1.0e-12,
            "fixed LQR returns system input");
}

void stateDependentLqr() {
    std::size_t generatorCalls = 0, inputCalls = 0;
    dynsyslab::core::StateSpaceLQRController controller(
        [&generatorCalls](const double* state) {
            ++generatorCalls; return std::vector<double>{-1.0 - 0.1 * state[0]};
        },
        [&inputCalls](const double* state) {
            ++inputCalls; return std::vector<double>{1.0 + 0.05 * state[0]};
        }, {1.0}, {1.0}, 1, 1,
        [](double) { return std::vector<double>{0.5}; });
    const double first[] = {1.0}, second[] = {2.0};
    require(controller.value(0.0, first).front() < 0.0,
            "state-dependent LQR restoring input");
    (void)controller.value(0.1, second);
    require(controller.stateDependent() && generatorCalls == 2 && inputCalls == 2,
            "state-dependent LQR refreshes A and B");
}

void finiteHorizonMpc() {
    const auto solution = dynsyslab::core::solveFiniteHorizonLQR(
        {0.0}, {1.0}, {1.0}, {1.0}, 1, 1, 5, 0.1);
    require(solution.gain.front() > 0.0 && solution.cost.front() > 1.0,
            "finite-horizon Riccati recursion");
    dynsyslab::core::StateSpaceMPCController controller(
        [](const double*) { return std::vector<double>{0.0}; },
        [](const double*) { return std::vector<double>{1.0}; },
        {1.0}, {1.0}, 1, 1, {}, 5, 0.1, {0.2});
    const double state[] = {10.0};
    require(std::abs(controller.value(0.0, state).front() + 0.2) < 1.0e-14,
            "finite-horizon MPC input limit");
}

void sampledPid() {
    dynsyslab::core::PIDController controller(
        {1.0}, {2.0}, {3.0}, 1, 1,
        [](double) { return std::vector<double>{1.0}; });
    const double zero[] = {0.0}, half[] = {0.5};
    require(controller.value(0.0, zero).front() == 1.0,
            "PID first sample is proportional only");
    require(controller.value(1.0, zero).front() == 3.0,
            "PID integral accumulation");
    require(controller.value(2.0, half).front() == 2.0,
            "PID derivative and integral terms");
    bool rejected = false;
    try { (void)controller.value(1.5, half); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "PID rejects a reversing sample clock");
    controller.reset();
    require(controller.value(0.0, zero).front() == 1.0, "PID reset");
}

}  // namespace

int main() {
    try {
        fixedLqr();
        stateDependentLqr();
        finiteHorizonMpc();
        sampledPid();
        std::cout << "dynsyslab_classical_control_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_classical_control_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
