#include <dynsyslab/core/ClassicalControl.hpp>

#include <algorithm>
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

void pidCompensatorOptions() {
    using dynsyslab::core::PIDCompensator;
    dynsyslab::core::PIDOptions options;
    options.integratorLimit = {0.5};
    PIDCompensator clamped({0.0}, {1.0}, {0.0}, 1, 1, options);
    (void)clamped.update(0.0, {1.0});
    require(std::abs(clamped.update(10.0, {1.0}).front() - 0.5) < 1e-12, "PID integrator limit");
    options = {};
    options.derivativeFilterTime = 1.0;
    PIDCompensator filtered({0.0}, {0.0}, {1.0}, 1, 1, options);
    (void)filtered.update(0.0, {0.0});
    const double raw = 1.0 / 0.1, expected = raw * (1.0 - std::exp(-0.1));
    require(std::abs(filtered.update(0.1, {1.0}).front() - expected) < 1e-12, "PID derivative filter");
    options = {};
    options.outputMaximum = 1.0; options.conditionalIntegration = true;
    PIDCompensator windup({0.5}, {1.0}, {0.0}, 1, 1, options);
    (void)windup.update(0.0, {4.0});
    require(windup.update(1.0, {4.0}).front() == 1.0 && windup.integral().front() == 0.0,
            "PID conditional integration holds the integral while saturated");
    require(windup.update(2.0, {-1.0}).front() < 0.0, "PID recovers immediately without windup");
}

void lqiRemovesOffset() {
    // x' = -x + u + d with an unmodelled constant d: LQR leaves an offset, LQI removes it.
    const double d = 0.5, target = 1.0, dt = 1e-3;
    const auto lqr = dynsyslab::core::solveContinuousLQR({-1.0}, {1.0}, {10.0}, {1.0}, 1, 1);
    const auto lqi = dynsyslab::core::solveContinuousLQI({-1.0}, {1.0}, {1.0}, {10.0}, {10.0}, {1.0}, 1, 1, 1);
    require(lqi.gain.size() == 2 && std::all_of(lqi.closedLoopPoles.begin(), lqi.closedLoopPoles.end(),
        [](auto pole) { return pole.real() < 0.0; }), "LQI gain and stable augmented closed loop");
    // Feedforward u_ff = target makes x = target for d = 0 (x' = -x + u).
    double x1 = 0.0, x2 = 0.0, q = 0.0;
    for (int k = 0; k < 40000; ++k) {
        const double u1 = target - lqr.gain[0] * (x1 - target);
        const double u2 = target - lqi.gain[0] * (x2 - target) - lqi.gain[1] * q;
        x1 += dt * (-x1 + u1 + d); x2 += dt * (-x2 + u2 + d); q += dt * (x2 - target);
    }
    require(std::abs(x1 - target) > 0.05, "LQR keeps a disturbance offset");
    require(std::abs(x2 - target) < 1e-6, "LQI removes the disturbance offset");
}

}  // namespace

int main() {
    try {
        fixedLqr();
        stateDependentLqr();
        finiteHorizonMpc();
        sampledPid();
        pidCompensatorOptions();
        lqiRemovesOffset();
        std::cout << "dynsyslab_classical_control_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_classical_control_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
