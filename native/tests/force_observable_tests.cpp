#include <dynsyslab/core/System.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace dynsyslab::core;

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, double tolerance = 1e-11) {
    check(std::isfinite(actual) && std::isfinite(expected)
        && std::abs(actual - expected) <= tolerance * std::max({1., std::abs(actual), std::abs(expected)}),
        "force/torque value mismatch");
}

void checkExpression(const ODESystem& system, const char* name, const double* state) {
    const auto observable = system.observableCatalog().at(name);
    check(observable.hasAnalyticJacobian() && observable.expressions().size() == 1,
          "force/torque must use standard analytic observables");
    const auto restored = Observable::fromExpressionDefinition(
        json::parse(json::dump(observable.expressionDefinition())));
    close(restored.evaluateReal(state)[0], observable.evaluateReal(state)[0]);
    check(restored.descriptor().units == observable.descriptor().units
        && restored.descriptor().stateCoordinateSystem == observable.descriptor().stateCoordinateSystem,
        "observable serialization lost unit/coordinate metadata");
    const auto gradient = observable.jacobian(state);
    const auto numerical = observable.jacobian(state, ObservableJacobianBackend::finiteDifference);
    for (std::size_t axis = 0; axis < 2; ++axis) {
        close(gradient[axis].real(), numerical[axis].real(), 1e-7);
        close(observable.partialDerivative(axis).evaluateReal(state)[0], gradient[axis].real());
    }
    const auto drift = system.driftExpressions();
    const auto inputs = system.inputExpressions();
    check(drift.size() == 2 && inputs.size() == 2, "missing analytic system fields");
    for (double input : {0., -.6, .9}) {
        double velocity[2];
        system.rhsWithInput(.3, state, input, velocity);
        const auto matrix = system.inputMatrix(state);
        for (std::size_t axis = 0; axis < 2; ++axis) {
            close(inputs[axis].evaluate(state), matrix[axis]);
            close(drift[axis].evaluate(state) + input * inputs[axis].evaluate(state), velocity[axis]);
        }
        close(observable.inputFlowDerivative(drift, inputs, {input}).evaluateReal(state)[0],
              gradient[0].real() * velocity[0] + gradient[1].real() * velocity[1]);
    }
    const auto rate = observable.flowDerivative(drift);
    const auto rateGradient = rate.jacobian(state);
    close(observable.flowDerivative(drift, 2).evaluateReal(state)[0],
          rateGradient[0].real() * drift[0].evaluate(state)
          + rateGradient[1].real() * drift[1].evaluate(state));
    const auto coordinates = system.observableCatalog().coordinates({name, "velocity"});
    close(coordinates.evaluateReal(state)[0], observable.evaluateReal(state)[0]);
    check(!coordinates.expressions().empty(), "stacked coordinates lost expression provenance");
}

void duffing() {
    for (const auto& descriptor : systemCatalog()) {
        if (descriptor.identifier != "duffing") continue;
        for (const auto& preset : descriptor.presets) {
            SystemConfiguration config;
            config.identifier = "duffing"; config.preset = preset.key;
            config.driver = {.8, .7, .2};
            const auto system = createSystem(config);
            const auto& parameters = system->configuration().parameters;
            const auto force = system->observableCatalog().at("force");
            check(force.descriptor().units == "dimensionless"
                && force.descriptor().stateCoordinateSystem == "dimensionless",
                "Duffing normalized force mislabeled as physical");
            for (const auto& state : std::vector<std::vector<double>>{{0,0},{.5,-.25},{-1.7,.9},{2.1,-1.2}}) {
                double derivative[2];
                system->rhsWithInput(.3, state.data(), 0., derivative);
                close(force.evaluateReal(state.data())[0], derivative[1]);
                close(force.evaluateReal(state.data())[0], parameters.at("alpha") * state[0]
                    - parameters.at("beta") * std::pow(state[0], 3) - parameters.at("delta") * state[1]);
                close(system->observableCatalog().at("force_dimensionless").evaluateReal(state.data())[0], derivative[1]);
                system->rhs(.3, state.data(), derivative);
                close(derivative[1], force.evaluateReal(state.data())[0]
                    + parameters.at("input_gain") * config.driver.value(.3));
                checkExpression(*system, "force", state.data());
            }
        }
    }
    // A physical calibration must recover m*d2X/dt2, not label x'' as newtons.
    const double mass = 2.3, lengthScale = .17, timeScale = .43, damping = .31;
    const double linear = 4.2, cubic = 28., forceScale = mass * lengthScale / (timeScale * timeScale);
    SystemConfiguration config; config.identifier = "duffing";
    config.parameters = {{"alpha", linear * timeScale * timeScale / mass},
        {"beta", cubic * lengthScale * lengthScale * timeScale * timeScale / mass},
        {"delta", damping * timeScale / mass}, {"input_gain", .4}};
    const auto system = createSystem(config);
    const double physicalPosition = -.23, physicalVelocity = .71;
    const double state[]{physicalPosition / lengthScale, physicalVelocity * timeScale / lengthScale};
    close(forceScale * system->observableCatalog().at("force").evaluateReal(state)[0],
          linear * physicalPosition - cubic * std::pow(physicalPosition, 3) - damping * physicalVelocity);
    checkExpression(*system, "force", state);
}

void pendulum() {
    SystemConfiguration config; config.identifier = "nonlinear_pendulum";
    config.parameters = {{"mass_kg", 2.3}, {"length_m", .7},
                         {"gravity_m_s2", 3.7}, {"damping_n_m_s", .19}};
    config.driver = {.8, .7, .2};
    for (double damping : {0., .19}) {
        config.parameters["damping_n_m_s"] = damping;
        const auto system = createSystem(config);
        const auto catalog = system->observableCatalog();
        const auto torque = catalog.at("torque"), normalized = catalog.at("torque_dimensionless");
        const double inertia = 2.3 * .7 * .7, scale = 2.3 * 3.7 * .7, T = std::sqrt(.7 / 3.7);
        check(system->configuration().parameterUnits.at("damping_n_m_s") == "N m s",
              "physical damping coefficient has incorrect units");
        check(torque.descriptor().units == "N m" && normalized.descriptor().units == "dimensionless",
              "torque units incorrect");
        check(normalized.descriptor().stateCoordinateSystem == "dimensional",
              "normalizing an output must not relabel its input coordinates");
        close(system->configuration().dimensionlessGroups.at("torque_scale_n_m"), scale);
        close(system->configuration().dimensionlessGroups.at("characteristic_time_s"), T);
        close(system->coordinates().secondsPerModelTime, 1.); // Native solver still advances seconds.
        for (const auto& state : std::vector<std::vector<double>>{{0,0},{.5,-.25},{-2.4,.9},{3.7,-1.2}}) {
            double derivative[2];
            system->rhsWithInput(.3, state.data(), 0., derivative);
            const double value = torque.evaluateReal(state.data())[0];
            close(value, inertia * derivative[1]);
            close(value, -scale * std::sin(state[0]) - damping * state[1]);
            close(normalized.evaluateReal(state.data())[0], value / scale);
            const double w = T * state[1]; // dtheta/d(t/T), not physical angular velocity.
            close(normalized.evaluateReal(state.data())[0], -std::sin(state[0]) - damping * T / inertia * w);
            for (double input : {-.6, .9}) {
                system->rhsWithInput(.3, state.data(), input, derivative);
                close(value + input, inertia * derivative[1]);
                close((value + input) / scale, T * T * derivative[1]);
            }
            system->rhs(.3, state.data(), derivative);
            close(value + config.driver.value(.3), inertia * derivative[1]);
            checkExpression(*system, "torque", state.data());
            checkExpression(*system, "torque_dimensionless", state.data());
        }
    }
    // Zero gravity has no mgL normalization, but its physical torque remains valid.
    config.parameters["gravity_m_s2"] = 0.;
    const auto zeroGravity = createSystem(config);
    check(zeroGravity->observableCatalog().contains("torque"), "zero-gravity torque missing");
    check(!zeroGravity->observableCatalog().contains("torque_dimensionless"), "zero normalization accepted");
    const double state[]{.2, .4};
    close(zeroGravity->observableCatalog().at("torque").evaluateReal(state)[0], -.19 * .4);
}
} // namespace

int main() {
    try {
        duffing(); pendulum();
        std::cout << "force observables: RHS, scaling, derivatives, input separation and serialization passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
