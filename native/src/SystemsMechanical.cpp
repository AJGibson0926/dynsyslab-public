#include "SystemSupport.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core::detail {
namespace {

std::vector<double> numericalJacobian(const RHSFunction& rhs, const SinusoidalDriver& driver,
                                      std::size_t dimension, double time,
                                      const double* state) {
    std::vector<double> matrix(dimension * dimension), plus(state, state + dimension),
        minus(plus), fPlus(dimension), fMinus(dimension);
    const double input = driver.value(time);
    for (std::size_t column = 0; column < dimension; ++column) {
        const double step = 1.0e-6 * std::max(1.0, std::abs(state[column]));
        plus[column] += step;
        minus[column] -= step;
        rhs(time, plus.data(), input, fPlus.data());
        rhs(time, minus.data(), input, fMinus.data());
        for (std::size_t row = 0; row < dimension; ++row)
            matrix[row * dimension + column] = (fPlus[row] - fMinus[row]) / (2.0 * step);
        plus[column] = state[column];
        minus[column] = state[column];
    }
    return matrix;
}

std::unique_ptr<ODESystem> nonlinearPendulum(SystemConfiguration configuration) {
    double length = 1.0, mass = 1.0, gravity = 9.80665, damping = 0.02;
    std::vector<double> initial{0.2, 0.0};
    if (configuration.preset == "large_angle_low_damping") damping = 0.003, initial = {2.4, 0.0};
    if (configuration.preset == "driven_damped") damping = 0.08, initial = {0.4, 0.0};
    if (configuration.preset == "holmes_chaotic") damping = 0.5, initial = {-3.141592653589793, 0.0};
    length = resolvedParameter(configuration, "length_m", length);
    mass = resolvedParameter(configuration, "mass_kg", mass);
    gravity = resolvedParameter(configuration, "gravity_m_s2", gravity);
    damping = resolvedParameter(configuration, "damping_n_m_s", damping, "N m s");
    const double inertia = mass * length * length;
    auto rhs = [length, gravity, damping, inertia](double, const double* x, double input,
                                                   double* out) {
        out[0] = x[1];
        out[1] = -gravity / length * std::sin(x[0]) - damping / inertia * x[1] + input / inertia;
    };
    auto jacobian = [length, gravity, damping, inertia](double, const double* x) {
        return std::vector<double>{0.0, 1.0, -gravity / length * std::cos(x[0]),
                                   -damping / inertia};
    };
    auto input = [inertia](const double*) { return std::vector<double>{0.0, 1.0 / inertia}; };
    const auto angle = AnalyticExpression::state(2, 0), omega = AnalyticExpression::state(2, 1);
    // Intrinsic torque = I * zero-input angular acceleration; not acceleration itself.
    const auto torque = -mass * gravity * length * angle.unary("sin") - damping * omega;
    std::vector<Observable> observables{Observable::fromExpressions(
        {"torque", "Intrinsic torque", "N m", "dimensional", {}, false}, {torque})};
    const double torqueScale = mass * gravity * length;
    if (mass > 0.0 && length > 0.0 && gravity > 0.0 && torqueScale > 0.0 && std::isfinite(torqueScale)) {
        // Normalize the output, not the input coordinates (theta, dtheta/dt).
        configuration.dimensionlessGroups["torque_scale_n_m"] = torqueScale;
        configuration.dimensionlessGroupUnits["torque_scale_n_m"] = "N m";
        configuration.dimensionlessGroups["characteristic_time_s"] = std::sqrt(length / gravity);
        configuration.dimensionlessGroupUnits["characteristic_time_s"] = "s";
        observables.push_back(Observable::fromExpressions(
            {"torque_dimensionless", "Normalized torque", "dimensionless", "dimensional", {}, false},
            {torque / torqueScale}));
    }
    return std::make_unique<FunctionalSystem>(descriptorFor("nonlinear_pendulum"), configuration,
        selectedInitial(configuration, std::move(initial)), rhs, jacobian, input,
        [](double, const double* x) { return genericInvalid(x, 2); }, std::vector<double>{},
        std::move(observables), std::map<std::string,std::string>{},
        std::map<std::string,std::vector<std::string>>{},
        std::vector<AnalyticExpression>{omega, torque / inertia},
        std::vector<AnalyticExpression>{AnalyticExpression::constant(0., 2),
                                      AnalyticExpression::constant(1.0 / inertia, 2)});
}

std::unique_ptr<ODESystem> doublePendulum(SystemConfiguration configuration) {
    double m1 = 1.0, m2 = 1.0, l1 = 1.0, l2 = 1.0, gravity = 9.80665;
    double damping1 = 0.01, damping2 = 0.01;
    std::vector<double> initial{0.5, 0.0, 0.7, 0.0};
    if (configuration.preset == "classic_chaotic_release")
        damping1 = damping2 = 0.001, initial = {2.1, 0.0, 1.2, 0.0};
    if (configuration.preset == "unequal_lab_rig") {
        m1 = 0.8; m2 = 1.2; l1 = 0.75; l2 = 1.1; damping1 = 0.015; damping2 = 0.02;
        initial = {0.9, 0.0, -0.6, 0.2};
    }
    m1 = resolvedParameter(configuration, "mass1_kg", m1);
    m2 = resolvedParameter(configuration, "mass2_kg", m2);
    l1 = resolvedParameter(configuration, "length1_m", l1);
    l2 = resolvedParameter(configuration, "length2_m", l2);
    gravity = resolvedParameter(configuration, "gravity_m_s2", gravity);
    damping1 = resolvedParameter(configuration, "damping1_n_m_s", damping1);
    damping2 = resolvedParameter(configuration, "damping2_n_m_s", damping2);
    auto rhs = [m1, m2, l1, l2, gravity, damping1, damping2]
        (double, const double* x, double input, double* out) {
        const double a = x[0], w1 = x[1], b = x[2], w2 = x[3], delta = a - b;
        const double denominator = 2.0 * m1 + m2 - m2 * std::cos(2.0 * delta);
        double dd1 = (-gravity * (2.0 * m1 + m2) * std::sin(a)
            - m2 * gravity * std::sin(a - 2.0 * b)
            - 2.0 * std::sin(delta) * m2
                * (w2 * w2 * l2 + w1 * w1 * l1 * std::cos(delta))) / (l1 * denominator);
        double dd2 = (2.0 * std::sin(delta)
            * (w1 * w1 * l1 * (m1 + m2) + gravity * (m1 + m2) * std::cos(a)
                + w2 * w2 * l2 * m2 * std::cos(delta))) / (l2 * denominator);
        dd1 += -damping1 * w1 / (m1 * l1 * l1) + input / (m1 * l1 * l1);
        dd2 -= damping2 * w2 / (m2 * l2 * l2);
        out[0] = w1; out[1] = dd1; out[2] = w2; out[3] = dd2;
    };
    const auto driver = configuration.driver;
    auto jacobian = [rhs, driver](double time, const double* x) {
        return numericalJacobian(rhs, driver, 4, time, x);
    };
    auto input = [m1, l1](const double*) {
        return std::vector<double>{0.0, 1.0 / (m1 * l1 * l1), 0.0, 0.0};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("double_pendulum"), configuration,
        selectedInitial(configuration, std::move(initial)), rhs, jacobian, input,
        [](double, const double* x) { return genericInvalid(x, 4); });
}

std::unique_ptr<ODESystem> springPendulums(SystemConfiguration configuration) {
    double m1 = 1.0, m2 = 1.0, l1 = 1.0, l2 = 1.0, pivot = 1.0, spring = 1.0,
        rest = 1.0, gravity = 9.80665, damping1 = 0.01, damping2 = 0.01;
    std::vector<double> initial{0.45, 0.0, -0.2, 0.0};
    if (configuration.preset == "symmetric_weak_coupling") spring = 0.5;
    if (configuration.preset == "strong_nonlinear_coupling")
        spring = 20.0, damping1 = damping2 = 0.005, initial = {1.1, 0.0, -0.9, 0.0};
    if (configuration.preset == "detuned_energy_exchange") {
        m1 = 0.8; m2 = 1.2; l1 = 0.85; l2 = 1.15; spring = 3.0;
        damping1 = 0.01; damping2 = 0.015; initial = {0.7, 0.0, 0.05, 0.0};
    }
    m1 = resolvedParameter(configuration, "mass1_kg", m1); m2 = resolvedParameter(configuration, "mass2_kg", m2);
    l1 = resolvedParameter(configuration, "length1_m", l1); l2 = resolvedParameter(configuration, "length2_m", l2);
    pivot = resolvedParameter(configuration, "pivot_separation_m", pivot);
    spring = resolvedParameter(configuration, "spring_constant_n_m", spring);
    rest = resolvedParameter(configuration, "spring_rest_length_m", rest);
    gravity = resolvedParameter(configuration, "gravity_m_s2", gravity);
    damping1 = resolvedParameter(configuration, "damping1_n_m_s", damping1);
    damping2 = resolvedParameter(configuration, "damping2_n_m_s", damping2);
    auto rhs = [m1, m2, l1, l2, pivot, spring, rest, gravity, damping1, damping2]
        (double, const double* x, double input, double* out) {
        const double a = x[0], w1 = x[1], b = x[2], w2 = x[3];
        const double dx = pivot + l2 * std::sin(b) - l1 * std::sin(a);
        const double dy = l1 * std::cos(a) - l2 * std::cos(b);
        const double distance = std::hypot(dx, dy);
        if (distance < 1.0e-12) throw std::runtime_error("spring bobs coincide");
        const double scale = spring * (distance - rest) / distance;
        const double grad1 = scale * (dx * -l1 * std::cos(a) + dy * -l1 * std::sin(a));
        const double grad2 = scale * (dx * l2 * std::cos(b) + dy * l2 * std::sin(b));
        out[0] = w1;
        out[1] = -gravity / l1 * std::sin(a) - (damping1 * w1 + grad1) / (m1 * l1 * l1)
            + input / (m1 * l1 * l1);
        out[2] = w2;
        out[3] = -gravity / l2 * std::sin(b) - (damping2 * w2 + grad2) / (m2 * l2 * l2);
    };
    const auto driver = configuration.driver;
    auto jacobian = [rhs, driver](double time, const double* x) {
        return numericalJacobian(rhs, driver, 4, time, x);
    };
    auto input = [m1, l1](const double*) {
        return std::vector<double>{0.0, 1.0 / (m1 * l1 * l1), 0.0, 0.0};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("spring_coupled_pendulums"),
        configuration, selectedInitial(configuration, std::move(initial)), rhs, jacobian, input,
        [](double, const double* x) { return genericInvalid(x, 4); });
}

std::unique_ptr<ODESystem> fput(SystemConfiguration configuration) {
    const bool preset32 = configuration.preset == "alpha_fput_32"
        || configuration.preset == "beta_fput_32";
    const std::size_t sites = resolvedIntegerParameter(configuration, "sites", preset32 ? 32 : 16);
    const double alpha = resolvedParameter(configuration, "alpha",
        configuration.preset == "beta_fput_32" ? 0.0 : 0.25);
    const double beta = resolvedParameter(configuration, "beta",
        configuration.preset == "alpha_fput_32" ? 0.0 : 0.25);
    auto rhs = [sites, alpha, beta](double, const double* x, double input, double* out) {
        std::copy(x + sites, x + 2 * sites, out);
        auto vPrime = [alpha, beta](double r) { return r + alpha*r*r + beta*r*r*r; };
        for (std::size_t index = 0; index < sites; ++index) {
            const double left = x[index] - (index == 0 ? 0.0 : x[index - 1]);
            const double right = (index + 1 == sites ? 0.0 : x[index + 1]) - x[index];
            out[sites + index] = vPrime(right) - vPrime(left);
        }
        out[sites] += input;
    };
    auto jacobian = [sites, alpha, beta](double, const double* x) {
        std::vector<double> matrix(4 * sites * sites, 0.0);
        for (std::size_t index = 0; index < sites; ++index)
            matrix[index * (2 * sites) + sites + index] = 1.0;
        auto vSecond = [alpha, beta](double r) { return 1.0 + 2.0*alpha*r + 3.0*beta*r*r; };
        for (std::size_t index = 0; index < sites; ++index) {
            const double left = x[index] - (index == 0 ? 0.0 : x[index - 1]);
            const double right = (index + 1 == sites ? 0.0 : x[index + 1]) - x[index];
            const double kLeft = vSecond(left), kRight = vSecond(right);
            const std::size_t row = sites + index;
            matrix[row * (2 * sites) + index] = -(kLeft + kRight);
            if (index) matrix[row * (2 * sites) + index - 1] = kLeft;
            if (index + 1 < sites) matrix[row * (2 * sites) + index + 1] = kRight;
        }
        return matrix;
    };
    auto input = [sites](const double*) {
        std::vector<double> matrix(2 * sites, 0.0); matrix[sites] = 1.0; return matrix;
    };
    std::vector<double> initial(2 * sites, 0.0); initial[0] = 0.1;
    return std::make_unique<FunctionalSystem>(descriptorFor("fput", 2 * sites), configuration,
        selectedInitial(configuration, std::move(initial)), rhs, jacobian, input,
        [sites](double, const double* x) { return genericInvalid(x, 2 * sites); });
}

std::unique_ptr<ODESystem> kuramoto(SystemConfiguration configuration) {
    const bool preset32 = configuration.preset == "weak_coupling"
        || configuration.preset == "synchronizing";
    const std::size_t count = resolvedIntegerParameter(configuration, "oscillators", preset32 ? 32 : 16);
    const double coupling = resolvedParameter(configuration, "coupling",
        configuration.preset == "weak_coupling" ? 0.1
            : configuration.preset == "synchronizing" ? 5.0 : 1.0);
    const double spread = resolvedParameter(configuration, "frequency_spread", 0.5);
    auto rhs = [count, coupling, spread](double, const double* x, double input, double* out) {
        for (std::size_t row = 0; row < count; ++row) {
            double interaction = 0.0;
            for (std::size_t column = 0; column < count; ++column)
                interaction += std::sin(x[column] - x[row]);
            const double fraction = count == 1 ? 0.0 : static_cast<double>(row) / (count - 1);
            out[row] = -spread + 2.0 * spread * fraction
                + coupling * interaction / static_cast<double>(count) + input;
        }
    };
    auto jacobian = [count, coupling](double, const double* x) {
        std::vector<double> matrix(count * count, 0.0);
        for (std::size_t row = 0; row < count; ++row) {
            double sum = 0.0;
            for (std::size_t column = 0; column < count; ++column) {
                if (row == column) continue;
                const double value = coupling * std::cos(x[column] - x[row]) / count;
                matrix[row * count + column] = value; sum += value;
            }
            matrix[row * count + row] = -sum;
        }
        return matrix;
    };
    auto input = [count](const double*) { return std::vector<double>(count, 1.0); };
    std::vector<double> initial(count);
    constexpr double twoPi = 6.283185307179586476925286766559;
    for (std::size_t index = 0; index < count; ++index)
        initial[index] = twoPi * static_cast<double>(index) / static_cast<double>(count);
    return std::make_unique<FunctionalSystem>(descriptorFor("kuramoto", count), configuration,
        selectedInitial(configuration, std::move(initial)), rhs, jacobian, input,
        [count](double, const double* x) { return genericInvalid(x, count, 1.0e12); });
}

}  // namespace

std::unique_ptr<ODESystem> makeMechanicalSystem(const SystemConfiguration& configuration) {
    if (configuration.identifier == "nonlinear_pendulum") return nonlinearPendulum(configuration);
    if (configuration.identifier == "double_pendulum") return doublePendulum(configuration);
    if (configuration.identifier == "spring_coupled_pendulums") return springPendulums(configuration);
    if (configuration.identifier == "fput") return fput(configuration);
    if (configuration.identifier == "kuramoto") return kuramoto(configuration);
    return {};
}

}  // namespace dynsyslab::core::detail
