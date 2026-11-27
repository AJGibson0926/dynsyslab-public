// Multi-input classical control and user state-dependent-coefficient factorisations.
//
// Multi-input: lqr, mpc and pid with an explicit n x m input matrix B (m > 1, given as
// settings.B or settings.input_matrix) act through additive-state actuation xdot = f(x) + B u.
// R is m x m (scalar, diagonal or full); control_limit is a scalar, an m-vector of symmetric
// bounds, or [[low, high], ...] per input.
//
// SDRE with a user factorisation (Cloutier 1997, "State-dependent Riccati equation
// techniques: an overview", Proc. American Control Conference, 932-936): settings.sdc_A is
// an n x n row-major list of expressions in the state names (and system parameters) with
// f(x) = A(x) x, and settings.sdc_B an n x m list for B(x); the gain solves the CARE at the
// current state with these factors instead of the Jacobian.
#include "ControlFactoryInternal.hpp"

#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/ClassicalControl.hpp>

#include <algorithm>
#include <map>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace dynsyslab::workflow::control_factory_detail {
namespace {

std::vector<std::pair<double, double>> inputCaps(const json::Value& limit, std::size_t m) {
    std::vector<std::pair<double, double>> caps;
    if (limit.isNull() || (std::holds_alternative<bool>(limit.storage()) && !limit.boolean())) return caps;
    if (!std::holds_alternative<json::Value::Array>(limit.storage())) {
        const auto cap = controlCap(limit);
        caps.assign(m, *cap);
        return caps;
    }
    for (const auto& item : limit.array()) caps.push_back(*controlCap(item));
    if (caps.size() == 1) caps.assign(m, caps.front());
    if (caps.size() != m) throw std::invalid_argument("control_limit needs one bound per input");
    return caps;
}

std::vector<double> window(const json::Value& settings) {
    auto active = numbers(settings.at("active_time_window"));
    if (!active.empty() && (active.size() != 2 || active[0] > active[1]))
        throw std::invalid_argument("active_time_window must be [start, stop]");
    return active;
}

}  // namespace

dynsyslab::core::StateMatrixFunction sdcFactor(const json::Value& entries, std::size_t rows, std::size_t columns,
                         const dynsyslab::core::ODESystem& system, const char* name) {
    const auto& descriptor = system.descriptor();
    std::vector<std::string> variables = descriptor.stateNames;
    if (variables.size() != descriptor.stateDimension) {
        variables.clear();
        for (std::size_t i = 0; i < descriptor.stateDimension; ++i) variables.push_back("x" + std::to_string(i));
    }
    const auto& parameters = system.configuration().parameters;
    std::map<std::string, double> constants(parameters.begin(), parameters.end());
    std::vector<dynsyslab::core::AnalyticExpression> expressions;
    for (const auto& item : entries.array())
        expressions.push_back(std::holds_alternative<std::string>(item.storage())
            ? dynsyslab::core::AnalyticExpression::parse(item.string(), variables, constants)
            : dynsyslab::core::AnalyticExpression::constant(item.number(), descriptor.stateDimension));
    if (expressions.size() != rows * columns)
        throw std::invalid_argument(std::string(name) + " needs " + std::to_string(rows * columns)
                                    + " row-major entries");
    return [expressions = std::move(expressions)](const double* state) {
        std::vector<double> value;
        value.reserve(expressions.size());
        for (const auto& expression : expressions) value.push_back(expression.evaluate(state));
        return value;
    };
}

bool buildMultiInputController(const std::string& kind, const json::Value& settings, const json::Value& solver,
                               const dynsyslab::core::ODESystem& system, std::size_t count,
                               double defaultSampleTime, const dynsyslab::core::StateTargetFunction& target,
                               ControllerSetup& result) {
    if (kind != "lqr" && kind != "mpc" && kind != "pid") return false;
    const auto& matrix = settings.at("input_matrix").isNull() ? settings.at("B") : settings.at("input_matrix");
    const std::size_t n = system.descriptor().stateDimension;
    const auto B = numbers(matrix);
    if (B.empty() || B.size() % n || B.size() / n < 2) return false;
    const std::size_t m = B.size() / n;
    if (settings.at("actuation").string("system_input") != "additive_state")
        throw std::invalid_argument("multi-input control needs actuation additive_state with an n x m B");
    const auto R = squareMatrix(settings.at("R"), m, "R");
    std::vector<double> A = numbers(settings.at("A"));
    auto point = numbers(settings.at("linearization_state"));
    if (point.empty()) point = system.configuration().initialState;
    point.resize(n, 0.0);
    if (A.empty()) A = system.jacobian(0.0, point.data());
    if (A.size() != n * n) throw std::invalid_argument("A must be n x n");
    using Vector = std::function<std::vector<double>(double, const double*)>;
    std::vector<Vector> bank;
    const auto clock = resolveClock(kind, solver, settings, defaultSampleTime, kind != "lqr");
    if (kind == "lqr") {
        const auto Q = squareMatrix(settings.at("Q"), n, "Q");
        for (std::size_t k = 0; k < count; ++k) {
            auto c = std::make_shared<dynsyslab::core::StateSpaceLQRController>(A, B, Q, R, n, m, target);
            bank.push_back([c](double t, const double* x) { return c->value(t, x); });
        }
    } else if (kind == "mpc") {
        const auto Q = squareMatrix(settings.at("Q"), n, "Q");
        const auto horizon = positiveCount(settings.at("horizon"), 20, "MPC horizon");
        for (std::size_t k = 0; k < count; ++k) {
            auto c = std::make_shared<dynsyslab::core::StateSpaceMPCController>(
                [A](const double*) { return A; }, [B](const double*) { return B; }, Q, R, n, m, target,
                horizon, clock.sampleTime);
            bank.push_back([c](double t, const double* x) { return c->value(t, x); });
        }
    } else {
        const auto gains = [&](const char* key, double fallback) {
            auto values = numbers(settings.at(key));
            if (values.empty()) {
                values.assign(m * n, 0.0);
                for (std::size_t i = 0; i < std::min(m, n); ++i) values[i * n + i] = settings.at(key).number(fallback);
            }
            if (values.size() != m * n) throw std::invalid_argument(std::string("PID ") + key + " must be m x n");
            return values;
        };
        const auto Kp = gains("Kp", 1.0), Ki = gains("Ki", 0.0), Kd = gains("Kd", 0.0);
        for (std::size_t k = 0; k < count; ++k) {
            auto c = std::make_shared<dynsyslab::core::PIDController>(Kp, Ki, Kd, n, m, target);
            bank.push_back([c](double t, const double* x) { return c->value(t, x); });
        }
    }
    const auto caps = inputCaps(settings.at("control_limit"), m);
    const auto active = window(settings);
    result.policy.vector = [bank = std::move(bank), caps, active, m](std::size_t trajectory, double time,
                                                                     const double* state, double* input) {
        if (!active.empty() && (time < active[0] || time > active[1])) {
            std::fill(input, input + m, 0.0); return; }
        const auto u = bank.at(trajectory)(time, state);
        for (std::size_t j = 0; j < m; ++j)
            input[j] = caps.empty() ? u.at(j) : std::clamp(u.at(j), caps[j].first, caps[j].second);
    };
    result.policy.inputDimension = m;
    result.policy.additiveMatrix = B;
    result.policy.evaluation = clock.mode;
    result.policy.sampleTime = clock.sampleTime;
    result.evaluation = clock.label;
    result.sampleTime = clock.sampleTime;
    result.kind = kind;
    result.family = "state_space";
    result.target = target;
    result.modelBindingEvidence = json::Value(json::Value::Object{{"input_dimension", json::Value(double(m))},
                                                                  {"actuation", json::Value("additive_state")}});
    return true;
}

}  // namespace dynsyslab::workflow::control_factory_detail
