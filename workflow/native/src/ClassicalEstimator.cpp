// Classical in-loop estimators for solve jobs (estimator.kind lqe | luenberger | static_map),
// over core::LinearObserver (Luenberger 1964; Kalman & Bucy 1961). Job keys:
//   "estimator": {"kind": "lqe", "A": [...], "B": [...] (default: the system's Jacobian and input
//     matrix at "linearization_state"), "C": [...] (p x n; default identity, or the numerical
//     Jacobian of the sensors), "measurement_observables": [...] (default: the state itself),
//     "process_covariance": W, "measurement_covariance": V   (lqe)
//     "L": [...] (luenberger), "M": [...], "offset": [...]   (static_map: x_hat = M y + c)
//     "initial_estimate": [...], "sample_time": dt}
#include "ClassicalEstimator.hpp"

#include <dynsyslab/core/LinearObserver.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>

namespace dynsyslab::workflow::detail {
namespace {

std::vector<double> numbers(const json::Value& v) {
    std::vector<double> result;
    for (const auto& item : v.array()) result.push_back(item.number());
    return result;
}

std::vector<double> square(const json::Value& v, std::size_t n, const char* label) {
    auto values = numbers(v);
    if (values.size() == n * n) return values;
    if (!values.empty() && values.size() != n)
        throw std::invalid_argument(std::string(label) + " must be scalar, diagonal or square");
    std::vector<double> result(n * n, 0.0);
    for (std::size_t k = 0; k < n; ++k) result[k * n + k] = values.empty() ? v.number(1.0) : values[k];
    return result;
}

std::vector<double> sized(const json::Value& v, std::size_t rows, std::size_t columns, const char* label) {
    auto values = numbers(v);
    if (values.size() != rows * columns)
        throw std::invalid_argument(std::string(label) + " needs " + std::to_string(rows * columns) + " row-major values");
    return values;
}

struct Sensors {
    std::vector<core::Observable> observables;
    std::size_t width{};
    std::vector<double> operator()(const double* state) const {
        std::vector<double> y;
        for (const auto& item : observables) {
            const auto values = item.evaluateReal(state);
            y.insert(y.end(), values.begin(), values.end());
        }
        return y;
    }
};

// Central-difference Jacobian of the sensors at x (used when C is not given).
std::vector<double> sensorJacobian(const Sensors& sensors, std::vector<double> x) {
    const std::size_t n = x.size(), p = sensors.width;
    std::vector<double> result(p * n);
    for (std::size_t j = 0; j < n; ++j) {
        const double h = 1e-6 * std::max(1.0, std::abs(x[j]));
        auto up = x, down = x;
        up[j] += h; down[j] -= h;
        const auto a = sensors(up.data()), b = sensors(down.data());
        for (std::size_t k = 0; k < p; ++k) result[k * n + j] = (a[k] - b[k]) / (2 * h);
    }
    return result;
}

}  // namespace

bool isClassicalEstimator(const std::string& kind) {
    return kind == "lqe" || kind == "luenberger" || kind == "static_map";
}

EstimatorSetup buildClassicalEstimator(const json::Value& spec, double sampleTime,
                                       const core::ODESystem& system, std::size_t count) {
    const auto kind = spec.at("kind").string();
    const std::size_t n = system.descriptor().stateDimension;
    auto xe = numbers(spec.at("linearization_state"));
    if (xe.empty()) xe = system.configuration().initialState;
    xe.resize(n, 0.0);
    Sensors sensors;
    const auto catalog = system.observableCatalog();
    for (const auto& name : spec.at("measurement_observables").array()) {
        sensors.observables.push_back(catalog.at(name.string()));
        sensors.width += sensors.observables.back().outputDimension();
    }
    const bool viaObservables = !sensors.observables.empty();
    const std::size_t p = viaObservables ? sensors.width
        : (spec.at("C").isNull() ? n : numbers(spec.at("C")).size() / n);
    std::vector<double> C;
    if (!spec.at("C").isNull()) C = sized(spec.at("C"), p, n, "estimator C");
    else if (viaObservables) C = sensorJacobian(sensors, xe);
    else { C.assign(n * n, 0.0); for (std::size_t k = 0; k < n; ++k) C[k * n + k] = 1.0; }
    // y from the plant: sensor observables, or C x for an explicit state map.
    auto measure = [sensors, viaObservables, C, n, p](const double* state) {
        if (viaObservables) return sensors(state);
        std::vector<double> y(p, 0.0);
        for (std::size_t k = 0; k < p; ++k) for (std::size_t j = 0; j < n; ++j) y[k] += C[k * n + j] * state[j];
        return y;
    };
    const std::vector<double> ye = viaObservables ? sensors(xe.data()) : std::vector<double>{};

    EstimatorSetup result;
    result.kind = kind; result.sampleTime = sampleTime;
    json::Value::Object evidence{{"estimator_kind", json::Value(kind)},
                                 {"estimator_sample_time", json::Value(sampleTime)}};
    if (kind == "static_map") {
        core::StaticMapEstimator map;
        map.states = n; map.outputs = p;
        map.map = sized(spec.at("M"), n, p, "static_map M");
        map.offset = numbers(spec.at("offset"));
        if (!map.offset.empty() && map.offset.size() != n) throw std::invalid_argument("static_map offset needs one value per state");
        result.observe = [map, measure](std::size_t, double time, const double* state, double) {
            core::KoopmanEstimate estimate;
            estimate.time = time;
            estimate.measurement = measure(state);
            estimate.output = map.estimate(estimate.measurement);
            estimate.lifted = estimate.output;
            estimate.innovation.assign(estimate.measurement.size(), 0.0);
            return estimate;
        };
        result.evidence = json::Value(std::move(evidence));
        return result;
    }
    const std::vector<double> A = spec.at("A").isNull() ? system.jacobian(0.0, xe.data())
                                                        : sized(spec.at("A"), n, n, "estimator A");
    std::vector<double> B = spec.at("B").isNull() ? system.inputMatrix(xe.data()) : numbers(spec.at("B"));
    if (B.size() % n) throw std::invalid_argument("estimator B must have n rows");
    const std::size_t m = B.size() / n;
    core::LinearObserverModel model = kind == "lqe"
        ? core::designLinearQuadraticEstimator(A, B, C, square(spec.at("process_covariance"), n, "process covariance"),
              square(spec.at("measurement_covariance"), p, "measurement covariance"), n, m, p)
        : core::luenbergerObserver(A, B, C, sized(spec.at("L"), n, p, "luenberger L"), n, m, p);
    if (kind != "lqe" && kind != "luenberger") throw std::invalid_argument("unknown classical estimator: " + kind);
    model.stateEquilibrium = xe;
    if (!ye.empty()) model.measurementEquilibrium = ye;
    else {
        model.measurementEquilibrium.assign(p, 0.0);
        for (std::size_t k = 0; k < p; ++k) for (std::size_t j = 0; j < n; ++j)
            model.measurementEquilibrium[k] += C[k * n + j] * xe[j];
    }
    auto initial = numbers(spec.at("initial_estimate"));
    struct Stream { std::unique_ptr<core::LinearObserver> observer; std::vector<double> time, estimate, innovation; };
    auto streams = std::make_shared<std::vector<Stream>>(count);
    for (auto& s : *streams) s.observer = std::make_unique<core::LinearObserver>(model, sampleTime, initial);
    result.observe = [streams, measure, m](std::size_t trajectory, double time, const double* state, double previous) {
        auto& s = streams->at(trajectory);
        core::KoopmanEstimate estimate;
        estimate.time = time;
        const auto y = measure(state);
        estimate.output = s.observer->observe(time, y, std::vector<double>(m, previous));
        estimate.lifted = estimate.output;
        estimate.measurement = s.observer->predictedMeasurement();
        estimate.innovation.resize(y.size());
        for (std::size_t k = 0; k < y.size(); ++k) estimate.innovation[k] = y[k] - estimate.measurement[k];
        s.time.push_back(time);
        s.estimate.insert(s.estimate.end(), estimate.output.begin(), estimate.output.end());
        s.innovation.insert(s.innovation.end(), estimate.innovation.begin(), estimate.innovation.end());
        return estimate;
    };
    result.appendDiagnostics = [streams, model, n, p](std::size_t trajectory, core::NpzRecord& archive) {
        const auto& s = streams->at(trajectory);
        if (s.time.empty()) return;
        archive.arrays["estimator_time"] = core::NpyArray::fromDoubles({s.time.size()}, s.time);
        archive.arrays["state_estimate"] = core::NpyArray::fromDoubles({s.time.size(), n}, s.estimate);
        archive.arrays["innovation"] = core::NpyArray::fromDoubles({s.time.size(), p}, s.innovation);
        archive.arrays["estimator_gain"] = core::NpyArray::fromDoubles({n, p}, model.L);
        archive.arrays["observer_poles"] = core::NpyArray::fromComplex({n}, model.observerPoles);
    };
    double slowest = -1e300;
    for (const auto& pole : model.observerPoles) slowest = std::max(slowest, pole.real());
    evidence["observer_slowest_real_part"] = slowest;
    evidence["observer_integration"] = "exact-zoh-previous-measurement-and-applied-input";
    result.evidence = json::Value(std::move(evidence));
    return result;
}

}  // namespace dynsyslab::workflow::detail
