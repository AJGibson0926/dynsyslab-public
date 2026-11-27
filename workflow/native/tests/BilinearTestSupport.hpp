#pragma once
// Helpers shared by the bilinear K-ROM tests (bilinear_certificate_tests.cpp,
// bilinear_workflow_tests.cpp): checks, synthetic forced Duffing records and job documents.
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/BilinearKrom.hpp>
#include <dynsyslab/core/EstimationRuntime.hpp>
#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/HistoryObservationArchive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanBilinearBridge.hpp>
#include <dynsyslab/core/LQR.hpp>
#include <dynsyslab/core/SupportedResidualSector.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/workflow/Postprocess.hpp>
#include <dynsyslab/workflow/Run.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <stdexcept>

namespace bilinear_test {
namespace c = dynsyslab::core;
namespace w = dynsyslab::workflow;
namespace json = c::json;
namespace fs = std::filesystem;
using Complex = std::complex<double>;

inline int checks = 0;
inline void check(bool ok, const std::string& why) {
    ++checks;
    if (!ok) throw std::runtime_error(why);
}

template <class F> void rejects(F f, const std::string& why) {
    ++checks;
    try { f(); } catch (const std::logic_error&) { return; }
    throw std::runtime_error("accepted: " + why);
}

inline std::string bytes(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

inline c::SystemConfiguration duffing(double beta) {
    c::SystemConfiguration config; config.identifier = "duffing";
    config.parameters = {{"alpha", 2}, {"beta", beta}, {"delta", 1}};
    return config;
}

inline json::Value postJob(const char* method, json::Value settings, const std::vector<fs::path>& sources) {
    json::Value::Array list;
    for (const auto& s : sources) list.emplace_back(s.generic_u8string());
    return json::Value::Object{{"schema", "dynsyslab.compute-job/v1"}, {"kind", "postprocess"},
        {"label", std::string("bilinear-") + method}, {"sources", std::move(list)},
        {"methods", json::Value::Array{json::Value::Object{{"name", method}, {"settings", std::move(settings)}}}}};
}

inline fs::path forcedDuffing(const fs::path& root, double beta, double amplitude) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"bilinear-forced",
      "system":{"name":"duffing","parameters":{"alpha":2,"delta":1}},
      "time":{"start":0,"stop":4,"output_step":0.02},
      "initial_conditions":{"mode":"explicit","values":[[0.04,0.02]],"count":1,"clusters":1},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12,"control_evaluation":"sampled","control_sample_time":0.02},
      "controller":{"kind":"random_pulse","settings":{"model_pulse_duration":0.1,"seed":7}},
      "archive":{"enabled":true}})");
    job.at("system").at("parameters").object()["beta"] = beta;
    job.at("controller").at("settings").object()["model_amplitude"] = amplitude;
    return w::executeJob(job, root).archives.front();
}

inline std::vector<fs::path> forcedWithHolds(const fs::path& root) {
    auto job = json::parse(R"({"schema":"dynsyslab.compute-job/v1","kind":"solve","label":"bilinear-holds",
      "system":{"name":"duffing","parameters":{"alpha":2,"beta":0.2,"delta":1}},
      "time":{"start":0,"stop":0.6,"output_step":0.02},
      "initial_conditions":{"mode":"explicit","count":128,"clusters":1,"values":[]},
      "solver":{"method":"RK45","rtol":1e-10,"atol":1e-12,"control_evaluation":"sampled","control_sample_time":0.02},
      "controller":{"kind":"random_pulse","settings":{"model_amplitude":0.75,"model_pulse_duration":0.06,"seed":7,
        "active_time_window":[0.3,0.6]}},
      "archive":{"enabled":true}})");
    json::Value::Array states;
    for (double radius : {0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.45, 0.5})
        for (int angle = 0; angle < 16; ++angle) {
            const double t = 2 * 3.141592653589793 * angle / 16, a = radius * std::cos(t), b = radius * std::sin(t);
            states.push_back(json::Value::Array{(a - b) / 3, (a + 2 * b) / 3});
        }
    job.at("initial_conditions").at("values") = std::move(states);
    return w::executeJob(job, root).archives;
}

inline std::pair<double, double> heldDuffingStep(double x, double v, double u, double dt) {
    const auto f = [u](double a, double b) { return 2 * a - 0.2 * a * a * a - b + u; };
    const int steps = 40; const double h = dt / steps;
    for (int i = 0; i < steps; ++i) {
        const double k1x = v, k1v = f(x, v);
        const double k2x = v + h / 2 * k1v, k2v = f(x + h / 2 * k1x, v + h / 2 * k1v);
        const double k3x = v + h / 2 * k2v, k3v = f(x + h / 2 * k2x, v + h / 2 * k2v);
        const double k4x = v + h * k3v, k4v = f(x + h * k3x, v + h * k3v);
        x += h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x);
        v += h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v);
    }
    return {x, v};
}

struct HeldOutResult { std::size_t points{}, violations{}; double worstRatio{}, worstZeroInputRatio{}; };

inline HeldOutResult heldOutSectorCheck(const c::NpzRecord& archive, const json::Value& record, double dt) {
    c::BilinearKrom model; model.stateDimension = 2; model.inputDimension = 1; model.sampleTime = dt;
    model.A = archive.at("model_A").doubles(); model.B0 = archive.at("model_B0").doubles();
    model.B = {archive.at("model_B").doubles()};
    const auto X = archive.at("certificate_X").doubles();
    const double l11 = std::sqrt(X[0]), l21 = X[2] / l11, l22 = std::sqrt(X[3] - l21 * l21);
    const double cx = record.at("sector").at("c_x").number(), cu = record.at("sector").at("c_u").number();
    const double ubar = record.at("input_bounds").array().front().number();
    HeldOutResult out;
    for (double scale : {0.2, 0.4, 0.6, 0.8, 1.0})
        for (int angle = 0; angle < 48; ++angle)
            for (double fraction : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
                const double t = 2 * 3.141592653589793 * (angle + 0.37) / 48, a = scale * std::cos(t), b = scale * std::sin(t);
                const std::vector<double> z{l11 * a, l21 * a + l22 * b}, u{fraction * ubar};
                const double x = (z[0] - z[1]) / 3, v = (z[0] + 2 * z[1]) / 3;
                const auto [xp, vp] = heldDuffingStep(x, v, u[0], dt);
                const auto predicted = c::stepBilinearKrom(model, z, u);
                const double r = std::hypot(2 * xp + vp - predicted[0], -xp + vp - predicted[1]);
                const double bound = cx * std::hypot(z[0], z[1]) + cu * std::abs(u[0]);
                const double ratio = r / bound;
                ++out.points;
                if (r > bound) ++out.violations;
                out.worstRatio = std::max(out.worstRatio, ratio);
                if (fraction == 0.0) out.worstZeroInputRatio = std::max(out.worstZeroInputRatio, ratio);
            }
    return out;
}
}  // namespace bilinear_test
