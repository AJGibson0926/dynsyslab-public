#pragma once
// JSON <-> core settings for noise-robust differentiation and measurement noise, shared by
// coordinate_transform (measured_derivatives), derivative_benchmark, coordinate SINDy, input_map
// and the KLQE estimator.
//
//   "derivative": {"method": "savitzky_golay", "cutoff": 20, "window": 0, "order": 3, "sigma": 0,
//                  "points": 0, "stride": 0, "lambda": 0, "penalty_order": 3, "process_noise": 0,
//                  "measurement_noise": 0, "weight": 0, "exponent": 4, "lipschitz": 0, "iterations": 0}
//   "noise": {"kind": "white", "level": 0.01, "relative": true, "correlation_time": 0,
//             "probability": 0, "bits": 0, "seed": 20261115, "parts": [...]}
// A bare string is shorthand for {"method": "<string>"} / {"kind": "<string>"}.
#include <dynsyslab/core/NumericalDifferentiation.hpp>
#include <dynsyslab/workflow/Json.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::workflow {

inline core::DerivativeSpec parseDerivativeSpec(const json::Value& v) {
    core::DerivativeSpec s;
    if (v.isNull()) return s;
    if (const auto* text = std::get_if<std::string>(&v.storage())) { s.method = *text; return s; }
    auto count = [&](const char* key, std::size_t fallback) {
        const double x = v.at(key).number(double(fallback));
        if (!(x >= 0) || x != std::floor(x)) throw std::invalid_argument(std::string("derivative.") + key + " must be a nonnegative integer");
        return static_cast<std::size_t>(x);
    };
    s.method = v.at("method").string(s.method);
    s.cutoff = v.at("cutoff").number(0);
    s.window = count("window", 0); s.order = count("order", 3); s.points = count("points", 0); s.stride = count("stride", 0);
    s.sigma = v.at("sigma").number(0); s.lambda = v.at("lambda").number(0);
    s.penaltyOrder = count("penalty_order", 3);
    s.processNoise = v.at("process_noise").number(0); s.measurementNoise = v.at("measurement_noise").number(0);
    s.weight = v.at("weight").number(0); s.exponent = v.at("exponent").number(4); s.lipschitz = v.at("lipschitz").number(0);
    s.iterations = count("iterations", 0);
    s.despike = count("despike", 0); s.despikeThreshold = v.at("despike_threshold").number(3);
    bool known = false; for (const auto& m : core::derivativeMethods()) known = known || m == s.method;
    if (!known) throw std::invalid_argument("unknown derivative method: " + s.method);
    return s;
}

inline core::NoiseSpec parseNoiseSpec(const json::Value& v) {
    core::NoiseSpec s;
    if (v.isNull()) return s;
    if (const auto* text = std::get_if<std::string>(&v.storage())) { s.kind = *text; return s; }
    s.kind = v.at("kind").string("none");
    s.level = v.at("level").number(0); s.relative = v.at("relative").boolean(true);
    s.correlationTime = v.at("correlation_time").number(0); s.probability = v.at("probability").number(0);
    s.bits = v.at("bits").number(0);
    const double seed = v.at("seed").number(20261115.0);
    if (!(seed >= 0)) throw std::invalid_argument("noise.seed must be nonnegative");
    s.seed = static_cast<std::uint64_t>(seed);
    for (const auto& p : v.at("parts").array()) { auto part = parseNoiseSpec(p); part.seed ^= s.seed + s.parts.size() * 7919u; s.parts.push_back(part); }
    if (s.kind == "sum" && s.parts.empty()) throw std::invalid_argument("noise kind sum needs parts");
    return s;
}

}  // namespace dynsyslab::workflow
