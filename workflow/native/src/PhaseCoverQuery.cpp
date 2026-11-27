// Analytic query operation "phase_cover": a finite prime-power orbit tower whose clocks are checked.
//
// {"schema":"dynsyslab.analytic-query/v1","operation":"phase_cover","system":{...},
//  "initial_states":[[...], ...],  (one orbit per tower level; level j should have period p^j T0)
//  "prime":2, "samples":1024, "relative_frequency_tolerance":1e-6}
// Each orbit is sampled uniformly in model time from its initial state (operation "periodic"), and
// its angular frequency must satisfy |omega_j * p^j / omega_0 - 1| <= tolerance, the clock ratio of
// the classical circle p-cover tower (math::phaseCoverTower; the inverse limit is the p-adic
// solenoid, Vietoris 1927, van Dantzig 1930). A failed check is an error, never a drawn cover.
// Returns dynsyslab.phase-cover-geometry/v1: per level the exact tower labels, frequency,
// period, clock error, closure error and samples. A finite sampled covering, not a conjugacy or an
// infinite-limit certificate.
#include <dynsyslab/workflow/AnalyticQuery.hpp>

#include <dynsyslab/math/Arithmetic.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

namespace dynsyslab::workflow {
namespace json = core::json;

json::Value evaluatePhaseCoverQuery(const json::Value& request) {
    const auto& states = request.at("initial_states").array();
    if (states.empty()) throw std::invalid_argument("phase_cover needs at least one initial state");
    const double prime = request.at("prime").number(0);
    if (prime != 2 && prime != 3 && prime != 5 && prime != 7)
        throw std::invalid_argument("phase_cover prime must be 2, 3, 5 or 7");
    const double tolerance = request.at("relative_frequency_tolerance").number(1e-6);
    if (!(tolerance >= 0 && tolerance <= 0.01)) throw std::invalid_argument("phase_cover tolerance must lie in [0, 0.01]");
    const auto levels = math::phaseCoverTower(static_cast<std::uint32_t>(prime),
                                              std::vector<std::uint32_t>(states.size() - 1, 0));
    json::Value::Array orbits;
    double base = 0;
    std::string fingerprint;
    for (std::size_t j = 0; j < states.size(); ++j) {
        auto sample = request;
        sample.at("operation") = json::Value("periodic");
        sample.at("initial_state") = states[j];
        sample.at("observables") = json::Value(json::Value::Array{});
        sample.at("samples") = json::Value(request.at("samples").number(1024));
        sample.object().erase("initial_states");
        const auto result = evaluateAnalyticQuery(sample);
        const double omega = result.at("angular_frequency").number();
        if (!(omega > 0)) throw std::invalid_argument("phase_cover: an orbit has no positive angular frequency");
        if (j == 0) base = omega;
        const double modulus = std::stod(levels[j].modulus);
        const double error = std::abs(omega * modulus / base - 1.0);
        if (error > tolerance)
            throw std::invalid_argument("phase_cover: orbit " + std::to_string(j) + " clock does not form the requested cover"
                                        " (relative error " + std::to_string(error) + ")");
        fingerprint = result.at("equations_fingerprint").string();
        orbits.emplace_back(json::Value::Object{{"modulus", json::Value(levels[j].modulus)},
            {"residue", json::Value(levels[j].residue)}, {"frequency_ratio", json::Value(levels[j].frequencyRatio)},
            {"phase_offset_turns", json::Value(levels[j].offsetTurns)}, {"angular_frequency", json::Value(omega)},
            {"period", result.at("period")}, {"relative_clock_error", json::Value(error)},
            {"state_dimension", result.at("state_dimension")}, {"values", result.at("values")},
            {"initial_state", states[j]}, {"closure_error", result.at("closure_error")},
            {"equations_fingerprint", result.at("equations_fingerprint")}});
    }
    return json::Value(json::Value::Object{{"schema", json::Value("dynsyslab.phase-cover-geometry/v1")},
        {"prime", json::Value(prime)}, {"samples", json::Value(request.at("samples").number(1024))},
        {"orbits", json::Value(std::move(orbits))}, {"equations_fingerprint", json::Value(fingerprint)},
        {"phase_convention", json::Value("Uniform physical time from each supplied initial state; endpoint omitted")},
        {"scope", json::Value("Finite sampled circle coverings; not a physical-state conjugacy or an infinite-limit certificate")}});
}

}  // namespace dynsyslab::workflow
