// mode_selection scoring on synthetic eigenfunction records with known answers.
#include "../src/ModeSelectionInternal.hpp"

#include <dynsyslab/core/KoopmanModel.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace ms = dynsyslab::workflow::postprocess_detail::mode_selection;
namespace json = dynsyslab::core::json;
using Complex = std::complex<double>;
int failures = 0;
void check(bool ok, const std::string& what) { if (!ok) { ++failures; std::cerr << "FAIL: " << what << '\n'; } }

const Complex lambda0(-0.1, 2.0), b0(0.5, 0.25), nu(0.4, -0.2);
const double lambda2 = -0.5, b2 = -0.3;

// psi(t) for a constant input u from psi(0): e^{lt} psi0 + B (e^{lt} - 1) / l.
Complex flow(Complex l, Complex start, Complex b, double u, double t) {
    return std::exp(l * t) * start + b * u * (std::exp(l * t) - 1.0) / l;
}

ms::Record record(bool forced) {
    ms::Record r; r.forced = forced; r.dt = 0.01; r.modes = 3; r.stateDimension = 2;
    const double u = forced ? 1.0 : 0.0;
    for (int k = 0; k <= 2000; ++k) {
        const double t = k * r.dt;
        const Complex p0 = flow(lambda0, 2.0, b0, u, t), p2 = flow(lambda2, 1.5, b2, u, t);
        r.time.push_back(t); r.input.push_back(u);
        r.psi.insert(r.psi.end(), {p0, std::conj(p0), p2});
        r.echo.push_back(0.3 + 2.0 * (nu * p0).real() + 0.7 * p2.real());
        r.states.insert(r.states.end(), {p0.real(), p2.real()});
    }
    r.samples = r.time.size();
    return r;
}
}  // namespace

int main() {
    try {
        dynsyslab::core::SystemConfiguration config; config.identifier = "lorenz";
        const auto model = dynsyslab::core::KoopmanModel::fromPolynomialCoefficients(config, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1},
            {{1, 0}, {0, 1}, {0, 0}, {1, 0}, {0, -1}, {0, 0}, {0, 0}, {0, 0}, {1, 0}},
            {lambda0, std::conj(lambda0), {lambda2, 0}});
        const std::vector<ms::Record> records{record(false), record(true)};
        ms::Options options;
        const auto units = ms::buildUnits(model, records, json::Value::Object{{"invariance_horizon", 0.25}}, options);
        check(units.size() == 2, "a conjugate pair and a real mode make two units");
        check(units[0].pair && units[0].modes[0] == 0 && units[0].modes[1] == 1, "pair unit keeps the upper mode first");
        check(units[0].invariance < 1e-10 && units[1].invariance < 1e-10, "exact exponentials are invariant");
        check(std::abs(units[0].inputMap - b0) < 2e-4 * std::abs(b0), "complex B_psi recovered from a constant input");
        check(std::abs(units[1].inputMap.real() - b2) < 2e-4 * std::abs(b2), "real B_psi recovered");
        check(units[0].inputMapResidual < 1e-3, "integral-form regression explains the forced response");
        const auto both = ms::fitReadouts(units, {0, 1}, records, 0.0);
        check(both.echoError < 1e-6 && both.stateError[0] < 1e-6 && both.stateError[1] < 1e-6, "exact readouts with both units");
        check(std::abs(both.echoOffset - 0.3) < 1e-8, "echo intercept");
        const auto one = ms::fitReadouts(units, {0}, records, 0.0);
        check(one.echoError > 1e-3, "dropping the real mode leaves an echo residual");
        options.ridge = 0.0;
        const auto score = ms::scoreSubset(units, {0, 1}, records, options);
        check(score.valid, "subset scores: " + score.failure);
        check(score.readouts.echoError < 1e-6 && score.inputMapResidual < 1e-3, "score carries readout and map quality");
        check(score.predictedStateError > 0 && std::isfinite(score.predictedStateError), "predicted KLQE state error finite");
        check(score.stateReach > 0 && score.observabilityBalance > 0, "controllable and observable");
    } catch (const std::exception& error) {
        std::cerr << "mode selection tests threw: " << error.what() << '\n';
        return 1;
    }
    if (failures) { std::cerr << failures << " mode selection checks failed\n"; return 1; }
    std::cout << "dynsyslab_mode_selection_tests: PASS\n";
    return 0;
}
