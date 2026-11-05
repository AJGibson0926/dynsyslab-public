#include "SystemSupport.hpp"
#include "HodgkinHuxleyRates.hpp"

#include <cmath>
#include <utility>

namespace dynsyslab::core::detail {
namespace {

std::unique_ptr<ODESystem> hindmarshRose(SystemConfiguration configuration) {
    double a = 1.0, b = 3.0, c = 1.0, d = 5.0, r = 0.006, s = 4.0,
        xRest = -1.6, current = 3.25;
    std::vector<double> initial{-0.5, -1.0, 3.0};
    if (configuration.preset == "standard_spiking") {
        b = 2.7; r = 0.01; current = 4.0; initial = {-1.0, -10.0, 1.0};
    }
    a = resolvedParameter(configuration, "a", a); b = resolvedParameter(configuration, "b", b);
    c = resolvedParameter(configuration, "c", c); d = resolvedParameter(configuration, "d", d);
    r = resolvedParameter(configuration, "r", r); s = resolvedParameter(configuration, "s", s);
    xRest = resolvedParameter(configuration, "x_r", xRest);
    current = resolvedParameter(configuration, "current", current);
    auto rhs = [a, b, c, d, r, s, xRest, current]
        (double, const double* x, double input, double* out) {
        out[0] = x[1] - a*x[0]*x[0]*x[0] + b*x[0]*x[0] - x[2] + current + input;
        out[1] = c - d*x[0]*x[0] - x[1];
        out[2] = r * (s * (x[0] - xRest) - x[2]);
    };
    auto jacobian = [a, b, d, r, s](double, const double* x) {
        return std::vector<double>{
            -3.0*a*x[0]*x[0] + 2.0*b*x[0], 1.0, -1.0,
            -2.0*d*x[0], -1.0, 0.0,
            r*s, 0.0, -r};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("hindmarsh_rose"), configuration,
        selectedInitial(configuration, std::move(initial)), rhs, jacobian,
        [](const double*) { return std::vector<double>{1.0, 0.0, 0.0}; },
        [](double, const double* x) { return genericInvalid(x, 3); });
}

HodgkinHuxleyRates rateDerivatives(double voltage) {
    constexpr double step = 1.0e-4;
    const auto plus = hodgkinHuxleyRates(voltage + step), minus = hodgkinHuxleyRates(voltage - step);
    return {
        (plus.am - minus.am) / (2.0 * step), (plus.bm - minus.bm) / (2.0 * step),
        (plus.ah - minus.ah) / (2.0 * step), (plus.bh - minus.bh) / (2.0 * step),
        (plus.an - minus.an) / (2.0 * step), (plus.bn - minus.bn) / (2.0 * step)};
}

std::unique_ptr<ODESystem> hodgkinHuxley(SystemConfiguration configuration) {
    double capacitance = 1.0, gNa = 120.0, gK = 36.0, gL = 0.3,
        eNa = 50.0, eK = -77.0, eL = -54.387, bias = 0.0;
    if (configuration.preset == "tonic_firing_10") bias = 10.0;
    capacitance = resolvedParameter(configuration, "capacitance", capacitance);
    gNa = resolvedParameter(configuration, "g_na", gNa); gK = resolvedParameter(configuration, "g_k", gK);
    gL = resolvedParameter(configuration, "g_l", gL); eNa = resolvedParameter(configuration, "e_na", eNa);
    eK = resolvedParameter(configuration, "e_k", eK); eL = resolvedParameter(configuration, "e_l", eL);
    bias = resolvedParameter(configuration, "bias_current", bias);
    auto rhs = [capacitance, gNa, gK, gL, eNa, eK, eL, bias]
        (double, const double* x, double input, double* out) {
        const auto rate = hodgkinHuxleyRates(x[0]);
        const double ionic = gNa*x[1]*x[1]*x[1]*x[2]*(x[0] - eNa)
            + gK*x[3]*x[3]*x[3]*x[3]*(x[0] - eK) + gL*(x[0] - eL);
        out[0] = 1000.0 * (bias + input - ionic) / capacitance;
        out[1] = 1000.0 * (rate.am * (1.0 - x[1]) - rate.bm * x[1]);
        out[2] = 1000.0 * (rate.ah * (1.0 - x[2]) - rate.bh * x[2]);
        out[3] = 1000.0 * (rate.an * (1.0 - x[3]) - rate.bn * x[3]);
    };
    auto jacobian = [capacitance, gNa, gK, gL, eNa, eK]
        (double, const double* x) {
        const auto rate = hodgkinHuxleyRates(x[0]), derivative = rateDerivatives(x[0]);
        const double scale = 1000.0 / capacitance;
        std::vector<double> matrix(16, 0.0);
        matrix[0] = -scale * (gNa*x[1]*x[1]*x[1]*x[2]
            + gK*x[3]*x[3]*x[3]*x[3] + gL);
        matrix[1] = -scale * 3.0*gNa*x[1]*x[1]*x[2]*(x[0] - eNa);
        matrix[2] = -scale * gNa*x[1]*x[1]*x[1]*(x[0] - eNa);
        matrix[3] = -scale * 4.0*gK*x[3]*x[3]*x[3]*(x[0] - eK);
        matrix[4] = 1000.0 * (derivative.am*(1.0 - x[1]) - derivative.bm*x[1]);
        matrix[5] = 1000.0 * (-rate.am - rate.bm);
        matrix[8] = 1000.0 * (derivative.ah*(1.0 - x[2]) - derivative.bh*x[2]);
        matrix[10] = 1000.0 * (-rate.ah - rate.bh);
        matrix[12] = 1000.0 * (derivative.an*(1.0 - x[3]) - derivative.bn*x[3]);
        matrix[15] = 1000.0 * (-rate.an - rate.bn);
        return matrix;
    };
    auto input = [capacitance](const double*) {
        return std::vector<double>{1000.0 / capacitance, 0.0, 0.0, 0.0};
    };
    auto invalid = [](double, const double* x) {
        return genericInvalid(x, 4) || std::abs(x[0]) > 1.0e5
            || x[1] < -1.0 || x[1] > 2.0 || x[2] < -1.0 || x[2] > 2.0
            || x[3] < -1.0 || x[3] > 2.0;
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("hodgkin_huxley"), configuration,
        selectedInitial(configuration, {-65.0, 0.0529, 0.596, 0.3177}),
        rhs, jacobian, input, invalid);
}

}  // namespace

std::unique_ptr<ODESystem> makeBiologicalSystem(const SystemConfiguration& configuration) {
    if (configuration.identifier == "hindmarsh_rose") return hindmarshRose(configuration);
    if (configuration.identifier == "hodgkin_huxley") return hodgkinHuxley(configuration);
    return {};
}

}  // namespace dynsyslab::core::detail
