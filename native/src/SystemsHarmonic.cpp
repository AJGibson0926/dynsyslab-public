#include "SystemSupport.hpp"
#include <cmath>
#include <stdexcept>
namespace dynsyslab::core::detail {
std::unique_ptr<ODESystem> makeHarmonicSystem(SystemConfiguration configuration) {
    const double mass = resolvedParameter(configuration, "mass_kg", 1);
    const double stiffness = resolvedParameter(configuration, "stiffness_n_m", 1);
    const double damping = resolvedParameter(configuration, "damping_kg_s",
                                             configuration.preset == "damped" ? 0.2 : 0.0);
    if (!std::isfinite(mass) || mass <= 0 || !std::isfinite(stiffness) || !std::isfinite(damping))
        throw std::invalid_argument(
            "Oscillator requires positive finite mass and finite coefficients");
    auto rhs = [=](double, const double *x, double u, double *f) {
        f[0] = x[1];
        f[1] = (u - damping * x[1] - stiffness * x[0]) / mass;
    };
    auto jacobian = [=](double, const double *) {
        return std::vector<double>{0, 1, -stiffness / mass, -damping / mass};
    };
    auto input = [=](const double *) { return std::vector<double>{0, 1 / mass}; };
    const auto q=AnalyticExpression::state(2,0), v=AnalyticExpression::state(2,1);
    return std::make_unique<FunctionalSystem>(
        descriptorFor("harmonic_oscillator"), configuration, selectedInitial(configuration, {1, 0}),
        rhs, jacobian, input, [](double, const double *x) { return genericInvalid(x, 2); },
        std::vector<double>{}, std::vector<Observable>{}, std::map<std::string,std::string>{},
        std::map<std::string,std::vector<std::string>>{},
        std::vector<AnalyticExpression>{v,(-damping*v-stiffness*q)/mass},
        std::vector<AnalyticExpression>{AnalyticExpression::constant(0,2),
                                        AnalyticExpression::constant(1/mass,2)});
}
} // namespace dynsyslab::core::detail
