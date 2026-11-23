#include "SystemSupport.hpp"
#include <dynsyslab/core/OscillatorObservables.hpp>

#include <algorithm>
#include <cmath>

namespace dynsyslab::core::detail {
namespace {

std::unique_ptr<ODESystem> duffing(SystemConfiguration configuration) {
    double delta = 0.2, alpha = 1.0, beta = 1.0, gain = 1.0;
    std::vector<double> initial{0.0, 0.0};
    if (configuration.preset == "double_well_reference") initial = {0.1, 0.0};
    if (configuration.preset == "hardening_single_well") initial = {0.2, 0.0}, alpha = -1.0;
    if (configuration.preset == "gibson_undamped_klqr_2022") delta = 0.0, initial = {0.0, -2.8};
    if (configuration.preset == "gibson_damped_stabilization_2022")
        delta = 0.5, beta = 4.0, initial = {0.5, 0.0};
    if (configuration.preset == "korda_mezic_optimal_2020") delta = 0.5, beta = 4.0, gain = 0.5;
    delta = resolvedParameter(configuration, "delta", delta);
    alpha = resolvedParameter(configuration, "alpha", alpha);
    beta = resolvedParameter(configuration, "beta", beta);
    gain = resolvedParameter(configuration, "input_gain", gain);
    auto rhs = [delta, alpha, beta, gain](double, const double* x, double input, double* out) {
        out[0] = x[1];
        out[1] = -delta * x[1] + alpha * x[0] - beta * x[0] * x[0] * x[0] + gain * input;
    };
    auto jacobian = [delta, alpha, beta](double, const double* x) {
        return std::vector<double>{0.0, 1.0, alpha - 3.0 * beta * x[0] * x[0], -delta};
    };
    auto input = [gain](const double*) { return std::vector<double>{0.0, gain}; };
    auto invalid = [](double, const double* x) { return genericInvalid(x, 2); };
    const auto scalar = [](const char* name, const char* label,
                           Observable::Evaluator evaluate,
                           Observable::Jacobian jacobian = {}, bool complex = false) {
        return Observable({name, label, "dimensionless", "dimensionless", {}, complex},
                          2, 1, std::move(evaluate), std::move(jacobian));
    };
    std::vector<Observable> observables;
    // Unit-mass, dimensionless Duffing approximation to a magnetic double well.
    // This state observable excludes the separately applied gain * input.
    const auto q = AnalyticExpression::state(2, 0), v = AnalyticExpression::state(2, 1);
    const auto force = alpha * q - beta * q.power(3) - delta * v;
    observables.push_back(Observable::fromExpressions(
        {"force", "Intrinsic force", "dimensionless", "dimensionless", {}, false}, {force}));
    observables.push_back(scalar("kinetic_energy", "Kinetic energy",
        [](const double* x) { return Observable::Values{0.5*x[1]*x[1]}; },
        [](const double* x) { return Observable::Values{0.0,x[1]}; }));
    observables.push_back(scalar("potential_energy", "Potential energy",
        [alpha,beta](const double* x) {
            return Observable::Values{-0.5*alpha*x[0]*x[0]+0.25*beta*std::pow(x[0],4)};
        }, [alpha,beta](const double* x) {
            return Observable::Values{-alpha*x[0]+beta*x[0]*x[0]*x[0],0.0};
        }));
    observables.push_back(scalar("total_energy", "Total energy",
        [alpha,beta](const double* x) {
            return Observable::Values{0.5*x[1]*x[1]-0.5*alpha*x[0]*x[0]
                +0.25*beta*std::pow(x[0],4)};
        }, [alpha,beta](const double* x) {
            return Observable::Values{-alpha*x[0]+beta*x[0]*x[0]*x[0],x[1]};
        }));
    observables.push_back(scalar("random_observation", "Random observation",
        [](const double* x) {
            const std::complex<double> value = std::exp(std::complex<double>(0.0,x[0]))
                / std::cosh(x[1]); return Observable::Values{value};
        }, [](const double* x) {
            const std::complex<double> value = std::exp(std::complex<double>(0.0,x[0]))
                / std::cosh(x[1]);
            return Observable::Values{std::complex<double>(0.0,1.0)*value,
                                      -value*std::tanh(x[1])};
        }, true));
    auto phase=phaseEnvelopeObservables("phase_sech",1,OscillatorEnvelope::sech);
    observables.insert(observables.end(),phase.begin(),phase.end());
    observables.push_back(gaussianSumObservable());
    const std::vector<double> controlledInitial{
        alpha > 0.0 && beta > 0.0 ? std::sqrt(alpha / beta) : 0.0, 0.0};
    return std::make_unique<FunctionalSystem>(descriptorFor("duffing"), configuration,
        selectedInitial(configuration, initial), rhs, jacobian, input, invalid, controlledInitial,
        std::move(observables), std::map<std::string,std::string>{{"x","position"},
        {"xdot","velocity"},{"force_dimensionless","force"}},
        std::map<std::string,std::vector<std::string>>{},
        std::vector<AnalyticExpression>{v, force},
        std::vector<AnalyticExpression>{AnalyticExpression::constant(0., 2),
                                      AnalyticExpression::constant(gain, 2)});
}

std::unique_ptr<ODESystem> vanDerPol(SystemConfiguration configuration) {
    double mu = configuration.preset == "weakly_nonlinear" ? 0.1 : 8.0;
    mu = resolvedParameter(configuration, "mu", mu);
    auto rhs = [mu](double, const double* x, double input, double* out) {
        out[0] = x[1]; out[1] = mu * (1.0 - x[0] * x[0]) * x[1] - x[0] + input;
    };
    auto jacobian = [mu](double, const double* x) {
        return std::vector<double>{0.0, 1.0, -2.0 * mu * x[0] * x[1] - 1.0,
                                   mu * (1.0 - x[0] * x[0])};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("van_der_pol"), configuration,
        selectedInitial(configuration, {0.1, 0.0}), rhs, jacobian,
        [](const double*) { return std::vector<double>{0.0, 1.0}; },
        [](double, const double* x) { return genericInvalid(x, 2); });
}

std::unique_ptr<ODESystem> quadraticSaddleNode(SystemConfiguration configuration) {
    const double mu = resolvedParameter(configuration, "mu", 0.0);
    const double eta = resolvedParameter(configuration, "eta", 0.0);
    auto rhs = [mu, eta](double, const double *x, double input, double *out) {
        out[0] = mu + eta * x[0] - x[0] * x[0] + input;
    };
    auto jacobian = [eta](double, const double *x) {
        return std::vector<double>{eta - 2 * x[0]};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("quadratic_saddle_node"),
        configuration, selectedInitial(configuration, {1.0}), rhs, jacobian,
        [](const double *) { return std::vector<double>{1}; },
        [](double, const double *x) { return genericInvalid(x, 1); });
}

std::unique_ptr<ODESystem> cubicBifurcation(SystemConfiguration configuration) {
    const double epsilon = resolvedParameter(configuration, "epsilon", 0.0);
    const double mu = resolvedParameter(configuration, "mu", 0.0);
    double quadratic = resolvedParameter(configuration, "quadratic", 0.0);
    double cubic = resolvedParameter(configuration, "cubic", -1.0);
    if (configuration.preset == "imperfect_pitchfork")
        quadratic = resolvedParameter(configuration, "quadratic", .2);
    auto rhs = [epsilon, mu, quadratic, cubic]
        (double, const double *x, double input, double *out) {
        out[0] = epsilon + mu * x[0] + quadratic * x[0] * x[0] +
                 cubic * x[0] * x[0] * x[0] + input;
    };
    auto jacobian = [mu, quadratic, cubic](double, const double *x) {
        return std::vector<double>{mu + 2 * quadratic * x[0] +
                                   3 * cubic * x[0] * x[0]};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("cubic_bifurcation"),
        configuration, selectedInitial(configuration, {0.0}), rhs, jacobian,
        [](const double *) { return std::vector<double>{1}; },
        [](double, const double *x) { return genericInvalid(x, 1); });
}

std::unique_ptr<ODESystem> radialLimitCycle(SystemConfiguration configuration) {
    const double mu = resolvedParameter(configuration, "mu", 1.0);
    const double alpha = resolvedParameter(configuration, "alpha", 1.0);
    const double omega = resolvedParameter(configuration, "omega", 1.0);
    auto rhs = [mu, alpha, omega](double, const double *x, double input, double *out) {
        const double radial = mu - alpha * (x[0] * x[0] + x[1] * x[1]);
        out[0] = radial * x[0] - omega * x[1] + input;
        out[1] = omega * x[0] + radial * x[1];
    };
    auto jacobian = [mu, alpha, omega](double, const double *x) {
        return std::vector<double>{
            mu - alpha * (3 * x[0] * x[0] + x[1] * x[1]),
            -omega - 2 * alpha * x[0] * x[1],
            omega - 2 * alpha * x[0] * x[1],
            mu - alpha * (x[0] * x[0] + 3 * x[1] * x[1])};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("radial_limit_cycle"),
        configuration, selectedInitial(configuration, {1.0, 0.0}), rhs, jacobian,
        [](const double *) { return std::vector<double>{1, 0}; },
        [](double, const double *x) { return genericInvalid(x, 2); });
}

std::unique_ptr<ODESystem> lorenz(SystemConfiguration configuration) {
    const double sigma = resolvedParameter(configuration, "sigma", 10.0);
    const double rho = resolvedParameter(configuration, "rho",
        configuration.preset == "near_classic_transition" ? 24.06 : 28.0);
    const double beta = resolvedParameter(configuration, "beta", 8.0 / 3.0);
    auto rhs = [sigma, rho, beta](double, const double* x, double input, double* out) {
        out[0] = sigma * (x[1] - x[0]) + input;
        out[1] = x[0] * (rho - x[2]) - x[1];
        out[2] = x[0] * x[1] - beta * x[2];
    };
    auto jacobian = [sigma, rho, beta](double, const double* x) {
        return std::vector<double>{-sigma, sigma, 0.0, rho - x[2], -1.0, -x[0],
                                   x[1], x[0], -beta};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("lorenz"), configuration,
        selectedInitial(configuration, {1.0, 1.0, 1.0}), rhs, jacobian,
        [](const double*) { return std::vector<double>{1.0, 0.0, 0.0}; },
        [](double, const double* x) { return genericInvalid(x, 3); });
}

std::unique_ptr<ODESystem> rossler(SystemConfiguration configuration) {
    double a = 0.2, b = 0.2, c = 5.7;
    std::vector<double> initial{0.1, 0.1, 0.1};
    if (configuration.preset == "common_alternative" || configuration.preset == "brunton_shavok_2021")
        a = 0.1, b = 0.1, c = 14.0;
    if (configuration.preset == "brunton_shavok_2021") initial = {1,1,1};
    if (configuration.preset == "funnel_chaos_2021") a = 0.2925, b = 0.1, c = 8.5, initial = {1,1,1};
    if (configuration.preset == "b045_chaos_2013") b = 0.45;
    a = resolvedParameter(configuration, "a", a); b = resolvedParameter(configuration, "b", b);
    c = resolvedParameter(configuration, "c", c);
    auto rhs = [a, b, c](double, const double* x, double input, double* out) {
        out[0] = -x[1] - x[2] + input; out[1] = x[0] + a * x[1];
        out[2] = b + x[2] * (x[0] - c);
    };
    auto jacobian = [a, c](double, const double* x) {
        return std::vector<double>{0,-1,-1, 1,a,0, x[2],0,x[0]-c};
    };
    return std::make_unique<FunctionalSystem>(descriptorFor("rossler"), configuration,
        selectedInitial(configuration, initial), rhs, jacobian,
        [](const double*) { return std::vector<double>{1,0,0}; },
        [](double, const double* x) { return genericInvalid(x, 3); });
}

std::unique_ptr<ODESystem> chua(SystemConfiguration configuration) {
    const double alpha = resolvedParameter(configuration, "alpha", 15.6);
    const double beta = resolvedParameter(configuration, "beta", 28.0);
    const double m0 = resolvedParameter(configuration, "m0", -1.143);
    const double m1 = resolvedParameter(configuration, "m1", -0.714);
    auto h = [m0, m1](double x) { return m1*x + 0.5*(m0-m1)*(std::abs(x+1)-std::abs(x-1)); };
    auto rhs = [alpha, beta, h](double, const double* x, double input, double* out) {
        out[0] = alpha * (x[1] - x[0] - h(x[0])) + input;
        out[1] = x[0] - x[1] + x[2]; out[2] = -beta * x[1];
    };
    auto jacobian = [alpha, beta, m0, m1](double, const double* x) {
        const double dh = std::abs(x[0]) < 1.0 ? m0 : m1;
        return std::vector<double>{-alpha*(1+dh),alpha,0, 1,-1,1, 0,-beta,0};
    };
    auto initial = std::vector<double>{0.1,0,0};
    if (configuration.preset == "canonical_double_scroll") initial = {0.7,0,0};
    if (configuration.preset == "double_scroll_seed_left") initial = {-0.7,0,0};
    return std::make_unique<FunctionalSystem>(descriptorFor("chua"), configuration,
        selectedInitial(configuration, initial), rhs, jacobian,
        [](const double*) { return std::vector<double>{1,0,0}; },
        [](double, const double* x) { return genericInvalid(x, 3); });
}

std::unique_ptr<ODESystem> henonHeiles(SystemConfiguration configuration) {
    const double lambda = resolvedParameter(configuration, "lam", 1.0);
    auto rhs = [lambda](double, const double* x, double input, double* out) {
        out[0]=x[2]; out[1]=x[3]; out[2]=-x[0]-2*lambda*x[0]*x[1]+input;
        out[3]=-x[1]-lambda*(x[0]*x[0]-x[1]*x[1]);
    };
    auto jacobian = [lambda](double,const double* x) {
        return std::vector<double>{0,0,1,0, 0,0,0,1,
            -1-2*lambda*x[1],-2*lambda*x[0],0,0,
            -2*lambda*x[0],-1+2*lambda*x[1],0,0};
    };
    auto initial = configuration.preset == "near_escape_energy"
        ? std::vector<double>{0,0,0,0.55} : std::vector<double>{0,0.1,0,0.1};
    return std::make_unique<FunctionalSystem>(descriptorFor("henon_heiles"), configuration,
        selectedInitial(configuration, initial), rhs, jacobian,
        [](const double*) { return std::vector<double>{0,0,1,0}; },
        [](double,const double* x) { return genericInvalid(x,4); });
}

}  // namespace

std::unique_ptr<ODESystem> makeBasicSystem(const SystemConfiguration& configuration) {
    if (configuration.identifier == "quadratic_triangular") return makeTriangularSystem(configuration);
    if (configuration.identifier == "harmonic_oscillator") return makeHarmonicSystem(configuration);
    if (configuration.identifier == "duffing") return duffing(configuration);
    if (configuration.identifier == "van_der_pol") return vanDerPol(configuration);
    if (configuration.identifier == "quadratic_saddle_node")
        return quadraticSaddleNode(configuration);
    if (configuration.identifier == "cubic_bifurcation")
        return cubicBifurcation(configuration);
    if (configuration.identifier == "radial_limit_cycle")
        return radialLimitCycle(configuration);
    if (configuration.identifier == "lorenz") return lorenz(configuration);
    if (configuration.identifier == "rossler") return rossler(configuration);
    if (configuration.identifier == "chua") return chua(configuration);
    if (configuration.identifier == "henon_heiles") return henonHeiles(configuration);
    return {};
}

}  // namespace dynsyslab::core::detail
