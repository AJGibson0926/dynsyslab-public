#include <algorithm>
#include <cmath>
#include <dynsyslab/core/Delay.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
namespace {
constexpr double pi = 3.14159265358979323846;
void check(bool ok, const char *text) {
    if (!ok)
        throw std::runtime_error(text);
}
void near(double a, double b, double tol, const char *text) {
    if (!std::isfinite(a) || !std::isfinite(b) || std::abs(a - b) > tol * (1 + std::abs(b)))
        throw std::runtime_error(std::string(text) + ": " + std::to_string(a) +
                                 " != " + std::to_string(b));
}
void transforms() {
    for (unsigned n : {1, 7, 16, 33}) {
        std::vector<std::complex<double>> input(n);
        for (unsigned i = 0; i < n; ++i)
            input[i] = {std::sin(i * .3), std::cos(i * .7)};
        auto frequency = discreteFourierTransform(input);
        auto back = discreteFourierTransform(frequency, true);
        for (unsigned k = 0; k < n; ++k) {
            std::complex<double> direct{};
            for (unsigned j = 0; j < n; ++j)
                direct += input[j] * std::polar(1., -2 * pi * j * k / n);
            near(std::abs(frequency[k] - direct), 0, 3e-13, "DFT versus defining sum");
            near(std::abs(back[k] - input[k]), 0, 3e-13, "Arbitrary-length inverse DFT");
        }
    }
}
void heat() {
    for (int boundary : {0, 1, 2}) {
        SystemConfiguration c;
        c.identifier = "heat_equation_1d";
        c.parameters = {{"points", 17},
                        {"length", 2},
                        {"diffusivity", .04},
                        {"decay_rate", .13},
                        {"boundary_code", double(boundary)}};
        c = resolveHeatEquation(c);
        auto coordinate = heatEquationCoordinates(c);
        AnalyticState initial(coordinate.size());
        const double wave = boundary == 0 ? 2 * pi : pi;
        for (unsigned j = 0; j < initial.size(); ++j)
            initial[j] = boundary == 1 ? std::cos(wave * coordinate[j] / 2)
                                       : std::sin(wave * coordinate[j] / 2);
        if (boundary == 2)
            initial.front() = initial.back() = 0;
        auto descriptor = heatEquationDescriptor(initial.size());
        check(bool(descriptor.analyticSolutions),
              "PDE analytic families are not a system property");
        auto families = descriptor.analyticSolutions(c);
        auto continuum = families[0].instantiate(initial, 2);
        auto discrete = families[1].instantiate(initial, 2);
        check(bool(continuum) && bool(discrete), "Heat family did not instantiate");
        const auto values = continuum.evaluate(2.7);
        for (unsigned j = 0; j < initial.size(); ++j)
            near(values[j], initial[j] * std::exp(-(.04 * wave * wave / 4 + .13) * .7), 2e-13,
                 "Continuum heat mode decay");
        auto rhs = heatEquationRightHandSide(c);
        auto mesh = heatEquationMesh(c);
        const double time = 2.23, epsilon = 1e-6;
        auto x = discrete.evaluate(time), left = discrete.evaluate(time - epsilon),
             right = discrete.evaluate(time + epsilon);
        std::vector<double> rate(x.size());
        rhs(time, mesh, x.data(), x.size(), rate.data());
        for (unsigned j = 0; j < x.size(); ++j)
            near((right[j] - left[j]) / (2 * epsilon), rate[j], 2e-9,
                 "Discrete semigroup solves actual boundary stencil");
        check(!descriptor.research.papers.empty(), "PDE system bibliography missing");
    }
    SystemConfiguration c;
    c.identifier = "heat_equation_1d";
    c.parameters = {{"points", 33}, {"boundary_code", 2}};
    c = resolveHeatEquation(c);
    check(!heatEquationAnalyticSolutions(c)[0].instantiate(AnalyticState(33, 1), 0),
          "Nonzero Dirichlet endpoints accepted");
    double priorError = 0;
    for (int n : {16, 32}) {
        c.parameters["points"] = n;
        c.parameters["boundary_code"] = 0;
        c.parameters["decay_rate"] = 0;
        auto initial = heatEquationInitialField(c);
        auto families = heatEquationAnalyticSolutions(c);
        auto continuum = families[0].instantiate(initial, 0).evaluate(.2);
        auto discrete = families[1].instantiate(initial, 0).evaluate(.2);
        double error = 0;
        for (unsigned i = 0; i < initial.size(); ++i)
            error = std::max(error, std::abs(continuum[i] - discrete[i]));
        if (priorError)
            check(priorError / error > 3.8 && priorError / error < 4.2,
                  "Spatial refinement is not second order");
        priorError = error;
    }
}
void delay() {
    for (double b : {-.8, .15, 1.2})
        for (double gamma : {-.2, 0., 1e-10, .35}) {
            const double tau = 1.1, h = 1.3;
            for (double t : {.2, 1.1, 1.9, 2.7, 3.8}) {
                const double eps = 2e-6;
                const double x = linearDelayConstantHistory(b, gamma, tau, h, t);
                const double delayed = linearDelayConstantHistory(b, gamma, tau, h, t - tau);
                const double derivative = (linearDelayConstantHistory(b, gamma, tau, h, t + eps) -
                                           linearDelayConstantHistory(b, gamma, tau, h, t - eps)) /
                                          (2 * eps);
                near(derivative, -gamma * x + b * delayed, 2e-6, "Delayed exponential residual");
                if (t <= tau) {
                    const double integral = gamma == 0 ? t : -std::expm1(-gamma * t) / gamma;
                    near(x, h * std::exp(-gamma * t) + b * h * integral, 2e-13,
                         "First method-of-steps interval");
                }
            }
        }
    const double t = 3.4, b = .3;
    const double polynomial = 1 + b * t + b * b * std::pow(t - 1, 2) / 2 +
                              std::pow(b, 3) * std::pow(t - 2, 3) / 6 +
                              std::pow(b, 4) * std::pow(t - 3, 4) / 24;
    near(linearDelayConstantHistory(b, 0, 1, 2, t), 2 * polynomial, 2e-14,
         "Pure-delay finite polynomial");
    near(linearDelayConstantHistory(-1, 0, 2, 1, 1), 0, 1e-25,
         "Exact zero is not a precision failure");
    near(linearDelayConstantHistory(0, .3, 2, 1, 100), std::exp(-30), 1e-20, "Feedback-free decay");
    SystemConfiguration c = resolveMackeyGlass({});
    c.parameters["delay"] = 1;
    auto descriptor = mackeyGlassDescriptor();
    check(bool(descriptor.analyticSolutions), "DDE family factory missing");
    auto families = descriptor.analyticSolutions(c);
    auto comparison =
        compareAnalyticSolution(families, {1.2}, 4, {4, 4.5, 5, 5.1}, std::vector<double>(4));
    check(comparison.status == "partial" && comparison.valid == std::vector<double>({1, 1, 1, 0}),
          "First-delay validity was extrapolated");
    check(std::isnan(comparison.values.back()), "Invalid delay interval has a fabricated value");
    c.parameters["exponent"] = 0;
    c.parameters["beta"] = .3;
    comparison = compareAnalyticSolution(mackeyGlassAnalyticSolutions(c), {1.2}, 4, {4, 5, 6, 7},
                                         std::vector<double>(4));
    check(comparison.status == "complete", "Global linear delay family was not selected");
    c.parameters["exponent"] = 10;
    c.parameters["beta"] = .2;
    comparison = compareAnalyticSolution(mackeyGlassAnalyticSolutions(c), {1}, 0, {0, 50, 100},
                                         std::vector<double>(3));
    check(comparison.status == "complete" &&
              comparison.descriptor.identifier == "mackey-glass-constant-equilibrium",
          "Exact equilibrium history not recognized");
    check(!descriptor.research.papers.empty(), "DDE system bibliography missing");
}
} // namespace
int main() {
    try {
        transforms();
        heat();
        delay();
        std::cout << "Non-ODE analytic modes, exact delay residuals, validity and system ownership "
                     "passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
