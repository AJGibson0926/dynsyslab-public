#include <algorithm>
#include <cmath>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include <dynsyslab/core/PDE.hpp>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
AnalyticTrajectory makeHeat(const SystemConfiguration &c, const AnalyticState &initial, double t0,
                            bool continuum) {
    const auto mesh = heatEquationMesh(c);
    const auto n = mesh.shape.front();
    const auto boundary = static_cast<int>(c.parameters.at("boundary_code"));
    if (initial.size() != n || !std::isfinite(t0) ||
        !std::all_of(initial.begin(), initial.end(), [](double x) { return std::isfinite(x); }))
        return {{}, "Initial field/time must be finite and match the configured mesh"};
    if (c.driver.amplitude != 0)
        return {{}, "Heat analytic families require no external driver"};
    if (boundary == 2 && (initial.front() != 0 || initial.back() != 0))
        return {{}, "Zero-value boundaries require zero initial endpoint values"};
    // Odd extension yields a sine interpolant; even extension yields a cosine
    // interpolant. The edge stencil's ghosts reflect about half-grid boundaries.
    const bool halfGrid = boundary == 1 && !continuum;
    const auto m = boundary == 0 ? n : halfGrid ? 2 * n : 2 * (n - 1);
    std::vector<std::complex<double>> extended(m);
    for (std::size_t j = 0; j < n; ++j)
        extended[j] = initial[j];
    if (halfGrid) {
        for (std::size_t j = 0; j < n; ++j)
            extended[m - 1 - j] = initial[j];
    } else if (boundary != 0) {
        for (std::size_t j = 1; j + 1 < n; ++j)
            extended[m - j] = (boundary == 2 ? -1.0 : 1.0) * initial[j];
    }
    auto coefficients = discreteFourierTransform(extended);
    std::vector<double> rates(m);
    const double length = c.parameters.at("length"), diffusivity = c.parameters.at("diffusivity"),
                 decay = c.parameters.at("decay_rate"), spacing = mesh.spacing.front();
    for (std::size_t k = 0; k < m; ++k) {
        const double mode = static_cast<double>(std::min(k, m - k));
        const double w = continuum ? (boundary == 0 ? 2 * pi : pi) * mode / length
                                   : 2 * std::sin(pi * mode / m) / spacing;
        rates[k] = diffusivity * w * w + decay;
    }
    return {[coefficients = std::move(coefficients), rates = std::move(rates), initial, n, boundary,
             t0](double time) {
                if (!std::isfinite(time) || time < t0)
                    throw std::domain_error(
                        "Heat analytic trajectory requires finite time at or after its IC");
                if (time == t0)
                    return initial;
                auto evolved = coefficients;
                for (std::size_t k = 0; k < evolved.size(); ++k)
                    evolved[k] *= std::exp(-rates[k] * (time - t0));
                evolved = discreteFourierTransform(evolved, true);
                AnalyticState result(n);
                for (std::size_t j = 0; j < n; ++j)
                    result[j] = evolved[j].real();
                if (boundary == 2)
                    result.front() = result.back() = 0;
                return result;
            },
            {}};
}
} // namespace
std::vector<AnalyticSolutionFamily>
heatEquationAnalyticSolutions(const SystemConfiguration &requested) {
    const auto c = resolveHeatEquation(requested);
    const std::string reference = "https://math.mit.edu/~dyatlov/18.03/lectures/Section10-1.pdf";
    AnalyticSolutionDescriptor continuum{
        "heat-continuum-fourier",
        "Analytic heat equation",
        R"(u(x,t)=\sum_k a_k\phi_k(x)e^{-(\kappa q_k^2+c)(t-t_0)})",
        "Constant diffusivity/decay; periodic, zero-flux or homogeneous zero-value boundaries. "
        "Initial field is the finite trigonometric interpolant of the supplied mesh samples.",
        "Exact continuum evolution of the specified interpolant, evaluated by FFT. Periodic "
        "modes use q=2*pi*k/L; Neumann/Dirichlet modes use q=pi*k/L. Errors versus the numerical "
        "trajectory include spatial and temporal discretization. Initial interpolation error "
        "relative to an unspecified continuous field is not measured.",
        {reference}};
    AnalyticSolutionDescriptor discrete{
        "heat-discrete-semigroup",
        "Analytic spatial discretization",
        R"(u_h(t)=e^{A_h(t-t_0)}u_h(t_0),\quad\lambda_k=-4\kappa h^{-2}\sin^2(\pi k/M)-c)",
        "The configured central-difference stencil; periodic, clamped-ghost zero-flux, or "
        "homogeneous zero-value boundaries; constant coefficients.",
        "Exact time evolution of the spatially discretized operator, evaluated by FFT. This "
        "is not a continuum PDE reference. The clamped-ghost flux stencil reflects at -h/2 "
        "and L+h/2, and has cosine modes with effective interval length N*h.",
        {reference}};
    return {
        {std::move(continuum), [c](const auto &x, double t0) { return makeHeat(c, x, t0, true); }},
        {std::move(discrete), [c](const auto &x, double t0) { return makeHeat(c, x, t0, false); }}};
}
} // namespace dynsyslab::core
