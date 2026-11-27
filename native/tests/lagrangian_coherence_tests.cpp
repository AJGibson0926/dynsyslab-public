// FTLE fields and LCS extraction (core/LagrangianCoherence.hpp) on analytic flows: a linear saddle
// (FTLE exactly 1) and the double gyre of Shadden, Lekien and Marsden 2005 (Physica D 212:271-304).
#include <dynsyslab/core/LagrangianCoherence.hpp>

#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;
constexpr double pi = 3.141592653589793;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

// A planar (possibly time-dependent) vector field as an ODESystem.
class PlanarFlow final : public c::ODESystem {
public:
    using Field = std::function<void(double, const double*, double*)>;
    explicit PlanarFlow(Field f) : f_(std::move(f)) {
        descriptor_.identifier = "test_flow";
        descriptor_.label = "test flow";
        descriptor_.stateDimension = 2;
        descriptor_.stateNames = {"x", "y"};
        descriptor_.stateLabels = {"x", "y"};
        configuration_.identifier = "test_flow";
    }
    const c::SystemDescriptor& descriptor() const override { return descriptor_; }
    const c::SystemConfiguration& configuration() const override { return configuration_; }
    std::vector<double> initialState() const override { return {0, 0}; }
    void rhs(double t, const double* x, double* dx) const override { f_(t, x, dx); }
    void rhsWithInput(double t, const double* x, double, double* dx) const override { f_(t, x, dx); }
    std::vector<double> jacobian(double, const double*) const override { return std::vector<double>(4, 0); }
    std::vector<double> inputMatrix(const double*) const override { return {0, 0}; }
    bool invalid(double, const double*) const override { return false; }

private:
    Field f_;
    c::SystemDescriptor descriptor_;
    c::SystemConfiguration configuration_;
};

c::DomainGrid box(double x0, double x1, double y0, double y1, std::size_t nx, std::size_t ny) {
    auto domain = c::FunctionDomain::physicalState(2);
    auto grid = c::DomainGrid::box(domain, {x0, y0}, {x1, y1}, nx);
    grid.axes[1].points = ny;
    grid.validate();
    return grid;
}

// The double gyre: psi = A sin(pi f(x, t)) sin(pi y), f = eps sin(w t) x^2 + (1 - 2 eps sin(w t)) x.
PlanarFlow doubleGyre(double eps) {
    return PlanarFlow([eps](double t, const double* p, double* v) {
        const double A = 0.1, w = 2 * pi / 10, a = eps * std::sin(w * t), b = 1 - 2 * a;
        const double f = a * p[0] * p[0] + b * p[0], df = 2 * a * p[0] + b;
        v[0] = -pi * A * std::sin(pi * f) * std::cos(pi * p[1]);
        v[1] = pi * A * std::cos(pi * f) * std::sin(pi * p[1]) * df;
    });
}

void symmetricEigenTest() {
    const std::vector<double> a{4, 1, -2, 1, 3, 0.5, -2, 0.5, 1};
    std::vector<double> values, vectors;
    c::symmetricEigen(a, 3, values, vectors);
    require(values[0] <= values[1] && values[1] <= values[2], "ascending eigenvalues");
    for (std::size_t j = 0; j < 3; ++j)
        for (std::size_t r = 0; r < 3; ++r) {
            double av = 0;
            for (std::size_t k = 0; k < 3; ++k) av += a[r * 3 + k] * vectors[k * 3 + j];
            require(std::abs(av - values[j] * vectors[r * 3 + j]) < 1e-12, "A v = lambda v");
        }
}

void saddleTest() {
    PlanarFlow saddle([](double, const double* p, double* v) { v[0] = p[0]; v[1] = -p[1]; });
    const auto grid = box(-1, 1, -1, 1, 11, 11);
    for (const double T : {2.0, -2.0})
        for (const double spacing : {0.0, 1e-4}) {
            const auto field = c::cauchyGreenField(saddle, grid, {0, T, 1e-3, spacing});
            const auto sigma = field.ftle();
            for (const auto s : sigma) require(std::abs(s - 1) < 1e-8, "linear saddle FTLE is exactly 1");
            // The stretching direction: x forward, y backward.
            const auto p = grid.size() / 2;
            const double major = field.eigenvectors[p * 4 + (T > 0 ? 1 : 3)];
            require(std::abs(std::abs(major) - 1) < 1e-8, "major Cauchy-Green eigenvector along the unstable axis of the flow direction");
        }
}

void steadyGyreTest() {
    // eps = 0: x = 1 joins the saddles (1, 1) -> (1, 0); its stable side is the repelling LCS. A short
    // horizon keeps the shear of the closed orbits from raising ridges of its own.
    const auto gyre = doubleGyre(0);
    const auto grid = box(0, 2, 0, 1, 121, 61);
    const auto field = c::cauchyGreenField(gyre, grid, {0, 5, 0.02, 0});
    const auto sigma = field.ftle();
    const auto ridges = c::ridgeNodes(grid, sigma, 0.3);
    for (std::size_t j = 12; j < 49; ++j) {
        bool found = false;
        for (std::size_t i = 0; i < 121; ++i) found = found || (ridges[i + 121 * j] && std::abs(grid.axes[0].value(i) - 1) < 0.035);
        require(found, "steady double-gyre ridge on x = 1 in every interior row");
    }
    const auto curves = c::strainlines(field);
    require(!curves.empty(), "steady double gyre has a strainline");
    double meanX = 0, lo = 1, hi = 0;
    const auto& best = curves.front();
    for (std::size_t s = 0; s < best.points.size() / 2; ++s) {
        meanX += best.points[2 * s];
        lo = std::min(lo, best.points[2 * s + 1]);
        hi = std::max(hi, best.points[2 * s + 1]);
    }
    meanX /= double(best.points.size() / 2);
    require(std::abs(meanX - 1) < 0.02 && hi - lo > 0.6, "the strongest steady strainline is the separatrix x = 1");
}

void periodicGyreTest() {
    // eps = 0.25, w = 2 pi / 10, t0 = 0, |T| = 10. The flow is symmetric under (x, y, t) -> (2 - x,
    // 1 - y, -t), so the backward FTLE is the forward one reflected through (1, 1/2).
    const auto gyre = doubleGyre(0.25);
    const auto grid = box(0, 2, 0, 1, 121, 61);
    std::vector<double> forward;
    for (const double T : {10.0, -10.0}) {
        const auto field = c::cauchyGreenField(gyre, grid, {0, T, 0.02, 0});
        const auto sigma = field.ftle();
        double top = 0;
        for (const auto s : sigma) top = std::max(top, s);
        require(top > 0.3 && top < 0.6, "double-gyre FTLE of the published magnitude");
        if (T > 0) forward = sigma;
        else
            for (std::size_t j = 0; j < 61; ++j)
                for (std::size_t i = 0; i < 121; ++i)
                    require(std::abs(sigma[i + 121 * j] - forward[(120 - i) + 121 * (60 - j)]) < 1e-6 * top,
                            "backward FTLE is the reflected forward FTLE");
        const auto curves = c::strainlines(field);
        require(!curves.empty(), "time-periodic double gyre has LCS");
        double lo = 1, hi = 0, meanX = 0;
        const auto& best = curves.front();
        for (std::size_t s = 0; s < best.points.size() / 2; ++s) {
            lo = std::min(lo, best.points[2 * s + 1]);
            hi = std::max(hi, best.points[2 * s + 1]);
            meanX += best.points[2 * s];
        }
        meanX /= double(best.points.size() / 2);
        if (T < 0) meanX = 2 - meanX;
        require(hi - lo > 0.8 && meanX > 0.9 && meanX < 1.4, "the strongest LCS crosses the domain beside the gyre midline");
        const auto mask = c::rasterizeCurves(grid, curves);
        std::size_t marked = 0;
        for (const auto m : mask) marked += m;
        require(marked > 60, "rasterized LCS");
    }
}

}  // namespace

int main() {
    try {
        symmetricEigenTest();
        saddleTest();
        steadyGyreTest();
        periodicGyreTest();
        std::cout << "lagrangian coherence tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
