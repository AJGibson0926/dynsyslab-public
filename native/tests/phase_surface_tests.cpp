#include "../src/PhaseSurfaceInternal.hpp"
#include <algorithm>
#include <cmath>
#include <dynsyslab/core/PhaseSpace.hpp>
#include <iostream>
#include <stdexcept>
using namespace dynsyslab::core;
void check(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
class CurvedSurface final : public ODESystem {
    SystemDescriptor d_;
    SystemConfiguration c_;
    double sign_;
    bool rotation_;

  public:
    CurvedSurface(double sign, bool rotation, std::size_t n) : sign_(sign), rotation_(rotation) {
        d_.identifier = "curved-surface-test";
        d_.stateDimension = n;
        c_.identifier = d_.identifier;
    }
    const SystemDescriptor &descriptor() const override { return d_; }
    const SystemConfiguration &configuration() const override { return c_; }
    std::vector<double> initialState() const override {
        return std::vector<double>(d_.stateDimension);
    }
    void rhs(double, const double *x, double *f) const override {
        const double a = rotation_ ? 1. : .7, b = rotation_ ? 1. : 1.3, w = rotation_ ? .8 : 0.;
        f[0] = sign_ * (a * x[0] - 2 * w * x[1]);
        f[1] = sign_ * (.5 * w * x[0] + b * x[1]);
        const double h = x[0] * x[0] + .5 * x[0] * x[1] + x[1] * x[1];
        f[2] = -sign_ * 1.7 * (x[2] - h) + (2 * x[0] + .5 * x[1]) * f[0] +
               (.5 * x[0] + 2 * x[1]) * f[1];
        if (d_.stateDimension == 4)
            f[3] = -sign_ * 2 * (x[3] - x[0] * x[1]) + x[1] * f[0] + x[0] * f[1];
    }
    void rhsWithInput(double t, const double *x, double, double *f) const override { rhs(t, x, f); }
    std::vector<double> jacobian(double t, const double *x) const override {
        const auto n = d_.stateDimension;
        std::vector<double> result(n * n), a(x, x + n), b = a, fa(n), fb(n);
        for (std::size_t j = 0; j < n; ++j) {
            a[j] += 1e-5;
            b[j] -= 1e-5;
            rhs(t, a.data(), fa.data());
            rhs(t, b.data(), fb.data());
            for (std::size_t i = 0; i < n; ++i)
                result[i * n + j] = (fa[i] - fb[i]) / 2e-5;
            a[j] = b[j] = x[j];
        }
        return result;
    }
    std::vector<double> inputMatrix(const double *) const override {
        return std::vector<double>(d_.stateDimension);
    }
    bool invalid(double, const double *x) const override {
        for (std::size_t j = 0; j < d_.stateDimension; ++j)
            if (!std::isfinite(x[j]) || std::abs(x[j]) > 20)
                return true;
        return false;
    }
};
PhaseSpaceOptions options(std::size_t n) {
    PhaseSpaceOptions o;
    o.lower.assign(n, -1.);
    o.upper.assign(n, 1.);
    o.equilibriumSeeds = {PhaseState(n)};
    o.periodicOrbits = false;
    o.vectorField = false;
    o.horizon = 10;
    o.outputStep = .1;
    o.surface.maximumDistance = .12;
    o.surface.maximumFronts = 48;
    o.surface.chordTolerance = .0005;
    return o;
}
void oracle(double sign, bool rotation, std::size_t n) {
    CurvedSurface system(sign, rotation, n);
    auto o = options(n);
    if (rotation) { o.surface.step = .01; o.surface.meshSize = .02; }
    const auto result = analyzePhaseSpace(system, o);
    check(result.manifoldSurfaces.size() >= 1, "Curved surface missing");
    const auto selected = std::find_if(
        result.manifoldSurfaces.begin(), result.manifoldSurfaces.end(), [&](const auto &m) {
            return m.kind == (sign > 0 ? "unstable_manifold_surface" : "stable_manifold_surface");
        });
    check(selected != result.manifoldSurfaces.end(), "Requested invariant surface missing");
    const auto &m = *selected;
    std::cout << "oracle " << sign << " " << rotation << " " << n << " " << m.status << " "
              << m.reason << " fronts=" << m.frontOffsets.size()
              << " vertices=" << m.vertices.size() / n << " error=" << m.refinementError
              << " tangent=" << m.tangentDefect << std::endl;
    check(m.status != "unresolved" && m.frontDistances.back() > .08,
          "Surface continuation failed validation or extent");
    check(!m.triangles.empty() && m.triangles.size() % 3 == 0, "Triangulated surface absent");
    double curvature = 0;
    for (std::size_t i = 0; i < m.vertices.size(); i += n) {
        const double u = m.vertices[i], v = m.vertices[i + 1];
        check(std::abs(m.vertices[i + 2] - u * u - .5 * u * v - v * v) < .004,
              "Surface violates exact invariant graph");
        if (n == 4)
            check(std::abs(m.vertices[i + 3] - u * v) < .004, "Full-state surface embedding lost");
        curvature = std::max(curvature, std::abs(m.vertices[i + 2]));
    }
    check(curvature > .01, "Surface is only a tangent plane");
    for (auto index : m.triangles)
        check(index < m.vertices.size() / n, "Invalid mesh index");
    o.surface.enabled = false;
    check(analyzePhaseSpace(system, o).manifoldSurfaces.empty(), "Surface disable flag ignored");
}
void lorenz() {
    SystemConfiguration c;
    c.identifier = "lorenz";
    auto system = createSystem(c);
    auto o = options(3);
    o.lower = {-10, -10, -5};
    o.upper = {10, 10, 10};
    o.surface.maximumDistance = .06;
    o.surface.step = .008;
    const auto result = analyzePhaseSpace(*system, o);
    check(!result.manifoldSurfaces.empty(), "Lorenz stable surface absent");
    const auto &m = result.manifoldSurfaces[0];
    std::cout << "Lorenz " << m.status << " " << m.reason << " fronts=" << m.frontOffsets.size()
              << " error=" << m.refinementError << " tangent=" << m.tangentDefect << std::endl;
    check(m.kind == "stable_manifold_surface" && m.status != "unresolved" &&
              m.frontDistances.back() > .04,
          "Lorenz stable surface failed");
}
int main() {
    try {
        const auto wrapped =
            phase_detail::surface::interpolate({{1., 0.}, {0., 1.}, {-1., 0.}}, -1e-18);
        check(wrapped.size() == 2 && wrapped[0] == 1.,
              "Negative roundoff ring parameter escaped the mesh");
        oracle(1, false, 3);
        oracle(-1, false, 3);
        oracle(1, true, 4);
        lorenz();
        CurvedSurface cancelledSystem(1, false, 3);
        auto cancelledOptions = options(3);
        ExecutionContext context;
        std::size_t probes = 0;
        context.cancellation = [&]() { return ++probes > 100; };
        bool cancelled = false;
        try {
            (void)phase_detail::surface::grow(cancelledSystem, {0., 0., 0.},
                                              {{1., 0., 0.}, {0., 1., 0.}}, 1., cancelledOptions,
                                              &context);
        } catch (const ExecutionCancelled &) {
            cancelled = true;
        }
        check(cancelled, "Surface continuation swallowed cancellation");
        std::cout << "Phase surface tests passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
