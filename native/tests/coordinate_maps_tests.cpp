// Domain grids (DomainGrid.hpp) and coordinate maps acting on functions (CoordinateMaps.hpp): grid
// ordering and status, the Koopman action g o Phi^t against the exact linear flow and a library
// eigenfunction, pull-back by composition, push-forward with conflict and no-preimage status, the
// round-trip error, and the fold set of a folding map.
#include <dynsyslab/core/CoordinateMaps.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;
using E = c::AnalyticExpression;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
c::Observable formula(const std::string& name, const std::vector<std::string>& texts, const std::vector<std::string>& variables) {
    std::vector<E> parts;
    for (const auto& t : texts) parts.push_back(E::parse(t, variables));
    return c::Observable::fromExpressions({name, name, "", "native", {}, false}, parts);
}
}  // namespace

int main() {
    try {
        // Grid ordering: the first axis varies fastest; a slice keeps the anchor.
        auto box = c::DomainGrid::box(c::FunctionDomain::physicalState(2), {0, 10}, {1, 12}, 3);
        double p[2];
        box.point(1, p);
        require(p[0] == 0.5 && p[1] == 10, "first axis fastest");
        box.point(3, p);
        require(p[0] == 0 && p[1] == 11, "second axis");
        require(box.locate(p) == 3, "locate inverts point");
        const auto field = c::evaluateOnGrid(formula("inv", {"1/x"}, {"x", "y"}), box);
        require(field.status[0] == c::GridStatus::undefined && field.count(c::GridStatus::undefined) == 3,
                "non-finite values are undefined");

        // Koopman action on the harmonic oscillator x'' = -2x: x o Phi^t = x cos(wt) + v sin(wt)/w.
        c::SystemConfiguration hc;
        hc.identifier = "harmonic_oscillator";
        hc.parameters = {{"mass_kg", 1.0}, {"stiffness_n_m", 2.0}, {"damping_kg_s", 0.0}};
        const auto harmonic = c::createSystem(hc);
        const auto grid = c::DomainGrid::box(c::FunctionDomain::physicalState(2), {-1, -1}, {1, 1}, 9);
        const auto position = formula("x", {"position"}, {"position", "velocity"});
        const auto frames = c::koopmanAction(position, *harmonic, grid, 0.25, 4, 1e-3);
        require(frames.size() == 5, "one frame per step and the start");
        const double w = std::sqrt(2.0), t = 1.0;
        for (std::size_t k = 0; k < grid.size(); ++k) {
            grid.point(k, p);
            const double exact = p[0] * std::cos(w * t) + p[1] * std::sin(w * t) / w;
            require(std::abs(frames[4].value(k).real() - exact) < 1e-9, "Koopman action matches the exact flow");
        }

        // A library eigenfunction under the Koopman action: psi o Phi^t = exp(lambda t) psi.
        const auto mode = c::libraryEigenfunctions(*harmonic, {"mode"}).front();
        const auto modeFrames = c::koopmanAction(mode, *harmonic, grid, 0.5, 2, 1e-3);
        const auto factor = mode.eigenvalue().transition(1.0);
        for (std::size_t k = 0; k < grid.size(); ++k)
            require(std::abs(modeFrames[2].value(k) - factor * modeFrames[0].value(k)) < 1e-9, "U^t psi = exp(lambda t) psi");

        // Pull-back by composition, with the exact chain-rule Jacobian.
        const auto h = formula("h", {"x^2", "y"}, {"x", "y"});
        const auto G = formula("G", {"u + 3*v"}, {"u", "v"});
        const auto pulled = c::pullBack(G, h);
        const double at[]{0.7, -0.2};
        require(std::abs(pulled.evaluate(at).front().real() - (0.49 - 0.6)) < 1e-12, "pull-back value");
        const auto J = pulled.jacobian(at, c::ObservableJacobianBackend::analytic);
        require(std::abs(J[0].real() - 1.4) < 1e-12 && std::abs(J[1].real() - 3) < 1e-12, "pull-back Jacobian");

        // Push-forward through h = (x^2, y) from [-1, 1]^2 onto [-1, 1] x [-1, 1].
        const auto source = c::DomainGrid::box(c::FunctionDomain::physicalState(2), {-1, -1}, {1, 1}, 41);
        const auto target = c::DomainGrid::box(c::FunctionDomain::supplied(2), {-1, -1}, {1, 1}, 21);
        // g = x^2 + y factors through h: consistent wherever reached; u < 0 has no preimage.
        const auto even = c::pushForward(formula("g", {"x^2 + y"}, {"x", "y"}), h, source, target, {0.06});
        const auto odd = c::pushForward(formula("g", {"x"}, {"x", "y"}), h, source, target, {0.06});
        std::size_t conflicts = 0, missing = 0;
        for (std::size_t k = 0; k < target.size(); ++k) {
            target.point(k, p);
            if (p[0] < -0.05) require(even.status[k] == c::GridStatus::noPreimage, "negative u has no preimage");
            if (even.status[k] == c::GridStatus::conflict) ++conflicts;
            if (odd.status[k] == c::GridStatus::noPreimage) ++missing;
        }
        require(conflicts == 0, "a function of h pushes forward consistently");
        require(odd.count(c::GridStatus::conflict) > target.size() / 4, "x is two-valued on the fibres of x^2");
        require(missing == even.count(c::GridStatus::noPreimage), "the image does not depend on g");

        // Round trip: zero where g factors through h, large where it does not.
        const auto back = c::roundTripError(formula("g", {"x"}, {"x", "y"}), h, source, target, {0.06});
        double worst = 0;
        for (std::size_t k = 0; k < source.size(); ++k)
            if (back.status[k] == c::GridStatus::conflict) worst = std::max(worst, back.value(k).real());
        require(worst > 0.5, "the round trip loses the sign of x");

        // Fold set: the volume factor of (x^2, y) is 2|x|, zero on x = 0.
        const auto fold = c::mapVolumeFactor(h, source);
        for (std::size_t k = 0; k < source.size(); ++k) {
            source.point(k, p);
            require(std::abs(fold.value(k).real() - 2 * std::abs(p[0])) < 1e-9, "volume factor 2|x|");
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "coordinate maps tests passed\n";
    return 0;
}
