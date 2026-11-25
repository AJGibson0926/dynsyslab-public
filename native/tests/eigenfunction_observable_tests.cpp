// Observable and its subclass Eigenfunction: the representation property (expansion or
// interpolation), the eigenvalue on the Eigenfunction, exact gradients and flow derivatives,
// projection of observables onto eigenfunctions through the shared dictionary, and a system's
// eigenfunction library (exact and linearized principal entries).
#include <dynsyslab/core/ExpansionProjection.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/KnownEigenfunctions.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;
using Complex = std::complex<double>;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
bool near(Complex a, Complex b, double tolerance = 1e-9) { return std::abs(a - b) <= tolerance * (1 + std::abs(b)); }

c::SystemConfiguration triangularConfiguration(double a, double b, double cc) {
    c::SystemConfiguration config;
    config.identifier = "quadratic_triangular";
    config.parameters = {{"a", a}, {"b", b}, {"c", cc}};
    return config;
}

c::FunctionDomain stateDomain() {
    auto d = c::FunctionDomain::physicalState(2);
    d.names = {"x_1", "x_2"};
    return d;
}

c::FunctionExpansion expressions(const std::vector<std::string>& texts) {
    std::vector<c::AnalyticExpression> parsed;
    for (const auto& t : texts) parsed.push_back(c::AnalyticExpression::parse(t, {"x_1", "x_2"}));
    return c::FunctionExpansion::fromExpressions(parsed, stateDomain());
}

}  // namespace

int main() {
    try {
        const double a = -1, b = -3, cc = 1, k = cc / (b - 2 * a);  // k = -1
        const auto config = triangularConfiguration(a, b, cc);
        const auto system = c::createSystem(config);
        const double x[]{0.7, -0.4};

        // Library: Eigenfunction objects (Observable + Eigenvalue) over one shared expansion.
        const auto library = c::libraryEigenfunctions(*system);
        require(library.size() == 4, "library entry count");
        const auto& eta = library[1];
        require(eta.descriptor().name == "eta", "library eigenfunctions are named by entry id");
        require(near(eta.eigenvalue().continuous, b), "eta's eigenvalue is b");
        require(eta.representationKind() == "expansion", "library eigenfunctions are expansions");
        require(near(eta.evaluate(x).front(), x[1] + k * x[0] * x[0]), "eta value");
        require(near(eta.propagate(2.0, 0.5), 2.0 * std::exp(b * 0.5)), "eigenfunction evolves by exp(lambda t)");

        // An Eigenfunction is an Observable: exact gradient through the base class.
        const c::Observable& asObservable = eta;
        require(asObservable.hasAnalyticJacobian(), "expansion eigenfunction has an exact Jacobian");
        const auto g = asObservable.jacobian(x, c::ObservableJacobianBackend::analytic);
        require(near(g[0], 2 * k * x[0]) && near(g[1], 1.0), "eta gradient");

        // Generator identity against the configured equations: exact entries vanish; the
        // linearized principal entry for b (x_2) is exact only at the origin.
        const std::vector<double> states{0.7, -0.4, -1.2, 0.9, 0.3, 0.2, 2.0, -1.5};
        const auto residual = c::eigenfunctionGeneratorResidual(*system, library, states);
        require(residual[0] < 1e-12 && residual[1] < 1e-12 && residual[2] < 1e-12, "exact entries satisfy L psi = lambda psi");
        require(residual[3] > 1e-3, "principal linearization x_2 is not exact away from the origin");
        const std::vector<double> nearOrigin{1e-4, 2e-4, -2e-4, 1e-4};
        require(c::eigenfunctionGeneratorResidual(*system, library, nearOrigin)[3] < 1e-3, "principal entry is accurate near x*");

        // The closed form as an expression Observable: its flow derivative along the drift is
        // exactly b eta.
        const auto bank = c::libraryEigenfunctionBank(*system, {"x_1", "eta"});
        const auto closed = c::expressionObservable(bank.functions, {"psi", "psi"});
        require(closed.has_value(), "polynomial library has closed forms");
        const auto drift = system->driftExpressions();
        if (!drift.empty()) {
            const auto lie = closed->flowDerivative(drift).evaluateReal(x);
            const auto value = closed->evaluateReal(x);
            require(near(lie[2], b * value[2], 1e-10), "symbolic flow derivative of eta");
        }

        // A K-ROM built from the bank hands out the same Eigenfunction objects.
        const c::KoopmanModel model(bank, config);
        const auto psi1 = model.eigenfunction(1);
        require(near(psi1.evaluate(x).front(), eta.evaluate(x).front()) && near(psi1.eigenvalue().continuous, b),
                "K-ROM eigenfunction matches the library entry");

        // Projection through the shared dictionary: 2 x_1 + 3 x_2 - 3 x_1^2 = 2 x_1 + 3 eta exactly.
        const auto target = expressions({"2*x_1 + 3*x_2 - 3*x_1^2"});
        const auto observable = c::Observable::fromRepresentation({"g", "g"}, std::make_shared<const c::FunctionExpansion>(target));
        const std::vector<c::Eigenfunction> pair{library[0], library[1]};
        const auto projected = c::projectObservable(observable, pair);
        require(projected.exact, "observable in the eigenfunction span projects exactly");
        require(near(projected.coefficients[0], 2.0) && near(projected.coefficients[1], 3.0), "analytic projection coefficients");

        // x_2 alone is not in span{x_1, eta}: the residual is the missing -k x_1^2 term.
        const auto partial = c::projectExpansion(expressions({"x_2"}), bank.functions);
        require(!partial.exact && partial.relativeResidual[0] > 0.1, "out-of-span observable reports its residual");
        require(near(partial.residual.evaluate(0, x) + partial.coefficients[0] * x[0] + partial.coefficients[1] * (x[1] + k * x[0] * x[0]),
                     x[1]), "residual completes the projection");

        // Normalized monomials canonicalize so that equal functions share terms.
        const c::FunctionExpansion scaled(c::PolynomialDictionary{2, {1, 0, 0, 1, 2, 0}, {0.5, 0.0}, {2.0, 1.0}},
                                          {Complex(2.0), 0.0, 0.0, 0.0, 1.0, Complex(-4.0)}, stateDomain());
        const auto viaScaled = c::projectExpansion(target, scaled);
        require(viaScaled.innerProduct == "coefficients", "coefficient inner product by default");
        require(near(viaScaled.residual.evaluate(0, x),
                     target.evaluate(0, x) - viaScaled.coefficients[0] * scaled.evaluate(0, x) - viaScaled.coefficients[1] * scaled.evaluate(1, x),
                     1e-9), "canonical residual is an identity of functions");

        // Sampled inner product agrees on an exact span.
        const auto sampled = c::projectExpansion(target, bank.functions, states);
        require(sampled.exact && near(sampled.coefficients[1], 3.0), "sampled projection");

        // The other representation: an interpolated eigenfunction is the same kind of object.
        c::InterpolationData samples;
        samples.trainingStates = states;
        for (std::size_t s = 0; s < 4; ++s) samples.trainingValues.push_back(eta.evaluate(states.data() + 2 * s).front());
        samples.samples = 4;
        samples.stateDimension = 2;
        samples.neighbors = 1;
        const c::KoopmanModel interpolated(samples, c::Spectrum::fromContinuousRates({b}), config);
        const auto sampledEta = interpolated.eigenfunction(0);
        require(sampledEta.representationKind() == "interpolation", "interpolation representation property");
        require(near(sampledEta.evaluate(states.data()).front(), eta.evaluate(states.data()).front()), "interpolated value");
        require(near(sampledEta.eigenvalue().continuous, b), "interpolated eigenfunction keeps its eigenvalue");

        // Resonant configuration: eta is not applicable and is rejected by name.
        const auto resonant = c::createSystem(triangularConfiguration(-1, -2, 1));
        bool rejected = false;
        try { (void)c::libraryEigenfunctions(*resonant, {"eta"}); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "non-applicable entry is rejected");
        require(c::libraryEigenfunctions(*resonant).size() == 3, "default library skips non-applicable entries");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "eigenfunction observable tests passed\n";
    return 0;
}
