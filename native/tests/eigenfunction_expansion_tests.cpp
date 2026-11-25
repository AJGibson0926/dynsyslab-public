// The expansion dictionary: exact conversion of every spectral family, gradients of a mixed
// dictionary, the expression parser, and a K-ROM archive round trip with a mixed dictionary.
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/EigenfunctionFit.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
bool near(double a, double b, double tolerance = 1e-9) { return std::abs(a - b) <= tolerance * (1 + std::abs(b)); }

std::vector<double> finiteGradient(const c::ExpansionDictionary& d, const double* q) {
    std::vector<double> out(d.size() * d.dimension()), point(q, q + d.dimension());
    for (std::size_t i = 0; i < d.dimension(); ++i) {
        const double h = 1e-6;
        point[i] = q[i] + h;
        const auto high = d.evaluate(point.data());
        point[i] = q[i] - h;
        const auto low = d.evaluate(point.data());
        point[i] = q[i];
        for (std::size_t k = 0; k < d.size(); ++k) out[k * d.dimension() + i] = (high[k] - low[k]) / (2 * h);
    }
    return out;
}

}  // namespace

int main() {
    try {
        const double q[]{1.3, 0.4};

        // Every spectral dictionary family converts to dictionary terms with identical values.
        for (const char* kind : {"linear", "polynomial", "chebyshev", "laguerre", "hermite", "fourier_hermite",
                                 "gaussian", "laplacian", "lorentzian", "matern", "partition"}) {
            c::SpectralDictionary s;
            s.kind = kind; s.dimension = 2; s.center = {0.2, -0.1}; s.scale = {1.5, 0.8}; s.bandwidth = 0.7;
            if (std::string(kind) == "linear") s.constant = true;
            else if (std::string(kind) == "polynomial" || std::string(kind) == "chebyshev" || std::string(kind) == "laguerre"
                     || std::string(kind) == "hermite" || std::string(kind) == "fourier_hermite")
                s.exponents = {0, 0, 1, 0, 0, 1, 2, 1, 1, 3};
            else s.nodes = {0, 0, 1, -1, 2, 0.5};
            if (std::string(kind) == "fourier_hermite") { s.center.clear(); s.scale.clear(); }  // raw angle by contract
            const auto expected = s.evaluate({q[0], q[1]});
            const auto converted = c::ExpansionDictionary::fromSpectral(s).evaluate(q);
            require(expected.size() == converted.size(), std::string(kind) + ": term count");
            for (std::size_t k = 0; k < expected.size(); ++k)
                require(expected[k] == converted[k], std::string(kind) + ": values differ");
        }

        // A mixed dictionary: monomial, exponential, cosine, gaussian and an expression.
        std::vector<c::DictionaryTerm> terms(5);
        terms[0].kind = "monomial"; terms[0].orders = {2, 1}; terms[0].center = {0.5, -1}; terms[0].scale = {2, 1};
        terms[1].kind = "exponential"; terms[1].vector = {-0.3, 0.7};
        terms[2].kind = "cosine"; terms[2].vector = {1.5, -2}; terms[2].phase = 0.2;
        terms[3].kind = "gaussian"; terms[3].vector = {1, 0}; terms[3].width = 0.6;
        terms[4].kind = "expression";
        terms[4].expression = std::make_shared<const c::AnalyticExpression>(
            c::AnalyticExpression::parse("r^2*cos(v) - 0.5*exp(-v) + sqrt(r)", {"r", "v"}));
        const auto mixed = std::make_shared<const c::ExpansionDictionary>(2, terms);
        require(mixed->exactGradient() && !mixed->monomialTable(), "mixed dictionary flags");
        const auto phi = mixed->evaluate(q);
        require(near(phi[1], std::exp(-0.3 * 1.3 + 0.7 * 0.4)), "exponential term");
        require(near(phi[4], 1.69 * std::cos(0.4) - 0.5 * std::exp(-0.4) + std::sqrt(1.3)), "expression term");
        const auto exact = mixed->gradient(q), numeric = finiteGradient(*mixed, q);
        for (std::size_t k = 0; k < exact.size(); ++k) require(near(exact[k], numeric[k], 1e-6), "mixed gradient");
        bool refused = false;
        try { (void)c::AnalyticExpression::parse("r + w", {"r", "v"}); } catch (const std::invalid_argument&) { refused = true; }
        require(refused, "unknown variables are refused");

        // A K-ROM with the mixed dictionary survives the archive.
        std::vector<std::complex<double>> coefficients{{1, 0.5}, {-0.2, 0.1}, {0.3, -0.4}, {0.1, 0}, {0, 0.2},
                                                       {1, -0.5}, {-0.2, -0.1}, {0.3, 0.4}, {0.1, 0}, {0, -0.2}};
        c::EigenfunctionFitResult bank{c::FunctionExpansion(mixed, coefficients, c::FunctionDomain::physicalState(2)), c::Spectrum::fromContinuousRates({{-0.1, 1.0}, {-0.1, -1.0}}), {}};
        c::SystemConfiguration duffing;
        duffing.identifier = "duffing"; duffing.preset = "double_well_reference";
        const c::KoopmanModel model(bank, duffing);
        require(model.eigenfunctionRepresentation() == "expansion", "representation name");
        const auto path = std::filesystem::temp_directory_path()
            / ("dynsyslab-expansion-" + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()) + ".npz");
        c::saveKoopmanModel(path, model, "{}");
        const auto loaded = c::loadKoopmanModel(path);
        std::filesystem::remove(path);
        const auto before = model.evaluate(q), after = loaded.evaluate(q);
        for (std::size_t j = 0; j < 2; ++j) require(std::abs(before[j] - after[j]) < 1e-12, "mixed K-ROM round trip");
        std::complex<double> expected{};
        for (std::size_t k = 0; k < 5; ++k) expected += coefficients[k] * phi[k];
        require(std::abs(before[0] - expected) < 1e-12, "expansion is the coefficient-weighted dictionary");
        std::cout << "eigenfunction expansion tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
