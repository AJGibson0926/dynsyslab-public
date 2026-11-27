// Model library families (core/ModelLibraries.hpp): sizes, the localized partition of unity, the
// harmonic phase basis against complex arithmetic, and a rational law recovered exactly by SINDy.
#include <dynsyslab/core/ModelLibraries.hpp>
#include <dynsyslab/core/SindyDynamics.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace c = dynsyslab::core;

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

c::TimeSeries training(std::size_t n) {
    std::vector<double> t, v;
    for (std::size_t k = 0; k < 400; ++k) {
        t.push_back(0.01 * double(k));
        for (std::size_t j = 0; j < n; ++j) v.push_back(2.5 * std::sin(0.03 * double(k) * double(j + 1) + double(j)));
    }
    return c::TimeSeries(t, v, n);
}

}  // namespace

int main() {
    try {
        const auto data = training(2);
        c::ModelLibraryOptions o;
        require(c::modelLibrary("polynomial", data, o).outputDimension() == 6, "1 + 5 monomials of degree <= 2 in two variables");
        require(c::modelLibrary("bounded", data, o).outputDimension() == 8, "1, x and the bounded monomials");
        require(c::modelLibrary("rational", data, o).outputDimension() == 8 + 12, "two terms per width and coordinate");
        require(c::modelLibrary("gated", data, o).outputDimension() == 10, "five gates times two coordinates");
        require(c::modelLibrary("paired_gated", data, o).outputDimension() == 12, "plus one pair gate");
        o.centers = 6;
        const auto localized = c::modelLibrary("localized", data, o);
        require(localized.outputDimension() == 8 + 6 * 3, "weights times [1, z] per center");
        const double q[2] = {0.4, -1.1};
        const auto values = localized.evaluateReal(q);
        double sum = 0;
        for (std::size_t a = 0; a < 6; ++a) sum += values[8 + 3 * a];
        require(std::abs(sum - 1) < 1e-12, "the localized weights are a partition of unity");
        o.harmonics = 3;
        const auto harmonic = c::modelLibrary("harmonic", data, o);
        require(harmonic.outputDimension() == 2 + 3 * 5 * 2, "x and gated Re/Im u^k");
        const auto h = harmonic.evaluateReal(q);
        const double p = q[0] * q[0] + q[1] * q[1], e = p / 2;
        const auto u = std::complex<double>(q[0], q[1]) / std::sqrt(p + 0.015);
        const auto u3 = u * u * u;
        const double amplitude = std::sqrt(e + 0.015);  // ungated (gate 1) column of k = 3
        require(std::abs(h[2 + 2 * 5 * 2] - amplitude * u3.real()) < 1e-12 && std::abs(h[2 + 2 * 5 * 2 + 1] - amplitude * u3.imag()) < 1e-12,
                "harmonic basis matches complex powers");

        // x' = -x / (1 + x^2) is the width-1 odd regularized reciprocal: SINDy finds it alone.
        std::vector<double> t, x, dx;
        for (int k = 0; k < 300; ++k) { const double s = -3 + 0.02 * k; t.push_back(0.01 * k); x.push_back(s); dx.push_back(-s / (1 + s * s)); }
        const c::TimeSeries qs(t, x, 1), rates(t, dx, 1);
        c::ModelLibraryOptions r;
        r.rationalWidths = {0.25, 1.0, 4.0};
        const auto library = c::modelLibrary("rational", qs, r);
        c::SindyDynamicsOptions options;
        options.threshold = 0.01;
        const auto model = c::fitSindyDynamics(qs, rates, library, options);
        std::size_t active = 0;
        for (std::size_t f = 0; f < library.outputDimension(); ++f)
            if (model.coefficients()[f] != 0) {
                ++active;
                require(library.descriptor().outputNames[f] == "odd_1.000000_x1" && std::abs(model.coefficients()[f] + 1) < 1e-8, "the rational term");
            }
        require(active == 1, "a single active rational term");
        std::cout << "model libraries tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
