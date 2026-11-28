// The expression language's special functions and series (ExpressionFunctions.cpp, AnalyticExpression.hpp):
// values against independent references (DLMF), defining identities, exact gradients against central
// differences, and JSON round trips.
#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/EllipticIntegrals.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace c = dynsyslab::core;
const double pi = std::acos(-1.0);

void require(bool ok, const std::string& what) { if (!ok) throw std::runtime_error(what); }
double value(const std::string& text, double x = 0, double y = 0) {
    const double p[]{x, y};
    return c::AnalyticExpression::parse(text, {"x", "y"}).evaluate(p);
}
void near(const std::string& text, double expected, double tolerance, double x = 0, double y = 0) {
    const double v = value(text, x, y);
    require(std::abs(v - expected) <= tolerance * std::max(1.0, std::abs(expected)),
            text + " = " + std::to_string(v) + ", expected " + std::to_string(expected));
}
// The exact gradient of `text` at (x, y) against central differences, with its JSON round trip.
void gradient(const std::string& text, double x, double y, double tolerance = 2e-7) {
    const auto e = c::AnalyticExpression::parse(text, {"x", "y"});
    const double q[]{x, y};
    for (std::size_t axis = 0; axis < 2; ++axis) {
        double a[]{x, y}, b[]{x, y};
        const double h = 1e-5 * std::max(1.0, std::abs(q[axis]));
        a[axis] += h; b[axis] -= h;
        const double numeric = (e.evaluate(a) - e.evaluate(b)) / (2 * h), exact = e.derivative(axis).evaluate(q);
        require(std::abs(exact - numeric) <= tolerance * std::max(1.0, std::abs(numeric)),
                "gradient of " + text + " in axis " + std::to_string(axis) + ": " + std::to_string(exact) + " vs " + std::to_string(numeric));
    }
    const auto back = c::AnalyticExpression::fromJson(e.toJson());
    require(std::abs(back.evaluate(q) - e.evaluate(q)) <= 1e-14 * std::max(1.0, std::abs(e.evaluate(q))), "JSON round trip of " + text);
}

}  // namespace

int main() {
    try {
        // References: DLMF tables and closed forms.
        near("besselj(0, 1)", 0.7651976865579666, 1e-14);
        near("bessely(0, 1)", 0.08825696421567696, 1e-13);
        near("besseli(1, 1)", 0.5651591039924851, 1e-14);
        near("besselk(0, 1)", 0.42102443824070834, 1e-14);
        near("sphbesselj(1, 1)", std::sin(1.0) - std::cos(1.0), 1e-14);
        near("airyai(0)", 0.3550280538878172, 1e-14);
        near("airybi(0)", 0.6149266274460007, 1e-14);
        near("gamma(0.5)", std::sqrt(pi), 1e-14);
        near("digamma(1)", -0.5772156649015329, 1e-14);
        near("zeta(2)", pi * pi / 6, 1e-14);
        near("erf(1)", 0.8427007929497149, 1e-14);
        near("expint(1)", 1.8951178163559368, 1e-14);
        near("lambertw(1)", 0.5671432904097838, 1e-14);
        near("beta(2, 3)", 1.0 / 12, 1e-14);
        near("gammap(1, 2)", 1 - std::exp(-2.0), 1e-14);
        near("hyp1f1(1, 1, 0.7)", std::exp(0.7), 1e-13);
        near("hyp2f1(1, 1, 2, 0.5)", -std::log(0.5) / 0.5, 1e-13);
        near("hyp0f1(0.5, -0.36)", std::cos(1.2), 1e-13);  // 0F1(; 1/2; -z^2/4) = cos z
        near("legendre(3, 0.4)", (5 * 0.064 - 3 * 0.4) / 2, 1e-15);
        near("hermite(3, 0.4)", 8 * 0.064 - 12 * 0.4, 1e-15);
        near("laguerre(2, 0.4)", (0.16 - 1.6 + 2) / 2, 1e-15);
        near("chebyshevt(4, 0.4)", std::cos(4 * std::acos(0.4)), 1e-14);
        near("tan(0.3) + asinh(0.3) + atanh(0.3)", std::tan(0.3) + std::asinh(0.3) + std::atanh(0.3), 1e-15);

        // Jacobi elliptic functions: sn^2 + cn^2 = 1, dn^2 + m sn^2 = 1, sn(K|m) = 1, sn(u|0) = sin u.
        near("sn(x, 0.6)^2 + cn(x, 0.6)^2", 1, 1e-14, 2.3);
        near("dn(x, 0.6)^2 + 0.6*sn(x, 0.6)^2", 1, 1e-14, 7.1);
        near("sn(x, 0.6)", 1, 1e-13, c::ellipticK(0.6));
        near("sn(x, 0)", std::sin(1.3), 1e-15, 1.3);
        near("jacobiam(x, 0.6)", pi / 2 + 3 * pi, 1e-13, 7 * c::ellipticK(0.6));
        // Weierstrass p: p'^2 = 4p^3 - g2 p - g3 (both signs of the discriminant), the derivative is p',
        // and p ~ 1/z^2 + g2 z^2/20 + g3 z^4/28 near 0.
        for (const auto& [g2, g3] : std::vector<std::pair<double, double>>{{4, 1}, {1, 2}}) {
            const auto args = std::to_string(g2) + ", " + std::to_string(g3);
            const double p = value("weierstrassp(" + args + ", x)", 0.37), dp = value("weierstrasspprime(" + args + ", x)", 0.37);
            require(std::abs(dp * dp - (4 * p * p * p - g2 * p - g3)) < 1e-10 * p * p * p, "Weierstrass differential equation");
            const auto wp = c::AnalyticExpression::parse("weierstrassp(" + args + ", x)", {"x"});
            const double z[]{0.37}, small[]{0.05};
            require(std::abs(wp.derivative(0).evaluate(z) - dp) < 1e-12 * std::abs(dp), "d p/dz = p'");
            require(std::abs(wp.evaluate(small) - (400 + g2 * 0.0025 / 20 + g3 * 6.25e-6 / 28)) < 1e-8, "Weierstrass Laurent series");
        }

        // Series: theta functions (Jacobi's identity, theta1'(0) = theta2 theta3 theta4), Kepler's equation
        // by its Kapteyn series E - M = sum 2/n J_n(n e) sin(n M) (Watson, Theory of Bessel Functions, 17.2),
        // a finite sum, a slowly convergent sum that the budget does not resolve (NaN), and gradients.
        const double q = 0.3;
        const double t2 = value("theta2(0, 0.3)"), t3 = value("theta3(0, 0.3)"), t4 = value("theta4(0, 0.3)");
        require(std::abs(std::pow(t3, 4) - std::pow(t2, 4) - std::pow(t4, 4)) < 1e-14, "theta3^4 = theta2^4 + theta4^4");
        const auto theta1 = c::AnalyticExpression::parse("theta1(x, 0.3)", {"x"});
        const double zero[]{0.0};
        require(std::abs(theta1.derivative(0).evaluate(zero) - t2 * t3 * t4) < 1e-14, "theta1'(0) = theta2 theta3 theta4");
        (void)q;
        double E = 1, e = 0.3, M = 1;
        for (int k = 0; k < 50; ++k) E -= (E - e * std::sin(E) - M) / (1 - e * std::cos(E));
        near("x + sum(n, 1, inf, 2/n*besselj(n, n*y)*sin(n*x))", E, 1e-14, M, e);
        near("sum(k, 0, 3, x^k)", 1 + 0.5 + 0.25 + 0.125, 1e-15, 0.5);
        near("sum(k, 0, inf, exp(k*log(x))/gamma(k + 1))", std::exp(0.5), 1e-14, 0.5);
        require(std::isnan(value("sum(k, 1, inf, 1/k^2)")), "a series the budget does not resolve is NaN, not truncated");

        // Distributions: H(0) = 1/2, H' = delta (0 off its support, no value on it), delta' = delta^(1).
        near("heaviside(x)", 0.5, 0, 0);
        near("heaviside(x) + sign(x)", 2, 0, 3);
        {
            const auto step = c::AnalyticExpression::parse("heaviside(x - 1)*x", {"x"});
            const double away[]{2.0}, on[]{1.0};
            const auto d = step.derivative(0);  // H(x-1) + x delta(x-1)
            require(d.evaluate(away) == 1 && std::isnan(d.evaluate(on)), "H' = delta: 0 off the support, no value on it");
            const auto dd = c::AnalyticExpression::parse("dirac(x - 1)", {"x"}).derivative(0);
            require(dd.evaluate(away) == 0 && dd.toJson().at("operation").string() == "diracderiv", "delta' = delta^(1)");
        }

        // Exact gradients of every function against central differences.
        for (const char* text : {"besselj(1.5, x*y + 1)", "bessely(2, x + 2)", "besseli(0.5, x*y)", "besselk(1, x + 1.5)",
                                 "sphbesselj(2, x + y)", "sphbessely(1, x + 2)", "airyai(x - y)", "airybi(x*y)",
                                 "gamma(x + 2)", "lgamma(x + y + 2)", "digamma(x + 3)", "polygamma(1, x + 2)", "beta(x + 1, y + 2)",
                                 "erf(x*y)", "erfc(x - y)", "expint(x + 1)", "expintn(2, x + 1)", "gammap(2.5, x + 1)",
                                 "gammaq(1.5, x + y + 1)", "lambertw(x*y + 1)", "lambertwm1(-0.2 + 0.01*x)", "hyp0f1(1.5, x*y)",
                                 "hyp1f1(0.5, 1.5, x - y)", "hyp2f1(0.5, 1, 1.5, 0.3*x*y)", "sn(x, 0.3 + 0.2*y)", "cn(x*y, 0.6)",
                                 "dn(x, 0.4*y)", "jacobizeta(x, 0.3 + 0.1*y)", "weierstrassp(4, 1, x + 0.3)",
                                 "weierstrasspprime(1, 2, 0.2*x + 0.3)", "theta3(x, 0.2 + 0.1*y)", "theta1(x*y, 0.4)",
                                 "legendre(4, x*y)", "hermite(3, x - y)", "tan(x*y)", "asin(0.3*x)", "acos(0.2*y)", "atan(x*y)",
                                 "sinh(x - y)", "cosh(x*y)", "asinh(x)", "acosh(x + 2)", "atanh(0.4*y)", "abs(x - 2)",
                                 "y + sum(n, 1, inf, 2/n*besselj(n, n*0.3)*sin(n*y))"})
            gradient(text, 0.7, 0.9);
        std::cout << "expression function tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "expression function tests failed: " << error.what() << '\n';
        return 1;
    }
}
