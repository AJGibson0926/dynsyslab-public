// The expression language's special functions (ExpressionFunctions.hpp). Values come from Boost.Math
// (Bessel, Airy, gamma family, error functions, exponential integrals, Lambert W, hypergeometric,
// zeta, incomplete gamma) and from the elliptic integrals here; a domain error is NaN, not an
// exception. Each derivative is an exact expression:
//   Bessel  J'_nu = (J_{nu-1} - J_{nu+1})/2, Y likewise, I'_nu = (I_{nu-1} + I_{nu+1})/2,
//           K'_nu = -(K_{nu-1} + K_{nu+1})/2, j'_n = (n j_{n-1} - (n+1) j_{n+1})/(2n+1) (DLMF 10.6, 10.51);
//   Airy    Ai'' = x Ai, Bi'' = x Bi (DLMF 9.2);  gamma: Gamma' = Gamma psi, psi^(n)' = psi^(n+1) (DLMF 5.15);
//   erf'    = 2 e^{-x^2}/sqrt(pi);  Ei' = e^x/x, E_n' = -E_{n-1} (DLMF 6.2, 8.19);  W' = W/(x(1 + W)) (DLMF 4.13);
//   1F1'    = (a/b) 1F1(a+1; b+1), 2F1' = (ab/c) 2F1(a+1, b+1; c+1), 0F1' = 0F1(; b+1)/b (DLMF 13.3, 15.5, 16.3);
//   P(a,x)' = x^{a-1} e^{-x}/Gamma(a) (DLMF 8.8);  beta: dB/da = B (psi(a) - psi(a+b)) (DLMF 5.12);
//   Jacobi  am(u|m): d am/du = dn = sqrt(1 - m sin^2 am), d am/dm = -dn dF/dm at phi = am (DLMF 22.16, 19.4);
//   Weierstrass p(z; g2, g3) through Jacobi functions (DLMF 23.6; Abramowitz & Stegun 18.9),
//           d p/dz = p', d p'/dz = 6 p^2 - g2/2 (DLMF 23.3.10).
#include "ExpressionFunctions.hpp"

#include <dynsyslab/core/EllipticIntegrals.hpp>
#include <dynsyslab/core/PolynomialRoots.hpp>

#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/zeta.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {
using E = AnalyticExpression;
using P = std::vector<double>;
using A = std::vector<double>;
using X = std::vector<AnalyticExpression>;
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
const double pi = std::acos(-1.0);

template <class F> double guarded(F f) {
    try { const double v = f(); return std::isfinite(v) ? v : nan; } catch (...) { return nan; }
}
E S(const std::string& name, X p, X a) { return E::special(name, std::move(p), std::move(a)); }
E C(double v, const E& like) { return E::constant(v, like.stateDimension()); }
[[noreturn]] void noPartial(const std::string& name, std::size_t i) {
    throw std::invalid_argument(name + " has no exact derivative in argument " + std::to_string(i + 1) + "; keep it constant");
}

// The Jacobi amplitude am(u|m), m < 1: F(am|m) = u, continued through every turn by
// am(u + 2K) = am(u) + pi, solved by Newton's method on F (dF/dphi = 1/Delta > 0).
double amplitude(double u, double m) {
    if (!std::isfinite(u) || !(m < 1)) return nan;
    if (m == 0) return u;
    const double K = ellipticK(m), turns = std::round(u / (2 * K)), r = u - 2 * K * turns;
    double phi = r * pi / (2 * K);
    for (int step = 0; step < 60; ++step) {
        const double delta = std::sqrt(1 - m * std::sin(phi) * std::sin(phi));
        const double change = (ellipticF(phi, m) - r) * delta;
        phi -= change;
        if (std::abs(change) < 1e-16 * std::max(1.0, std::abs(phi))) break;
    }
    return phi + turns * pi;
}

// Weierstrass p and p' for real invariants (DLMF 23.6; A&S 18.9): with three real roots e1 > e2 > e3
// of 4t^3 - g2 t - g3, p = e3 + (e1 - e3)/sn^2(sqrt(e1 - e3) z | (e2 - e3)/(e1 - e3)); with one real root
// e2, H = sqrt(3 e2^2 - g2/4), p = e2 + H (1 + cn(2 sqrt(H) z | m))/(1 - cn), m = 1/2 - 3 e2/(4H).
void weierstrass(double z, double g2, double g3, double* p, double* dp) {
    const auto roots = realPolynomialRoots({4, 0, -g2, -g3});
    if (roots.size() == 3) {
        const double e1 = roots[0], e2 = roots[1], e3 = roots[2], lam = std::sqrt(e1 - e3), m = (e2 - e3) / (e1 - e3);
        const double am = amplitude(lam * z, m), sn = std::sin(am), cn = std::cos(am), dn = std::sqrt(1 - m * sn * sn);
        *p = e3 + (e1 - e3) / (sn * sn);
        *dp = -2 * (e1 - e3) * lam * cn * dn / (sn * sn * sn);
    } else if (roots.size() == 1) {
        const double e2 = roots[0], H = std::sqrt(3 * e2 * e2 - g2 / 4), m = 0.5 - 3 * e2 / (4 * H), s = 2 * std::sqrt(H);
        const double am = amplitude(s * z, m), sn = std::sin(am), cn = std::cos(am), dn = std::sqrt(1 - m * sn * sn);
        *p = e2 + H * (1 + cn) / (1 - cn);
        *dp = -4 * H * std::sqrt(H) * sn * dn / ((1 - cn) * (1 - cn));
    } else {
        *p = *dp = nan;
    }
}

std::map<std::string, SpecialFunction> registry() {
    namespace bm = boost::math;
    std::map<std::string, SpecialFunction> f;
    const auto unary = [&](const std::string& name, double (*value)(double), std::function<E(const E&, const E&)> derivative) {
        f[name] = {0, 1, [value](const P&, const A& a) { return guarded([&] { return value(a[0]); }); },
                   [derivative](const X&, const X& a, const E& self, std::size_t) { return derivative(a[0], self); }};
    };
    // Elementary functions beyond exp, log, sin, cos, tanh, sech.
    unary("tan", [](double x) { return std::tan(x); }, [](const E&, const E& t) { return 1 + t.power(2); });
    unary("asin", [](double x) { return std::asin(x); }, [](const E& x, const E&) { return (1 - x.power(2)).power(-0.5); });
    unary("acos", [](double x) { return std::acos(x); }, [](const E& x, const E&) { return -(1 - x.power(2)).power(-0.5); });
    unary("atan", [](double x) { return std::atan(x); }, [](const E& x, const E&) { return 1 / (1 + x.power(2)); });
    unary("sinh", [](double x) { return std::sinh(x); }, [](const E& x, const E&) { return S("cosh", {}, {x}); });
    unary("cosh", [](double x) { return std::cosh(x); }, [](const E& x, const E&) { return S("sinh", {}, {x}); });
    unary("asinh", [](double x) { return std::asinh(x); }, [](const E& x, const E&) { return (x.power(2) + 1).power(-0.5); });
    unary("acosh", [](double x) { return std::acosh(x); }, [](const E& x, const E&) { return (x.power(2) - 1).power(-0.5); });
    unary("atanh", [](double x) { return std::atanh(x); }, [](const E& x, const E&) { return 1 / (1 - x.power(2)); });
    unary("abs", [](double x) { return std::abs(x); }, [](const E& x, const E& self) { return x / self; });
    // Error functions, gamma family, exponential integral, Lambert W.
    unary("erf", [](double x) { return bm::erf(x); }, [](const E& x, const E&) { return 2 / std::sqrt(pi) * (-x.power(2)).unary("exp"); });
    unary("erfc", [](double x) { return bm::erfc(x); }, [](const E& x, const E&) { return -2 / std::sqrt(pi) * (-x.power(2)).unary("exp"); });
    unary("gamma", [](double x) { return bm::tgamma(x); }, [](const E& x, const E& self) { return self * S("digamma", {}, {x}); });
    unary("lgamma", [](double x) { return bm::lgamma(x); }, [](const E& x, const E&) { return S("digamma", {}, {x}); });
    unary("digamma", [](double x) { return bm::digamma(x); }, [](const E& x, const E&) { return S("polygamma", {C(1, x)}, {x}); });
    unary("expint", [](double x) { return bm::expint(x); }, [](const E& x, const E&) { return x.unary("exp") / x; });
    unary("lambertw", [](double x) { return bm::lambert_w0(x); }, [](const E& x, const E& w) { return w / (x * (1 + w)); });
    unary("lambertwm1", [](double x) { return bm::lambert_wm1(x); }, [](const E& x, const E& w) { return w / (x * (1 + w)); });
    unary("airyai", [](double x) { return bm::airy_ai(x); }, [](const E& x, const E&) { return S("airyaiprime", {}, {x}); });
    unary("airybi", [](double x) { return bm::airy_bi(x); }, [](const E& x, const E&) { return S("airybiprime", {}, {x}); });
    unary("airyaiprime", [](double x) { return bm::airy_ai_prime(x); }, [](const E& x, const E&) { return x * S("airyai", {}, {x}); });
    unary("airybiprime", [](double x) { return bm::airy_bi_prime(x); }, [](const E& x, const E&) { return x * S("airybi", {}, {x}); });
    unary("zeta", [](double s) { return bm::zeta(s); }, [](const E&, const E&) -> E { noPartial("zeta", 0); });
    f["polygamma"] = {1, 1, [](const P& p, const A& a) { return guarded([&] { return bm::polygamma(int(p[0]), a[0]); }); },
                      [](const X& p, const X& a, const E&, std::size_t) { return S("polygamma", {p[0] + 1}, {a[0]}); }};
    f["expintn"] = {1, 1, [](const P& p, const A& a) { return guarded([&] { return bm::expint(unsigned(p[0]), a[0]); }); },
                    [](const X& p, const X& a, const E&, std::size_t) {  // E_n' = -E_{n-1}; E_0' = -e^{-x}(1 + 1/x)/x
                        const double n = p[0].evaluate(std::vector<double>(p[0].stateDimension(), nan).data());
                        return n == 0 ? -(-a[0]).unary("exp") * (1 + 1 / a[0]) / a[0] : -S("expintn", {p[0] - 1}, {a[0]}); }};
    f["beta"] = {0, 2, [](const P&, const A& a) { return guarded([&] { return bm::beta(a[0], a[1]); }); },
                 [](const X&, const X& a, const E& self, std::size_t i) {
                     return self * (S("digamma", {}, {a[i]}) - S("digamma", {}, {a[0] + a[1]})); }};
    f["gammap"] = {1, 1, [](const P& p, const A& a) { return guarded([&] { return bm::gamma_p(p[0], a[0]); }); },
                   [](const X& p, const X& a, const E&, std::size_t) {
                       return ((p[0] - 1) * a[0].unary("log") - a[0]).unary("exp") / S("gamma", {}, {p[0]}); }};
    f["gammaq"] = {1, 1, [](const P& p, const A& a) { return guarded([&] { return bm::gamma_q(p[0], a[0]); }); },
                   [](const X& p, const X& a, const E&, std::size_t) {
                       return -((p[0] - 1) * a[0].unary("log") - a[0]).unary("exp") / S("gamma", {}, {p[0]}); }};
    // Distributions (Gel'fand & Shilov, Generalized Functions 1, 1964; DLMF 1.16): the Heaviside step
    // H (1/2 at 0), sign, the Dirac delta and its derivatives delta^(n). Their derivatives are exact as
    // distributions, H' = delta, sign' = 2 delta, (delta^(n))' = delta^(n+1); a term with a delta is 0
    // off its support and NaN on it (it has no value there, only pairings with test functions).
    unary("heaviside", [](double x) { return x > 0 ? 1.0 : x < 0 ? 0.0 : 0.5; }, [](const E& x, const E&) { return S("dirac", {}, {x}); });
    unary("sign", [](double x) { return x > 0 ? 1.0 : x < 0 ? -1.0 : 0.0; }, [](const E& x, const E&) { return 2 * S("dirac", {}, {x}); });
    unary("dirac", [](double x) { return x == 0 ? nan : 0.0; }, [](const E& x, const E&) { return S("diracderiv", {C(1, x)}, {x}); });
    f["diracderiv"] = {1, 1, [](const P&, const A& a) { return a[0] == 0 ? nan : 0.0; },
                       [](const X& p, const X& a, const E&, std::size_t) { return S("diracderiv", {p[0] + 1}, {a[0]}); }};
    // Powers with an exponent free of the state (a summation variable): d x^p/dx = p x^(p-1).
    f["pow"] = {1, 1, [](const P& p, const A& a) { return guarded([&] { return std::pow(a[0], p[0]); }); },
                [](const X& p, const X& a, const E&, std::size_t) { return p[0] * S("pow", {p[0] - 1}, {a[0]}); }};
    // Bessel functions of constant real order nu.
    const auto bessel = [&](const std::string& name, double (*value)(double, double), double sign, double scale) {
        f[name] = {1, 1, [value](const P& p, const A& a) { return guarded([&] { return value(p[0], a[0]); }); },
                   [name, sign, scale](const X& p, const X& a, const E&, std::size_t) {
                       return scale * (S(name, {p[0] - 1}, {a[0]}) + sign * S(name, {p[0] + 1}, {a[0]})); }};
    };
    bessel("besselj", [](double n, double x) { return bm::cyl_bessel_j(n, x); }, -1, 0.5);
    bessel("bessely", [](double n, double x) { return bm::cyl_neumann(n, x); }, -1, 0.5);
    bessel("besseli", [](double n, double x) { return bm::cyl_bessel_i(n, x); }, 1, 0.5);
    bessel("besselk", [](double n, double x) { return bm::cyl_bessel_k(n, x); }, 1, -0.5);
    const auto spherical = [&](const std::string& name, double (*value)(unsigned, double)) {
        f[name] = {1, 1, [value](const P& p, const A& a) { return guarded([&] { return value(unsigned(p[0]), a[0]); }); },
                   [name](const X& p, const X& a, const E&, std::size_t) {
                       // j'_0 = -j_1; else (n j_{n-1} - (n+1) j_{n+1})/(2n+1) (orders n >= 1 in a series)
                       const auto& n = p[0];
                       if (n.evaluate(std::vector<double>(n.stateDimension(), nan).data()) == 0) return -S(name, {C(1, n)}, {a[0]});
                       return (n * S(name, {n - 1}, {a[0]}) - (n + 1) * S(name, {n + 1}, {a[0]})) / (2 * n + 1); }};
    };
    spherical("sphbesselj", [](unsigned n, double x) { return bm::sph_bessel(n, x); });
    spherical("sphbessely", [](unsigned n, double x) { return bm::sph_neumann(n, x); });
    // Hypergeometric functions of constant parameters.
    f["hyp0f1"] = {1, 1, [](const P& p, const A& a) { return guarded([&] { return bm::hypergeometric_0F1(p[0], a[0]); }); },
                   [](const X& p, const X& a, const E&, std::size_t) { return S("hyp0f1", {p[0] + 1}, {a[0]}) / p[0]; }};
    f["hyp1f1"] = {2, 1, [](const P& p, const A& a) { return guarded([&] { return bm::hypergeometric_1F1(p[0], p[1], a[0]); }); },
                   [](const X& p, const X& a, const E&, std::size_t) { return p[0] / p[1] * S("hyp1f1", {p[0] + 1, p[1] + 1}, {a[0]}); }};
    f["hyp2f1"] = {3, 1, [](const P& p, const A& a) {
                       return guarded([&] { return bm::hypergeometric_pFq({p[0], p[1]}, {p[2]}, a[0]); }); },
                   [](const X& p, const X& a, const E&, std::size_t) {
                       return p[0] * p[1] / p[2] * S("hyp2f1", {p[0] + 1, p[1] + 1, p[2] + 1}, {a[0]}); }};
    // Jacobi amplitude; sn, cn, dn and the rest are built from it (expandSugar).
    f["jacobiam"] = {0, 2, [](const P&, const A& a) { return amplitude(a[0], a[1]); },
                     [](const X&, const X& a, const E& am, std::size_t i) {
                         const auto& m = a[1];
                         const auto s = am.unary("sin"), delta = (1 - m * s.power(2)).power(0.5);
                         if (i == 0) return delta;
                         // dF/dm at phi = am (DLMF 19.4.1-2), with F(am|m) = u.
                         const auto dF = am.function("ellipeinc", m) / (2 * m * (1 - m)) - a[0] / (2 * m)
                                         - (2 * am).unary("sin") / (4 * (1 - m) * delta);
                         return -delta * dF; }};
    // Weierstrass p and p' for constant real invariants g2, g3.
    f["weierstrassp"] = {2, 1, [](const P& p, const A& a) { double v, d; weierstrass(a[0], p[0], p[1], &v, &d); return v; },
                         [](const X& p, const X& a, const E&, std::size_t) { return S("weierstrasspprime", p, {a[0]}); }};
    f["weierstrasspprime"] = {2, 1, [](const P& p, const A& a) { double v, d; weierstrass(a[0], p[0], p[1], &v, &d); return d; },
                              [](const X& p, const X& a, const E&, std::size_t) {
                                  return 6 * S("weierstrassp", p, {a[0]}).power(2) - p[0] / 2; }};
    return f;
}

const std::map<std::string, SpecialFunction>& table() {
    static const auto t = registry();
    return t;
}

// The orthogonal polynomial of degree n in x by its three-term recurrence (DLMF 18.9.1), so it is an
// ordinary polynomial expression with an exact derivative.
E polynomial(const std::string& name, int n, const E& x) {
    if (n < 0 || n > 200) throw std::invalid_argument(name + " needs a degree from 0 to 200");
    E previous = C(1, x);
    if (n == 0) return previous;
    E current = name == "hermite" ? 2 * x : name == "laguerre" ? 1 - x : name == "chebyshevu" ? 2 * x : x;
    for (int k = 1; k < n; ++k) {
        E next = C(0, x);
        if (name == "legendre") next = ((2 * k + 1) * x * current - k * previous) / (k + 1);
        else if (name == "hermite") next = 2 * x * current - 2 * k * previous;
        else if (name == "laguerre") next = ((2 * k + 1 - x) * current - k * previous) / (k + 1);
        else next = 2 * x * current - previous;  // Chebyshev T and U
        previous = current;
        current = next;
    }
    return current;
}

}  // namespace

const SpecialFunction* specialFunction(const std::string& name) {
    const auto& t = table();
    const auto found = t.find(name);
    return found == t.end() ? nullptr : &found->second;
}

std::vector<std::string> specialFunctionNames() {
    std::vector<std::string> out;
    for (const auto& [name, f] : table()) out.push_back(name);
    return out;
}

bool expandSugar(const std::string& name, const std::vector<AnalyticExpression>& a, std::size_t slot, AnalyticExpression& result) {
    const auto require = [&](std::size_t count) {
        if (a.size() != count) throw std::invalid_argument(name + " takes " + std::to_string(count) + " arguments");
    };
    const auto integerConstant = [&](const E& e) {
        const double v = e.evaluate(std::vector<double>(e.stateDimension(), std::nan("")).data());
        if (!std::isfinite(v) || v != std::floor(v)) throw std::invalid_argument(name + " needs a constant integer degree");
        return int(v);
    };
    if (name == "sn" || name == "cn" || name == "dn") {  // Jacobi elliptic functions of u and m
        require(2);
        const auto am = S("jacobiam", {}, {a[0], a[1]});
        result = name == "sn" ? am.unary("sin") : name == "cn" ? am.unary("cos") : (1 - a[1] * am.unary("sin").power(2)).power(0.5);
        return true;
    }
    if (name == "jacobizeta") {  // Z(phi|m) = E(phi|m) - E(m) F(phi|m)/K(m) (DLMF 22.16.32)
        require(2);
        result = a[0].function("ellipeinc", a[1]) - a[1].unary("ellipe") * a[0].function("ellipf", a[1]) / a[1].unary("ellipk");
        return true;
    }
    if (name == "legendre" || name == "hermite" || name == "laguerre" || name == "chebyshevt" || name == "chebyshevu") {
        require(2);
        result = polynomial(name, integerConstant(a[0]), a[1]);
        return true;
    }
    if (name == "theta1" || name == "theta2" || name == "theta3" || name == "theta4") {
        // Jacobi theta functions of z and nome q, 0 < q < 1, as their Fourier series (DLMF 20.2.1-4),
        // q^s = exp(s log q), (-1)^n = cos(pi n).
        require(2);
        const auto& z = a[0];
        const auto logq = a[1].unary("log");
        const auto n = E::index(z.stateDimension(), slot);
        if (name == "theta1")
            result = 2 * E::sum(slot, 0, INFINITY, (pi * n).unary("cos") * (logq * (n + 0.5).power(2)).unary("exp") * ((2 * n + 1) * z).unary("sin"));
        else if (name == "theta2")
            result = 2 * E::sum(slot, 0, INFINITY, (logq * (n + 0.5).power(2)).unary("exp") * ((2 * n + 1) * z).unary("cos"));
        else if (name == "theta3")
            result = 1 + 2 * E::sum(slot, 1, INFINITY, (logq * n.power(2)).unary("exp") * (2 * n * z).unary("cos"));
        else
            result = 1 + 2 * E::sum(slot, 1, INFINITY, (pi * n).unary("cos") * (logq * n.power(2)).unary("exp") * (2 * n * z).unary("cos"));
        return true;
    }
    return false;
}

}  // namespace dynsyslab::core::detail
