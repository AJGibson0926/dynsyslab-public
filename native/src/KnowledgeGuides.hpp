#pragma once
#include <algorithm>
#include <cmath>
#include <dynsyslab/core/DynamicalKnowledge.hpp>
#include <limits>
namespace dynsyslab::core::detail {
inline bool guideInside(const std::vector<double> &x, const std::vector<double> &lo,
                        const std::vector<double> &hi) {
    if (x.size() != lo.size() || x.size() != hi.size())
        return false;
    for (std::size_t j = 0; j < x.size(); ++j)
        if (!std::isfinite(x[j]) || !(lo[j] <= x[j] && x[j] <= hi[j]))
            return false;
    return true;
}
// Bounded sign-change search supplies seeds, never an exhaustive root certificate.
// Callers split at branch boundaries; the sidecar checks the complete native RHS.
template <class F> std::vector<double> bracketedGuideRoots(F f, double lo, double hi) {
    std::vector<double> roots;
    if (!std::isfinite(lo) || !std::isfinite(hi) || lo > hi)
        return roots;
    auto add = [&](double x) {
        if (roots.empty() || std::abs(x - roots.back()) > 1e-10 * (1 + std::abs(x)))
            roots.push_back(x);
    };
    double a = lo, fa = f(a);
    if (fa == 0)
        add(a);
    for (int j = 1; j <= 512; ++j) {
        const double weight = double(j) / 512, b = (1 - weight) * lo + weight * hi, fb = f(b);
        if (std::isfinite(fa) && std::isfinite(fb) && fa != 0 && fb != 0 &&
            std::signbit(fa) != std::signbit(fb)) {
            double l = a, r = b, fl = fa;
            for (int k = 0; k < 80; ++k) {
                double m = .5 * l + .5 * r, fm = f(m);
                if (!std::isfinite(fm))
                    break;
                if (fm == 0 || m == l || m == r) {
                    l = r = m;
                    break;
                }
                if (std::signbit(fl) != std::signbit(fm))
                    r = m;
                else {
                    l = m;
                    fl = fm;
                }
            }
            add(.5 * l + .5 * r);
        }
        if (fb == 0)
            add(b);
        a = b;
        fa = fb;
    }
    return roots;
}
inline void knowledgeFact(DynamicalKnowledge &k, const char *id, const char *title,
                          const char *kind, const char *statement, const char *formula,
                          const char *conditions, const char *coordinates,
                          std::vector<std::string> references, bool applies) {
    k.facts.push_back({id, title, kind, statement, formula, conditions, coordinates,
                       std::move(references), applies,
                       applies ? "Native parameters satisfy the stated conditions"
                               : "Configured parameters exclude this relation",
                       "Native equation reduction with primary model/method context; independently "
                       "validate field residuals, branch membership and stability"});
}
} // namespace dynsyslab::core::detail
