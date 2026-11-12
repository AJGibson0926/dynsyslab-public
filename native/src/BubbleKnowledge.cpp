#include "DynamicalKnowledgeInternal.hpp"
#include "KnowledgeGuides.hpp"
namespace dynsyslab::core::detail {
DynamicalKnowledge marmottantKnowledge(const SystemConfiguration &c) {
    DynamicalKnowledge k;
    k.researchStatus = "reviewed_shell_branches_and_pressure_balance";
    const auto p = c.parameters, g = c.dimensionlessGroups;
    const bool free = c.driver.amplitude == 0;
    const std::vector<std::string> refs{"https://doi.org/10.1121/1.2109427"};
    knowledgeFact(
        k, "marmottant-pressure-balance", "Coated-bubble stationary pressure balance", "equilibria",
        "Zero wall velocity removes liquid, shell and acoustic damping terms. The remaining gas "
        "pressure balances ambient and branch-specific Laplace pressure. The supplied reference "
        "radius is not assumed to solve the balance for every shell configuration.",
        R"(v=0,\quad G r^{-3\kappa}=1+2\sigma(r)/(p_0R_0r))",
        "Zero driver; positive radius; use the actual shell branch",
        "Native dimensionless radius/velocity", refs, free);
    knowledgeFact(
        k, "marmottant-shell-history", "Shell history and smooth-phase-space limits", "domain",
        "Between rupture and breaking radii, surface tension depends on hidden rupture history. A "
        "two-coordinate autonomous phase portrait cannot identify that state. Numerical manifolds "
        "require a domain wholly within one history-independent smooth branch.",
        R"(r_b=R_b/R_0,\quad r_r=r_b\sqrt{1+\sigma_w/\chi},\quad r_{break}=r_b\sqrt{1+\sigma_{crit}/\chi})",
        "Avoid shell boundaries and the history-dependent interval",
        "Native model including reset shell history", refs, true);
    knowledgeFact(k, "marmottant-analytic-scope",
                  "Exact stationary and approximate oscillatory references", "analytic-limit",
                  "The stationary family is exact when the native pressure residual vanishes. "
                  "Published damped harmonic oscillations are a local linearization; they are not "
                  "registered as exact full nonlinear trajectories.",
                  R"(r(t)=r_*,\quad v(t)=0)", "Zero driver and verified native pressure balance",
                  "Native coordinates", refs, free);
    k.equilibriumGuides = [=](const auto &lo, const auto &hi) {
        std::vector<EquilibriumGuide> result;
        if (!free || lo.size() != 2 || hi.size() != 2 || lo[0] <= 0 || lo[1] > 0 || hi[1] < 0)
            return result;
        const double b = g.at("buckling_radius_ratio"), r = g.at("rupture_radius_ratio"),
                     br = g.at("breaking_radius_ratio");
        int branch = -1;
        if (hi[0] < std::min(b, r))
            branch = 0;
        else if (lo[0] > b && hi[0] < r)
            branch = 1;
        else if (lo[0] > std::max(br, r))
            branch = 2;
        if (branch < 0)
            return result;
        auto balance = [=](double x) {
            double sigma = branch == 0   ? 0
                           : branch == 1 ? p.at("shell_elasticity_n_m") * ((x / b) * (x / b) - 1)
                                         : p.at("sigma_n_m");
            return g.at("gas_pressure_ratio") * std::pow(x, -3 * p.at("kappa")) - 1 -
                   2 * sigma / (p.at("ambient_pressure_pa") * p.at("equilibrium_radius_m") * x);
        };
        for (double r0 : bracketedGuideRoots(balance, lo[0], hi[0]))
            result.push_back({"marmottant-pressure-balance", {r0, 0}});
        return result;
    };
    return k;
}
DynamicalKnowledge reddySzeriKnowledge(const SystemConfiguration &c) {
    DynamicalKnowledge k;
    k.researchStatus = "reviewed_native_radial_shape_reduction";
    const auto p = c.parameters, g = c.dimensionlessGroups;
    const bool quiet =
        c.driver.amplitude == 0 &&
        (p.at("translation_enabled") == 0 || g.at("angular_frequency") == 0 ||
         (g.at("pressure_amplitude") == 0 && g.at("pressure_derivative_amplitude") == 0));
    const std::vector<std::string> refs{"https://doi.org/10.1121/1.1502899",
                                        "https://doi.org/10.1063/1.1483840"};
    knowledgeFact(
        k, "reddy-szeri-equilibria", "Stationary-radius and zero-shape equilibria", "equilibria",
        "The native radial equation is independent of shape coordinates. Stationarity reduces to a "
        "pressure balance, with the configured exponent below unit radius and exponent one above "
        "it. Axial position is free, so stationary states form a continuum.",
        R"(v=w=a_n=b_n=0,\quad G r^{-3\kappa(r)}=1+C/r)",
        "Zero effective forcing; positive radius and nonsingular native denominator; n=2,3,4,5",
        "Native dimensionless 12-state closure", refs, quiet);
    knowledgeFact(
        k, "reddy-szeri-shape-limit",
        "Exact shape oscillations in the stationary-radius native closure", "analytic-limit",
        "With stationary radius and zero translation, shape modes decouple into linear "
        "oscillators. This is an exact limit of the implemented perturbative closure, not a claim "
        "of exact finite-amplitude bubble hydrodynamics. Its shape equations have no separate "
        "viscous damping term at zero radial velocity.",
        R"(a_n''+\Omega_n^2a_n=0,\quad \Omega_n^2=\tfrac12 Eu\,C(n-1)(n+1)(n+2)/r_*^3)",
        "Verified radial equilibrium, zero radial/translational velocity and effective forcing",
        "Native shape amplitudes and model time", refs, quiet);
    knowledgeFact(
        k, "reddy-szeri-phase-domain", "Radial branch and translation configuration", "domain",
        "At unit radius the native gas exponent changes branch; derivatives need one-sided "
        "interpretation unless the exponents agree. Disabling translation freezes both axial "
        "coordinates. Full coupled manifolds must use the complete 12-state equations.",
        R"(\kappa(r)=\kappa_0\ (r\le1),\quad \kappa(r)=1\ (r>1))",
        "No global smoothness claim across branch boundaries", "Native closure", refs, true);
    k.equilibriumGuides = [=](const auto &lo, const auto &hi) {
        std::vector<EquilibriumGuide> result;
        if (!quiet || lo.size() != 12 || hi.size() != 12 || lo[0] <= 0)
            return result;
        auto balance = [=](double r) {
            const double exponent = r > 1 ? 1 : p.at("kappa");
            return g.at("gas_pressure_ratio") * std::pow(r, -3 * exponent) - 1 -
                   g.at("capillarity") / r;
        };
        // The pressure balance is continuous at one; split the search to retain branch identity.
        std::vector<double> roots;
        if (lo[0] <= 1)
            roots = bracketedGuideRoots(balance, lo[0], std::min(hi[0], 1.));
        if (hi[0] > 1) {
            auto above = bracketedGuideRoots(balance, std::max(lo[0], 1.), hi[0]);
            roots.insert(roots.end(), above.begin(), above.end());
        }
        for (double r : roots) {
            std::vector<double> x(12, 0);
            x[0] = r;
            x[2] = std::clamp(0., lo[2], hi[2]);
            if (guideInside(x, lo, hi))
                result.push_back({"reddy-szeri-equilibria", std::move(x)});
        }
        return result;
    };
    return k;
}
} // namespace dynsyslab::core::detail
