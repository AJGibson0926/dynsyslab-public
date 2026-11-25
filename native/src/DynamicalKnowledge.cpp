#include "DynamicalKnowledgeInternal.hpp"
#include <algorithm>
#include <cmath>
namespace dynsyslab::core::detail {
namespace {
constexpr double pi = 3.14159265358979323846;
double p(const SystemConfiguration &c, const char *key) { return c.parameters.at(key); }
DynamicalKnowledge builtInFacts(const SystemDescriptor &d, const SystemConfiguration &c);
} // namespace
DynamicalKnowledge builtInDynamicalKnowledge(const SystemDescriptor &d,
                                             const SystemConfiguration &c) {
    auto result = builtInFacts(d, c);
    auto library = builtInEigenfunctions(d, c);
    result.eigenfunctions.insert(result.eigenfunctions.end(), library.begin(), library.end());
    return result;
}
namespace {
DynamicalKnowledge builtInFacts(const SystemDescriptor &d, const SystemConfiguration &c) {
    if (d.identifier == "quadratic_triangular") {
        auto result = triangularKnowledge(c);
        result.papers = d.research.papers;
        return result;
    }
    if (d.identifier == "hodgkin_huxley" || d.identifier == "kuramoto" ||
        d.identifier == "spring_coupled_pendulums" || d.identifier == "marmottant" ||
        d.identifier == "reddy_szeri") {
        auto result = d.identifier == "hodgkin_huxley" ? hodgkinHuxleyKnowledge(c) :
            d.identifier == "kuramoto" ? kuramotoKnowledge(c) :
            d.identifier == "spring_coupled_pendulums" ? springPendulumKnowledge(c) :
            d.identifier == "marmottant" ? marmottantKnowledge(c) : reddySzeriKnowledge(c);
        result.papers = d.research.papers;
        return result;
    }
    if (d.identifier == "double_pendulum") {
        auto result = doublePendulumKnowledge(c);
        result.papers = d.research.papers;
        return result;
    }
    if (d.identifier == "rpe") {
        auto result = rpeDynamicalKnowledge(c);
        result.papers = d.research.papers;
        return result;
    }
    DynamicalKnowledge result;
    result.papers = d.research.papers;
    const bool free = c.driver.amplitude == 0;
    const std::string id = d.identifier + "-equilibria";
    std::vector<std::vector<double>> candidates;
    auto fact = [&](const char *title, const char *statement, const char *formula,
                    const char *conditions, const char *reference) {
        result.facts.push_back({id,
                                title,
                                "equilibria",
                                statement,
                                formula,
                                conditions,
                                d.preferredCoordinates,
                                {reference},
                                free,
                                free ? "Zero-driver equations selected"
                                     : "Driver excludes autonomous equilibrium claims"});
        result.facts.back().evidence =
            "Algebraic candidate derived from native equations; cited literature provides model or "
            "method context. Validate against the configured field.";
        result.researchStatus = "reviewed_exact_reduction";
    };
    const char *teschl = "https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf";
    if (d.identifier == "duffing") {
        fact(
            "Duffing critical points and energy",
            "Stationary positions follow the exact force polynomial in DynSysLab's sign "
            "convention. The conservative saddle energy guides separatrix checks; damping changes "
            "the manifolds.",
            R"(v=0,\ x(\alpha-\beta x^2)=0,\quad H=\tfrac12v^2-\tfrac12\alpha x^2+\tfrac14\beta x^4)",
            "No forcing; energy conservation additionally requires delta=0",
            "https://dlmf.nist.gov/22.19");
        const double a = p(c, "alpha"), b = p(c, "beta");
        candidates = {{0, 0}};
        if (b != 0 && a / b > 0) {
            double r = std::sqrt(a / b);
            candidates.push_back({r, 0});
            candidates.push_back({-r, 0});
        }
        if (free && p(c, "delta") == 0)
            result.invariantGuides.push_back({id, "Duffing saddle energy", [=](const auto &x) {
                                                  return .5 * x[1] * x[1] - .5 * a * x[0] * x[0] +
                                                         .25 * b * std::pow(x[0], 4);
                                              }});
    } else if (d.identifier == "nonlinear_pendulum") {
        fact("Pendulum equilibria and separatrices",
             "Stationary angles repeat by integer multiples of pi. Conservative separatrices lie "
             "at the energy of the upright saddle; angles retain their native winding.",
             R"(\theta=n\pi,\ \omega=0,\quad H=\tfrac12\omega^2+(g/l)(1-\cos\theta))",
             "No forcing; energy conservation additionally requires zero damping", teschl);
        const double g = p(c, "gravity_m_s2") / p(c, "length_m");
        if (free && p(c, "damping_n_m_s") == 0)
            result.invariantGuides.push_back({id, "Pendulum saddle energy", [=](const auto &x) {
                                                  return .5 * x[1] * x[1] +
                                                         g * (1 - std::cos(x[0]));
                                              }});
        result.equilibriumGuides = [id, free](const auto &lo, const auto &hi) {
            std::vector<EquilibriumGuide> seeds;
            if (free && lo.size() == 2 && hi.size() == 2 && lo[1] <= 0 && hi[1] >= 0)
                for (double k = std::ceil(lo[0] / pi);
                     k <= std::floor(hi[0] / pi) && seeds.size() < 4096; ++k)
                    seeds.push_back({id, {k * pi, 0}});
            return seeds;
        };
    } else if (d.identifier == "harmonic_oscillator") {
        fact("Linear oscillator equilibrium",
             "The origin is stationary. If stiffness vanishes there is an equilibrium line, so an "
             "isolated-point interpretation is invalid.",
             R"(v=0,\ kx=0)", "No forcing; positive mass", teschl);
        candidates = {{0, 0}};
    } else if (d.identifier == "van_der_pol") {
        fact("Van der Pol stationary point and cycle guide",
             "The origin is the unique equilibrium. The standard positive-mu oscillator has a "
             "stable self-sustained cycle; its numerical period and Floquet multiplier still "
             "require computation.",
             R"(x=v=0)", "Unforced standard Van der Pol model; cycle statement requires mu>0",
             teschl);
        candidates = {{0, 0}};
    } else if (d.identifier == "lorenz") {
        fact("Lorenz equilibrium branches",
             "Use the origin and real nonzero symmetry pair as seeds. Their stability changes with "
             "parameters and must be computed from the actual Jacobian.",
             R"((0,0,0),\ (\pm\sqrt{\beta(\rho-1)},\pm\sqrt{\beta(\rho-1)},\rho-1))",
             "No forcing; nonzero branch requires beta(rho-1)>0 and nondegenerate parameters",
             "https://journals.ametsoc.org/view/journals/atsc/20/2/"
             "1520-0469_1963_020_0130_dnf_2_0_co_2.xml");
        candidates = {{0, 0, 0}};
        double a = p(c, "beta") * (p(c, "rho") - 1);
        if (a > 0) {
            double x = std::sqrt(a);
            candidates.push_back({x, x, p(c, "rho") - 1});
            candidates.push_back({-x, -x, p(c, "rho") - 1});
        }
    } else if (d.identifier == "rossler") {
        fact("Rossler stationary branches",
             "The equilibrium equations reduce to one quadratic. Degenerate parameter limits may "
             "contain additional sets; the formula supplies candidates only.",
             R"(x=az,\ y=-z,\ az^2-cz+b=0)", "No forcing",
             "https://www.cfm.brown.edu/people/dobrush/am34/tutorial/"
             "AN-EQUATION-FOR-CONTINUOUS-CHAOS_ROSSLER.pdf");
        double a = p(c, "a"), b = p(c, "b"), cc = p(c, "c"), disc = cc * cc - 4 * a * b;
        if (a != 0 && disc >= 0) {
            for (double sign : {-1., 1.}) {
                double z = (cc + sign * std::sqrt(disc)) / (2 * a);
                candidates.push_back({a * z, -z, z});
            }
        } else if (a == 0 && cc != 0)
            candidates.push_back({0, -b / cc, b / cc});
    } else if (d.identifier == "chua") {
        fact("Chua piecewise-linear equilibria",
             "Solve the exact piecewise force balance; branch membership must be checked. The kink "
             "in the characteristic makes smooth-manifold assumptions local to a branch.",
             R"(y=0,\ z=-x,\ x+h(x)=0)", "No forcing; nondegenerate alpha and beta",
             "https://www2.eecs.berkeley.edu/Pubs/TechRpts/1984/307.html");
        candidates = {{0, 0, 0}};
        if (p(c, "m1") != -1) {
            double x = (p(c, "m1") - p(c, "m0")) / (1 + p(c, "m1"));
            if (x >= 1) {
                candidates.push_back({x, 0, -x});
                candidates.push_back({-x, 0, x});
            }
        }
    } else if (d.identifier == "henon_heiles") {
        fact("Exact symmetry-axis family", "Each symmetry axis reduces to a cubic-potential oscillator. Its bounded orbits and homoclinic boundary have exact elliptic or elementary solutions.", R"(q=\lambda s,\ q''+q-q^2=0,\ H_q=\dot q^2/2+q^2/2-q^3/3)", "No forcing; position and momentum on the same symmetry axis; bounded well", "https://epub.uni-regensburg.de/1660/1/brack11.pdf");
        result.facts.back().id = "henon-heiles-straight-line";
        result.facts.back().kind = "analytic_family";
        result.facts.back().applicable = free && p(c, "lam") != 0;
        result.facts.back().evidence = "Primary equations 27-33 inspected; exact symmetry identity and native RHS verified";
        fact("Henon-Heiles critical points",
             "The origin and three escape saddles follow from the gradient of the potential. These "
             "positions are exact candidate identities for the native lambda convention.",
             R"(p_x=p_y=0,\ x(1+2\lambda y)=0,\ y+\lambda(x^2-y^2)=0)", "No forcing",
             "https://doi.org/10.1086/109234");
        candidates = {{0, 0, 0, 0}};
        double l = p(c, "lam");
        if (l != 0) {
            candidates.push_back({0, 1 / l, 0, 0});
            for (double sign : {-1., 1.})
                candidates.push_back({sign * std::sqrt(3.) / (2 * l), -1 / (2 * l), 0, 0});
        }
    } else if (d.identifier == "fput") {
        fact("FPUT zero-displacement equilibrium",
             "The zero position/momentum configuration is stationary for the fixed-end chain. "
             "Nonconvex spring potentials can admit further equilibria; do not stop the numerical "
             "search at this seed.",
             R"(q_j=p_j=0)", "Fixed ends and no forcing", "https://arxiv.org/abs/nlin/0504036");
        candidates = {std::vector<double>(d.stateDimension, 0)};
    } else if (d.identifier == "hindmarsh_rose") {
        fact("Hindmarsh-Rose equilibrium constraint",
             "Stationarity reduces to a scalar cubic followed by reconstruction of the recovery "
             "variables. The native search checks the complete equations.",
             R"(-ax^3+(b-d)x^2-sx+c+sx_r+I=0,\ y=c-dx^2,\ z=s(x-x_r))",
             "No forcing; r nonzero for the z reconstruction",
             "https://pubmed.ncbi.nlm.nih.gov/6144106/");
    } else {
        result.researchStatus = "catalogue_review_pending";
        result.facts.push_back({d.identifier + "-research-status",
                                "Dynamical knowledge coverage",
                                "coverage",
                                "No verified parameter-specific equilibrium or separatrix guide "
                                "has yet been registered for this system. Numerical analysis must "
                                "retain its search-domain and convergence limitations.",
                                "",
                                "No completeness claim",
                                d.preferredCoordinates,
                                {},
                                false,
                                "Further primary-source review required",
                                "coverage record"});
    }
    if (!result.equilibriumGuides)
        result.equilibriumGuides = [=](const auto &lo, const auto &hi) {
            std::vector<EquilibriumGuide> guides;
            if (!free || lo.size() != d.stateDimension || hi.size() != d.stateDimension)
                return guides;
            for (const auto &x : candidates) {
                bool inside = true;
                for (std::size_t j = 0; j < x.size(); ++j)
                    inside = inside && x[j] >= lo[j] && x[j] <= hi[j];
                if (inside)
                    guides.push_back({id, x});
            }
            return guides;
        };
    return result;
}
} // namespace
} // namespace dynsyslab::core::detail
