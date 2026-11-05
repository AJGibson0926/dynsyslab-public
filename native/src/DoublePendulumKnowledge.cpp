#include "DynamicalKnowledgeInternal.hpp"
#include <algorithm>
#include <cmath>
namespace dynsyslab::core::detail {
std::vector<ResearchPaper> doublePendulumResearchPapers() {
    ResearchPaper older;
    older.id = "enolskii-pronine-richter-2003";
    older.title = "Double Pendulum and theta-Divisor";
    older.authors = {"V. Z. Enolskii", "M. Pronine", "P. H. Richter"};
    older.year = 2003;
    older.venue = "Journal of Nonlinear Science 13, 157-174";
    older.doi = "10.1007/s00332-002-0514-0";
    older.urls = {"https://link.springer.com/article/10.1007/s00332-002-0514-0"};
    older.findings = "Gravity-free planar double-pendulum motion reduces to a genus-two inversion "
                     "on a theta divisor.";
    older.applicability =
        "Conservative zero-gravity motion. The native real-quadrature evaluator uses the same "
        "energy/momentum reduction; it does not implement complex Kleinian sigma functions.";
    older.verification = "Publisher abstract inspected; native absolute-angle convention "
                         "independently derived and tested.";
    older.topics = {"analytic solutions", "hyperelliptic functions", "zero gravity"};
    older.relatedFactIds = {"double-pendulum-integrable-reduction"};
    ResearchPaper newer;
    newer.id = "bocci-mingari-scarpello-2022";
    newer.title = "The Differential Equations of Gravity-free Double Pendulum: Lauricella "
                  "Hypergeometric Solutions and Their Inversion";
    newer.authors = {"Alessio Bocci", "Giovanni Mingari Scarpello"};
    newer.year = 2022;
    newer.venue = "arXiv:2104.13211v2";
    newer.urls = {"https://arxiv.org/abs/2104.13211"};
    newer.findings = "Energy and angular momentum yield real hyperelliptic time and angle "
                     "integrals with a Lauricella representation; the relative angle librates or "
                     "rotates according to the invariant level.";
    newer.applicability = "Positive point masses and lengths, zero gravity, no dissipation or "
                          "external input. Convert the paper's relative second angle/velocity into "
                          "DynSysLab's absolute second angle/velocity.";
    newer.coordinates = "theta_native_1=theta_paper; theta_native_2=theta_paper+psi_paper; same "
                        "derivative relation";
    newer.verification =
        "Primary full PDF inspected, sections 3-7 and appendix 9. Native equations, independent "
        "invariant identities and numerical comparisons checked.";
    newer.topics = {"analytic solutions", "Lauricella functions", "turning points",
                    "relative separatrix"};
    newer.relatedFactIds = {"double-pendulum-integrable-reduction", "double-pendulum-momentum"};
    return {older, newer};
}
DynamicalKnowledge doublePendulumKnowledge(const SystemConfiguration &c) {
    const auto &p = c.parameters;
    const double m1 = p.at("mass1_kg"), m2 = p.at("mass2_kg"), l1 = p.at("length1_m"),
                 l2 = p.at("length2_m"), g = p.at("gravity_m_s2");
    const double a = (m1 + m2) * l1 * l1, b = m2 * l1 * l2, cc = m2 * l2 * l2;
    const bool positive = m1 > 0 && m2 > 0 && l1 > 0 && l2 > 0;
    const bool free = c.driver.amplitude == 0;
    const bool conservative = free && p.at("damping1_n_m_s") == 0 && p.at("damping2_n_m_s") == 0;
    DynamicalKnowledge k;
    k.researchStatus = "reviewed_native_invariants_and_integrable_limit";
    auto fact = [&](const char *id, const char *title, const char *kind, const char *statement,
                    const char *formula, const char *conditions, bool applies) {
        k.facts.push_back({id,
                           title,
                           kind,
                           statement,
                           formula,
                           conditions,
                           "dimensional",
                           {"https://arxiv.org/abs/2104.13211"},
                           applies,
                           applies ? "Native parameters satisfy the stated conditions"
                                   : "Configured parameters exclude this relation",
                           "Derived in native absolute angles; the sidecar checks actual field "
                           "residuals and invariant defects"});
    };
    fact("double-pendulum-equilibria", "Double-pendulum stationary configurations", "equilibria",
         "With nonzero gravity both rods are vertical. At zero gravity every zero-velocity "
         "configuration is stationary; these form a continuum, not isolated fixed points.",
         R"(\omega_1=\omega_2=0;\quad g\ne0:\ \theta_1=n\pi,\ \theta_2=m\pi)",
         "Positive masses/lengths and zero external input", positive && free);
    fact(
        "double-pendulum-energy", "Double-pendulum mechanical energy", "first-integral",
        "Mechanical energy includes both point-mass kinetic energies and the gravitational "
        "potential.",
        R"(H=\tfrac12(A\omega_1^2+2B\cos\delta\,\omega_1\omega_2+C\omega_2^2)-g[(m_1+m_2)l_1\cos\theta_1+m_2l_2\cos\theta_2])",
        "Zero damping and external input", positive && conservative);
    fact("double-pendulum-momentum", "Double-pendulum angular momentum", "first-integral",
         "Rotational symmetry supplies a second invariant only when gravity is absent.",
         R"(L=(A+B\cos\delta)\omega_1+(C+B\cos\delta)\omega_2)",
         "Zero gravity, damping and external input", positive && conservative && g == 0);
    fact(
        "double-pendulum-integrable-reduction", "Relative motion and separatrix",
        "separatrix-guide",
        "For positive kinetic energy define epsilon=(L^2-2H(A+C))/(4HB). Relative libration has "
        "-1<epsilon<1; rotation has epsilon<-1. The epsilon=-1 level approaches folded rigid "
        "rotation. That relative separatrix is not a manifold of an isolated equilibrium in the "
        "unreduced four-state space.",
        R"(\delta=\theta_1-\theta_2,\quad \dot\delta^2=\frac{2H(A+C+2B\cos\delta)-L^2}{AC-B^2\cos^2\delta})",
        "Positive energy, zero gravity/damping/input; check the actual invariant level",
        positive && conservative && g == 0);
    if (positive && conservative) {
        k.invariantGuides.push_back(
            {"double-pendulum-energy", "Mechanical energy", [=](const auto &x) {
                 return .5 * (a * x[1] * x[1] + 2 * b * std::cos(x[0] - x[2]) * x[1] * x[3] +
                              cc * x[3] * x[3]) -
                        g * ((m1 + m2) * l1 * std::cos(x[0]) + m2 * l2 * std::cos(x[2]));
             }});
        if (g == 0)
            k.invariantGuides.push_back(
                {"double-pendulum-momentum", "Angular momentum", [=](const auto &x) {
                     const double co = std::cos(x[0] - x[2]);
                     return (a + b * co) * x[1] + (cc + b * co) * x[3];
                 }});
    }
    k.equilibriumGuides = [=](const auto &lo, const auto &hi) {
        std::vector<EquilibriumGuide> result;
        if (!positive || !free || lo.size() != 4 || hi.size() != 4 || lo[1] > 0 || hi[1] < 0 ||
            lo[3] > 0 || hi[3] < 0)
            return result;
        if (g == 0) {
            result.push_back(
                {"double-pendulum-equilibria", {.5 * (lo[0] + hi[0]), 0, .5 * (lo[2] + hi[2]), 0}});
            return result;
        }
        constexpr double pi = 3.14159265358979323846;
        if (std::ceil(lo[2] / pi) > std::floor(hi[2] / pi))
            return result;
        for (double i = std::ceil(lo[0] / pi); i <= std::floor(hi[0] / pi) && result.size() < 4096;
             ++i)
            for (double j = std::ceil(lo[2] / pi);
                 j <= std::floor(hi[2] / pi) && result.size() < 4096; ++j)
                result.push_back({"double-pendulum-equilibria", {i * pi, 0, j * pi, 0}});
        return result;
    };
    return k;
}
} // namespace dynsyslab::core::detail
