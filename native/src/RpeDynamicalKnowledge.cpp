#include "DynamicalKnowledgeInternal.hpp"
#include <algorithm>
#include <cmath>
namespace dynsyslab::core::detail {
DynamicalKnowledge rpeDynamicalKnowledge(const SystemConfiguration &c) {
    DynamicalKnowledge result;
    result.researchStatus = "reviewed_primary_sources";
    const auto &g = c.dimensionlessGroups;
    const double G = g.at("gas_pressure_ratio"), V = g.at("vapor_pressure_ratio"),
                 P = g.at("ambient_pressure_ratio"), S = g.at("capillarity"),
                 nu = g.at("viscosity"), k = c.parameters.at("kappa");
    const bool autonomous = c.driver.amplitude == 0;
    const std::string modern = "https://pp.bme.hu/me/article/download/1243/647#page=2";
    result.facts.push_back(
        {"rpe-equilibria",
         "Radial equilibria and Blake turning point",
         "equilibria",
         "Positive roots of the static pressure balance are candidates. Below vapor pressure two "
         "roots may occur; the larger root can be a saddle. The current model chooses gas content "
         "so r=1 balances pressure. Recompute the residual and eigenvalues for every parameter "
         "set.",
         R"(F(r)=Gr^{-3\kappa}+V-P-S/r=0,\quad r_B=(3\kappa G/S)^{1/(3\kappa-1)})",
         "Autonomous zero-input polytropic RPE; r>0. Blake radius formula requires S>0, G>0, "
         "kappa>1/3.",
         "native (r,dr/dtau); pressure groups from resolved configuration",
         {modern},
         autonomous,
         autonomous ? "Static pressure balance applies" : "Driver is time dependent"});
    result.facts.push_back(
        {"rpe-energy",
         "Conservative energy and saddle separatrix",
         "invariant",
         "At zero viscosity the saddle-energy level guides a candidate homoclinic separatrix when "
         "a saddle exists. Viscosity makes energy decrease; a conservative contour is then not the "
         "manifold.",
         R"(E=\tfrac12r^3v^2-\tfrac{G r^{3-3\kappa}}{3-3\kappa}-\tfrac{V-P}{3}r^3+\tfrac12Sr^2,\quad E'=-\nu r v^2)",
         "r>0; autonomous zero-input RPE. At kappa=1 replace the gas term by -G log(r). "
         "Conservative separatrix requires viscosity=0 and a verified saddle.",
         "native radius and radial velocity",
         {"https://arxiv.org/abs/1409.6699", modern},
         autonomous && nu == 0,
         autonomous && nu == 0
             ? "Conservative energy can validate manifold candidates"
             : "Energy conservation is not applicable with viscosity or driving"});
    result.facts.push_back(
        {"rpe-prosperetti-1975",
         "Subharmonic amplitude-plane separatrices (1975)",
         "reduced_phase_portrait",
         "Prosperetti's equations 9 and 12-13 and Figures 2 and 4 describe critical points and "
         "separatrices in the slowly varying subharmonic amplitude plane. These are approximate "
         "averaged coordinates, not the physical radius-velocity plane; do not copy their curves "
         "into native phase plots.",
         R"(u=C\cos\phi,\quad v=C\sin\phi)",
         "Weakly nonlinear, near-resonant acoustically forced bubble; specific averaging and "
         "damping assumptions. Requires a separately validated reduction and coordinate "
         "reconstruction.",
         "slow subharmonic amplitude/phase (u,v)",
         {"https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/"
          "JASA_NonlinOscilTrans.pdf#page=3",
          "https://thesis.caltech.edu/4443/"},
         false,
         "Advisory literature only; no native-coordinate guide is emitted"});
    result.facts.push_back(
        {"rpe-periodic-manifolds",
         "Invariant manifolds of forced bubble orbits",
         "periodic_orbit_manifolds",
         "Hegedus and Kullmann compute stable and unstable manifolds of saddle periodic orbits on "
         "a stroboscopic section in a low-ambient-pressure regime. Their full forced-system "
         "coordinates and period must be matched before reuse.",
         R"(P^{N}(x)=x,\quad P=\Phi^{2\pi/\Omega})",
         "Harmonically forced RPE; example uses ambient pressure below vapor pressure, isothermal "
         "gas, and finite viscosity. Published parameter-specific curves are not universal.",
         "fixed forcing-phase Poincare section",
         {"https://pp.bme.hu/me/article/download/1243/647#page=5"},
         false,
         "Current sidecar is autonomous; forced stroboscopic manifold continuation is not "
         "implemented"});
    result.equilibriumGuides = [=](const auto &lower, const auto &upper) {
        std::vector<EquilibriumGuide> guides;
        if (!autonomous || lower.size() != 2 || upper.size() != 2 || lower[1] > 0 || upper[1] < 0)
            return guides;
        auto add = [&](double r) {
            if (r >= lower[0] && r <= upper[0] && r > .01)
                guides.push_back({"rpe-equilibria", {r, 0}});
        };
        add(1);
        auto f = [=](double r) { return G * std::pow(r, -3 * k) + V - P - S / r; };
        std::vector<double> breaks{std::max(.010000001, lower[0]), upper[0]};
        if (S > 0 && G > 0 && k > 1. / 3.) {
            double rb = std::pow(3 * k * G / S, 1 / (3 * k - 1));
            if (rb > breaks[0] && rb < breaks[1])
                breaks.push_back(rb);
        }
        std::sort(breaks.begin(), breaks.end());
        for (std::size_t i = 1; i < breaks.size(); ++i) {
            double a = breaks[i - 1], b = breaks[i], fa = f(a), fb = f(b);
            if (!std::isfinite(fa) || !std::isfinite(fb) || fa * fb > 0)
                continue;
            for (int j = 0; j < 70; ++j) {
                double m = (a + b) / 2, fm = f(m);
                if (fa * fm <= 0)
                    b = m;
                else {
                    a = m;
                    fa = fm;
                }
            }
            add((a + b) / 2);
        }
        return guides;
    };
    if (autonomous && nu == 0)
        result.invariantGuides.push_back(
            {"rpe-energy", "RPE saddle energy", [=](const auto &x) {
                 double r = x[0], v = x[1];
                 return .5 * r * r * r * v * v -
                        (k == 1 ? G * std::log(r) : G * std::pow(r, 3 - 3 * k) / (3 - 3 * k)) -
                        (V - P) * r * r * r / 3 + .5 * S * r * r;
             }});
    return result;
}
} // namespace dynsyslab::core::detail
