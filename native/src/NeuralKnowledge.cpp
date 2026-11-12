#include "DynamicalKnowledgeInternal.hpp"
#include "HodgkinHuxleyRates.hpp"
#include "KnowledgeGuides.hpp"
namespace dynsyslab::core::detail {
DynamicalKnowledge hodgkinHuxleyKnowledge(const SystemConfiguration &c) {
    DynamicalKnowledge k;
    k.researchStatus = "reviewed_native_stationarity_and_gating_limit";
    const auto p = c.parameters;
    const bool free = c.driver.amplitude == 0, valid = p.at("capacitance") > 0;
    bool constantPossible = free && valid;
    if (p.at("g_na") != 0 && p.at("g_k") != 0 && p.at("e_na") != p.at("e_k"))
        constantPossible = false;
    if (p.at("g_na") != 0 || p.at("g_k") != 0) {
        const double voltage = p.at("g_na") != 0 ? p.at("e_na") : p.at("e_k");
        constantPossible =
            constantPossible && p.at("bias_current") == p.at("g_l") * (voltage - p.at("e_l"));
    } else if (p.at("g_l") == 0)
        constantPossible = constantPossible && p.at("bias_current") == 0;
    const std::vector<std::string> refs{"https://doi.org/10.1113/jphysiol.1952.sp004764"};
    knowledgeFact(
        k, "hh-stationary-current", "Hodgkin-Huxley stationary current balance", "equilibria",
        "At stationarity each gate equals its voltage-dependent steady value. A scalar current "
        "balance supplies candidate voltages; firing cycles and their stability require numerical "
        "analysis.",
        R"(q_\infty(V)=\alpha_q/(\alpha_q+\beta_q),\quad I_b=g_{Na}m_\infty^3h_\infty(V-E_{Na})+g_Kn_\infty^4(V-E_K)+g_L(V-E_L))",
        "Positive capacitance; zero external driver; roots searched only within the requested "
        "voltage interval",
        "Native voltage in mV; native time in seconds (rates multiplied by 1000)", refs,
        free && valid);
    knowledgeFact(k, "hh-fixed-voltage-gates", "Exact gating on a constant-voltage invariant set",
                  "analytic-limit",
                  "The original voltage-clamp exponentials apply to the free-voltage native model "
                  "only when its current balance keeps voltage constant throughout gate "
                  "relaxation. The registered family enforces that condition.",
                  R"(q(t)=q_\infty+(q_0-q_\infty)e^{-1000(\alpha_q+\beta_q)(t-t_0)})",
                  "No driver; positive capacitance; each voltage-dependent channel has zero "
                  "conductance or zero driving voltage; bias balances leak",
                  "Native four-state model, no imposed voltage clamp", refs, constantPossible);
    k.equilibriumGuides = [p, free, valid](const auto &lo, const auto &hi) {
        std::vector<EquilibriumGuide> result;
        if (!free || !valid || lo.size() != 4 || hi.size() != 4)
            return result;
        auto state = [](double v) {
            const auto r = hodgkinHuxleyRates(v);
            return std::vector<double>{v, r.am / (r.am + r.bm), r.ah / (r.ah + r.bh),
                                       r.an / (r.an + r.bn)};
        };
        auto balance = [&](double v) {
            const auto x = state(v);
            return p.at("bias_current") -
                   p.at("g_na") * std::pow(x[1], 3) * x[2] * (v - p.at("e_na")) -
                   p.at("g_k") * std::pow(x[3], 4) * (v - p.at("e_k")) -
                   p.at("g_l") * (v - p.at("e_l"));
        };
        for (double v : bracketedGuideRoots(balance, lo[0], hi[0])) {
            auto x = state(v);
            if (guideInside(x, lo, hi))
                result.push_back({"hh-stationary-current", std::move(x)});
        }
        return result;
    };
    return k;
}
} // namespace dynsyslab::core::detail
