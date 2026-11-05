#include "DynamicalKnowledgeInternal.hpp"
#include "KnowledgeGuides.hpp"
#include <numeric>
namespace dynsyslab::core::detail {
namespace {
constexpr double pi = 3.14159265358979323846;
}
DynamicalKnowledge kuramotoKnowledge(const SystemConfiguration &c) {
    DynamicalKnowledge k;
    k.researchStatus = "reviewed_finite_array_symmetry_and_two_oscillator_reduction";
    const auto p = c.parameters;
    const std::size_t n = static_cast<std::size_t>(p.at("oscillators"));
    const double coupling = p.at("coupling"), spread = p.at("frequency_spread");
    const bool free = c.driver.amplitude == 0;
    const std::vector<std::string> refs{"https://doi.org/10.4171/EMSS/17",
                                        "https://arxiv.org/abs/1403.2083"};
    knowledgeFact(k, "kuramoto-mean-phase", "Finite-array phase symmetry", "first-integral",
                  "Pairwise sine coupling cancels in the phase sum. For two or more native evenly "
                  "spaced frequencies the sum is conserved without input. For one oscillator its "
                  "frequency is minus the supplied spread.",
                  R"(\frac{d}{dt}\sum_j\theta_j=\sum_j\omega_j+Nu;\quad N>1:\ \sum_j\omega_j=0)",
                  "Zero external driver; phases kept unwrapped", "Native unwrapped phases", refs,
                  free);
    if (free && (n > 1 || spread == 0))
        k.invariantGuides.push_back(
            {"kuramoto-mean-phase", "Sum of unwrapped phases",
             [](const auto &x) { return std::accumulate(x.begin(), x.end(), 0.); }});
    knowledgeFact(
        k, "kuramoto-relative-locking", "Two-oscillator locking and phase-slip threshold",
        "relative-equilibria",
        "Two oscillators reduce to the Adler equation. Locked phase differences exist when the "
        "detuning magnitude does not exceed the coupling magnitude. Global phase-shift symmetry "
        "leaves a neutral direction; these are not isolated hyperbolic fixed points.",
        R"(\delta'=2s-K\sin\delta,\quad \sin\delta_*=2s/K,\quad \lambda_{relative}=-K\cos\delta_*)",
        "N=2 and zero external driver; K=0 requires separate detuning-zero interpretation",
        "delta=theta2-theta1", refs, free && n == 2);
    knowledgeFact(
        k, "kuramoto-identical-synchrony", "Identical-frequency synchronous manifold", "equilibria",
        "At zero frequency spread every equal-phase configuration is stationary. Positive coupling "
        "contracts phase differences locally; negative coupling reverses transverse stability. "
        "This guide makes no population-limit threshold claim.",
        R"(s=0:\ \theta_1=\cdots=\theta_N,\quad \lambda_{transverse}=-K)",
        "Zero frequency spread and external driver; N>1 for transverse directions", "Native phases",
        refs, free && spread == 0);
    k.equilibriumGuides = [=](const auto &lo, const auto &hi) {
        std::vector<EquilibriumGuide> result;
        if (!free || lo.size() != n || hi.size() != n)
            return result;
        if (spread == 0) {
            double l = *std::max_element(lo.begin(), lo.end()),
                   h = *std::min_element(hi.begin(), hi.end());
            if (l <= h)
                result.push_back(
                    {"kuramoto-identical-synchrony", std::vector<double>(n, .5 * l + .5 * h)});
        }
        if (n == 2 && coupling != 0 && std::abs(2 * spread / coupling) <= 1) {
            const double angle = std::asin(2 * spread / coupling);
            for (double base : {angle, pi - angle}) {
                double a = std::ceil((lo[1] - hi[0] - base) / (2 * pi)),
                       b = std::floor((hi[1] - lo[0] - base) / (2 * pi));
                for (std::size_t j = 0; j < 4096 && result.size() < 4096 && a + double(j) <= b;
                     ++j) {
                    double delta = base + 2 * pi * (a + double(j));
                    double l = std::max(lo[0], lo[1] - delta), h = std::min(hi[0], hi[1] - delta);
                    if (l <= h) {
                        double phase = .5 * l + .5 * h;
                        result.push_back({"kuramoto-relative-locking", {phase, phase + delta}});
                    }
                }
            }
        }
        return result;
    };
    return k;
}
DynamicalKnowledge springPendulumKnowledge(const SystemConfiguration &c) {
    DynamicalKnowledge k;
    k.researchStatus = "reviewed_native_geometry_energy_and_synchronous_limit";
    const auto p = c.parameters;
    const double m1 = p.at("mass1_kg"), m2 = p.at("mass2_kg"), l1 = p.at("length1_m"),
                 l2 = p.at("length2_m"), pivot = p.at("pivot_separation_m"),
                 rest = p.at("spring_rest_length_m"), spring = p.at("spring_constant_n_m"),
                 g = p.at("gravity_m_s2");
    const bool free = c.driver.amplitude == 0, valid = m1 > 0 && m2 > 0 && l1 > 0 && l2 > 0;
    const bool conservative = free && p.at("damping1_n_m_s") == 0 && p.at("damping2_n_m_s") == 0;
    const std::vector<std::string> refs{"https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf",
                                        "https://dlmf.nist.gov/22.19"};
    knowledgeFact(k, "spring-pendulum-force-balance", "Spring-pendulum stationary geometry",
                  "equilibria",
                  "Equilibria solve both angular force balances including spring prestress. "
                  "Vertical rods supply seeds only when the connecting spring exerts no angular "
                  "force. Coincident bobs are outside the native smooth domain.",
                  R"(\omega_1=\omega_2=0,\quad \partial_{\theta_1}U=\partial_{\theta_2}U=0)",
                  "Positive masses/lengths, noncoincident bobs and zero driver",
                  "Native absolute angles", refs, free && valid);
    knowledgeFact(
        k, "spring-pendulum-energy", "Spring-pendulum total mechanical energy", "first-integral",
        "Use the actual bob-to-bob spring length in the elastic energy; linear normal-mode "
        "superposition applies only after a small-angle approximation.",
        R"(H=\tfrac12\sum_jm_jl_j^2\omega_j^2-\sum_jm_jgl_j\cos\theta_j+\tfrac12k(d-d_0)^2)",
        "Zero damping and external driver; positive masses/lengths", "Native dimensional energy",
        refs, conservative && valid);
    if (conservative && valid)
        k.invariantGuides.push_back(
            {"spring-pendulum-energy", "Total mechanical energy", [=](const auto &x) {
                 double distance = std::hypot(pivot + l2 * std::sin(x[2]) - l1 * std::sin(x[0]),
                                              l1 * std::cos(x[0]) - l2 * std::cos(x[2]));
                 return .5 * (m1 * l1 * l1 * x[1] * x[1] + m2 * l2 * l2 * x[3] * x[3]) -
                        g * (m1 * l1 * std::cos(x[0]) + m2 * l2 * std::cos(x[2])) +
                        .5 * spring * (distance - rest) * (distance - rest);
             }});
    knowledgeFact(k, "spring-pendulum-synchrony", "Exact common pendulum motion", "analytic-limit",
                  "Equal angles and velocities keep an unstretched spring at the fixed pivot "
                  "separation when rod lengths agree. Both bobs then follow the same full "
                  "nonlinear pendulum, even at nonzero spring stiffness and unequal masses.",
                  R"(\theta_1=\theta_2=\theta,\quad \theta''+(g/l)\sin\theta=0)",
                  "Equal positive lengths; rest length equals positive pivot separation; zero "
                  "damping/input; identical initial angles/velocities",
                  "Native absolute angles", refs,
                  conservative && valid && l1 == l2 && pivot > 0 && rest == pivot);
    k.equilibriumGuides = [=](const auto &lo, const auto &hi) {
        std::vector<EquilibriumGuide> result;
        if (!free || !valid || lo.size() != 4 || hi.size() != 4 || lo[1] > 0 || hi[1] < 0 ||
            lo[3] > 0 || hi[3] < 0)
            return result;
        double a = std::ceil(lo[0] / pi), b = std::floor(hi[0] / pi), u = std::ceil(lo[2] / pi),
               v = std::floor(hi[2] / pi);
        for (std::size_t i = 0, count = 0; i < 4096 && a + double(i) <= b; ++i)
            for (std::size_t j = 0; j < 4096 && u + double(j) <= v && count++ < 4096; ++j) {
                double x = (a + double(i)) * pi, y = (u + double(j)) * pi;
                double distance = std::hypot(pivot, l1 * std::cos(x) - l2 * std::cos(y));
                if (distance > 1e-12 && (spring == 0 || pivot == 0 ||
                                         std::abs(distance - rest) <= 1e-13 * (1 + distance)))
                    result.push_back({"spring-pendulum-force-balance", {x, 0, y, 0}});
            }
        return result;
    };
    return k;
}
} // namespace dynsyslab::core::detail
