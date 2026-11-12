#include "AnalyticQuadrature.hpp"
#include "AnalyticSystems.hpp"
#include <algorithm>
#include <cmath>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <memory>
#include <stdexcept>

namespace dynsyslab::core::detail {
namespace {
class ParametricBubble {
    special::QuarticOrbit orbit_;
    double power_, clockPower_, period_;
    std::unique_ptr<PeriodicPrimitive> clock_;
    double clockRate(double s) const {
        const auto q = orbit_.evaluate(s);
        if (!(q[0] > 0))
            throw std::domain_error("RPE elliptic branch left positive radius domain");
        return std::pow(q[0], power_ * clockPower_);
    }

  public:
    ParametricBubble(std::array<double, 5> polynomial, double q0, double dq0, double power,
                     double clockPower)
        : orbit_(polynomial, q0, dq0), power_(power), clockPower_(clockPower),
          period_(orbit_.realPeriod()) {
        if (!std::isfinite(period_) || period_ <= 0)
            throw std::domain_error("Degenerate RPE elliptic orbit requires its equilibrium limit");
        clock_ = std::make_unique<PeriodicPrimitive>([this](double s) { return clockRate(s); },
                                                     period_, true);
    }
    AnalyticState evaluate(double t) const {
        double s = std::fmod(clock_->inverse(t), period_);
        if (s < 0)
            s += period_;
        const auto q = orbit_.evaluate(s);
        const double r = std::pow(q[0], power_);
        const double v = power_ * std::pow(q[0], power_ - 1) * q[1] / std::pow(r, clockPower_);
        if (!(r > .01))
            throw std::domain_error("RPE analytic radius reached the model collapse boundary");
        return {r, v};
    }
};
} // namespace
std::vector<AnalyticSolutionFamily> rpeAnalyticSolutions(const SystemConfiguration &c) {
    if (c.identifier != "rpe")
        return {};
    AnalyticSolutionDescriptor d{
        "rpe-inviscid-elliptic",
        "Analytic RPE",
        R"(\kappa=4/3:\ r=q^{-1},\ q_s^2=-2Gq^4+Cq^3-Sq+\tfrac23(V-P),\ dt=r^2ds;\quad \kappa=3/2,S=0:\ r=q^{-2/3},\ q_s^2=-3Gq^3+\tfrac94 Cq^2+\tfrac32(V-P),\ dt=r^{5/2}ds)",
        "viscosity=0 and no forcing; kappa=4/3 with surface tension, or kappa=3/2 with surface "
        "tension=0; positive bounded gas-filled motion",
        "elliptic-parametric; adaptive quadrature and monotone inversion of the analytic clock",
        {"https://arxiv.org/abs/1409.6699", "https://arxiv.org/abs/1508.01157",
         "https://doi.org/10.1088/1361-6382/ac95f2"}};
    return {
        {d, [c](const AnalyticState &x, double t0) -> AnalyticTrajectory {
             if (x.size() != 2 || x[0] <= .01)
                 return {{}, "Positive radius and radial velocity required"};
             const auto &g = c.dimensionlessGroups;
             const double gas = g.at("gas_pressure_ratio"), vapor = g.at("vapor_pressure_ratio"),
                          external = g.at("ambient_pressure_ratio"), cap = g.at("capillarity");
             if (g.at("viscosity") != 0 || c.driver.amplitude != 0)
                 return {{}, "Zero viscosity and forcing required"};
             const double k = c.parameters.at("kappa"), a = 2 * (vapor - external) / 3;
             const bool four = std::abs(k - 4. / 3.) < 4e-16, three = k == 1.5 && cap == 0;
             if (!four && !three)
                 return {{},
                         "Supported polytropic cases are kappa=4/3, or kappa=3/2 without surface "
                         "tension"};
             if (!(gas > 0 && external > vapor && cap >= 0))
                 return {{}, "Bounded positive gas-filled branch required"};
             const double balance = gas * std::pow(x[0], -3 * k) + vapor - external - cap / x[0];
             if (x[1] == 0 && std::abs(balance) < 2e-15 * (1 + gas + external + cap))
                 return {[x](double) { return x; }, {}};
             const double c1 =
                 std::pow(x[0], 3) *
                 (x[1] * x[1] - 2 * gas / (3 * (1 - k)) * std::pow(x[0], -3 * k) - a + cap / x[0]);
             std::shared_ptr<ParametricBubble> orbit;
             if (four)
                 orbit = std::make_shared<ParametricBubble>(
                     std::array<double, 5>{a, -cap, 0, c1, -2 * gas}, 1 / x[0], -x[1], -1, 2);
             else
                 orbit = std::make_shared<ParametricBubble>(
                     std::array<double, 5>{2.25 * a, 0, 2.25 * c1, -3 * gas, 0},
                     std::pow(x[0], -1.5), -1.5 * x[1], -2. / 3., 2.5);
             return {[orbit, t0](double t) { return orbit->evaluate(t - t0); }, {}};
         }}};
}
} // namespace dynsyslab::core::detail
