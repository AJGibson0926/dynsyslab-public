#include "AnalyticQuadrature.hpp"
#include "AnalyticSystems.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
namespace dynsyslab::core::detail {
namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
class FreeDoublePendulum {
    double a_, b_, c_, energy_, momentum_, k2_, direction_, u0_, offset_, t0_ = 0, angle0_ = 0;
    AnalyticState initial_;
    std::unique_ptr<PeriodicPrimitive> clock_, angle_;
    double q(double cosine) const { return a_ * c_ - b_ * b_ * cosine * cosine; }
    double d(double cosine) const { return a_ + c_ + 2 * b_ * cosine; }
    double delta(double u) const {
        if (k2_ < 1)
            return 2 * std::atan2(std::sqrt(k2_) * std::sin(u),
                                  std::sqrt((1 - k2_) + k2_ * std::pow(std::cos(u), 2)));
        if (k2_ == 1)
            return direction_ *
                   (u >= 0 ? pi - 4 * std::atan(std::exp(-u)) : -pi + 4 * std::atan(std::exp(u)));
        return direction_ * u;
    }
    double cosine(double u) const {
        if (k2_ < 1)
            return 1 - 2 * k2_ * std::pow(std::sin(u), 2);
        if (k2_ == 1)
            return 1 - 2 * std::pow(std::tanh(u), 2);
        return std::cos(u);
    }
    double deltaRate(double u) const {
        if (k2_ < 1)
            return 2 * std::sqrt(k2_) * std::cos(u) /
                   std::sqrt((1 - k2_) + k2_ * std::pow(std::cos(u), 2));
        if (k2_ == 1) {
            const double z = std::exp(-std::abs(u));
            return direction_ * 4 * z / (1 + z * z);
        }
        return direction_;
    }
    double clockRate(double u) const {
        const double co = cosine(u), base = q(co) / (2 * energy_ * b_);
        if (k2_ < 1)
            return std::sqrt(base / ((1 - k2_) + k2_ * std::pow(std::cos(u), 2)));
        if (k2_ == 1)
            return std::sqrt(base);
        return std::sqrt(base / (4 * ((k2_ - 1) + std::pow(std::cos(u / 2), 2))));
    }
    double angleRate(double u) const {
        const double co = cosine(u);
        return (momentum_ * clockRate(u) - (a_ + b_ * co) * deltaRate(u)) / d(co);
    }
    double separatrixParameter(double elapsed) const {
        const double minimum = std::sqrt((a_ * c_ - b_ * b_) / (2 * energy_ * b_));
        const double maximum = std::sqrt(a_ * c_ / (2 * energy_ * b_));
        double lo = u0_ + elapsed / maximum, hi = u0_ + elapsed / minimum;
        if (lo > hi)
            std::swap(lo, hi);
        double u = (lo + hi) / 2;
        for (int i = 0; i < 60; ++i) {
            const double residual =
                analyticIntegral([this](double s) { return clockRate(s); }, u0_, u) - elapsed;
            if (std::abs(residual) < 2e-12 * (1 + std::abs(elapsed)))
                return u;
            if (residual > 0)
                hi = u;
            else
                lo = u;
            const double next = u - residual / clockRate(u);
            u = next > lo && next < hi ? next : (lo + hi) / 2;
        }
        throw std::runtime_error("Double-pendulum separatrix clock inversion did not converge");
    }

  public:
    FreeDoublePendulum(double a, double b, double c, AnalyticState initial)
        : a_(a), b_(b), c_(c), initial_(std::move(initial)) {
        const double dd = initial_[0] - initial_[2], v = initial_[1] - initial_[3],
                     co = std::cos(dd);
        energy_ = .5 * (a_ * initial_[1] * initial_[1] + 2 * b_ * co * initial_[1] * initial_[3] +
                        c_ * initial_[3] * initial_[3]);
        momentum_ = (a_ + b_ * co) * initial_[1] + (c_ + b_ * co) * initial_[3];
        if (!(energy_ > 0) || !std::isfinite(energy_) || !std::isfinite(momentum_))
            throw std::domain_error(
                "Double-pendulum invariants exceed the finite positive-energy domain");
        // This form avoids cancellation in 1-epsilon for small librations.
        k2_ = std::pow(std::sin(dd / 2), 2) + q(co) * v * v / (8 * energy_ * b_);
        if (!(k2_ > 0) || !std::isfinite(k2_))
            throw std::domain_error("Unresolved double-pendulum relative-motion scale");
        direction_ = v < 0 ? -1 : 1;
        const double principal = std::remainder(dd, 2 * pi);
        if (k2_ < 1) {
            const double argument = std::sin(principal / 2) / std::sqrt(k2_);
            if (std::abs(argument) > 1 + 2e-12)
                throw std::domain_error("Inconsistent double-pendulum turning point");
            u0_ = std::asin(std::clamp(argument, -1.0, 1.0));
            if (v < 0)
                u0_ = pi - u0_;
        } else if (k2_ == 1) {
            u0_ = std::asinh(std::tan(direction_ * principal / 2));
        } else
            u0_ = direction_ * principal;
        if (k2_ != 1) {
            if (u0_ < 0)
                u0_ += 2 * pi;
            clock_ = std::make_unique<PeriodicPrimitive>([this](double u) { return clockRate(u); },
                                                         2 * pi, true);
            angle_ = std::make_unique<PeriodicPrimitive>([this](double u) { return angleRate(u); },
                                                         2 * pi);
            t0_ = clock_->value(u0_);
            angle0_ = angle_->value(u0_);
        }
        offset_ = dd - delta(u0_);
    }
    AnalyticState evaluate(double elapsed) const {
        if (!std::isfinite(elapsed))
            throw std::invalid_argument("Nonfinite double-pendulum time");
        if (elapsed == 0)
            return initial_;
        const double u = clock_ ? clock_->inverse(t0_ + elapsed) : separatrixParameter(elapsed);
        const double co = cosine(u), v = deltaRate(u) / clockRate(u);
        const double theta2 =
            initial_[2] +
            (angle_ ? angle_->value(u) - angle0_
                    : analyticIntegral([this](double s) { return angleRate(s); }, u0_, u));
        const double w2 = (momentum_ - (a_ + b_ * co) * v) / d(co);
        return {theta2 + delta(u) + offset_, w2 + v, theta2, w2};
    }
};
} // namespace
std::vector<AnalyticSolutionFamily> doublePendulumAnalyticSolutions(const SystemConfiguration &c) {
    if (c.identifier != "double_pendulum")
        return {};
    AnalyticSolutionDescriptor descriptor{
        "double-pendulum-gravity-free",
        "Analytic gravity-free double pendulum",
        R"(\dot\delta^2=(2HD-L^2)/Q,\quad \dot\theta_2=L/D-(E/D)\dot\delta,\quad D=A+C+2B\cos\delta,\ Q=AC-B^2\cos^2\delta,\ E=A+B\cos\delta)",
        "Zero gravity, both damping coefficients zero, no forcing; positive masses and lengths. "
        "Includes relative libration, rotation, separatrix, rest and collinear rigid rotation.",
        "Exact energy/angular-momentum reduction; real hyperelliptic quadrature and monotone "
        "clock inversion, with turning points regularized. Unwrapped angles are retained. "
        "Lauricella FD integral representation is supported separately by the native "
        "special-functions API; "
        "this evaluator does not approximate gravity-free motion with an ODE time stepper.",
        {"https://doi.org/10.1007/s00332-002-0514-0", "https://arxiv.org/abs/2104.13211"}};
    return {
        {descriptor, [c](const AnalyticState &x, double t0) -> AnalyticTrajectory {
             const auto &p = c.parameters;
             if (p.at("gravity_m_s2") != 0 || p.at("damping1_n_m_s") != 0 ||
                 p.at("damping2_n_m_s") != 0 || c.driver.amplitude != 0)
                 return {{}, "Requires zero gravity, zero damping and no forcing"};
             const double m1 = p.at("mass1_kg"), m2 = p.at("mass2_kg"), l1 = p.at("length1_m"),
                          l2 = p.at("length2_m");
             const double a = (m1 + m2) * l1 * l1, b = m2 * l1 * l2, cc = m2 * l2 * l2;
             if (!(m1 > 0 && m2 > 0 && l1 > 0 && l2 > 0 && a * cc - b * b > 0) ||
                 !std::isfinite(a * cc) || !std::isfinite(b * b))
                 return {
                     {},
                     "Requires finite positive masses, lengths and resolvable inertia determinant"};
             if (x[1] == 0 && x[3] == 0)
                 return {[x](double) { return x; }, {}};
             if (x[1] == x[3] && std::remainder(x[0] - x[2], pi) == 0)
                 return {[x, t0](double t) {
                             return AnalyticState{x[0] + x[1] * (t - t0), x[1],
                                                  x[2] + x[3] * (t - t0), x[3]};
                         },
                         {}};
             auto orbit = std::make_shared<FreeDoublePendulum>(a, b, cc, x);
             return {[orbit, t0](double t) { return orbit->evaluate(t - t0); }, {}};
         }}};
}
} // namespace dynsyslab::core::detail
