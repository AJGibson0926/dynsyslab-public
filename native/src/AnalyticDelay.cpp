#include <algorithm>
#include <cmath>
#include <dynsyslab/core/Delay.hpp>
#include <stdexcept>
namespace dynsyslab::core {
namespace {
AnalyticSolutionFamily
family(const char *id, const char *label, const char *formula, const char *conditions,
       std::function<AnalyticTrajectory(const AnalyticState &, double)> build) {
    return {
        {id,
         label,
         formula,
         conditions,
         "Closed-form evaluation for a constant prehistory; native requested delay, not the "
         "rounded numerical grid delay.",
         {"https://doi.org/10.1126/science.267326", "https://doi.org/10.1186/s13662-015-0538-z"}},
        [build](const AnalyticState &x, double t0) -> AnalyticTrajectory {
            if (x.size() != 1 || !std::isfinite(x[0]) || !std::isfinite(t0))
                return {{}, "Requires one finite constant history value and initial time"};
            return build(x, t0);
        }};
}
} // namespace
std::vector<AnalyticSolutionFamily>
mackeyGlassAnalyticSolutions(const SystemConfiguration &requested) {
    const auto c = resolveMackeyGlass(requested);
    const double beta = c.parameters.at("beta"), gamma = c.parameters.at("gamma"),
                 delay = c.parameters.at("delay"), exponent = c.parameters.at("exponent");
    if (!std::isfinite(beta) || !std::isfinite(gamma) || !std::isfinite(delay) || delay <= 0 ||
        !std::isfinite(exponent) || c.driver.amplitude != 0)
        throw std::invalid_argument("Mackey-Glass analytic families require finite parameters, "
                                    "positive delay and no driver");
    std::vector<AnalyticSolutionFamily> result;
    result.push_back(family(
        "mackey-glass-linear-delay", "Analytic linear delay limit",
        R"(x'=-\gamma x+(\beta/2)x(t-\tau),\quad x(t)=h[R(t)+(\beta/2)\int_{\max(0,t-\tau)}^t R(s)ds])",
        "exponent=0; constant prehistory; finite delayed-exponential sum (50-digit working "
        "precision)",
        [=](const auto &x, double t0) -> AnalyticTrajectory {
            if (exponent != 0)
                return {{}, "Requires exponent=0"};
            return {[=](double t) {
                        return AnalyticState{
                            linearDelayConstantHistory(beta / 2, gamma, delay, x[0], t - t0)};
                    },
                    {}};
        }));
    result.push_back(family(
        "mackey-glass-no-feedback", "Analytic feedback-free decay", R"(x(t)=h e^{-\gamma(t-t_0)})",
        "beta=0; real nonsingular feedback domain; constant prehistory",
        [=](const auto &x, double t0) -> AnalyticTrajectory {
            if (beta != 0 ||
                (x[0] < 0 && (std::floor(exponent) != exponent || std::fmod(exponent, 2) != 0)))
                return {{}, "Requires beta=0 and a nonsingular real feedback domain"};
            return {[=](double t) {
                        return AnalyticState{t <= t0 ? x[0] : x[0] * std::exp(-gamma * (t - t0))};
                    },
                    {}};
        }));
    result.push_back(family(
        "mackey-glass-constant-equilibrium", "Analytic constant equilibrium",
        R"(x(t)=h,\quad\beta h/(1+h^n)-\gamma h=0)",
        "The supplied constant history satisfies the configured equilibrium relation exactly",
        [=](const auto &x, double) -> AnalyticTrajectory {
            const double forcing = beta * x[0] / (1 + std::pow(x[0], exponent));
            if (!std::isfinite(forcing) || forcing != gamma * x[0])
                return {{}, "History is not an exact equilibrium of the configured feedback"};
            return {[=](double) { return x; }, {}};
        }));
    result.push_back(family(
        "mackey-glass-first-delay", "Analytic first delay interval",
        R"(x(t)=h e^{-\gamma s}+\frac{\beta h}{1+h^n}\frac{1-e^{-\gamma s}}{\gamma},\quad s=t-t_0\in[0,\tau])",
        "Constant prehistory; finite feedback; valid only through the first requested delay "
        "interval; continuous gamma=0 limit",
        [=](const auto &x, double t0) -> AnalyticTrajectory {
            const double forcing = beta * x[0] / (1 + std::pow(x[0], exponent));
            if (!std::isfinite(forcing))
                return {{}, "History is outside the real finite feedback domain"};
            return {
                [=](double t) {
                    const double s = t - t0;
                    if (!std::isfinite(t) || s > delay)
                        throw std::domain_error(
                            "Analytic first-delay formula ends at initial_time + requested delay");
                    if (s <= 0)
                        return x;
                    const double integral = gamma == 0 ? s : -std::expm1(-gamma * s) / gamma;
                    return AnalyticState{x[0] * std::exp(-gamma * s) + forcing * integral};
                },
                {}};
        }));
    return result;
}
} // namespace dynsyslab::core
