#include "IntegratorBackend.hpp"
#include "IntegratorTableaus.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Explicit steppers:
//   euler / rk4 take one step across each interval.
//   dormandPrince45: Dormand & Prince, J. Comput. Appl. Math. 6 (1980) 19-26, with the
//     historical DynSysLab max-norm controller (kept bit-for-bit for archived parity).
//   bogackiShampine23 and dormandPrince853 use the standard RMS-norm step controller of
//     Hairer, Norsett & Wanner, Solving ODEs I (1993), Sec. II.4, with FSAL reuse.

namespace dynsyslab::core::detail {
namespace {

class FixedStepper final : public Stepper {
public:
    FixedStepper(IntegrationMethod method, std::size_t n)
        : method_(method), k1_(n), k2_(n), k3_(n), k4_(n), work_(n) {}
    void restart(RhsContext&, double, const std::vector<double>&) override {}
    void advance(RhsContext& context, double time, double end, std::vector<double>& state) override {
        const std::size_t n = state.size();
        const double dt = end - time;
        context.countStep(true);
        context.rhs(time, state.data(), k1_.data());
        if (method_ == IntegrationMethod::euler) {
            for (std::size_t i = 0; i < n; ++i) state[i] += dt * k1_[i];
        } else {
            for (std::size_t i = 0; i < n; ++i) work_[i] = state[i] + 0.5 * dt * k1_[i];
            context.rhs(time + 0.5 * dt, work_.data(), k2_.data());
            for (std::size_t i = 0; i < n; ++i) work_[i] = state[i] + 0.5 * dt * k2_[i];
            context.rhs(time + 0.5 * dt, work_.data(), k3_.data());
            for (std::size_t i = 0; i < n; ++i) work_[i] = state[i] + dt * k3_[i];
            context.rhs(time + dt, work_.data(), k4_.data());
            for (std::size_t i = 0; i < n; ++i)
                state[i] += dt * (k1_[i] + 2.0 * k2_[i] + 2.0 * k3_[i] + k4_[i]) / 6.0;
        }
        context.accept(end, state.data());
    }

private:
    IntegrationMethod method_;
    std::vector<double> k1_, k2_, k3_, k4_, work_;
};

double roundoff(double time, double end) {
    return 16.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(time), std::abs(end)});
}

double firstStep(const IntegrationOptions& options) {
    if (options.initialStep > 0.0) return options.initialStep;
    const double span = options.finalTime - options.initialTime;
    return options.outputTimes.empty() ? std::min(options.outputStep, span / 100.0) : span / 100.0;
}

class DormandPrince45 final : public Stepper {
public:
    DormandPrince45(const IntegrationOptions& options, std::size_t n)
        : options_(options), k_(7, std::vector<double>(n)), work_(n), fifth_(n), fourth_(n) {}
    void restart(RhsContext&, double time, const std::vector<double>&) override {
        if (time == options_.initialTime) candidate_ = firstStep(options_);
    }
    void advance(RhsContext& context, double time, double end, std::vector<double>& state) override {
        const double maximum = options_.maximumStep > 0.0 ? options_.maximumStep
            : (options_.outputTimes.empty() ? options_.outputStep
                                            : (std::numeric_limits<double>::infinity)());
        while (time < end) {
            const double remaining = end - time;
            const double tiny = roundoff(time, end);
            if (remaining <= tiny) break;
            const double dt = std::min({candidate_, maximum, remaining});
            if (dt <= std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(time)))
                throw std::runtime_error("native adaptive integrator step underflow");
            bool accepted = false;
            const double suggested = step(context, time, dt, state, accepted);
            candidate_ = std::min(maximum, suggested);
            context.countStep(accepted);
            if (accepted) {
                time += dt;
                if (end - time <= tiny) time = end;
                context.accept(time, state.data());
            }
        }
    }

private:
    double step(RhsContext& c, double t, double dt, std::vector<double>& y, bool& accepted) {
        const std::size_t n = y.size();
        auto &k1 = k_[0], &k2 = k_[1], &k3 = k_[2], &k4 = k_[3], &k5 = k_[4], &k6 = k_[5], &k7 = k_[6];
        auto& w = work_;
        c.rhs(t, y.data(), k1.data());
        for (std::size_t i=0;i<n;++i) w[i]=y[i]+dt*(1.0/5.0)*k1[i];
        c.rhs(t+dt*(1.0/5.0),w.data(),k2.data());
        for (std::size_t i=0;i<n;++i) w[i]=y[i]+dt*(3.0/40.0*k1[i]+9.0/40.0*k2[i]);
        c.rhs(t+dt*(3.0/10.0),w.data(),k3.data());
        for (std::size_t i=0;i<n;++i) w[i]=y[i]+dt*(44.0/45.0*k1[i]-56.0/15.0*k2[i]+32.0/9.0*k3[i]);
        c.rhs(t+dt*(4.0/5.0),w.data(),k4.data());
        for (std::size_t i=0;i<n;++i) w[i]=y[i]+dt*(19372.0/6561.0*k1[i]-25360.0/2187.0*k2[i]+64448.0/6561.0*k3[i]-212.0/729.0*k4[i]);
        c.rhs(t+dt*(8.0/9.0),w.data(),k5.data());
        for (std::size_t i=0;i<n;++i) w[i]=y[i]+dt*(9017.0/3168.0*k1[i]-355.0/33.0*k2[i]+46732.0/5247.0*k3[i]+49.0/176.0*k4[i]-5103.0/18656.0*k5[i]);
        c.rhs(t+dt,w.data(),k6.data());
        for (std::size_t i=0;i<n;++i) fifth_[i]=y[i]+dt*(35.0/384.0*k1[i]+500.0/1113.0*k3[i]+125.0/192.0*k4[i]-2187.0/6784.0*k5[i]+11.0/84.0*k6[i]);
        c.rhs(t+dt,fifth_.data(),k7.data());
        double error = 0.0;
        for (std::size_t i=0;i<n;++i) {
            fourth_[i]=y[i]+dt*(5179.0/57600.0*k1[i]+7571.0/16695.0*k3[i]+393.0/640.0*k4[i]-92097.0/339200.0*k5[i]+187.0/2100.0*k6[i]+1.0/40.0*k7[i]);
            const double scale=options_.absoluteTolerance+options_.relativeTolerance*std::max(std::abs(y[i]),std::abs(fifth_[i]));
            error=std::max(error,std::abs(fifth_[i]-fourth_[i])/scale);
        }
        accepted = error <= 1.0;
        const double factor = error == 0.0 ? 5.0 : std::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
        if (accepted) y.swap(fifth_);
        return dt * factor;
    }

    const IntegrationOptions& options_;
    std::vector<std::vector<double>> k_;
    std::vector<double> work_, fifth_, fourth_;
    double candidate_{};
};

// Generic FSAL embedded pair: s stages, y_new = y + h sum b_i k_i, k_s = f(t + h, y_new).
class EmbeddedRungeKutta final : public Stepper {
public:
    EmbeddedRungeKutta(const IntegrationOptions& options, std::size_t n, bool dop853)
        : options_(options), dop853_(dop853), stages_(dop853 ? 12 : 3),
          k_(stages_ + 1, std::vector<double>(n)), work_(n), next_(n),
          exponent_(dop853 ? -1.0 / 8.0 : -1.0 / 3.0) {}
    void restart(RhsContext&, double time, const std::vector<double>&) override {
        if (time == options_.initialTime) candidate_ = firstStep(options_);
        fsal_ = false;
    }
    void advance(RhsContext& context, double time, double end, std::vector<double>& state) override {
        const double maximum = options_.maximumStep > 0.0 ? options_.maximumStep
                                                          : (std::numeric_limits<double>::infinity)();
        bool rejected = false;
        while (end - time > roundoff(time, end)) {
            const double h = std::min({candidate_, maximum, end - time});
            if (h <= 10.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(time)))
                throw std::runtime_error("native adaptive integrator step underflow");
            const double error = step(context, time, h, state);
            const bool accepted = error < 1.0;
            context.countStep(accepted);
            if (accepted) {
                double factor = error == 0.0 ? 10.0 : std::min(10.0, 0.9 * std::pow(error, exponent_));
                if (rejected) factor = std::min(1.0, factor);
                candidate_ = h * factor;
                time = end - time - h <= roundoff(time, end) ? end : time + h;
                state.swap(next_);
                k_[0].swap(k_[stages_]);
                fsal_ = true;
                rejected = false;
                context.accept(time, state.data());
            } else {
                candidate_ = h * std::max(0.2, 0.9 * std::pow(error, exponent_));
                rejected = true;
            }
        }
    }

private:
    double a(std::size_t i, std::size_t j) const {
        return dop853_ ? tableau::dop853A[i][j] : tableau::bs23A[i][j];
    }
    double b(std::size_t i) const { return dop853_ ? tableau::dop853B[i] : tableau::bs23B[i]; }
    double c(std::size_t i) const { return dop853_ ? tableau::dop853C[i] : tableau::bs23C[i]; }

    double step(RhsContext& context, double t, double h, const std::vector<double>& y) {
        const std::size_t n = y.size();
        if (!fsal_) context.rhs(t, y.data(), k_[0].data());
        for (std::size_t s = 1; s < stages_; ++s) {
            for (std::size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (std::size_t j = 0; j < s; ++j) sum += a(s, j) * k_[j][i];
                work_[i] = y[i] + h * sum;
            }
            context.rhs(t + c(s) * h, work_.data(), k_[s].data());
        }
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t j = 0; j < stages_; ++j) sum += b(j) * k_[j][i];
            next_[i] = y[i] + h * sum;
        }
        context.rhs(t + h, next_.data(), k_[stages_].data());
        double first = 0.0, second = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double scale = toleranceScale(options_, y[i], next_[i]);
            double e5 = 0.0, e3 = 0.0;
            for (std::size_t j = 0; j <= stages_; ++j) {
                if (dop853_) { e5 += tableau::dop853E5[j] * k_[j][i]; e3 += tableau::dop853E3[j] * k_[j][i]; }
                else e5 += tableau::bs23E[j] * k_[j][i];
            }
            e5 /= scale; e3 /= scale;
            first += e5 * e5; second += e3 * e3;
        }
        if (!std::isfinite(first) || !std::isfinite(second))
            return (std::numeric_limits<double>::infinity)();
        if (!dop853_) return h * std::sqrt(first / static_cast<double>(n));
        if (first == 0.0 && second == 0.0) return 0.0;
        // DOP853 error: |h| e5^2 / sqrt((e5^2 + 0.01 e3^2) n)  (Hairer et al. 1993, II.10).
        return h * first / std::sqrt((first + 0.01 * second) * static_cast<double>(n));
    }

    const IntegrationOptions& options_;
    bool dop853_;
    std::size_t stages_;
    std::vector<std::vector<double>> k_;
    std::vector<double> work_, next_;
    double exponent_;
    double candidate_{};
    bool fsal_{false};
};

}  // namespace

std::unique_ptr<Stepper> makeExplicitStepper(const IntegrationOptions& options, std::size_t n) {
    switch (options.method) {
    case IntegrationMethod::euler: case IntegrationMethod::rk4:
        return std::make_unique<FixedStepper>(options.method, n);
    case IntegrationMethod::dormandPrince45: return std::make_unique<DormandPrince45>(options, n);
    case IntegrationMethod::bogackiShampine23: return std::make_unique<EmbeddedRungeKutta>(options, n, false);
    case IntegrationMethod::dormandPrince853: return std::make_unique<EmbeddedRungeKutta>(options, n, true);
    default: throw std::logic_error("not an explicit ODE method");
    }
}

}  // namespace dynsyslab::core::detail
