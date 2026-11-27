#include "IntegratorBackend.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

// Radau IIA of order 5 (three-stage implicit Runge-Kutta collocation) with simplified Newton
// iterations on the transformed real/complex systems, the embedded error estimate and the
// predictive step controller of
//   E. Hairer and G. Wanner, Solving Ordinary Differential Equations II: Stiff and
//   Differential-Algebraic Problems, 2nd ed. (Springer 1996), Sec. IV.8 (code RADAU5).
// Constants and the iteration policy follow that section; the Jacobian is a forward
// difference of the controlled right-hand side, refreshed when Newton converges slowly.

namespace dynsyslab::core::detail {
namespace {

const double s6 = std::sqrt(6.0);
const double nodes[3] = {(4.0 - s6) / 10.0, (4.0 + s6) / 10.0, 1.0};
const double errorWeights[3] = {(-13.0 - 7.0 * s6) / 3.0, (-13.0 + 7.0 * s6) / 3.0, -1.0 / 3.0};
const double muReal = 3.0 + std::cbrt(9.0) - std::cbrt(3.0);
const std::complex<double> muComplex{3.0 + 0.5 * (std::cbrt(3.0) - std::cbrt(9.0)),
                                     -0.5 * (std::pow(3.0, 5.0 / 6.0) + std::pow(3.0, 7.0 / 6.0))};
const double T[3][3] = {{0.09443876248897524, -0.14125529502095421, 0.03002919410514742},
                        {0.25021312296533332, 0.20412935229379994, -0.38294211275726192},
                        {1.0, 1.0, 0.0}};
const double TI[3][3] = {{4.17871859155190428, 0.32768282076106237, 0.52337644549944951},
                         {-4.17871859155190428, -0.32768282076106237, 0.47662355450055044},
                         {0.50287263494578682, -2.57192694985560522, 0.59603920482822492}};
// Dense-output collocation polynomial (used as the Newton predictor of the next step).
const double P[3][3] = {{13.0 / 3.0 + 7.0 * s6 / 3.0, -23.0 / 3.0 - 22.0 * s6 / 3.0, 10.0 / 3.0 + 5.0 * s6},
                        {13.0 / 3.0 - 7.0 * s6 / 3.0, -23.0 / 3.0 + 22.0 * s6 / 3.0, 10.0 / 3.0 - 5.0 * s6},
                        {1.0 / 3.0, -8.0 / 3.0, 10.0 / 3.0}};
constexpr int newtonIterations = 6;

template <class Scalar>
class DenseLU {
public:
    void factor(std::vector<Scalar> matrix, std::size_t n) {
        a_ = std::move(matrix); n_ = n; pivot_.resize(n);
        for (std::size_t k = 0; k < n; ++k) {
            std::size_t best = k;
            for (std::size_t r = k + 1; r < n; ++r)
                if (std::abs(a_[r * n + k]) > std::abs(a_[best * n + k])) best = r;
            pivot_[k] = best;
            if (best != k) for (std::size_t c = 0; c < n; ++c) std::swap(a_[k * n + c], a_[best * n + c]);
            if (a_[k * n + k] == Scalar{}) throw std::runtime_error("Radau Newton matrix is singular");
            for (std::size_t r = k + 1; r < n; ++r) {
                const Scalar f = a_[r * n + k] / a_[k * n + k];
                a_[r * n + k] = f;
                for (std::size_t c = k + 1; c < n; ++c) a_[r * n + c] -= f * a_[k * n + c];
            }
        }
    }
    void solve(std::vector<Scalar>& b) const {
        for (std::size_t k = 0; k < n_; ++k) if (pivot_[k] != k) std::swap(b[k], b[pivot_[k]]);
        for (std::size_t r = 0; r < n_; ++r) for (std::size_t c = 0; c < r; ++c) b[r] -= a_[r * n_ + c] * b[c];
        for (std::size_t r = n_; r-- > 0;) {
            for (std::size_t c = r + 1; c < n_; ++c) b[r] -= a_[r * n_ + c] * b[c];
            b[r] /= a_[r * n_ + r];
        }
    }

private:
    std::vector<Scalar> a_;
    std::vector<std::size_t> pivot_;
    std::size_t n_{};
};

class RadauStepper final : public Stepper {
public:
    RadauStepper(const IntegrationOptions& options, std::size_t n)
        : options_(options), n_(n), f_(n), jacobian_(n * n), z_(3 * n), q_(3 * n),
          newtonTolerance_(std::max(10.0 * std::numeric_limits<double>::epsilon()
                                    / options.relativeTolerance,
                                    std::min(0.03, std::sqrt(options.relativeTolerance)))) {}

    void restart(RhsContext& context, double time, const std::vector<double>& state) override {
        if (time == options_.initialTime) {
            hAbs_ = options_.initialStep > 0.0 ? options_.initialStep
                                               : 0.01 * (options_.finalTime - options_.initialTime);
            hOld_ = 0.0; errorOld_ = 0.0;
        }
        context.rhs(time, state.data(), f_.data());
        computeJacobian(context, time, state);
        haveDense_ = false;
        factoredFor_ = 0.0;
    }

    void advance(RhsContext& context, double time, double end, std::vector<double>& y) override {
        const double maximum = options_.maximumStep > 0.0 ? options_.maximumStep
                                                          : (std::numeric_limits<double>::infinity)();
        std::vector<double> next(n_), error(n_), work(n_);
        while (end - time > 16.0 * std::numeric_limits<double>::epsilon()
                                * std::max({1.0, std::abs(time), std::abs(end)})) {
            const double minimum = 10.0 * std::abs(std::nextafter(time, end) - time);
            double h = hAbs_;
            if (h > maximum) { h = maximum; hOld_ = 0.0; errorOld_ = 0.0; }
            bool accepted = false, rejected = false;
            int iterations = 0; double rate = 0.0, errorNorm = 0.0;
            while (!accepted) {
                if (h < minimum) throw std::runtime_error("native Radau step underflow");
                if (time + h > end) h = end - time;
                const double tNew = (end - time - h) <= 0.0 ? end : time + h;
                h = tNew - time;
                predictor(h);
                bool converged = false;
                while (!converged) {
                    if (factoredFor_ != h) factorize(h);
                    converged = collocation(context, time, y, h, iterations, rate);
                    if (!converged) {
                        if (jacobianCurrent_) break;
                        computeJacobian(context, time, y);
                        factoredFor_ = 0.0;
                    }
                }
                if (!converged) { h *= 0.5; factoredFor_ = 0.0; context.countStep(false); continue; }
                for (std::size_t i = 0; i < n_; ++i) next[i] = y[i] + z_[2 * n_ + i];
                std::vector<double> rhsError(n_);
                for (std::size_t i = 0; i < n_; ++i) {
                    double ze = 0.0;
                    for (int s = 0; s < 3; ++s) ze += errorWeights[s] * z_[s * n_ + i];
                    rhsError[i] = f_[i] + ze / h;
                }
                error = rhsError; luReal_.solve(error);
                errorNorm = norm(error, y, next);
                const double safety = 0.9 * (2.0 * newtonIterations + 1.0) / (2.0 * newtonIterations + iterations);
                if (rejected && errorNorm > 1.0) {
                    for (std::size_t i = 0; i < n_; ++i) work[i] = y[i] + error[i];
                    std::vector<double> f2(n_);
                    context.rhs(time, work.data(), f2.data());
                    for (std::size_t i = 0; i < n_; ++i) {
                        double ze = 0.0;
                        for (int s = 0; s < 3; ++s) ze += errorWeights[s] * z_[s * n_ + i];
                        error[i] = f2[i] + ze / h;
                    }
                    luReal_.solve(error);
                    errorNorm = norm(error, y, next);
                }
                if (errorNorm > 1.0) {
                    h *= std::max(0.2, safety * predictFactor(h, errorNorm));
                    factoredFor_ = 0.0; rejected = true;
                    context.countStep(false);
                } else {
                    accepted = true;
                    const bool refresh = iterations > 2 && rate > 1.0e-3;
                    double factor = std::min(10.0, safety * predictFactor(h, errorNorm));
                    if (!refresh && factor < 1.2) factor = 1.0;
                    else factoredFor_ = 0.0;
                    hOld_ = h; errorOld_ = errorNorm;
                    hAbs_ = h * factor;
                    for (std::size_t i = 0; i < n_; ++i)
                        for (int k = 0; k < 3; ++k) {
                            double sum = 0.0;
                            for (int s = 0; s < 3; ++s) sum += z_[s * n_ + i] * P[s][k];
                            q_[i * 3 + k] = sum;
                        }
                    haveDense_ = true; denseStep_ = h;
                    time = tNew;
                    y = next;
                    context.countStep(true);
                    context.accept(time, y.data());
                    context.rhs(time, y.data(), f_.data());
                    if (refresh) computeJacobian(context, time, y);
                    else jacobianCurrent_ = false;
                }
            }
        }
    }

private:
    double norm(const std::vector<double>& v, const std::vector<double>& a,
                const std::vector<double>& b) const {
        double sum = 0.0;
        for (std::size_t i = 0; i < n_; ++i) {
            const double scaled = v[i] / toleranceScale(options_, a[i], b[i]);
            sum += scaled * scaled;
        }
        return std::sqrt(sum / static_cast<double>(n_));
    }
    double predictFactor(double h, double errorNorm) const {
        double multiplier = 1.0;
        if (hOld_ > 0.0 && errorOld_ > 0.0 && errorNorm > 0.0)
            multiplier = h / hOld_ * std::pow(errorOld_ / errorNorm, 0.25);
        return errorNorm == 0.0 ? 10.0 : std::min(1.0, multiplier) * std::pow(errorNorm, -0.25);
    }
    void computeJacobian(RhsContext& context, double time, const std::vector<double>& y) {
        std::vector<double> shifted(y), f1(n_);
        for (std::size_t c = 0; c < n_; ++c) {
            const double step = std::sqrt(std::numeric_limits<double>::epsilon())
                * std::max(std::abs(y[c]), options_.absoluteTolerance / options_.relativeTolerance);
            const double delta = step > 0.0 ? step : 1.0e-8;
            shifted[c] = y[c] + delta;
            context.rhs(time, shifted.data(), f1.data());
            for (std::size_t r = 0; r < n_; ++r) jacobian_[r * n_ + c] = (f1[r] - f_[r]) / delta;
            shifted[c] = y[c];
        }
        jacobianCurrent_ = true;
    }
    void factorize(double h) {
        std::vector<double> real(n_ * n_);
        std::vector<std::complex<double>> complex(n_ * n_);
        for (std::size_t r = 0; r < n_; ++r)
            for (std::size_t c = 0; c < n_; ++c) {
                const double identity = r == c ? 1.0 : 0.0;
                real[r * n_ + c] = muReal / h * identity - jacobian_[r * n_ + c];
                complex[r * n_ + c] = muComplex / h * identity - jacobian_[r * n_ + c];
            }
        luReal_.factor(std::move(real), n_);
        luComplex_.factor(std::move(complex), n_);
        factoredFor_ = h;
    }
    void predictor(double h) {
        if (!haveDense_) { std::fill(z_.begin(), z_.end(), 0.0); return; }
        // Extrapolate the previous collocation polynomial to the new nodes.
        for (int s = 0; s < 3; ++s) {
            const double x = 1.0 + nodes[s] * h / denseStep_;
            for (std::size_t i = 0; i < n_; ++i) {
                double power = x, value = 0.0;
                for (int k = 0; k < 3; ++k) { value += q_[i * 3 + k] * power; power *= x; }
                // Relative to the new base point y_new = y_old + Z_old[2].
                z_[s * n_ + i] = value - q_[i * 3 + 0] - q_[i * 3 + 1] - q_[i * 3 + 2];
            }
        }
    }
    bool collocation(RhsContext& context, double t, const std::vector<double>& y, double h,
                     int& iterations, double& rate) {
        std::vector<double> w(3 * n_), stage(n_), f(3 * n_), dReal(n_);
        std::vector<std::complex<double>> dComplex(n_);
        for (int s = 0; s < 3; ++s)
            for (std::size_t i = 0; i < n_; ++i) {
                double sum = 0.0;
                for (int k = 0; k < 3; ++k) sum += TI[s][k] * z_[k * n_ + i];
                w[s * n_ + i] = sum;
            }
        std::vector<double> scale(n_);
        for (std::size_t i = 0; i < n_; ++i) scale[i] = toleranceScale(options_, y[i], y[i]);
        double previous = -1.0; rate = 0.0;
        for (int k = 0; k < newtonIterations; ++k) {
            iterations = k + 1;
            for (int s = 0; s < 3; ++s) {
                for (std::size_t i = 0; i < n_; ++i) stage[i] = y[i] + z_[s * n_ + i];
                context.rhs(t + nodes[s] * h, stage.data(), f.data() + s * n_);
            }
            if (!std::all_of(f.begin(), f.end(), [](double v) { return std::isfinite(v); }))
                return false;
            for (std::size_t i = 0; i < n_; ++i) {
                double fr = 0.0; std::complex<double> fc{};
                for (int s = 0; s < 3; ++s) {
                    fr += TI[0][s] * f[s * n_ + i];
                    fc += std::complex<double>(TI[1][s], TI[2][s]) * f[s * n_ + i];
                }
                dReal[i] = fr - muReal / h * w[i];
                dComplex[i] = fc - muComplex / h * std::complex<double>(w[n_ + i], w[2 * n_ + i]);
            }
            luReal_.solve(dReal); luComplex_.solve(dComplex);
            double sum = 0.0;
            for (std::size_t i = 0; i < n_; ++i) {
                const double a = dReal[i] / scale[i], b = dComplex[i].real() / scale[i],
                             c = dComplex[i].imag() / scale[i];
                sum += a * a + b * b + c * c;
            }
            const double dNorm = std::sqrt(sum / (3.0 * static_cast<double>(n_)));
            if (previous >= 0.0) rate = previous > 0.0 ? dNorm / previous : 0.0;
            if (previous >= 0.0 && (rate >= 1.0
                || std::pow(rate, newtonIterations - k) / (1.0 - rate) * dNorm > newtonTolerance_))
                return false;
            for (std::size_t i = 0; i < n_; ++i) {
                w[i] += dReal[i]; w[n_ + i] += dComplex[i].real(); w[2 * n_ + i] += dComplex[i].imag();
            }
            for (int s = 0; s < 3; ++s)
                for (std::size_t i = 0; i < n_; ++i) {
                    double sumZ = 0.0;
                    for (int q = 0; q < 3; ++q) sumZ += T[s][q] * w[q * n_ + i];
                    z_[s * n_ + i] = sumZ;
                }
            if (dNorm == 0.0 || (previous >= 0.0 && rate / (1.0 - rate) * dNorm < newtonTolerance_))
                return true;
            previous = dNorm;
        }
        return false;
    }

    const IntegrationOptions& options_;
    std::size_t n_;
    std::vector<double> f_, jacobian_, z_, q_;
    DenseLU<double> luReal_;
    DenseLU<std::complex<double>> luComplex_;
    double newtonTolerance_;
    double hAbs_{}, hOld_{}, errorOld_{}, factoredFor_{}, denseStep_{};
    bool jacobianCurrent_{false}, haveDense_{false};
};

}  // namespace

std::unique_ptr<Stepper> makeRadauStepper(const IntegrationOptions& options, std::size_t n) {
    return std::make_unique<RadauStepper>(options, n);
}

}  // namespace dynsyslab::core::detail
