#include <dynsyslab/core/SignalMetrics.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

// Solves the symmetric positive system M x = r (dimension d) by Gaussian elimination with
// partial pivoting.
std::vector<double> solve(std::vector<double> M, std::vector<double> r, std::size_t d) {
    for (std::size_t c = 0; c < d; ++c) {
        std::size_t pivot = c;
        for (std::size_t i = c + 1; i < d; ++i) if (std::abs(M[i * d + c]) > std::abs(M[pivot * d + c])) pivot = i;
        if (std::abs(M[pivot * d + c]) < 1e-300) throw std::invalid_argument("tone fit is singular (too few samples)");
        if (pivot != c) { for (std::size_t j = 0; j < d; ++j) std::swap(M[c * d + j], M[pivot * d + j]); std::swap(r[c], r[pivot]); }
        for (std::size_t i = c + 1; i < d; ++i) {
            const double f = M[i * d + c] / M[c * d + c];
            for (std::size_t j = c; j < d; ++j) M[i * d + j] -= f * M[c * d + j];
            r[i] -= f * r[c];
        }
    }
    std::vector<double> x(d);
    for (std::size_t i = d; i-- > 0;) {
        double sum = r[i];
        for (std::size_t j = i + 1; j < d; ++j) sum -= M[i * d + j] * x[j];
        x[i] = sum / M[i * d + i];
    }
    return x;
}

}  // namespace

ToneFit fitTone(const std::vector<double>& time, const std::vector<double>& values, double frequency,
                std::size_t harmonics) {
    if (time.size() != values.size()) throw std::invalid_argument("tone fit needs matching time and values");
    if (!(frequency > 0.0) || !harmonics) throw std::invalid_argument("tone fit needs a positive frequency and harmonic count");
    const std::size_t d = 1 + 2 * harmonics;
    if (time.size() < d) throw std::invalid_argument("tone fit needs more samples than parameters");
    constexpr double twoPi = 6.283185307179586476925286766559;
    std::vector<double> M(d * d, 0.0), r(d, 0.0), basis(d);
    const auto fill = [&](double t) {
        basis[0] = 1.0;
        for (std::size_t k = 1; k <= harmonics; ++k) {
            basis[2 * k - 1] = std::cos(twoPi * frequency * double(k) * t);
            basis[2 * k] = std::sin(twoPi * frequency * double(k) * t);
        }
    };
    for (std::size_t s = 0; s < time.size(); ++s) {
        fill(time[s]);
        for (std::size_t i = 0; i < d; ++i) {
            r[i] += basis[i] * values[s];
            for (std::size_t j = 0; j < d; ++j) M[i * d + j] += basis[i] * basis[j];
        }
    }
    const auto x = solve(M, r, d);
    ToneFit fit;
    fit.offset = x[0];
    double higher = 0.0;
    for (std::size_t k = 1; k <= harmonics; ++k) {
        const double a = x[2 * k - 1], b = x[2 * k];  // a cos + b sin = A sin(. + phi)
        fit.amplitudes.push_back(std::hypot(a, b));
        fit.phases.push_back(std::atan2(a, b));
        if (k > 1) higher += a * a + b * b;
    }
    double squares = 0.0;
    for (std::size_t s = 0; s < time.size(); ++s) {
        fill(time[s]);
        double model = 0.0;
        for (std::size_t i = 0; i < d; ++i) model += basis[i] * x[i];
        squares += (values[s] - model) * (values[s] - model);
    }
    fit.residualRms = std::sqrt(squares / double(time.size()));
    fit.harmonicContent = fit.amplitudes[0] > 0.0 ? std::sqrt(higher) / fit.amplitudes[0]
                                                  : std::numeric_limits<double>::quiet_NaN();
    return fit;
}

double settlingTime(const std::vector<double>& time, const std::vector<double>& error, double band) {
    if (time.size() != error.size() || time.empty()) return std::numeric_limits<double>::quiet_NaN();
    if (!(band > 0.0)) throw std::invalid_argument("settling band must be positive");
    std::size_t last = time.size();
    for (std::size_t k = 0; k < time.size(); ++k) if (!(std::abs(error[k]) <= band)) last = k;
    if (last == time.size()) return time.front();
    if (last + 1 >= time.size()) return std::numeric_limits<double>::quiet_NaN();
    return time[last + 1];
}

std::vector<double> lowPass(const std::vector<double>& time, const std::vector<double>& values, double cutoff) {
    if (time.size() != values.size()) throw std::invalid_argument("low-pass needs matching time and values");
    if (!(cutoff > 0.0)) throw std::invalid_argument("low-pass cutoff must be positive");
    std::vector<double> result(values.size());
    constexpr double twoPi = 6.283185307179586476925286766559;
    for (std::size_t k = 0; k < values.size(); ++k) {
        if (!k) { result[k] = values[k]; continue; }
        const double alpha = 1.0 - std::exp(-twoPi * cutoff * (time[k] - time[k - 1]));
        result[k] = result[k - 1] + alpha * (values[k] - result[k - 1]);
    }
    return result;
}

double autocorrelation(const std::vector<double>& values, std::size_t lag) {
    if (values.size() <= lag + 1) return std::numeric_limits<double>::quiet_NaN();
    double mean = 0.0;
    for (double v : values) mean += v;
    mean /= double(values.size());
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t k = 0; k < values.size(); ++k) {
        denominator += (values[k] - mean) * (values[k] - mean);
        if (k + lag < values.size()) numerator += (values[k] - mean) * (values[k + lag] - mean);
    }
    return denominator > 0.0 ? numerator / denominator : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace dynsyslab::core
