#include <dynsyslab/core/NumericalDifferentiation.hpp>
#include <dynsyslab/core/DiscreteFourier.hpp>
#include "NumericalDifferentiationSmoothers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
namespace differentiation_detail {
double rmsAboutMean(const std::vector<double>& y) {
    if (y.empty()) return 0;
    const double mean = std::accumulate(y.begin(), y.end(), 0.0) / double(y.size());
    double s = 0; for (double v : y) s += (v - mean) * (v - mean);
    return std::sqrt(s / double(y.size()));
}
std::vector<double> centralDifference(const std::vector<double>& y, double dt) {
    const std::size_t n = y.size(); std::vector<double> d(n, 0.0);
    if (n < 2) return d;
    if (n == 2) { d[0] = d[1] = (y[1] - y[0]) / dt; return d; }
    for (std::size_t i = 1; i + 1 < n; ++i) d[i] = (y[i + 1] - y[i - 1]) / (2 * dt);
    d[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * dt);
    d[n - 1] = (3 * y[n - 1] - 4 * y[n - 2] + y[n - 3]) / (2 * dt);
    return d;
}
}  // namespace differentiation_detail
using namespace differentiation_detail;
namespace {
constexpr double kPi = 3.14159265358979323846;

// Odd (point) reflection about each end: preserves value and slope at the boundary.
std::vector<double> oddExtend(const std::vector<double>& y, std::size_t pad) {
    const std::size_t n = y.size(); pad = std::min(pad, n - 1);
    std::vector<double> e(n + 2 * pad);
    for (std::size_t k = 0; k < pad; ++k) e[pad - 1 - k] = 2 * y.front() - y[k + 1];
    std::copy(y.begin(), y.end(), e.begin() + std::ptrdiff_t(pad));
    for (std::size_t k = 0; k < pad; ++k) e[pad + n + k] = 2 * y.back() - y[n - 2 - k];
    return e;
}

// out_i = sum_k w[k] * y_{i + k - half}, edges by odd reflection (pads further by repetition if needed).
std::vector<double> correlate(const std::vector<double>& y, const std::vector<double>& w) {
    const std::size_t half = w.size() / 2, n = y.size();
    auto e = oddExtend(y, half);
    const std::size_t pad = (e.size() - n) / 2;
    std::vector<double> out(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0;
        for (std::size_t k = 0; k < w.size(); ++k) {
            const std::ptrdiff_t j = std::ptrdiff_t(i + pad + k) - std::ptrdiff_t(half);
            const std::ptrdiff_t c = std::clamp<std::ptrdiff_t>(j, 0, std::ptrdiff_t(e.size()) - 1);
            s += w[k] * e[std::size_t(c)];
        }
        out[i] = s;
    }
    return out;
}


// Weights w_j (one per offset) such that sum w_j y(s_j) is the `deriv`-th derivative at 0 of the
// least-squares polynomial of degree p through (s_j, y_j). Offsets in samples.
std::vector<double> polynomialWeights(const std::vector<double>& s, std::size_t p, std::size_t deriv) {
    const std::size_t m = p + 1, n = s.size();
    if (n < m) throw std::invalid_argument("polynomial window has fewer samples than coefficients");
    std::vector<double> G(m * m, 0.0);  // V^T V
    for (std::size_t a = 0; a < m; ++a) for (std::size_t b = 0; b < m; ++b)
        for (double sj : s) G[a * m + b] += std::pow(sj, double(a + b));
    // solve G c = e_deriv * deriv!
    std::vector<double> rhs(m, 0.0); double fact = 1; for (std::size_t k = 2; k <= deriv; ++k) fact *= double(k);
    rhs[deriv] = fact;
    for (std::size_t c = 0; c < m; ++c) {  // Gaussian elimination with partial pivoting
        std::size_t piv = c; for (std::size_t r = c + 1; r < m; ++r) if (std::abs(G[r * m + c]) > std::abs(G[piv * m + c])) piv = r;
        if (piv != c) { for (std::size_t k = 0; k < m; ++k) std::swap(G[c * m + k], G[piv * m + k]); std::swap(rhs[c], rhs[piv]); }
        for (std::size_t r = c + 1; r < m; ++r) {
            const double f = G[r * m + c] / G[c * m + c];
            for (std::size_t k = c; k < m; ++k) G[r * m + k] -= f * G[c * m + k];
            rhs[r] -= f * rhs[c];
        }
    }
    std::vector<double> coef(m);
    for (std::size_t c = m; c-- > 0;) {
        double v = rhs[c]; for (std::size_t k = c + 1; k < m; ++k) v -= G[c * m + k] * coef[k];
        coef[c] = v / G[c * m + c];
    }
    std::vector<double> w(n, 0.0);  // w_j = sum_a coef_a s_j^a  (G symmetric, so this is e^T G^-1 V^T)
    for (std::size_t j = 0; j < n; ++j) for (std::size_t a = 0; a < m; ++a) w[j] += coef[a] * std::pow(s[j], double(a));
    return w;
}

std::size_t oddAtLeast(double v, std::size_t minimum) {
    auto k = std::size_t(std::max(1.0, std::round(v)));
    k = std::max(k, minimum); if (k % 2 == 0) ++k; return k;
}

std::size_t sgWindow(const DerivativeSpec& spec, double dt, std::size_t order) {
    if (spec.window) return spec.window;
    if (!(spec.cutoff > 0)) throw std::invalid_argument("savitzky_golay needs window or cutoff");
    // Schafer (2011): -3 dB at f/f_nyquist ~ (N+1)/(3.2 M - 4.6), W = 2M + 1.
    const double fn = 2 * spec.cutoff * dt;
    const double M = ((double(order) + 1) / fn + 4.6) / 3.2;
    return oddAtLeast(2 * M + 1, order + 2);
}

// ---- Butterworth (bilinear transform) as cascaded sections: {b0,b1,b2,a1,a2}
std::vector<std::array<double, 5>> butterworthSections(std::size_t order, double cutoff, double dt) {
    if (!(cutoff > 0) || cutoff * dt >= 0.5) throw std::invalid_argument("butterworth cutoff must be in (0, Nyquist)");
    if (order < 1) order = 1;
    const double K = std::tan(kPi * cutoff * dt);
    std::vector<std::array<double, 5>> s;
    for (std::size_t k = 0; k < order / 2; ++k) {
        const double q = 1.0 / (2 * std::sin(kPi * (2.0 * double(k) + 1) / (2.0 * double(order))));
        const double norm = 1 / (1 + K / q + K * K), b0 = K * K * norm;
        s.push_back({b0, 2 * b0, b0, 2 * (K * K - 1) * norm, (1 - K / q + K * K) * norm});
    }
    if (order % 2) { const double norm = 1 / (1 + K); s.push_back({K * norm, K * norm, 0, (K - 1) * norm, 0}); }
    return s;
}
void filterInPlace(std::vector<double>& x, const std::vector<std::array<double, 5>>& sections) {
    for (const auto& c : sections) {
        double x1 = x.empty() ? 0 : x[0], x2 = x1, y1 = x1, y2 = x1;  // steady state at the first value
        for (auto& v : x) {
            const double y = c[0] * v + c[1] * x1 + c[2] * x2 - c[3] * y1 - c[4] * y2;
            x2 = x1; x1 = v; y2 = y1; y1 = y; v = y;
        }
    }
}

// ---- banded symmetric positive definite solve (I + lambda D_d^T D_d) z = y
std::vector<double> whittaker(const std::vector<double>& y, double lambda, std::size_t d) {
    const std::size_t n = y.size();
    if (n <= d) return y;
    std::vector<double> coeff(d + 1);  // difference stencil of order d: (-1)^(d-k) C(d,k)
    for (std::size_t k = 0; k <= d; ++k) {
        double c = 1; for (std::size_t j = 1; j <= k; ++j) c = c * double(d - k + j) / double(j);
        coeff[k] = ((d - k) % 2 ? -1.0 : 1.0) * c;
    }
    const std::size_t bw = d;  // band storage: A[i][j-i], j = i..i+bw
    std::vector<double> A(n * (bw + 1), 0.0);
    for (std::size_t r = 0; r + d < n; ++r)  // row r of D touches columns r..r+d
        for (std::size_t a = 0; a <= d; ++a) for (std::size_t b = a; b <= d; ++b)
            A[(r + a) * (bw + 1) + (b - a)] += lambda * coeff[a] * coeff[b];
    for (std::size_t i = 0; i < n; ++i) A[i * (bw + 1)] += 1.0;
    // banded Cholesky A = L L^T, stored in place (L[i][i-k] at A[(i-k)*(bw+1)+k])
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = std::min(i, bw); k >= 1; --k) {  // L[i][i-k]
            const std::size_t j = i - k;
            double s = A[j * (bw + 1) + k];
            for (std::size_t m = 1; m + k <= bw && m <= j; ++m)
                s -= A[(j - m) * (bw + 1) + m] * A[(j - m) * (bw + 1) + m + k];
            A[j * (bw + 1) + k] = s / A[j * (bw + 1)];
        }
        double s = A[i * (bw + 1)];
        for (std::size_t k = 1; k <= std::min(i, bw); ++k) s -= A[(i - k) * (bw + 1) + k] * A[(i - k) * (bw + 1) + k];
        if (!(s > 0)) throw std::runtime_error("whittaker system is not positive definite");
        A[i * (bw + 1)] = std::sqrt(s);
    }
    std::vector<double> z(y);
    for (std::size_t i = 0; i < n; ++i) {  // L w = y
        double s = z[i];
        for (std::size_t k = 1; k <= std::min(i, bw); ++k) s -= A[(i - k) * (bw + 1) + k] * z[i - k];
        z[i] = s / A[i * (bw + 1)];
    }
    for (std::size_t i = n; i-- > 0;) {  // L^T z = w
        double s = z[i];
        for (std::size_t k = 1; k <= bw && i + k < n; ++k) s -= A[i * (bw + 1) + k] * z[i + k];
        z[i] = s / A[i * (bw + 1)];
    }
    return z;
}


DerivativeResult spectral(const std::vector<double>& y, double dt, const DerivativeSpec& spec) {
    if (!(spec.cutoff > 0)) throw std::invalid_argument("spectral needs cutoff");
    const std::size_t n = y.size(), N = 2 * n;
    std::vector<std::complex<double>> z(N);
    for (std::size_t i = 0; i < n; ++i) { z[i] = y[i]; z[N - 1 - i] = y[i]; }  // even extension: no jump
    auto Z = discreteFourierTransform(z);
    auto D = Z;
    for (std::size_t k = 0; k < N; ++k) {
        const double kk = k <= N / 2 ? double(k) : double(k) - double(N);
        const double f = kk / (double(N) * dt), af = std::abs(f);
        double h = 1;  // raised-cosine taper from 0.75 fc to 1.25 fc (gain 1/2 at fc)
        if (af >= 1.25 * spec.cutoff) h = 0;
        else if (af > 0.75 * spec.cutoff) h = 0.5 * (1 + std::cos(kPi * (af - 0.75 * spec.cutoff) / (0.5 * spec.cutoff)));
        if (2 * k == N) h = 0;
        Z[k] *= h; D[k] = Z[k] * std::complex<double>(0, 2 * kPi * f);
    }
    const auto zs = discreteFourierTransform(Z, true), ds = discreteFourierTransform(D, true);
    DerivativeResult out; out.values.resize(n); out.derivative.resize(n);
    for (std::size_t i = 0; i < n; ++i) { out.values[i] = zs[i].real(); out.derivative[i] = ds[i].real(); }
    return out;
}

std::vector<double> holoborodkoWeights(std::size_t N) {  // c_k for k = 1..M, f' ~ sum c_k (f_k - f_-k) / h
    if (N < 5 || N % 2 == 0) throw std::invalid_argument("holoborodko points must be odd and >= 5");
    const std::size_t m = (N - 3) / 2, M = (N - 1) / 2;
    auto binom = [](std::ptrdiff_t a, std::ptrdiff_t b) { if (b < 0 || b > a) return 0.0; double c = 1; for (std::ptrdiff_t j = 1; j <= b; ++j) c = c * double(a - b + j) / double(j); return c; };
    std::vector<double> c(M);
    for (std::size_t k = 1; k <= M; ++k)
        c[k - 1] = (binom(std::ptrdiff_t(2 * m), std::ptrdiff_t(m) - std::ptrdiff_t(k) + 1) - binom(std::ptrdiff_t(2 * m), std::ptrdiff_t(m) - std::ptrdiff_t(k) - 1))
                   / std::pow(2.0, double(2 * m + 1));
    return c;
}
}  // namespace

std::vector<std::string> derivativeMethods() {
    return {"central", "central_high_order", "savitzky_golay", "gaussian", "holoborodko", "weak", "butterworth", "whittaker", "kalman",
            "spectral", "tv_denoise", "tikhonov_integral", "tvr",
            "backward", "backward2", "causal_savitzky_golay", "causal_butterworth", "causal_kalman", "levant"};
}
bool isCausalDerivativeMethod(const std::string& m) {
    return m == "backward" || m == "backward2" || m == "causal_savitzky_golay" || m == "causal_butterworth" || m == "causal_kalman" || m == "levant";
}

double estimateWhiteNoiseStd(const std::vector<double>& y) {
    if (y.size() < 8) return 0;
    std::vector<double> d; d.reserve(y.size() - 3);
    for (std::size_t i = 3; i < y.size(); ++i) d.push_back(y[i] - 3 * y[i - 1] + 3 * y[i - 2] - y[i - 3]);
    auto med = [](std::vector<double> v) { const auto mid = v.begin() + std::ptrdiff_t(v.size() / 2); std::nth_element(v.begin(), mid, v.end()); return *mid; };
    const double m = med(d); for (auto& v : d) v = std::abs(v - m);
    return 1.4826 * med(d) / std::sqrt(20.0);
}

std::vector<double> hampel(const std::vector<double>& y, std::size_t half, double threshold) {
    const std::size_t n = y.size(); std::vector<double> out(y), window;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t a = i >= half ? i - half : 0, b = std::min(n, i + half + 1);
        window.assign(y.begin() + std::ptrdiff_t(a), y.begin() + std::ptrdiff_t(b));
        const auto mid = window.begin() + std::ptrdiff_t(window.size() / 2);
        std::nth_element(window.begin(), mid, window.end()); const double med = *mid;
        for (auto& v : window) v = std::abs(v - med);
        std::nth_element(window.begin(), mid, window.end()); const double mad = 1.4826 * *mid;
        if (std::abs(y[i] - med) > threshold * mad && mad > 0) out[i] = med;
    }
    return out;
}

DerivativeResult differentiate(const std::vector<double>& given, double dt, const DerivativeSpec& spec) {
    const auto y = spec.despike ? hampel(given, spec.despike, spec.despikeThreshold) : given;
    const std::size_t n = y.size();
    if (!(dt > 0) || !std::isfinite(dt)) throw std::invalid_argument("differentiation needs a positive sample time");
    if (n < 3) throw std::invalid_argument("differentiation needs at least 3 samples");
    const auto& m = spec.method;
    DerivativeResult out; out.values = y;
    if (m == "central") { out.derivative = centralDifference(y, dt); return out; }
    if (m == "central_high_order") { out.derivative = highOrderCentralDifference(y, dt, spec.points ? spec.points : 9); return out; }
    if (m == "savitzky_golay") {
        const std::size_t p = std::max<std::size_t>(spec.order, 1), W = sgWindow(spec, dt, p);
        if (W % 2 == 0) throw std::invalid_argument("savitzky_golay window must be odd");
        std::vector<double> s(W); for (std::size_t j = 0; j < W; ++j) s[j] = double(j) - double(W / 2);
        auto w1 = polynomialWeights(s, p, 1), w0 = polynomialWeights(s, p, 0);
        for (auto& v : w1) v /= dt;
        out.derivative = correlate(y, w1); out.values = correlate(y, w0); return out;
    }
    if (m == "gaussian") {
        const double sigma = spec.sigma > 0 ? spec.sigma : (spec.cutoff > 0 ? std::sqrt(2 * std::log(2.0)) / (2 * kPi * spec.cutoff) : 0);
        if (!(sigma > 0)) throw std::invalid_argument("gaussian needs sigma or cutoff");
        const double ss = sigma / dt; const auto half = std::size_t(std::ceil(4 * ss));
        std::vector<double> g(2 * half + 1), w(2 * half + 1); double sg = 0, sk = 0;
        for (std::size_t j = 0; j < g.size(); ++j) { const double k = double(j) - double(half); g[j] = std::exp(-0.5 * k * k / (ss * ss)); sg += g[j]; sk += k * k * g[j]; }
        for (std::size_t j = 0; j < g.size(); ++j) { const double k = double(j) - double(half); w[j] = k * g[j] / (dt * sk); g[j] /= sg; }
        out.derivative = correlate(y, w); out.values = correlate(y, g); return out;
    }
    if (m == "holoborodko") {
        const std::size_t N = spec.points ? spec.points : 7, mm = (N - 3) / 2;
        std::size_t s = spec.stride;
        if (!s) {
            if (!(spec.cutoff > 0)) s = 1;
            else { const double x = std::acos(std::pow(2.0, -1.0 / (2.0 * double(std::max<std::size_t>(mm, 1)))));
                   s = std::size_t(std::max(1.0, std::round(2 * x / (2 * kPi * spec.cutoff) / dt))); }
        }
        const auto c = holoborodkoWeights(N); const std::size_t M = c.size();
        std::vector<double> w(2 * M * s + 1, 0.0);
        for (std::size_t k = 1; k <= M; ++k) { w[M * s + k * s] = c[k - 1] / (double(s) * dt); w[M * s - k * s] = -c[k - 1] / (double(s) * dt); }
        out.derivative = correlate(y, w); return out;
    }
    if (m == "weak") {
        const double p = spec.exponent;
        std::size_t W = spec.window;
        if (!W) {
            if (!(spec.cutoff > 0)) throw std::invalid_argument("weak needs window or cutoff");
            const double sigmaT = std::sqrt(2 * std::log(2.0)) / (2 * kPi * spec.cutoff);
            W = oddAtLeast(2 * sigmaT * std::sqrt(2 * p + 3) / dt + 1, 5);
        }
        const double half = double(W / 2);
        std::vector<double> w(W), phi(W); double norm = 0, sphi = 0;
        for (std::size_t j = 0; j < W; ++j) {
            const double s = (double(j) - half) / (half + 1);  // phi vanishes just outside the window
            phi[j] = std::pow(1 - s * s, p); sphi += phi[j];
            w[j] = 2 * p * s * std::pow(1 - s * s, p - 1);  // -phi'(s)
            norm += w[j] * (double(j) - half) * dt;
        }
        for (auto& v : w) v /= norm; for (auto& v : phi) v /= sphi;
        out.derivative = correlate(y, w); out.values = correlate(y, phi); return out;
    }
    if (m == "butterworth") {
        const auto sec = butterworthSections(spec.order ? spec.order : 4, spec.cutoff, dt);
        const std::size_t pad = std::min<std::size_t>(n - 1, std::size_t(std::ceil(6.0 / (spec.cutoff * dt))));
        auto e = oddExtend(y, pad); const std::size_t padUsed = (e.size() - n) / 2;
        filterInPlace(e, sec); std::reverse(e.begin(), e.end()); filterInPlace(e, sec); std::reverse(e.begin(), e.end());
        out.values.assign(e.begin() + std::ptrdiff_t(padUsed), e.begin() + std::ptrdiff_t(padUsed + n));
        out.derivative = centralDifference(out.values, dt); return out;
    }
    if (m == "whittaker") {
        const std::size_t d = std::max<std::size_t>(spec.penaltyOrder, 1);
        const double lambda = spec.lambda > 0 ? spec.lambda : (spec.cutoff > 0 ? std::pow(2 * std::sin(kPi * spec.cutoff * dt), -2.0 * double(d)) : 0);
        if (!(lambda > 0)) throw std::invalid_argument("whittaker needs lambda or cutoff");
        out.values = whittaker(y, lambda, d); out.derivative = centralDifference(out.values, dt); return out;
    }
    if (m == "kalman") return kalmanSmoother(y, dt, spec);
    if (m == "spectral") return spectral(y, dt, spec);
    if (m == "tv_denoise") {
        const double lambda = spec.weight > 0 ? spec.weight : 3 * estimateWhiteNoiseStd(y);
        out.values = lambda > 0 ? tvDenoise(y, lambda) : y; out.derivative = centralDifference(out.values, dt); return out;
    }
    if (m == "tikhonov_integral") return integralRegularized(y, dt, spec, false);
    if (m == "tvr") return integralRegularized(y, dt, spec, true);
    if (isCausalDerivativeMethod(m)) {
        CausalDifferentiator c(spec, dt, rmsAboutMean(y));
        out.derivative.resize(n); for (std::size_t i = 0; i < n; ++i) out.derivative[i] = c.step(y[i]);
        return out;
    }
    throw std::invalid_argument("unknown derivative method: " + m);
}

std::vector<double> despike(const std::vector<double>& y, std::size_t half, double threshold) { return hampel(y, half, threshold); }

CausalDifferentiator::CausalDifferentiator(const DerivativeSpec& spec, double dt, double signalScale) : method_(spec.method), dt_(dt) {
    if (!(dt > 0)) throw std::invalid_argument("causal differentiator needs a positive sample time");
    if (!isCausalDerivativeMethod(method_)) throw std::invalid_argument("not a causal derivative method: " + method_);
    const double scale = signalScale > 0 ? signalScale : 1.0;
    if (method_ == "causal_savitzky_golay") {
        const std::size_t p = std::max<std::size_t>(spec.order, 1);
        std::size_t W = spec.window ? spec.window : sgWindow(spec, dt, p);
        std::vector<double> s(W); for (std::size_t j = 0; j < W; ++j) s[j] = double(j) - double(W - 1);
        coefficients_ = polynomialWeights(s, p, 1); for (auto& v : coefficients_) v /= dt;
    } else if (method_ == "causal_butterworth") {
        for (const auto& c : butterworthSections(spec.order ? spec.order : 2, spec.cutoff, dt)) { b_.insert(b_.end(), c.begin(), c.end()); }
    } else if (method_ == "causal_kalman") {
        d_ = std::clamp<std::size_t>(spec.order, 2, 5);
        r_ = spec.measurementNoise > 0 ? spec.measurementNoise : std::pow(1e-3 * scale, 2);
        integratedModel(d_, dt, kalmanQ(spec, r_, dt, d_), F_, Q_);
        x_.assign(d_, 0.0); P_.assign(d_ * d_, 0.0);
        const double wc = 2 * kPi * std::max(spec.cutoff, 1e-6);
        for (std::size_t k = 0; k < d_; ++k) P_[k * d_ + k] = 1e2 * scale * scale * std::pow(wc, 2.0 * double(k));
    } else if (method_ == "levant") {
        L_ = spec.lipschitz > 0 ? spec.lipschitz : (spec.cutoff > 0 ? scale * std::pow(2 * kPi * spec.cutoff, 3) : 0);
        if (!(L_ > 0)) throw std::invalid_argument("levant needs lipschitz or cutoff");
        z_.assign(3, 0.0);
    }
}

double CausalDifferentiator::step(double y) {
    if (method_.empty()) throw std::logic_error("unconfigured differentiator");
    ++count_;
    if (method_ == "causal_kalman") {
        const std::size_t d = d_;
        if (count_ == 1) { x_[0] = y; P_[0] = r_; }
        else {
            std::vector<double> xn(d, 0.0); for (std::size_t a = 0; a < d; ++a) for (std::size_t b = 0; b < d; ++b) xn[a] += F_[a * d + b] * x_[b];
            x_ = xn; P_ = matmul(matmul(F_, P_, d), transpose(F_, d), d); for (std::size_t k = 0; k < d * d; ++k) P_[k] += Q_[k];
            const double S = P_[0] + r_, innov = y - x_[0];
            std::vector<double> K(d); for (std::size_t a = 0; a < d; ++a) K[a] = P_[a * d] / S;
            for (std::size_t a = 0; a < d; ++a) x_[a] += K[a] * innov;
            std::vector<double> Pn(P_); for (std::size_t a = 0; a < d; ++a) for (std::size_t b = 0; b < d; ++b) Pn[a * d + b] -= K[a] * P_[b];
            P_ = Pn;
        }
        value_ = x_[0]; return derivative_ = x_[1];
    }
    if (method_ == "levant") {
        if (count_ == 1) { z_ = {y, 0, 0}; value_ = y; return derivative_ = 0; }
        const int sub = 20; const double h = dt_ / sub;
        auto sgn = [](double v) { return double((v > 0) - (v < 0)); };
        for (int k = 0; k < sub; ++k) {  // Levant (2003), n = 2, lambda = (1.1, 1.5, 2)
            const double e0 = z_[0] - y;
            const double v0 = -2.0 * std::cbrt(L_) * std::pow(std::abs(e0), 2.0 / 3.0) * sgn(e0) + z_[1];
            const double e1 = z_[1] - v0;
            const double v1 = -1.5 * std::sqrt(L_) * std::sqrt(std::abs(e1)) * sgn(e1) + z_[2];
            const double v2 = -1.1 * L_ * sgn(z_[2] - v1);
            z_[0] += h * v0; z_[1] += h * v1; z_[2] += h * v2;
        }
        value_ = z_[0]; return derivative_ = z_[1];
    }
    double v = y;
    if (method_ == "causal_butterworth") {
        const std::size_t ns = b_.size() / 5;
        if (count_ == 1) { xs_.assign(2 * ns, y); ys_.assign(2 * ns, y); }
        for (std::size_t s = 0; s < ns; ++s) {
            const double* c = &b_[5 * s];
            const double out = c[0] * v + c[1] * xs_[2 * s] + c[2] * xs_[2 * s + 1] - c[3] * ys_[2 * s] - c[4] * ys_[2 * s + 1];
            xs_[2 * s + 1] = xs_[2 * s]; xs_[2 * s] = v; ys_[2 * s + 1] = ys_[2 * s]; ys_[2 * s] = out; v = out;
        }
    }
    history_.push_back(v);
    const std::size_t keep = std::max<std::size_t>(coefficients_.size(), 3);
    if (history_.size() > keep) history_.erase(history_.begin());
    const std::size_t h = history_.size();
    double d = 0;
    if (method_ == "causal_savitzky_golay" && h == coefficients_.size()) {
        for (std::size_t j = 0; j < h; ++j) d += coefficients_[j] * history_[j];
    } else if ((method_ == "backward2" || method_ == "causal_savitzky_golay") && h >= 3) {
        d = (3 * history_[h - 1] - 4 * history_[h - 2] + history_[h - 3]) / (2 * dt_);
    } else if (h >= 2) d = (history_[h - 1] - history_[h - 2]) / dt_;
    value_ = v; return derivative_ = d;
}

// ---------------- noise ----------------
StreamingNoise::StreamingNoise(const NoiseSpec& spec, double dt, double scale, double range) : dt_(dt) {
    std::vector<NoiseSpec> flat;
    if (spec.kind == "sum") flat = spec.parts; else if (spec.kind != "none" && !spec.kind.empty()) flat.push_back(spec);
    rng_ ^= spec.seed * 0xBF58476D1CE4E5B9ull;
    for (const auto& s : flat) {
        if (s.kind == "sum") throw std::invalid_argument("noise parts cannot nest");
        Part p; p.spec = s;
        p.sigma = s.relative ? s.level * scale : s.level;
        if (s.kind == "ar1") {
            if (!(s.correlationTime > 0)) throw std::invalid_argument("ar1 noise needs correlation_time");
            p.a = std::exp(-dt / s.correlationTime);
        } else if (s.kind == "quantization") {
            if (s.bits > 0) { if (!(range > 0)) throw std::invalid_argument("quantization by bits needs the signal range"); p.step = range / std::pow(2.0, s.bits); }
            else p.step = s.relative ? s.level * scale : s.level;
            if (!(p.step > 0)) throw std::invalid_argument("quantization needs bits or a step");
        } else if (s.kind == "random_walk") {
            p.a = p.sigma * std::sqrt(dt);  // level = RMS drift accumulated over one time unit
        } else if (s.kind == "impulse") {
            if (!(s.probability >= 0 && s.probability <= 1)) throw std::invalid_argument("impulse probability must be in [0, 1]");
        } else if (s.kind != "white" && s.kind != "uniform" && s.kind != "pink")
            throw std::invalid_argument("unknown noise kind: " + s.kind);
        parts_.push_back(p);
    }
    for (auto& p : parts_) {
        if (p.spec.kind == "ar1") p.state = p.sigma * normal();  // stationary start
        if (p.spec.kind == "pink") {  // calibrate Kellet's filter gain so the output RMS is sigma
            Part q = p; double s2 = 0; const int N = 1 << 15;
            for (int k = 0; k < N + 2000; ++k) {
                const double w = normal(); double* b = q.pink;
                b[0] = 0.99886 * b[0] + w * 0.0555179; b[1] = 0.99332 * b[1] + w * 0.0750759; b[2] = 0.96900 * b[2] + w * 0.1538520;
                b[3] = 0.86650 * b[3] + w * 0.3104856; b[4] = 0.55000 * b[4] + w * 0.5329522; b[5] = -0.7616 * b[5] - w * 0.0168980;
                const double o = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6] + w * 0.5362; b[6] = w * 0.115926;
                if (k >= 2000) s2 += o * o;
            }
            p.a = p.sigma / std::sqrt(s2 / N);
            std::copy(std::begin(q.pink), std::end(q.pink), std::begin(p.pink));
        }
    }
}
double StreamingNoise::uniform() {
    rng_ += 0x9E3779B97F4A7C15ull; std::uint64_t z = rng_;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; z ^= z >> 31;
    return (double(z >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}
double StreamingNoise::normal() {
    const double u1 = uniform(), u2 = uniform();
    return std::sqrt(-2 * std::log(u1)) * std::cos(2 * kPi * u2);
}
double StreamingNoise::apply(double y) {
    for (auto& p : parts_) {
        const auto& k = p.spec.kind;
        if (k == "white") y += p.sigma * normal();
        else if (k == "uniform") y += p.sigma * std::sqrt(3.0) * (2 * uniform() - 1);
        else if (k == "ar1") { p.state = p.a * p.state + std::sqrt(1 - p.a * p.a) * p.sigma * normal(); y += p.state; }
        else if (k == "random_walk") { p.state += p.a * normal(); y += p.state; }
        else if (k == "pink") {
            const double w = normal(); double* b = p.pink;
            b[0] = 0.99886 * b[0] + w * 0.0555179; b[1] = 0.99332 * b[1] + w * 0.0750759; b[2] = 0.96900 * b[2] + w * 0.1538520;
            b[3] = 0.86650 * b[3] + w * 0.3104856; b[4] = 0.55000 * b[4] + w * 0.5329522; b[5] = -0.7616 * b[5] - w * 0.0168980;
            const double o = b[0] + b[1] + b[2] + b[3] + b[4] + b[5] + b[6] + w * 0.5362; b[6] = w * 0.115926;
            y += p.a * o;
        } else if (k == "impulse") { if (uniform() < p.spec.probability) y += p.sigma * (uniform() < 0.5 ? -1 : 1) * (0.5 + 0.5 * uniform()); }
        else if (k == "quantization") y = p.step * std::round(y / p.step);
    }
    return y;
}

std::vector<double> addMeasurementNoise(const std::vector<double>& y, double dt, const NoiseSpec& spec, double scale) {
    if (spec.kind.empty() || spec.kind == "none") return y;
    const double s = scale > 0 ? scale : rmsAboutMean(y);
    const auto [lo, hi] = std::minmax_element(y.begin(), y.end());
    StreamingNoise noise(spec, dt, s, y.empty() ? 0 : *hi - *lo);
    std::vector<double> out(y.size());
    for (std::size_t i = 0; i < y.size(); ++i) out[i] = noise.apply(y[i]);
    return out;
}

}  // namespace dynsyslab::core
