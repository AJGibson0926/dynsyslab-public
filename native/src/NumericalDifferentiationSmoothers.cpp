// Model-based derivative estimators; see NumericalDifferentiationSmoothers.hpp for the papers.
#include "NumericalDifferentiationSmoothers.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core::differentiation_detail {
namespace {
constexpr double kPi = 3.14159265358979323846;
}  // namespace

// ---- integrated white-noise model of order d: state (y, y', ..., y^(d-1))
void integratedModel(std::size_t d, double dt, double q, std::vector<double>& F, std::vector<double>& Q) {
    auto fact = [](std::size_t k) { double f = 1; for (std::size_t j = 2; j <= k; ++j) f *= double(j); return f; };
    F.assign(d * d, 0.0); Q.assign(d * d, 0.0);
    for (std::size_t i = 0; i < d; ++i) for (std::size_t j = i; j < d; ++j) F[i * d + j] = std::pow(dt, double(j - i)) / fact(j - i);
    for (std::size_t i = 0; i < d; ++i) for (std::size_t j = 0; j < d; ++j) {
        const double p = double(2 * d - 1 - i - j);
        Q[i * d + j] = q * std::pow(dt, p) / (p * fact(d - 1 - i) * fact(d - 1 - j));
    }
}
std::vector<double> matmul(const std::vector<double>& A, const std::vector<double>& B, std::size_t d) {
    std::vector<double> C(d * d, 0.0);
    for (std::size_t i = 0; i < d; ++i) for (std::size_t k = 0; k < d; ++k) { const double a = A[i * d + k]; if (a == 0) continue;
        for (std::size_t j = 0; j < d; ++j) C[i * d + j] += a * B[k * d + j]; }
    return C;
}
std::vector<double> transpose(const std::vector<double>& A, std::size_t d) {
    std::vector<double> T(d * d); for (std::size_t i = 0; i < d; ++i) for (std::size_t j = 0; j < d; ++j) T[j * d + i] = A[i * d + j]; return T;
}
std::vector<double> inverse(std::vector<double> A, std::size_t d) {
    std::vector<double> I(d * d, 0.0); for (std::size_t i = 0; i < d; ++i) I[i * d + i] = 1;
    for (std::size_t c = 0; c < d; ++c) {
        std::size_t p = c; for (std::size_t r = c + 1; r < d; ++r) if (std::abs(A[r * d + c]) > std::abs(A[p * d + c])) p = r;
        for (std::size_t k = 0; k < d; ++k) { std::swap(A[c * d + k], A[p * d + k]); std::swap(I[c * d + k], I[p * d + k]); }
        const double piv = A[c * d + c];
        for (std::size_t k = 0; k < d; ++k) { A[c * d + k] /= piv; I[c * d + k] /= piv; }
        for (std::size_t r = 0; r < d; ++r) if (r != c) { const double f = A[r * d + c];
            for (std::size_t k = 0; k < d; ++k) { A[r * d + k] -= f * A[c * d + k]; I[r * d + k] -= f * I[c * d + k]; } }
    }
    return I;
}

double kalmanQ(const DerivativeSpec& spec, double r, double dt, std::size_t d) {
    if (spec.processNoise > 0) return spec.processNoise;
    if (!(spec.cutoff > 0)) throw std::invalid_argument("kalman needs cutoff or process_noise");
    // RTS smoother of the d-integrated model == smoothing spline: gain 1/(1 + (w/wc)^(2d)), wc^(2d) = q / (r dt).
    return r * dt * std::pow(2 * kPi * spec.cutoff, 2.0 * double(d));
}


DerivativeResult kalmanSmoother(const std::vector<double>& y, double dt, const DerivativeSpec& spec) {
    const std::size_t n = y.size(), d = std::clamp<std::size_t>(spec.order, 2, 5);
    const double scale = std::max(rmsAboutMean(y), 1e-300);
    double r = spec.measurementNoise > 0 ? spec.measurementNoise : std::pow(estimateWhiteNoiseStd(y), 2);
    r = std::max(r, std::pow(1e-7 * scale, 2));
    std::vector<double> F, Q; integratedModel(d, dt, kalmanQ(spec, r, dt, d), F, Q);
    const auto Ft = transpose(F, d);
    std::vector<double> xp(n * d), Pp(n * d * d), xf(n * d), Pf(n * d * d);
    std::vector<double> x(d, 0.0), P(d * d, 0.0); x[0] = y[0];
    const double wc = 2 * kPi * std::max(spec.cutoff, 1.0 / (double(n) * dt));
    for (std::size_t k = 0; k < d; ++k) P[k * d + k] = 1e2 * scale * scale * std::pow(wc, 2.0 * double(k));
    P[0] = r;
    for (std::size_t i = 0; i < n; ++i) {
        if (i) {  // predict
            std::vector<double> xn(d, 0.0); for (std::size_t a = 0; a < d; ++a) for (std::size_t b = 0; b < d; ++b) xn[a] += F[a * d + b] * x[b];
            x = xn; P = matmul(matmul(F, P, d), Ft, d); for (std::size_t k = 0; k < d * d; ++k) P[k] += Q[k];
        }
        std::copy(x.begin(), x.end(), xp.begin() + std::ptrdiff_t(i * d)); std::copy(P.begin(), P.end(), Pp.begin() + std::ptrdiff_t(i * d * d));
        const double S = P[0] + r, innov = y[i] - x[0];  // update with H = e1
        std::vector<double> K(d); for (std::size_t a = 0; a < d; ++a) K[a] = P[a * d] / S;
        for (std::size_t a = 0; a < d; ++a) x[a] += K[a] * innov;
        std::vector<double> Pn(P); for (std::size_t a = 0; a < d; ++a) for (std::size_t b = 0; b < d; ++b) Pn[a * d + b] -= K[a] * P[b];
        P = Pn;
        std::copy(x.begin(), x.end(), xf.begin() + std::ptrdiff_t(i * d)); std::copy(P.begin(), P.end(), Pf.begin() + std::ptrdiff_t(i * d * d));
    }
    std::vector<double> xs(xf);  // RTS backward pass
    std::vector<double> Ps(Pf.end() - std::ptrdiff_t(d * d), Pf.end());
    for (std::size_t i = n - 1; i-- > 0;) {
        std::vector<double> Pfi(Pf.begin() + std::ptrdiff_t(i * d * d), Pf.begin() + std::ptrdiff_t((i + 1) * d * d));
        std::vector<double> Ppn(Pp.begin() + std::ptrdiff_t((i + 1) * d * d), Pp.begin() + std::ptrdiff_t((i + 2) * d * d));
        const auto C = matmul(matmul(Pfi, Ft, d), inverse(Ppn, d), d);
        std::vector<double> dx(d); for (std::size_t a = 0; a < d; ++a) dx[a] = xs[(i + 1) * d + a] - xp[(i + 1) * d + a];
        for (std::size_t a = 0; a < d; ++a) { double s = 0; for (std::size_t b = 0; b < d; ++b) s += C[a * d + b] * dx[b]; xs[i * d + a] = xf[i * d + a] + s; }
        std::vector<double> dP(d * d); for (std::size_t k = 0; k < d * d; ++k) dP[k] = Ps[k] - Ppn[k];
        const auto corr = matmul(matmul(C, dP, d), transpose(C, d), d);
        for (std::size_t k = 0; k < d * d; ++k) Ps[k] = Pfi[k] + corr[k];
    }
    DerivativeResult out; out.values.resize(n); out.derivative.resize(n);
    for (std::size_t i = 0; i < n; ++i) { out.values[i] = xs[i * d]; out.derivative[i] = xs[i * d + 1]; }
    return out;
}

// Condat (2013), direct 1-D TV denoising: argmin 1/2|x - y|^2 + lambda sum |x_{k+1} - x_k|.
std::vector<double> tvDenoise(const std::vector<double>& input, double lambda) {
    const std::ptrdiff_t width = std::ptrdiff_t(input.size());
    std::vector<double> output(input.size());
    if (width <= 0) return output;
    std::ptrdiff_t k = 0, k0 = 0, kplus = 0, kminus = 0;
    double umin = lambda, umax = -lambda, vmin = input[0] - lambda, vmax = input[0] + lambda;
    const double twolambda = 2.0 * lambda, minlambda = -lambda;
    for (;;) {
        while (k == width - 1) {
            if (umin < 0.0) {
                do output[std::size_t(k0++)] = vmin; while (k0 <= kminus);
                umax = (vmin = input[std::size_t(kminus = k = k0)]) + (umin = lambda) - vmax;
            } else if (umax > 0.0) {
                do output[std::size_t(k0++)] = vmax; while (k0 <= kplus);
                umin = (vmax = input[std::size_t(kplus = k = k0)]) + (umax = minlambda) - vmin;
            } else {
                vmin += umin / double(k - k0 + 1);
                do output[std::size_t(k0++)] = vmin; while (k0 <= k);
                return output;
            }
        }
        if ((umin += input[std::size_t(k + 1)] - vmin) < minlambda) {
            do output[std::size_t(k0++)] = vmin; while (k0 <= kminus);
            vmax = (vmin = input[std::size_t(kplus = kminus = k = k0)]) + twolambda;
            umin = lambda; umax = minlambda;
        } else if ((umax += input[std::size_t(k + 1)] - vmax) > lambda) {
            do output[std::size_t(k0++)] = vmax; while (k0 <= kplus);
            vmin = (vmax = input[std::size_t(kplus = kminus = k = k0)]) - twolambda;
            umin = lambda; umax = minlambda;
        } else {
            ++k;
            if (umin >= lambda) { vmin += (umin - lambda) / double((kminus = k) - k0 + 1); umin = lambda; }
            if (umax <= minlambda) { vmax += (umax + lambda) / double((kplus = k) - k0 + 1); umax = minlambda; }
        }
    }
}

// Integral formulation (Chartrand 2011): unknown u = y'. A u = trapezoid cumulative integral,
// data term |P(A u - y)|^2 with P removing the mean (the unknown y(0)); penalty alpha |D u|^2 or
// alpha |D u|_1 (lagged diffusivity). Solved by matrix-free conjugate gradients.
DerivativeResult integralRegularized(const std::vector<double>& y, double dt, const DerivativeSpec& spec, bool tv) {
    const std::size_t n = y.size();
    auto A = [&](const std::vector<double>& u) {
        std::vector<double> v(n, 0.0); for (std::size_t i = 1; i < n; ++i) v[i] = v[i - 1] + 0.5 * dt * (u[i - 1] + u[i]); return v; };
    auto At = [&](const std::vector<double>& v) {  // adjoint of the trapezoid integral
        std::vector<double> u(n, 0.0); double tail = 0;  // tail = sum_{i>=k+1} v_i
        for (std::size_t k = n; k-- > 0;) {
            const double after = tail;  // sum_{i>k} v_i
            u[k] = 0.5 * dt * ((k + 1 < n ? after : 0.0) + (k >= 1 ? after + v[k] : 0.0));
            tail += v[k];
        }
        return u; };
    auto center = [&](std::vector<double> v) { const double m = std::accumulate(v.begin(), v.end(), 0.0) / double(n); for (auto& x : v) x -= m; return v; };
    if (!(spec.cutoff > 0) && !(spec.weight > 0)) throw std::invalid_argument("integral differentiation needs cutoff or weight");
    const double alphaL2 = spec.weight > 0 && !tv ? spec.weight : std::pow(2 * kPi * std::max(spec.cutoff, 1e-300), -4.0);
    std::vector<double> u = centralDifference(y, dt);
    std::vector<double> w(n > 1 ? n - 1 : 0, 1.0);
    double alpha = alphaL2;
    if (tv) {
        if (spec.cutoff > 0) {  // warm start from the L2 solution; its slope variation sets the TV weight
            DerivativeSpec l2; l2.method = "tikhonov_integral"; l2.cutoff = spec.cutoff;
            u = differentiate(y, dt, l2).derivative;
        }
        double s = 0; for (std::size_t i = 0; i + 1 < n; ++i) s += std::pow((u[i + 1] - u[i]) / dt, 2);
        const double rmsDu = std::sqrt(s / double(std::max<std::size_t>(1, n - 1)));
        alpha = spec.weight > 0 ? spec.weight : alphaL2 * rmsDu;
    }
    const auto rhs = At(center(y));
    auto op = [&](const std::vector<double>& x) {
        auto v = At(center(A(x)));
        for (std::size_t i = 0; i + 1 < n; ++i) {  // + alpha D^T W D with D u = (u_{i+1}-u_i)/dt, sum weighted by dt
            // penalty alpha * sum w (Du)^2 dt against the data term sum (Au - y)^2 dt: dividing by dt
            // keeps the continuous weighting, so cutoff means the same as for the other smoothers
            const double g = alpha * w[i] * (x[i + 1] - x[i]) / (dt * dt);
            v[i] -= g; v[i + 1] += g;
        }
        return v; };
    const std::size_t outer = tv ? (spec.iterations ? spec.iterations : 20) : 1;
    for (std::size_t it = 0; it < outer; ++it) {
        if (tv) {
            double s = 0; for (std::size_t i = 0; i + 1 < n; ++i) s += std::abs(u[i + 1] - u[i]) / dt;
            const double eps = 1e-3 * s / double(std::max<std::size_t>(1, n - 1)) + 1e-300;
            for (std::size_t i = 0; i + 1 < n; ++i) w[i] = 1.0 / std::sqrt(std::pow((u[i + 1] - u[i]) / dt, 2) + eps * eps);
        }
        auto r = op(u); for (std::size_t i = 0; i < n; ++i) r[i] = rhs[i] - r[i];
        auto p = r; double rr = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
        const double r0 = std::sqrt(std::inner_product(rhs.begin(), rhs.end(), rhs.begin(), 0.0)) + 1e-300;
        const std::size_t maxCg = std::min<std::size_t>(4000, 4 * n);
        for (std::size_t k = 0; k < maxCg && std::sqrt(rr) > 1e-10 * r0; ++k) {
            const auto Ap = op(p); const double pAp = std::inner_product(p.begin(), p.end(), Ap.begin(), 0.0);
            if (!(pAp > 0)) break;
            const double a = rr / pAp;
            for (std::size_t i = 0; i < n; ++i) { u[i] += a * p[i]; r[i] -= a * Ap[i]; }
            const double rn = std::inner_product(r.begin(), r.end(), r.begin(), 0.0);
            for (std::size_t i = 0; i < n; ++i) p[i] = r[i] + (rn / rr) * p[i];
            rr = rn;
        }
    }
    DerivativeResult out; out.derivative = u;
    auto v = A(u); const double shift = std::accumulate(y.begin(), y.end(), 0.0) / double(n) - std::accumulate(v.begin(), v.end(), 0.0) / double(n);
    for (auto& x : v) x += shift; out.values = std::move(v);
    return out;
}
}  // namespace dynsyslab::core::differentiation_detail
