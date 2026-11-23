#include <algorithm>
#include <cmath>
#include <dynsyslab/core/SpecialFunctions.hpp>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core::special {
QuarticOrbit::QuarticOrbit(std::array<double, 5> c, double x, double v)
    : coefficients_(c), x_(x), v_(v) {
    for (const auto value : c)
        if (!std::isfinite(value))
            throw std::invalid_argument("Nonfinite quartic coefficient");
    const double f = ((((c[4] * x + c[3]) * x + c[2]) * x + c[1]) * x + c[0]);
    if (!std::isfinite(x) || !std::isfinite(v) ||
        std::abs(f - v * v) > 1e-10 * (1 + std::abs(f) + v * v))
        throw std::invalid_argument(
            "Quartic initial velocity is incompatible with the first integral");
    f1_ = ((4 * c[4] * x + 3 * c[3]) * x + 2 * c[2]) * x + c[1];
    f2_ = (12 * c[4] * x + 6 * c[3]) * x + 2 * c[2];
    f3_ = 24 * c[4] * x + 6 * c[3];
    f4_ = 24 * c[4];
    const double a = c[4], b = c[3] / 4, d = c[2] / 6, e = c[1] / 4, h = c[0];
    g2_ = a * h - 4 * b * e + 3 * d * d;
    g3_ = a * d * h + 2 * b * d * e - d * d * d - a * e * e - b * b * h;
    if(c[1]==0&&c[3]==0&&c[0]>0&&c[4]<0){
        // x=A cn(nu*t+u0 | m), from (cn')^2=(1-cn^2)(1-m+m cn^2).
        // DLMF 22.13. The state and its derivative share one stable chart.
        const double quartic=-c[4],disc=std::sqrt(c[2]*c[2]+4*quartic*c[0]);
        if(!std::isfinite(disc))throw std::overflow_error("Even quartic scaling overflow");
        amplitude_=std::sqrt(c[2]>0?(disc+c[2])/(2*quartic):2*c[0]/(disc-c[2]));frequency_=std::sqrt(disc);
        parameter_=quartic*amplitude_*amplitude_/disc;
        const double cn=std::clamp(x/amplitude_,-1.,1.),dn=std::sqrt(1-parameter_+parameter_*cn*cn);
        const double sn=-v/(amplitude_*frequency_*dn);
        phase_=ellipticF(std::atan2(sn,cn),parameter_);period_=4*ellipticK(parameter_)/frequency_;
        evenJacobi_=true;return;
    }
    const double scale =
        std::max({1., std::pow(std::abs(g2_), .25), std::pow(std::abs(g3_), 1. / 6.)});
    period_ = weierstrass(.137 / scale, g2_, g3_).realPeriod;
}
double QuarticOrbit::realPeriod() const { return period_; }
std::array<double, 2> QuarticOrbit::evaluate(double elapsed) const {
    if (!std::isfinite(elapsed))
        throw std::invalid_argument("Nonfinite analytic time");
    if(evenJacobi_){const double t=std::remainder(elapsed,period_);const auto j=jacobi(phase_+frequency_*t,parameter_);
        return {amplitude_*j.cn,-amplitude_*frequency_*j.sn*j.dn};}
    if (v_ == 0 && f1_ == 0)
        return {x_, v_};
    double t = std::isfinite(period_) ? std::remainder(elapsed, period_) : elapsed;
    const double scale = std::max({1., std::pow(std::abs(g2_), .25),
                                   std::pow(std::abs(g3_), 1. / 6.), std::sqrt(std::abs(f2_))});
    if (std::abs(t) * scale < 1e-5) {
        const double a = f1_ / 2, j = f2_ * v_ / 2, k = f3_ * v_ * v_ / 2 + f2_ * f1_ / 4;
        return {x_ + t * (v_ + t * (a / 2 + t * (j / 6 + t * k / 24))),
                v_ + t * (a + t * (j / 2 + t * k / 6))};
    }
    const auto w = weierstrass(t, g2_, g3_);
    const double z = w.value - f2_ / 24;
    const double n = f1_ * z / 2 + v_ * v_ * f3_ / 24 - v_ * w.derivative;
    const double d = 2 * z * z - v_ * v_ * f4_ / 48;
    if (d == 0)
        throw std::domain_error("Quartic analytic solution has a pole");
    const double nd = f1_ * w.derivative / 2 - v_ * (6 * w.value * w.value - g2_ / 2),
                 dd = 4 * z * w.derivative;
    const double x = x_ + n / d, v = (nd * d - n * dd) / (d * d);
    if (!std::isfinite(x) || !std::isfinite(v))
        throw std::domain_error("Quartic analytic solution is nonfinite");
    return {x, v};
}
std::vector<double> affineExponential(const std::vector<double> &a, const std::vector<double> &b,
                                      const std::vector<double> &x, double t) {
    const auto n = x.size(), m = n + 1;
    if (n == 0 || a.size() != n * n || b.size() != n || !std::isfinite(t))
        throw std::invalid_argument("Affine exponential shape/time mismatch");
    std::vector<double> z(m * m), sum(m * m), term(m * m);
    double norm = 0;
    for (std::size_t i = 0; i < n; ++i) {
        double row = std::abs(b[i] * t);
        z[i * m + n] = b[i] * t;
        for (std::size_t j = 0; j < n; ++j) {
            z[i * m + j] = a[i * n + j] * t;
            row += std::abs(z[i * m + j]);
        }
        norm = std::max(norm, row);
    }
    if (!std::isfinite(norm))
        throw std::domain_error("Affine exponential overflow");
    const int squarings = norm > .5 ? static_cast<int>(std::ceil(std::log2(norm / .5))) : 0;
    for (auto &q : z)
        q = std::ldexp(q, -squarings);
    for (std::size_t i = 0; i < m; ++i)
        sum[i * m + i] = term[i * m + i] = 1;
    auto multiply = [m](const auto &p, const auto &q) {
        std::vector<double> r(m * m);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t k = 0; k < m; ++k)
                for (std::size_t j = 0; j < m; ++j)
                    r[i * m + j] += p[i * m + k] * q[k * m + j];
        return r;
    };
    for (int k = 1; k <= 24; ++k) {
        term = multiply(term, z);
        double size = 0;
        for (std::size_t i = 0; i < term.size(); ++i) {
            term[i] /= k;
            sum[i] += term[i];
            size = std::max(size, std::abs(term[i]));
        }
        if (size < 2e-17)
            break;
    }
    for (int k = 0; k < squarings; ++k)
        sum = multiply(sum, sum);
    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = sum[i * m + n];
        for (std::size_t j = 0; j < n; ++j)
            out[i] += sum[i * m + j] * x[j];
        if (!std::isfinite(out[i]))
            throw std::domain_error("Affine exponential result overflow");
    }
    return out;
}
} // namespace dynsyslab::core::special
