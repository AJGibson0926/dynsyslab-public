#include <dynsyslab/core/PortableNormal.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

// Natural logarithm from frexp and the atanh series: ln(m) = 2 sum z^(2k+1)/(2k+1), z=(m-1)/(m+1),
// with m in [sqrt(1/2), sqrt(2)), so |z| < 0.1716 and 13 terms reach double precision.
double portableLog(double x) {
    int exponent = 0;
    double m = std::frexp(x, &exponent);
    if (m < 0.70710678118654752440) { m *= 2.0; --exponent; }
    const double z = (m - 1.0) / (m + 1.0), z2 = z * z;
    double sum = 0.0;
    for (int k = 12; k >= 0; --k) sum = sum * z2 + 1.0 / (2.0 * k + 1.0);
    return exponent * 0.69314718055994530942 + 2.0 * z * sum;
}

template <int N>
double polynomial(const double (&coefficients)[N], double x) {
    double value = coefficients[N - 1];
    for (int k = N - 2; k >= 0; --k) value = value * x + coefficients[k];
    return value;
}

}  // namespace

double portableNormalQuantile(double p) {
    if (!(p > 0.0 && p < 1.0)) throw std::invalid_argument("normal quantile needs p in (0, 1)");
    static const double a[8]{3.3871328727963666080e0, 1.3314166789178437745e+2, 1.9715909503065514427e+3,
                             1.3731693765509461125e+4, 4.5921953931549871457e+4, 6.7265770927008700853e+4,
                             3.3430575583588128105e+4, 2.5090809287301226727e+3};
    static const double b[8]{1.0, 4.2313330701600911252e+1, 6.8718700749205790830e+2, 5.3941960214247511077e+3,
                             2.1213794301586595867e+4, 3.9307895800092710610e+4, 2.8729085735721942674e+4,
                             5.2264952788528545610e+3};
    static const double c[8]{1.42343711074968357734e0, 4.63033784615654529590e0, 5.76949722146069140550e0,
                             3.64784832476320460504e0, 1.27045825245236838258e0, 2.41780725177450611770e-1,
                             2.27238449892691845833e-2, 7.74545014278341407640e-4};
    static const double d[8]{1.0, 2.05319162663775882187e0, 1.67638483018380384940e0, 6.89767334985100004550e-1,
                             1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4,
                             1.05075007164441684324e-9};
    static const double e[8]{6.65790464350110377720e0, 5.46378491116411436990e0, 1.78482653991729133580e0,
                             2.96560571828504891230e-1, 2.65321895265761230930e-2, 1.24266094738807843860e-3,
                             2.71155556874348757815e-5, 2.01033439929228813265e-7};
    static const double f[8]{1.0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2,
                             7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7,
                             2.04426310338993978564e-15};
    const double q = p - 0.5;
    if (q >= -0.425 && q <= 0.425) {
        const double r = 0.180625 - q * q;
        return q * polynomial(a, r) / polynomial(b, r);
    }
    double r = std::sqrt(-portableLog(q < 0.0 ? p : 1.0 - p));
    double value;
    if (r <= 5.0) { r -= 1.6; value = polynomial(c, r) / polynomial(d, r); }
    else { r -= 5.0; value = polynomial(e, r) / polynomial(f, r); }
    return q < 0.0 ? -value : value;
}

double portableNormal(std::mt19937_64& generator) {
    const std::uint64_t word = generator() >> 11;  // 53 bits
    return portableNormalQuantile((static_cast<double>(word) + 0.5) * 0x1p-53);
}

}  // namespace dynsyslab::core
