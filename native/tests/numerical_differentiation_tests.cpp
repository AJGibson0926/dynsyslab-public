#include <dynsyslab/core/NumericalDifferentiation.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
void require(bool ok, const std::string& what) { if (!ok) throw std::runtime_error(what); }

double relativeError(const std::vector<double>& est, const std::vector<double>& ref, std::size_t skip) {
    double e = 0, r = 0;
    for (std::size_t i = skip; i + skip < ref.size(); ++i) { e += (est[i] - ref[i]) * (est[i] - ref[i]); r += ref[i] * ref[i]; }
    return std::sqrt(e / r);
}

void everyMethodOnACleanSignal() {
    const double dt = 1e-3; const std::size_t n = 4000;
    std::vector<double> y(n), dy(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double t = double(i) * dt;
        y[i] = std::sin(2 * pi * t) + 0.3 * std::cos(2 * pi * 3 * t);
        dy[i] = 2 * pi * std::cos(2 * pi * t) - 0.3 * 2 * pi * 3 * std::sin(2 * pi * 3 * t);
    }
    for (const auto& m : dynsyslab::core::derivativeMethods()) {
        dynsyslab::core::DerivativeSpec spec; spec.method = m; spec.cutoff = 25;
        if (m == "levant") spec.lipschitz = 5000;  // a bound on |y'''| of the signal (here ~2400)
        const auto out = dynsyslab::core::differentiate(y, dt, spec);
        const bool causal = dynsyslab::core::isCausalDerivativeMethod(m);
        const double err = relativeError(out.derivative, dy, causal ? 400 : 300);
        // Causal methods carry a lag (and levant a sliding-mode ripple); smoothers must be accurate.
        const double tol = (m == "levant" || m == "causal_butterworth") ? 0.25 : causal ? 0.08 : (m == "tv_denoise" || m == "tvr") ? 0.03 : 0.01;
        std::cout << m << " clean error " << err << "\n";
        require(err < tol, m + " is inaccurate on a clean two-tone signal");
    }
}

void smoothersBeatCentralDifferencesUnderNoise() {
    const double dt = 1e-3; const std::size_t n = 4000;
    std::vector<double> y(n), dy(n);
    for (std::size_t i = 0; i < n; ++i) { const double t = double(i) * dt; y[i] = std::sin(2 * pi * t); dy[i] = 2 * pi * std::cos(2 * pi * t); }
    dynsyslab::core::NoiseSpec noise; noise.kind = "white"; noise.level = 0.01;
    const auto noisy = dynsyslab::core::addMeasurementNoise(y, dt, noise);
    double rms = 0; for (std::size_t i = 0; i < n; ++i) rms += (noisy[i] - y[i]) * (noisy[i] - y[i]);
    rms = std::sqrt(rms / double(n));
    require(std::abs(rms / (1 / std::sqrt(2.0)) - 0.01) < 0.002, "white noise level is not relative to the signal RMS");
    require(std::abs(dynsyslab::core::estimateWhiteNoiseStd(noisy) - rms) < 0.2 * rms, "noise estimate is off");
    dynsyslab::core::DerivativeSpec central; central.method = "central";
    const double base = relativeError(dynsyslab::core::differentiate(noisy, dt, central).derivative, dy, 200);
    for (const char* m : {"savitzky_golay", "gaussian", "holoborodko", "weak", "butterworth", "whittaker", "kalman", "spectral", "tikhonov_integral"}) {
        dynsyslab::core::DerivativeSpec spec; spec.method = m; spec.cutoff = 5;
        const double err = relativeError(dynsyslab::core::differentiate(noisy, dt, spec).derivative, dy, 200);
        std::cout << m << " noisy error " << err << " (central " << base << ")\n";
        require(err < 0.1 * base, std::string(m) + " does not suppress white noise");
    }
}

void noiseKinds() {
    std::vector<double> y(20000); for (std::size_t i = 0; i < y.size(); ++i) y[i] = std::sin(2 * pi * 1e-3 * double(i));
    for (const char* k : {"white", "uniform", "ar1", "pink"}) {
        dynsyslab::core::NoiseSpec s; s.kind = k; s.level = 0.05; s.correlationTime = 0.01;
        const auto z = dynsyslab::core::addMeasurementNoise(y, 1e-3, s);
        double e = 0; for (std::size_t i = 0; i < y.size(); ++i) e += (z[i] - y[i]) * (z[i] - y[i]);
        const double rel = std::sqrt(e / double(y.size())) / (1 / std::sqrt(2.0));
        std::cout << k << " relative noise " << rel << "\n";
        require(std::abs(rel - 0.05) < 0.015, std::string(k) + " noise level is wrong");
    }
    dynsyslab::core::NoiseSpec q; q.kind = "quantization"; q.bits = 4;
    const auto z = dynsyslab::core::addMeasurementNoise(y, 1e-3, q);
    for (double v : z) require(std::abs(v / (2.0 / 16) - std::round(v / (2.0 / 16))) < 1e-9, "quantization step is wrong");
    dynsyslab::core::NoiseSpec a; a.kind = "white"; a.level = 0.1; a.seed = 7;
    require(dynsyslab::core::addMeasurementNoise(y, 1e-3, a) == dynsyslab::core::addMeasurementNoise(y, 1e-3, a), "noise is not reproducible");
}
void despikingRemovesImpulses() {
    const double dt = 1e-3; const std::size_t n = 4000;
    std::vector<double> y(n), dy(n);
    for (std::size_t i = 0; i < n; ++i) { const double t = double(i) * dt; y[i] = std::sin(2 * pi * t); dy[i] = 2 * pi * std::cos(2 * pi * t); }
    dynsyslab::core::NoiseSpec spikes; spikes.kind = "impulse"; spikes.probability = 0.005; spikes.level = 1.0;
    const auto noisy = dynsyslab::core::addMeasurementNoise(y, dt, spikes);
    dynsyslab::core::DerivativeSpec plain; plain.method = "kalman"; plain.cutoff = 10;
    auto cleaned = plain; cleaned.despike = 5;
    const double a = relativeError(dynsyslab::core::differentiate(noisy, dt, plain).derivative, dy, 200);
    const double b = relativeError(dynsyslab::core::differentiate(noisy, dt, cleaned).derivative, dy, 200);
    std::cout << "impulses: kalman " << a << ", with despike " << b << "\n";
    require(b < 0.2 * a, "despiking does not remove impulses");
}
}  // namespace

int main() {
    try {
        everyMethodOnACleanSignal();
        smoothersBeatCentralDifferencesUnderNoise();
        noiseKinds();
        despikingRemovesImpulses();
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << "\n"; return 1;
    }
    std::cout << "numerical differentiation tests passed\n";
    return 0;
}
