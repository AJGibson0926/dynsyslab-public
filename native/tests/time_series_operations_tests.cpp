// TimeSeries operations (TimeSeriesOperations.hpp): crop by index and time, resampling with its
// boundary modes, the Fourier round trip, and the window-history transformation built from crops
// (too short, exactly one window, several windows, histories of histories), checked against the
// fast FFT history projection it defines.
#include <dynsyslab/core/HistoryObservation.hpp>
#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
namespace c = dynsyslab::core;
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
template <class F> bool rejects(F f) {
    try { (void)f(); } catch (const std::invalid_argument&) { return true; }
    return false;
}
c::TimeSeries series(std::size_t n, double dt = 0.1) {  // two components: t^2 and sin(3t)
    std::vector<double> t(n), v;
    for (std::size_t i = 0; i < n; ++i) {
        t[i] = double(i) * dt;
        v.push_back(t[i] * t[i]);
        v.push_back(std::sin(3 * t[i]));
    }
    return c::TimeSeries(std::move(t), std::move(v), 2);
}
}  // namespace

int main() {
    try {
        const auto s = series(20);
        // Crop by index and by time.
        const auto cropped = s.crop(5, 4);
        require(cropped.sampleCount() == 4 && cropped.time(0) == s.time(5) && cropped.value(3, 1) == s.value(8, 1), "crop");
        require(rejects([&] { return s.crop(18, 3); }) && rejects([&] { return s.crop(0, 0); }), "crop bounds");
        const auto window = s.cropTime(0.45, 0.85);
        require(window.sampleCount() == 4 && std::abs(window.time(0) - 0.5) < 1e-12, "crop time");

        // Resample: linear is exact on lines, cubic Hermite exact on the interior of a quadratic.
        const auto line = c::TimeSeries({0, 1, 2}, {0, 2, 4}, 1).resample({0.5, 1.75});
        require(std::abs(line.value(0, 0) - 1) < 1e-14 && std::abs(line.value(1, 0) - 3.5) < 1e-14, "linear resample");
        const auto cubic = s.resample({0.55, 1.23}, {c::Interpolation::cubic, c::Extrapolation::error});
        require(std::abs(cubic.value(0, 0) - 0.55 * 0.55) < 1e-12 && std::abs(cubic.value(1, 0) - 1.23 * 1.23) < 1e-12,
                "cubic resample of a quadratic");
        require(rejects([&] { return s.resample({-1.0}); }), "outside the span is an error by default");
        const auto held = s.resample({-1.0, 5.0}, {c::Interpolation::linear, c::Extrapolation::hold});
        require(held.value(0, 0) == s.value(0, 0) && held.value(1, 0) == s.value(19, 0), "hold the ends");
        require(std::isnan(s.resample({9.0}, {c::Interpolation::linear, c::Extrapolation::nan}).value(0, 1)), "NaN outside");
        const auto onto = s.resampleOnto(cropped);
        require(onto.sampleCount() == 4 && std::abs(onto.value(2, 1) - cropped.value(2, 1)) < 1e-14, "resample onto another grid");

        // Fourier: to the frequency domain and back.
        for (const bool oneSided : {true, false}) {
            c::FourierOptions options;
            options.oneSided = oneSided;
            const auto back = c::TimeSeries::fromSpectrum(*s.spectrum(options));
            for (std::size_t i = 0; i < 20; ++i)
                for (std::size_t k = 0; k < 2; ++k)
                    require(std::abs(back.value(i, k) - s.value(i, k)) < 1e-10, "Fourier round trip");
        }

        // Window history, fixed by its window: too short, exactly one window, several windows.
        const c::WindowHistory h(6, 2);
        require(h.windowCount(5) == 0 && rejects([&] { return h.apply(s.crop(0, 5)); }), "shorter than the window");
        const auto one = h.apply(s.crop(0, 6));
        require(one.sampleCount() == 1 && one.componentCount() == 12 && one.time(0) == s.time(5), "exactly one window");
        const auto many = h.apply(s);
        require(many.sampleCount() == 8 && many.time(1) == s.time(7), "(n - window) / shift + 1 windows, labelled at the end");
        require(many.value(1, 0) == s.value(2, 0) && many.value(1, 11) == s.value(7, 1), "oldest sample first");
        // A history is a time series: a history of histories is just the transform applied again.
        const auto twice = s.history(3).history(2);
        require(twice.componentCount() == 2 * 3 * 2 && twice.sampleCount() == 17, "history of histories");

        // The FFT history projection is the window history followed by (U^T h)/sigma.
        c::PooledHankelBasis basis;
        basis.delays = 4;
        basis.rank = 2;
        basis.leftSingularVectors = {0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5};  // delay x rank
        basis.singularValues = {2.0, 1.0};
        std::vector<double> scalar;
        for (std::size_t i = 0; i < 20; ++i) scalar.push_back(s.value(i, 1));
        const auto fast = c::projectHistory(scalar, basis);
        const auto windows = c::TimeSeries({s.timeData(), s.timeData() + 20}, scalar, 1).history(4);
        require(fast.windows == windows.sampleCount(), "same window count");
        for (std::size_t w = 0; w < windows.sampleCount(); ++w)
            for (std::size_t j = 0; j < 2; ++j) {
                double q = 0;
                for (std::size_t d = 0; d < 4; ++d) q += basis.leftSingularVectors[d * 2 + j] * windows.value(w, d);
                require(std::abs(fast.coordinates[w * 2 + j] - q / basis.singularValues[j]) < 1e-12, "projection of the cropped windows");
            }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "time series operations tests passed\n";
    return 0;
}
