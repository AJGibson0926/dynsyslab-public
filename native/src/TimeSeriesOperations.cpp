#include <dynsyslab/core/DiscreteFourier.hpp>
#include <dynsyslab/core/TimeSeriesOperations.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace dynsyslab::core {

TimeSeries cropSamples(const TimeSeries& s, std::size_t first, std::size_t count) {
    if (!count || first > s.sampleCount() || count > s.sampleCount() - first)
        throw std::invalid_argument("crop [" + std::to_string(first) + ", +" + std::to_string(count) + ") is outside the "
                                    + std::to_string(s.sampleCount()) + " samples");
    const auto d = s.componentCount();
    return TimeSeries({s.timeData() + first, s.timeData() + first + count},
                      {s.valuesData() + first * d, s.valuesData() + (first + count) * d}, d);
}

TimeSeries cropTime(const TimeSeries& s, double start, double stop) {
    if (!(stop >= start)) throw std::invalid_argument("crop time needs start <= stop");
    const auto* t = s.timeData();
    const auto first = std::size_t(std::lower_bound(t, t + s.sampleCount(), start) - t);
    const auto last = std::size_t(std::upper_bound(t, t + s.sampleCount(), stop) - t);
    if (last <= first) throw std::invalid_argument("no samples between the crop times");
    return cropSamples(s, first, last - first);
}

namespace {

// PCHIP slopes (Fritsch & Carlson 1980; Fritsch & Butland 1984 weighted harmonic mean): zero at
// extrema and flat segments, so the interpolant is monotone wherever the data are.
std::vector<double> pchipSlopes(const TimeSeries& s) {
    const std::size_t n = s.sampleCount(), d = s.componentCount();
    const auto* t = s.timeData();
    std::vector<double> slopes(n * d, 0.0);
    if (n < 2) return slopes;
    const auto secant = [&](std::size_t i, std::size_t c) {
        const double h = t[i + 1] - t[i];
        return h > 0 ? (s.value(i + 1, c) - s.value(i, c)) / h : 0.0;
    };
    for (std::size_t c = 0; c < d; ++c) {
        if (n == 2) { slopes[c] = slopes[d + c] = secant(0, c); continue; }
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const double h0 = t[i] - t[i - 1], h1 = t[i + 1] - t[i];
            const double m0 = secant(i - 1, c), m1 = secant(i, c);
            if (m0 * m1 <= 0.0) continue;
            const double w1 = 2.0 * h1 + h0, w2 = h1 + 2.0 * h0;
            slopes[i * d + c] = (w1 + w2) / (w1 / m0 + w2 / m1);
        }
        const auto end = [&](double h0, double h1, double m0, double m1) {
            double slope = ((2.0 * h0 + h1) * m0 - h0 * m1) / (h0 + h1);
            if (slope * m0 <= 0.0) slope = 0.0;
            else if (m0 * m1 <= 0.0 && std::abs(slope) > std::abs(3.0 * m0)) slope = 3.0 * m0;
            return slope;
        };
        slopes[c] = end(t[1] - t[0], t[2] - t[1], secant(0, c), secant(1, c));
        slopes[(n - 1) * d + c] = end(t[n - 1] - t[n - 2], t[n - 2] - t[n - 3],
                                      secant(n - 2, c), secant(n - 3, c));
    }
    return slopes;
}

}  // namespace

Interpolation interpolationFromName(const std::string& name) {
    if (name == "linear") return Interpolation::linear;
    if (name == "cubic" || name == "hermite") return Interpolation::cubic;
    if (name == "pchip") return Interpolation::pchip;
    if (name == "nearest") return Interpolation::nearest;
    throw std::invalid_argument("interpolation method must be linear, cubic, pchip or nearest");
}

const char* interpolationName(Interpolation method) {
    switch (method) {
    case Interpolation::linear: return "linear";
    case Interpolation::cubic: return "cubic";
    case Interpolation::pchip: return "pchip";
    case Interpolation::nearest: return "nearest";
    }
    return "linear";
}

TimeSeries resample(const TimeSeries& s, const std::vector<double>& times, ResampleOptions options) {
    const std::size_t n = s.sampleCount(), d = s.componentCount();
    if (!n) throw std::invalid_argument("resampling needs at least one sample");
    if (n < 2 && options.method != Interpolation::linear && options.method != Interpolation::nearest)
        throw std::invalid_argument("cubic resampling needs two samples");
    const auto* t = s.timeData();
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> out(times.size() * d);
    const auto slope = [&](std::size_t i, std::size_t c) {  // Catmull-Rom slope at sample i
        const auto a = i ? i - 1 : i, b = i + 1 < n ? i + 1 : i;
        return b == a ? 0.0 : (s.value(b, c) - s.value(a, c)) / (t[b] - t[a]);
    };
    const auto monotone = options.method == Interpolation::pchip ? pchipSlopes(s) : std::vector<double>{};
    for (std::size_t k = 0; k < times.size(); ++k) {
        const double x = times[k];
        if (!std::isfinite(x)) throw std::invalid_argument("resample times must be finite");
        const bool below = x < t[0], above = x > t[n - 1];
        if (below || above) {
            if (options.outside == Extrapolation::error) throw std::invalid_argument("resample time outside the series span");
            for (std::size_t c = 0; c < d; ++c)
                out[k * d + c] = options.outside == Extrapolation::nan ? nan : s.value(below ? 0 : n - 1, c);
            continue;
        }
        if (n == 1) { for (std::size_t c = 0; c < d; ++c) out[k * d + c] = s.value(0, c); continue; }
        auto j = std::size_t(std::upper_bound(t, t + n, x) - t);
        j = std::clamp<std::size_t>(j, 1, n - 1);  // interval [j - 1, j]
        const double h = t[j] - t[j - 1], u = h > 0 ? (x - t[j - 1]) / h : 0.0;
        if (options.method == Interpolation::nearest) {
            const std::size_t pick = x - t[j - 1] <= t[j] - x ? j - 1 : j;
            for (std::size_t c = 0; c < d; ++c) out[k * d + c] = s.value(pick, c);
            continue;
        }
        for (std::size_t c = 0; c < d; ++c) {
            const double y0 = s.value(j - 1, c), y1 = s.value(j, c);
            const bool linear = options.method == Interpolation::linear
                || (options.method == Interpolation::cubic && n < 3);
            if (linear) { out[k * d + c] = y0 + u * (y1 - y0); continue; }
            const double h00 = (1 + 2 * u) * (1 - u) * (1 - u), h10 = u * (1 - u) * (1 - u);
            const double h01 = u * u * (3 - 2 * u), h11 = u * u * (u - 1);
            const double s0 = monotone.empty() ? slope(j - 1, c) : monotone[(j - 1) * d + c];
            const double s1 = monotone.empty() ? slope(j, c) : monotone[j * d + c];
            out[k * d + c] = h00 * y0 + h10 * h * s0 + h01 * y1 + h11 * h * s1;
        }
    }
    return TimeSeries(times, std::move(out), d);
}

TimeSeries resampleOnto(const TimeSeries& s, const TimeSeries& reference, ResampleOptions options) {
    return resample(s, {reference.timeData(), reference.timeData() + reference.sampleCount()}, options);
}

TimeSeries timeSeriesFromSpectrum(const FourierTransform& f, double startTime) {
    if (f.window != FourierWindow::rectangular) throw std::invalid_argument("only a rectangular-window transform inverts exactly");
    const std::size_t n = f.sourceSampleCount, d = f.componentCount;
    if (n < 2 || !(f.deltaTime > 0)) throw std::invalid_argument("the transform has no source clock");
    std::vector<double> values(n * d), time(n);
    for (std::size_t i = 0; i < n; ++i) time[i] = startTime + double(i) * f.deltaTime;
    for (std::size_t c = 0; c < d; ++c) {
        std::vector<std::complex<double>> full(n);  // X_k = N c_k; one-sided spectra are Hermitian
        for (std::size_t k = 0; k < n; ++k) {
            if (!f.oneSided) full[k] = f.coefficient(k, c);
            else full[k] = k < f.binCount ? f.coefficient(k, c) : std::conj(f.coefficient(n - k, c));
            full[k] *= double(n);
        }
        const auto x = discreteFourierTransform(full, true);
        for (std::size_t i = 0; i < n; ++i) values[i * d + c] = x[i].real();
    }
    return TimeSeries(std::move(time), std::move(values), d);
}

WindowHistory::WindowHistory(std::size_t window, std::size_t shift, Label label) : window_(window), shift_(shift), label_(label) {
    if (!window || !shift) throw std::invalid_argument("a window history needs a positive window and shift");
}

std::size_t WindowHistory::windowCount(std::size_t samples) const noexcept {
    return samples < window_ ? 0 : (samples - window_) / shift_ + 1;
}

TimeSeries WindowHistory::apply(const TimeSeries& series) const {
    const auto windows = windowCount(series.sampleCount());
    if (!windows)
        throw std::invalid_argument("the series (" + std::to_string(series.sampleCount()) + " samples) is shorter than the "
                                    + std::to_string(window_) + "-sample window");
    const auto width = window_ * series.componentCount();
    std::vector<double> time, values;
    time.reserve(windows);
    values.reserve(windows * width);
    for (std::size_t k = 0; k < windows; ++k) {
        const auto piece = cropSamples(series, k * shift_, window_);  // the window is a crop
        time.push_back(piece.time(label_ == Label::end ? window_ - 1 : 0));
        values.insert(values.end(), piece.valuesData(), piece.valuesData() + width);
    }
    return TimeSeries(std::move(time), std::move(values), width);
}

// Member shortcuts (TimeSeries.hpp).
TimeSeries TimeSeries::crop(std::size_t first, std::size_t count) const { return cropSamples(*this, first, count); }
TimeSeries TimeSeries::cropTime(double start, double stop) const { return core::cropTime(*this, start, stop); }
TimeSeries TimeSeries::resample(const std::vector<double>& times, ResampleOptions options) const {
    return core::resample(*this, times, options);
}
TimeSeries TimeSeries::resampleOnto(const TimeSeries& reference, ResampleOptions options) const {
    return core::resampleOnto(*this, reference, options);
}
std::shared_ptr<const FourierTransform> TimeSeries::spectrum(FourierOptions options) const { return fourier().transform(options); }
TimeSeries TimeSeries::fromSpectrum(const FourierTransform& transform, double startTime) {
    return timeSeriesFromSpectrum(transform, startTime);
}
TimeSeries TimeSeries::history(std::size_t window, std::size_t shift) const { return WindowHistory(window, shift).apply(*this); }

}  // namespace dynsyslab::core
