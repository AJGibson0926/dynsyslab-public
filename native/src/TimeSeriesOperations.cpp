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

TimeSeries resample(const TimeSeries& s, const std::vector<double>& times, ResampleOptions options) {
    const std::size_t n = s.sampleCount(), d = s.componentCount();
    if (n < 2 && options.method != Interpolation::linear) throw std::invalid_argument("cubic resampling needs two samples");
    const auto* t = s.timeData();
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> out(times.size() * d);
    const auto slope = [&](std::size_t i, std::size_t c) {  // Catmull-Rom slope at sample i
        const auto a = i ? i - 1 : i, b = i + 1 < n ? i + 1 : i;
        return b == a ? 0.0 : (s.value(b, c) - s.value(a, c)) / (t[b] - t[a]);
    };
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
        auto j = std::size_t(std::upper_bound(t, t + n, x) - t);
        j = std::clamp<std::size_t>(j, 1, n - 1);  // interval [j - 1, j]
        const double h = t[j] - t[j - 1], u = h > 0 ? (x - t[j - 1]) / h : 0.0;
        for (std::size_t c = 0; c < d; ++c) {
            const double y0 = s.value(j - 1, c), y1 = s.value(j, c);
            if (options.method == Interpolation::linear || n < 3) { out[k * d + c] = y0 + u * (y1 - y0); continue; }
            const double h00 = (1 + 2 * u) * (1 - u) * (1 - u), h10 = u * (1 - u) * (1 - u);
            const double h01 = u * u * (3 - 2 * u), h11 = u * u * (u - 1);
            out[k * d + c] = h00 * y0 + h10 * h * slope(j - 1, c) + h01 * y1 + h11 * h * slope(j, c);
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
