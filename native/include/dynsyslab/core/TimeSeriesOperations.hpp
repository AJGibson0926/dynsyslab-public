#pragma once
// Operations on TimeSeries that return new TimeSeries (TimeSeries has member shortcuts to each):
//
//   crop(first, count)          samples [first, first + count): the one windowing primitive
//   cropTime(start, stop)       samples with start <= t <= stop
//   resample(times)             values at new times (linear, cubic Hermite, monotone PCHIP or
//   resampleOnto(reference)     nearest; outside the source span: error, hold the end value, or NaN)
//   interpolationFromName(name) "linear" | "cubic" (alias "hermite") | "pchip" | "nearest"
//   spectrum(options)           the Fourier transform (Fourier.cpp, cached on the series)
//   fromSpectrum(transform)     back to the time domain (rectangular window, mean kept)
//   WindowHistory               the window-history coordinate transformation, fixed by its window
//                               length and shift: output sample k is crop(k * shift, window)
//                               flattened oldest first, labelled at the window's last time.
//                               A series shorter than the window has no windows; exactly the window
//                               gives one; longer gives (n - window) / shift + 1. The output is a
//                               TimeSeries again, so histories of histories work the same way.
//
// Cubic resampling uses Hermite cubics with Catmull-Rom (centred-difference) slopes (Fritsch and
// Carlson 1980, SIAM J. Numer. Anal. 17:238, for the Hermite form). PCHIP is the monotone
// variant of Fritsch and Carlson 1980 with the weighted-harmonic-mean interior slopes of Fritsch
// and Butland 1984 (SIAM J. Sci. Stat. Comput. 5:300) and shape-preserving one-sided end slopes
// (Moler, Numerical Computing with MATLAB, 2004, Sec. 3.4); it never overshoots the data. Delay embedding: Takens 1981,
// Lecture Notes in Mathematics 898:366.
#include <dynsyslab/core/TimeSeries.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

[[nodiscard]] TimeSeries cropSamples(const TimeSeries& source, std::size_t first, std::size_t count);
[[nodiscard]] TimeSeries cropTime(const TimeSeries& source, double start, double stop);
[[nodiscard]] TimeSeries resample(const TimeSeries& source, const std::vector<double>& times, ResampleOptions options = {});
[[nodiscard]] TimeSeries resampleOnto(const TimeSeries& source, const TimeSeries& reference, ResampleOptions options = {});
[[nodiscard]] Interpolation interpolationFromName(const std::string& name);
[[nodiscard]] const char* interpolationName(Interpolation method);
// The time series a transform came from: exact for a rectangular window (mean removed or not; a
// removed mean is not restored). Starts at `startTime`.
[[nodiscard]] TimeSeries timeSeriesFromSpectrum(const FourierTransform& transform, double startTime = 0);

class WindowHistory {
public:
    enum class Label { end, start };
    explicit WindowHistory(std::size_t window, std::size_t shift = 1, Label label = Label::end);
    [[nodiscard]] std::size_t window() const noexcept { return window_; }
    [[nodiscard]] std::size_t shift() const noexcept { return shift_; }
    // Number of windows in a series of `samples` samples (0 when it is shorter than the window).
    [[nodiscard]] std::size_t windowCount(std::size_t samples) const noexcept;
    // One sample per window; window * components coordinates, oldest sample first. Throws when
    // the series is shorter than the window.
    [[nodiscard]] TimeSeries apply(const TimeSeries& series) const;

private:
    std::size_t window_{}, shift_{};
    Label label_{Label::end};
};

}  // namespace dynsyslab::core
