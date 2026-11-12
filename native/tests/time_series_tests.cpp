#include <dynsyslab/core/TimeSeries.hpp>
#include <dynsyslab/core/TimeSeriesArchive.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void transformAndCache() {
    constexpr std::size_t samples = 1000;
    constexpr double dt = 0.001;
    std::vector<double> time(samples), values(samples * 2);
    for (std::size_t sample = 0; sample < samples; ++sample) {
        time[sample] = sample * dt;
        values[sample * 2] = 3.0 + 2.0 * std::sin(2.0 * pi * 50.0 * time[sample]);
        values[sample * 2 + 1] = 0.5 * std::cos(2.0 * pi * 120.0 * time[sample]);
    }
    dynsyslab::core::TimeSeries series(time, values, 2);
    require(series.uniformTimeSteps() && series.sampleCount() == samples,
            "uniform TimeSeries metadata");
    dynsyslab::core::FourierOptions options;
    options.removeMean = true;
    const auto sidecar = series.fourier();
    require(!sidecar.cached(options), "Fourier analysis ran eagerly");
    const auto first = sidecar.transform(options);
    require(sidecar.cached(options) && sidecar.transform(options) == first,
            "Fourier transform cache was not reused");
    require(first->binCount == samples / 2 + 1 && first->componentCount == 2,
            "one-sided Fourier dimensions");
    require(std::abs(first->dominantFrequency(0) - 50.0) < 1.0e-10
            && std::abs(first->amplitude(50, 0) - 2.0) < 1.0e-10,
            "first component frequency or amplitude");
    require(std::abs(first->dominantFrequency(1) - 120.0) < 1.0e-10
            && std::abs(first->amplitude(120, 1) - 0.5) < 1.0e-10,
            "second component frequency or amplitude");
    options.window = dynsyslab::core::FourierWindow::hann;
    const auto windowed = sidecar.transform(options);
    require(windowed != first && sidecar.cached(options),
            "Fourier option change did not replace the active cache");
    require(std::abs(windowed->dominantFrequency(0) - 50.0) < 1.0e-10,
            "windowed Fourier frequency");
    sidecar.clear(); require(!sidecar.cached(options), "Fourier cache clear");
}

void viewAndArchiveAdapter() {
    std::vector<double> time{0.0, 0.1, 0.2};
    std::vector<double> values{1.0, 0.0, 1.1, 0.2, 1.2, 0.1};
    const auto view = dynsyslab::core::TimeSeries::view(time, values.data(), 3, 2);
    require(view.valuesData() == values.data(), "TimeSeries view copied its source data");
    dynsyslab::core::NpzRecord archive;
    dynsyslab::core::storeTimeSeries(archive, view);
    require(archive.contains("time") && archive.contains("values")
            && !archive.contains("frequencies") && !archive.contains("fourier_amplitudes"),
            "archive adapter persisted derived Fourier data");
    const auto loaded = dynsyslab::core::loadTimeSeries(archive);
    require(loaded.sampleCount() == 3 && loaded.componentCount() == 2
            && loaded.value(2, 0) == 1.2, "archive TimeSeries round trip");
}

void rejectNonuniformTransform() {
    dynsyslab::core::TimeSeries series(
        {0.0, 0.1, 0.21}, {1.0, 2.0, 3.0}, 1);
    require(!series.uniformTimeSteps(), "nonuniform grid classification");
    bool rejected = false;
    try { static_cast<void>(series.fourier().transform()); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "nonuniform Fourier input was accepted");
}

}  // namespace

int main() {
    try {
        transformAndCache();
        viewAndArchiveAdapter();
        rejectNonuniformTransform();
        std::cout << "dynsyslab_time_series_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_time_series_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
