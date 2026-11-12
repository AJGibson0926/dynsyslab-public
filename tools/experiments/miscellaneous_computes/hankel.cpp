#include "hankel.hpp"
#include <unsupported/Eigen/FFT>
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

extern "C" int mc_hankel_gram_v1(const double* x, std::size_t samples,
    std::size_t channels, std::size_t taps, double* gram) noexcept {
    const auto limit = std::numeric_limits<std::size_t>::max();
    if (!x || !gram || !channels || !taps || taps > samples ||
        samples > limit / channels || samples > (1u << 29)) return 1;
    const auto dim = taps * channels;
    if (dim > limit / dim || !std::all_of(x, x + samples * channels,
        [](double value) { return std::isfinite(value); })) return 1;
    try {
        using Complex = std::complex<double>;
        const auto windows = samples - taps + 1;
        std::size_t fftSize = 2;
        while (fftSize < samples) fftSize *= 2;
        Eigen::FFT<double> fft;
        std::vector<std::vector<Complex>> full(channels), prefix(channels);
        std::vector<double> signal(fftSize);
        for (std::size_t c = 0; c < channels; ++c) {
            std::fill(signal.begin(), signal.end(), 0.);
            for (std::size_t t = 0; t < samples; ++t) signal[t] = x[t*channels+c];
            fft.fwd(full[c], signal);
            std::fill(signal.begin() + windows, signal.end(), 0.);
            fft.fwd(prefix[c], signal);
        }
        std::vector<Complex> product(fftSize);
        std::vector<std::vector<double>> edge(channels*channels);
        for (std::size_t a = 0; a < channels; ++a) {
            for (std::size_t b = 0; b < channels; ++b) {
                for (std::size_t k = 0; k < fftSize; ++k)
                    product[k] = std::conj(prefix[a][k]) * full[b][k];
                fft.inv(signal, product);
                edge[a*channels+b].assign(signal.begin(), signal.begin()+taps);
            }
        }
        for (std::size_t a = 0; a < channels; ++a) {
            for (std::size_t b = 0; b < channels; ++b) {
                auto previous = edge[a*channels+b];
                std::vector<double> current(taps);
                for (std::size_t i = 0; i < taps; ++i) {
                    if (i) {
                        current[0] = edge[b*channels+a][i];
                        const auto enter = x[(windows+i-1)*channels+a];
                        const auto leave = x[(i-1)*channels+a];
                        for (std::size_t j = 1; j < taps; ++j)
                            current[j] = previous[j-1] + enter*x[(windows+j-1)*channels+b] - leave*x[(j-1)*channels+b];
                        previous.swap(current);
                    }
                    for (std::size_t j = 0; j < taps; ++j)
                        gram[(i*channels+a)*dim+j*channels+b] += previous[j];
                }
            }
        }
    } catch (...) { return 2; }
    return 0;
}
