#pragma once
#include <cstddef>
#include <vector>

namespace dynsyslab::core {

// Tap lags and both strides count ticks on the decimated observation grid.
struct HankelChannel {
    std::size_t input{};
    std::size_t lag{};
    double scale{1.0};
};
struct HankelWindow {
    std::size_t delays{};
    std::size_t sampleStride{1};
    std::size_t delayStride{1};
    std::size_t snapshotStride{1};
    std::size_t columns{}; // zero uses all complete windows
};
// Adapt the original H-DMD convention: keep the requested number of pairs
// and shorten the window if the source has fewer complete delays.
[[nodiscard]] HankelWindow fixedPairHankelWindow(std::size_t samples,std::size_t maxLag,
    std::size_t requestedDelays,std::size_t pairs,std::size_t sampleStride,std::size_t snapshotStride);

struct HankelMatrix {
    HankelWindow window;
    std::vector<HankelChannel> channels;
    std::size_t sourceSamples{};
    std::size_t rows{};
    std::size_t columns{};
    std::vector<double> x, y; // feature by snapshot, row-major; y is optional
    std::vector<std::size_t> starts, ends; // original source ticks at augmented-window start/end
};

// Each block row contains every augmented channel, oldest block first.
// X[d*C+c,j] = values[(starts[j]+s*(d*delayStride-lag[c])),input[c]] / scale[c].
// Y evaluates those same entries one analysis tick s later.
[[nodiscard]] HankelMatrix buildHankelMatrix(const std::vector<double>& values,
    std::size_t samples, std::size_t inputs, const std::vector<HankelChannel>& channels,
    const HankelWindow& window, bool shiftedPair);

} // namespace dynsyslab::core
