#include <dynsyslab/core/PathSimplification.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {

std::vector<std::size_t> boundedPathVertices(const std::vector<double>& v, std::size_t n, std::size_t width,
                                             double fraction) {
    if (!n || !width || v.size() != n * width || !(fraction >= 0))
        throw std::invalid_argument("path simplification needs a rows x width matrix and a nonnegative tolerance");
    if (n <= 2) { std::vector<std::size_t> all(n); for (std::size_t k = 0; k < n; ++k) all[k] = k; return all; }
    std::vector<std::size_t> anchors{0, n - 1};
    std::vector<double> tolerance(width);
    for (std::size_t c = 0; c < width; ++c) {
        std::size_t lo = 0, hi = 0;
        for (std::size_t k = 1; k < n; ++k) {
            if (v[k * width + c] < v[lo * width + c]) lo = k;
            if (v[k * width + c] > v[hi * width + c]) hi = k;
        }
        anchors.push_back(lo); anchors.push_back(hi);
        tolerance[c] = fraction * (v[hi * width + c] - v[lo * width + c]);
    }
    std::sort(anchors.begin(), anchors.end());
    anchors.erase(std::unique(anchors.begin(), anchors.end()), anchors.end());
    std::vector<std::size_t> out{0};
    std::vector<std::pair<std::size_t, std::size_t>> stack;
    for (std::size_t p = 1; p < anchors.size(); ++p)
        for (std::size_t start = anchors[p - 1]; start < anchors[p];) {
            const auto stop = std::min(anchors[p], start + 8192);
            stack.emplace_back(start, stop);
            while (!stack.empty()) {  // depth-first, so indices come out in order
                const auto [a, b] = stack.back();
                stack.pop_back();
                bool fits = true;
                for (std::size_t i = a + 1; i < b && fits; ++i) {
                    const double f = double(i - a) / double(b - a);
                    for (std::size_t c = 0; c < width; ++c) {
                        const double x0 = v[a * width + c], x1 = v[b * width + c];
                        if (std::abs(v[i * width + c] - (x0 + f * (x1 - x0))) > tolerance[c]) { fits = false; break; }
                    }
                }
                if (fits) { out.push_back(b); continue; }
                const auto middle = a + (b - a) / 2;
                stack.emplace_back(middle, b);
                stack.emplace_back(a, middle);
            }
            start = stop;
        }
    return out;
}

std::vector<std::size_t> simplifiedPathIndices(const std::vector<double>& v, std::size_t n, std::size_t width,
                                               std::size_t maximum, double fraction) {
    auto even = [&] {
        std::vector<std::size_t> out;
        const std::size_t count = std::min(n, std::max<std::size_t>(maximum, 2));
        for (std::size_t k = 0; k < count; ++k) out.push_back(count > 1 ? k * (n - 1) / (count - 1) : 0);
        return out;
    };
    if (n <= maximum) { std::vector<std::size_t> all(n); for (std::size_t k = 0; k < n; ++k) all[k] = k; return all; }
    for (const double x : v) if (!std::isfinite(x)) return even();
    for (double f = fraction; f < 1; f *= 2) {
        auto kept = boundedPathVertices(v, n, width, f);
        if (kept.size() <= maximum) return kept;
    }
    return even();
}

} // namespace dynsyslab::core
