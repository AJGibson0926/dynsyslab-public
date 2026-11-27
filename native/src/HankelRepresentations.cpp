#include <dynsyslab/core/Execution.hpp>
#include <dynsyslab/core/HankelRepresentations.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <thread>

namespace dynsyslab::core {
namespace {

// G = L L^H (Hermitian positive definite).
std::vector<std::complex<double>> cholesky(const std::vector<std::complex<double>>& gram, std::size_t n) {
    std::vector<std::complex<double>> L(n * n);
    for (std::size_t j = 0; j < n; ++j) {
        std::complex<double> diagonal = gram[j * n + j];
        for (std::size_t k = 0; k < j; ++k) diagonal -= L[j * n + k] * std::conj(L[j * n + k]);
        if (!(diagonal.real() > 0.0)) throw std::runtime_error("window-map Gram matrix is not positive definite");
        const double d = std::sqrt(diagonal.real());
        L[j * n + j] = d;
        for (std::size_t i = j + 1; i < n; ++i) {
            std::complex<double> value = gram[i * n + j];
            for (std::size_t k = 0; k < j; ++k) value -= L[i * n + k] * std::conj(L[j * n + k]);
            L[i * n + j] = value / d;
        }
    }
    return L;
}

std::vector<std::complex<double>> solve(const std::vector<std::complex<double>>& L, std::vector<std::complex<double>> b,
                                        std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) {  // L y = b
        for (std::size_t k = 0; k < i; ++k) b[i] -= L[i * n + k] * b[k];
        b[i] /= L[i * n + i];
    }
    for (std::size_t i = n; i-- > 0;) {  // L^H x = y
        for (std::size_t k = i + 1; k < n; ++k) b[i] -= std::conj(L[k * n + i]) * b[k];
        b[i] /= std::conj(L[i * n + i]);
    }
    return b;
}

}  // namespace

std::shared_ptr<const FunctionExpansion> windowMapFromMap(const HankelWindowLayout& layout,
                                                          const std::vector<std::complex<double>>& map, std::size_t modes) {
    const std::size_t rows = layout.rows(), channels = layout.channels.size();
    if (!modes || map.size() != modes * rows) throw std::invalid_argument("window map must be modes x window rows");
    // psi_j(h) = sum_r T_jr h_r / scale_r: one coordinate term per window row.
    std::vector<DictionaryTerm> terms(rows);
    for (std::size_t r = 0; r < rows; ++r) { terms[r].kind = "coordinate"; terms[r].index = r; }
    std::vector<std::complex<double>> c(modes * rows);
    for (std::size_t j = 0; j < modes; ++j)
        for (std::size_t r = 0; r < rows; ++r) c[j * rows + r] = map[j * rows + r] / layout.channels[r % channels].scale;
    return std::make_shared<const FunctionExpansion>(std::make_shared<const ExpansionDictionary>(rows, std::move(terms)),
                                                     std::move(c), FunctionDomain::seriesWindow(layout));
}

std::shared_ptr<const FunctionExpansion> windowMapFromModes(const HankelWindowLayout& layout,
    const std::vector<std::complex<double>>& modes, const std::vector<std::complex<double>>& gram, std::size_t rank) {
    const std::size_t rows = layout.rows();
    if (!rank || modes.size() != rows * rank || gram.size() != rank * rank)
        throw std::invalid_argument("window map needs rows x rank modes and a rank x rank Gram matrix");
    const auto L = cholesky(gram, rank);
    // Row s of T = G^-1 Phi^H is conj(Phi g_s)^T with g_s = G^-1 e_s (G Hermitian). Rows are
    // independent, so large maps are split across threads.
    std::vector<std::complex<double>> T(rank * rows);
    auto fill = [&](std::size_t begin, std::size_t end) {
        for (std::size_t s = begin; s < end; ++s) {
            std::vector<std::complex<double>> unit(rank);
            unit[s] = 1.0;
            const auto g = solve(L, std::move(unit), rank);
            for (std::size_t r = 0; r < rows; ++r) {
                std::complex<double> sum{};
                for (std::size_t a = 0; a < rank; ++a) sum += g[a] * modes[r * rank + a];
                T[s * rows + r] = std::conj(sum);
            }
        }
    };
    const std::size_t workers = rank * rows * rank < (std::size_t(1) << 24) ? 1
        : dynsyslab::core::workerThreads(rank);  // honours --max-threads
    if (workers == 1) fill(0, rank);
    else {
        std::vector<std::thread> threads;
        std::vector<std::exception_ptr> errors(workers);
        for (std::size_t w = 0; w < workers; ++w)
            threads.emplace_back([&, w] {
                try { fill(rank * w / workers, rank * (w + 1) / workers); }
                catch (...) { errors[w] = std::current_exception(); }
            });
        for (auto& thread : threads) thread.join();
        for (const auto& error : errors) if (error) std::rethrow_exception(error);
    }
    return windowMapFromMap(layout, T, rank);
}

std::vector<std::complex<double>> windowMapMatrix(const FunctionExpansion& windowMap) {
    const auto& window = windowMap.coordinates().window;
    if (!window) throw std::invalid_argument("not a window map: its domain is not a series window");
    const std::size_t rows = window->rows(), channels = window->channels.size(), modes = windowMap.outputCount();
    const auto& c = windowMap.coefficients();
    if (c.size() != modes * rows) throw std::invalid_argument("not a window map: one coordinate term per window row is expected");
    std::vector<std::complex<double>> T(c.size());
    for (std::size_t j = 0; j < modes; ++j)
        for (std::size_t r = 0; r < rows; ++r) T[j * rows + r] = c[j * rows + r] * window->channels[r % channels].scale;
    return T;
}

json::Value windowMapDescription(const HankelWindowLayout& layout) {
    json::Value::Array channels;
    for (std::size_t c = 0; c < layout.channels.size(); ++c)
        channels.emplace_back(json::Value::Object{{"column", double(layout.channels[c].input)},
            {"lag", double(layout.channels[c].lag)}, {"scale", layout.channels[c].scale},
            {"observable", c < layout.channelObservables.size() ? layout.channelObservables[c] : ""}});
    return json::Value::Object{{"kind", "expansion"}, {"name", "time-delay"}, {"family", "linear-window"},
        {"form", "psi = W^-1 U^T h"}, {"causal", layout.alignEnd}, {"label", layout.alignEnd ? "newest sample" : "oldest sample"},
        {"delays", double(layout.window.delays)}, {"sample_stride", double(layout.window.sampleStride)},
        {"delay_stride", double(layout.window.delayStride)}, {"label_offset", double(layout.labelOffset)},
        {"channels", channels}, {"rows", double(layout.rows())}, {"materialized", true}};
}

}  // namespace dynsyslab::core
