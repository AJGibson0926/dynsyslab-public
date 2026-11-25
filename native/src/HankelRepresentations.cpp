#include <dynsyslab/core/HankelRepresentations.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <thread>

namespace dynsyslab::core {

std::shared_ptr<WindowMapRepresentation> WindowMapRepresentation::fromModes(HankelWindowLayout layout,
    std::shared_ptr<const std::vector<std::complex<double>>> modes,
    std::shared_ptr<const std::vector<std::complex<double>>> gram, std::size_t rank) {
    std::shared_ptr<WindowMapRepresentation> out(new WindowMapRepresentation);
    if (!modes || !gram || !rank || modes->size() != layout.rows() * rank || gram->size() != rank * rank)
        throw std::invalid_argument("window map needs rows x rank modes and a rank x rank Gram matrix");
    out->layout_ = std::move(layout);
    out->modeCount_ = rank;
    out->modes_ = std::move(modes);
    out->gram_ = std::move(gram);
    return out;
}

std::shared_ptr<WindowMapRepresentation> WindowMapRepresentation::fromMap(HankelWindowLayout layout,
    std::vector<std::complex<double>> map, std::size_t modes) {
    std::shared_ptr<WindowMapRepresentation> out(new WindowMapRepresentation);
    if (!modes || map.size() != modes * layout.rows())
        throw std::invalid_argument("window map must be modes x window rows");
    out->layout_ = std::move(layout);
    out->modeCount_ = modes;
    out->map_ = std::move(map);
    return out;
}

const std::vector<std::complex<double>>& WindowMapRepresentation::factor() const {
    std::call_once(choleskyOnce_, [&] {  // Hermitian positive definite: G = L L^H
        const std::size_t n = modeCount_;
        auto& L = cholesky_;
        L.assign(n * n, {});
        for (std::size_t j = 0; j < n; ++j) {
            std::complex<double> diagonal = (*gram_)[j * n + j];
            for (std::size_t k = 0; k < j; ++k) diagonal -= L[j * n + k] * std::conj(L[j * n + k]);
            if (!(diagonal.real() > 0.0)) throw std::runtime_error("window-map Gram matrix is not positive definite");
            const double d = std::sqrt(diagonal.real());
            L[j * n + j] = d;
            for (std::size_t i = j + 1; i < n; ++i) {
                std::complex<double> value = (*gram_)[i * n + j];
                for (std::size_t k = 0; k < j; ++k) value -= L[i * n + k] * std::conj(L[j * n + k]);
                L[i * n + j] = value / d;
            }
        }
    });
    return cholesky_;
}

std::vector<std::complex<double>> WindowMapRepresentation::solveGram(std::vector<std::complex<double>> b) const {
    const auto& L = factor();
    const std::size_t n = modeCount_;
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

std::vector<std::complex<double>> WindowMapRepresentation::mapRows(const std::vector<std::size_t>& modes) const {
    const std::size_t rows = layout_.rows();
    std::vector<std::complex<double>> out(modes.size() * rows);
    if (materialized()) {
        for (std::size_t i = 0; i < modes.size(); ++i)
            std::copy_n(map_.begin() + std::ptrdiff_t(modes.at(i) * rows), rows, out.begin() + std::ptrdiff_t(i * rows));
        return out;
    }
    // Row s of T = G^-1 Phi^H is conj(Phi g_s)^T with g_s = G^-1 e_s (G Hermitian). Rows are
    // independent, so large requests are split across threads.
    const std::size_t rank = modeCount_;
    (void)factor();  // factor once, before the threads share it
    auto fill = [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            std::vector<std::complex<double>> unit(rank);
            unit.at(modes.at(i)) = 1.0;
            const auto column = solveGram(std::move(unit));
            for (std::size_t r = 0; r < rows; ++r) {
                std::complex<double> sum{};
                for (std::size_t a = 0; a < rank; ++a) sum += column[a] * (*modes_)[r * rank + a];
                out[i * rows + r] = std::conj(sum);
            }
        }
    };
    const std::size_t workers = modes.size() * rows * rank < (std::size_t(1) << 24) ? 1
        : std::max<std::size_t>(1, std::min<std::size_t>(std::thread::hardware_concurrency(), modes.size()));
    if (workers == 1) { fill(0, modes.size()); return out; }
    std::vector<std::thread> threads;
    std::vector<std::exception_ptr> errors(workers);
    for (std::size_t w = 0; w < workers; ++w)
        threads.emplace_back([&, w] {
            try { fill(modes.size() * w / workers, modes.size() * (w + 1) / workers); }
            catch (...) { errors[w] = std::current_exception(); }
        });
    for (auto& thread : threads) thread.join();
    for (const auto& error : errors) if (error) std::rethrow_exception(error);
    return out;
}

std::vector<std::complex<double>> WindowMapRepresentation::evaluate(const double* point, const std::string& domain) const {
    if (domain != "window") throw std::invalid_argument("the window expansion has only the window domain");
    const std::size_t rows = layout_.rows(), channels = layout_.channels.size();
    std::vector<double> scaled(point, point + rows);
    for (std::size_t r = 0; r < rows; ++r) scaled[r] /= layout_.channels[r % channels].scale;
    std::vector<std::complex<double>> psi(modeCount_);
    if (materialized()) {
        for (std::size_t j = 0; j < modeCount_; ++j)
            for (std::size_t r = 0; r < rows; ++r) psi[j] += map_[j * rows + r] * scaled[r];
        return psi;
    }
    for (std::size_t r = 0; r < rows; ++r)  // b = Phi^H h, then psi = G^-1 b
        for (std::size_t a = 0; a < modeCount_; ++a) psi[a] += std::conj((*modes_)[r * modeCount_ + a]) * scaled[r];
    return solveGram(std::move(psi));
}

std::vector<std::complex<double>> WindowMapRepresentation::evaluateTraining(const EigenfunctionSamples& samples,
    const std::string& domain, const TimeSeriesRecord* input) const {
    if (!input) throw std::invalid_argument("window-map errors need the input time series the set was fitted on");
    std::vector<std::complex<double>> out;
    out.reserve(samples.samples() * modeCount_);
    for (const auto label : samples.sampleIndices) {
        const auto h = layout_.window_at(input->series, label);
        const auto psi = evaluate(h.data(), domain);
        out.insert(out.end(), psi.begin(), psi.end());
    }
    return out;
}

std::shared_ptr<const FunctionRepresentation> WindowMapRepresentation::function(const std::string& domain) const {
    if (domain != "window") throw std::invalid_argument("the window expansion has only the window domain");
    // psi_j(h) = sum_r T_jr h_r / scale_r: one coordinate term per window row, on the window domain,
    // so the function is evaluated along any series by cutting its window at each sample.
    const std::size_t rows = layout_.rows(), channels = layout_.channels.size();
    std::vector<std::size_t> modes(modeCount_);
    for (std::size_t j = 0; j < modeCount_; ++j) modes[j] = j;
    const auto T = materialized() ? map_ : mapRows(modes);
    std::vector<DictionaryTerm> terms(rows);
    for (std::size_t r = 0; r < rows; ++r) { terms[r].kind = "coordinate"; terms[r].index = r; }
    std::vector<std::complex<double>> c(modeCount_ * rows);
    for (std::size_t j = 0; j < modeCount_; ++j)
        for (std::size_t r = 0; r < rows; ++r) c[j * rows + r] = T[j * rows + r] / layout_.channels[r % channels].scale;
    return std::make_shared<const FunctionExpansion>(std::make_shared<const ExpansionDictionary>(rows, std::move(terms)),
                                                     std::move(c), FunctionDomain::seriesWindow(layout_));
}

std::shared_ptr<const EigenfunctionRepresentation> WindowMapRepresentation::extract(
    const std::vector<std::size_t>& modes) const {
    return fromMap(layout_, mapRows(modes), modes.size());
}

json::Value WindowMapRepresentation::describe() const {
    json::Value::Array channels;
    for (std::size_t c = 0; c < layout_.channels.size(); ++c)
        channels.emplace_back(json::Value::Object{{"column", double(layout_.channels[c].input)},
            {"lag", double(layout_.channels[c].lag)}, {"scale", layout_.channels[c].scale},
            {"observable", c < layout_.channelObservables.size() ? layout_.channelObservables[c] : ""}});
    return json::Value::Object{{"kind", "expansion"}, {"family", "linear-window"}, {"form", "psi = W^-1 U^T h"},
        {"causal", layout_.alignEnd}, {"label", layout_.alignEnd ? "newest sample" : "oldest sample"},
        {"delays", double(layout_.window.delays)}, {"sample_stride", double(layout_.window.sampleStride)},
        {"delay_stride", double(layout_.window.delayStride)}, {"label_offset", double(layout_.labelOffset)},
        {"channels", channels}, {"rows", double(layout_.rows())}, {"materialized", materialized()}};
}

}  // namespace dynsyslab::core
