// Candidate library families frozen on training data (ModelLibraries.hpp).
#include <dynsyslab/core/ModelLibraries.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using E = AnalyticExpression;

// Index lists of the monomials of degree 1..degree in n variables (combinations with replacement).
std::vector<std::vector<std::size_t>> monomials(std::size_t n, std::size_t degree) {
    std::vector<std::vector<std::size_t>> out;
    std::vector<std::size_t> current;
    std::function<void(std::size_t, std::size_t)> grow = [&](std::size_t start, std::size_t left) {
        if (!left) { out.push_back(current); return; }
        for (std::size_t j = start; j < n; ++j) { current.push_back(j); grow(j, left - 1); current.pop_back(); }
    };
    for (std::size_t d = 1; d <= degree; ++d) grow(0, d);
    return out;
}

E product(const std::vector<E>& factors, const std::vector<std::size_t>& indices, std::size_t n) {
    auto term = E::constant(1, n);
    for (const auto i : indices) term = term * factors[i];
    return term;
}

std::string monomialName(const std::string& symbol, const std::vector<std::size_t>& indices) {
    std::string name;
    for (const auto i : indices) name += (name.empty() ? "" : "*") + symbol + std::to_string(i + 1);
    return name;
}

}  // namespace

std::vector<std::string> modelLibraryFamilies() {
    return {"polynomial", "bounded", "rational", "localized", "gated", "paired_gated", "harmonic"};
}

Observable modelLibrary(const std::string& family, const TimeSeries& training, const ModelLibraryOptions& o) {
    const auto n = training.componentCount(), rows = training.sampleCount();
    if (!n || !rows || !o.degree || !o.centers || !o.harmonics) throw std::invalid_argument("model library: empty training data or zero sizes");
    std::vector<E> x, z;
    for (std::size_t j = 0; j < n; ++j) {
        x.push_back(E::state(n, j));
        z.push_back(x.back() / (1 + (x.back() / 3).power(2)).power(0.5));
    }
    std::vector<E> terms;
    std::vector<std::string> names;
    const auto add = [&](E term, std::string name) { terms.push_back(std::move(term)); names.push_back(std::move(name)); };
    const auto energy = [&] { auto e = E::constant(0, n); for (const auto& v : x) e = e + v.power(2); return e / double(n); };
    if (family == "polynomial" || family == "bounded" || family == "rational" || family == "localized") {
        add(E::constant(1, n), "1");
        if (family != "polynomial") for (std::size_t j = 0; j < n; ++j) add(x[j], "x" + std::to_string(j + 1));
        const auto& factors = family == "polynomial" ? x : z;
        for (const auto& m : monomials(n, o.degree)) add(product(factors, m, n), monomialName(family == "polynomial" ? "x" : "z", m));
        if (family == "rational")
            for (const double w : o.rationalWidths) {
                if (!(w > 0) || !std::isfinite(w)) throw std::invalid_argument("rational widths must be positive");
                for (std::size_t j = 0; j < n; ++j) {
                    const auto denominator = w * w + x[j].power(2);
                    add(E::constant(w * w, n) / denominator, "even_" + std::to_string(w) + "_x" + std::to_string(j + 1));
                    add(w * x[j] / denominator, "odd_" + std::to_string(w) + "_x" + std::to_string(j + 1));
                }
            }
        if (family == "localized") {
            // Centers: farthest-point sampling of the bounded training coordinates.
            const auto stride = std::max<std::size_t>(1, rows / 20000);
            std::vector<double> data;
            for (std::size_t k = 0; k < rows; k += stride)
                for (std::size_t j = 0; j < n; ++j) { const double v = training.value(k, j); data.push_back(v / std::hypot(1.0, v / 3)); }
            const auto count = data.size() / n;
            const auto norm2 = [&](std::size_t a, std::size_t b) {
                double s = 0; for (std::size_t j = 0; j < n; ++j) s += std::pow(data[a * n + j] - (b == count ? 0.0 : data[b * n + j]), 2); return s;
            };
            std::vector<std::size_t> chosen;
            std::size_t first = 0;
            for (std::size_t k = 1; k < count; ++k) if (norm2(k, count) < norm2(first, count)) first = k;
            chosen.push_back(first);
            std::vector<double> distance(count, std::numeric_limits<double>::infinity());
            while (chosen.size() < std::min(o.centers, count)) {
                for (std::size_t k = 0; k < count; ++k) distance[k] = std::min(distance[k], norm2(k, chosen.back()));
                const auto far = std::size_t(std::max_element(distance.begin(), distance.end()) - distance.begin());
                if (distance[far] <= 1e-14) break;
                chosen.push_back(far);
            }
            const auto C = chosen.size(), kth = std::min<std::size_t>(3, C - 1);
            std::vector<E> weights;
            auto total = E::constant(0, n);
            for (std::size_t a = 0; a < C; ++a) {
                std::vector<double> d2;
                for (std::size_t b = 0; b < C; ++b) d2.push_back(norm2(chosen[a], chosen[b]));
                std::nth_element(d2.begin(), d2.begin() + std::ptrdiff_t(kth), d2.end());
                const double width2 = std::max(0.05, d2[kth]);
                auto d = E::constant(0, n);
                for (std::size_t j = 0; j < n; ++j) d = d + (z[j] - data[chosen[a] * n + j]).power(2);
                weights.push_back((1 + d / width2).power(-2));
                total = total + weights.back();
            }
            for (std::size_t a = 0; a < C; ++a) {
                const auto w = weights[a] / total;
                add(w, "w" + std::to_string(a + 1));
                for (std::size_t j = 0; j < n; ++j) add(w * z[j], "w" + std::to_string(a + 1) + "*z" + std::to_string(j + 1));
            }
        }
    } else if (family == "gated" || family == "paired_gated") {
        const auto e = energy();
        std::vector<std::pair<E, std::string>> gates{{E::constant(1, n), "1"}};
        for (const double s : {0.04, 0.4, 4.0, 16.0}) gates.emplace_back(e / (e + s), "E/(E+" + std::to_string(s) + ")");
        if (family == "paired_gated")
            for (std::size_t j = 0; j + 1 < n; j += 2) {
                const auto p = (x[j].power(2) + x[j + 1].power(2)) / 2;
                gates.emplace_back(p / (1 + p), "p" + std::to_string(j / 2 + 1) + "/(1+p)");
            }
        for (const auto& [g, label] : gates)
            for (std::size_t j = 0; j < n; ++j) add(g * x[j], label + "*x" + std::to_string(j + 1));
    } else if (family == "harmonic") {
        if (n < 2) throw std::invalid_argument("harmonic libraries need a modal pair");
        const auto p = x[0].power(2) + x[1].power(2), e = energy();
        const std::vector<E> gates{E::constant(1, n), e / (e + 0.1), e / (e + 1), e / (e + 6), p / (p + 1)};
        const auto amplitude = (e + 0.015).power(0.5);
        for (std::size_t j = 0; j < n; ++j) add(x[j], "x" + std::to_string(j + 1));
        auto re = E::constant(1, n), im = E::constant(0, n);  // (x1 + i x2)^k by recurrence
        for (std::size_t k = 1; k <= o.harmonics; ++k) {
            const auto nextRe = re * x[0] - im * x[1], nextIm = re * x[1] + im * x[0];
            re = nextRe; im = nextIm;
            const auto scale = amplitude / (p + 0.015).power(double(k) / 2);
            for (std::size_t g = 0; g < gates.size(); ++g) {
                add(gates[g] * scale * re, "g" + std::to_string(g) + "*Re u^" + std::to_string(k));
                add(gates[g] * scale * im, "g" + std::to_string(g) + "*Im u^" + std::to_string(k));
            }
        }
    } else {
        throw std::invalid_argument("unknown model library family " + family);
    }
    if (terms.size() > 4096) throw std::invalid_argument("model library exceeds 4096 terms");
    return Observable::fromExpressions({"model_library", family + " model library", "", "model", names, false}, terms);
}

}  // namespace dynsyslab::core
