#include <dynsyslab/core/ExpansionProjection.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Complex = std::complex<double>;

// Solves (M) X = R, M n x n Hermitian positive semidefinite (normal equations), R n x m, with a tiny ridge.
std::vector<Complex> solveNormal(std::vector<Complex> M, std::size_t n, std::vector<Complex> R, std::size_t m) {
    double scale = 0;
    for (std::size_t i = 0; i < n; ++i) scale = std::max(scale, std::abs(M[i * n + i]));
    for (std::size_t i = 0; i < n; ++i) M[i * n + i] += 1e-14 * (scale > 0 ? scale : 1.0);
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r) if (std::abs(M[r * n + c]) > std::abs(M[p * n + c])) p = r;
        if (!(std::abs(M[p * n + c]) > 0)) throw std::invalid_argument("projection basis is rank deficient");
        if (p != c) {
            for (std::size_t j = 0; j < n; ++j) std::swap(M[c * n + j], M[p * n + j]);
            for (std::size_t j = 0; j < m; ++j) std::swap(R[c * m + j], R[p * m + j]);
        }
        for (std::size_t r = c + 1; r < n; ++r) {
            const auto f = M[r * n + c] / M[c * n + c];
            for (std::size_t j = c; j < n; ++j) M[r * n + j] -= f * M[c * n + j];
            for (std::size_t j = 0; j < m; ++j) R[r * m + j] -= f * R[c * m + j];
        }
    }
    for (std::size_t r = n; r-- > 0;)
        for (std::size_t j = 0; j < m; ++j) {
            auto v = R[r * m + j];
            for (std::size_t k = r + 1; k < n; ++k) v -= M[r * n + k] * R[k * m + j];
            R[r * m + j] = v / M[r * n + r];
        }
    return R;
}

}  // namespace

std::optional<Observable> expressionObservable(const FunctionExpansion& expansion, ObservableDescriptor descriptor) {
    auto expressions = expansion.realExpressions();
    if (!expressions) return std::nullopt;
    descriptor.complexValued = false;
    descriptor.outputNames.clear();
    for (std::size_t r = 0; r < expansion.outputCount(); ++r) {
        const auto base = r < expansion.outputNames().size() ? expansion.outputNames()[r] : descriptor.name + "[" + std::to_string(r) + "]";
        descriptor.outputNames.push_back(base + ".re");
        descriptor.outputNames.push_back(base + ".im");
    }
    return Observable::fromExpressions(std::move(descriptor), std::move(*expressions));
}

ExpansionProjection projectExpansion(const FunctionExpansion& observable, const FunctionExpansion& onto,
                                     const std::vector<double>& points, double tolerance) {
    if (observable.coordinateDimension() != onto.coordinateDimension())
        throw std::invalid_argument("observable and eigenfunctions live on different domains");
    const auto g0 = observable.canonical(), p0 = onto.canonical();
    const auto merged = mergeDictionaries(g0.dictionary(), p0.dictionary());
    const auto g = g0.inDictionary(merged), p = p0.inDictionary(merged);
    const std::size_t K = merged->size(), outputs = g.outputCount(), modes = p.outputCount(), d = merged->dimension();
    const auto& G = g.coefficients();  // outputs x K
    const auto& C = p.coefficients();  // modes x K
    ExpansionProjection out{points.empty() ? "coefficients" : "samples", {}, {}, false, g};
    // Rows of the design (in the chosen inner product): coefficient space, or sampled values.
    std::vector<Complex> A, Y;  // A: rows x modes, Y: rows x outputs
    std::size_t rows = 0;
    if (points.empty()) {
        rows = K;
        A.resize(K * modes);
        Y.resize(K * outputs);
        for (std::size_t k = 0; k < K; ++k) {
            for (std::size_t j = 0; j < modes; ++j) A[k * modes + j] = C[j * K + k];
            for (std::size_t r = 0; r < outputs; ++r) Y[k * outputs + r] = G[r * K + k];
        }
    } else {
        if (points.size() % d) throw std::invalid_argument("projection sample points must be rows of the domain coordinates");
        rows = points.size() / d;
        A.resize(rows * modes);
        Y.resize(rows * outputs);
        for (std::size_t s = 0; s < rows; ++s) {
            const auto phi = merged->evaluate(points.data() + s * d);
            for (std::size_t j = 0; j < modes; ++j) for (std::size_t k = 0; k < K; ++k) A[s * modes + j] += C[j * K + k] * phi[k];
            for (std::size_t r = 0; r < outputs; ++r) for (std::size_t k = 0; k < K; ++k) Y[s * outputs + r] += G[r * K + k] * phi[k];
        }
    }
    std::vector<Complex> M(modes * modes), R(modes * outputs);
    for (std::size_t s = 0; s < rows; ++s)
        for (std::size_t a = 0; a < modes; ++a) {
            for (std::size_t b = 0; b < modes; ++b) M[a * modes + b] += std::conj(A[s * modes + a]) * A[s * modes + b];
            for (std::size_t r = 0; r < outputs; ++r) R[a * outputs + r] += std::conj(A[s * modes + a]) * Y[s * outputs + r];
        }
    const auto X = solveNormal(std::move(M), modes, std::move(R), outputs);  // modes x outputs
    out.coefficients.resize(outputs * modes);
    for (std::size_t r = 0; r < outputs; ++r) for (std::size_t j = 0; j < modes; ++j) out.coefficients[r * modes + j] = X[j * outputs + r];
    out.exact = true;
    for (std::size_t r = 0; r < outputs; ++r) {
        double residual = 0, norm = 0;
        for (std::size_t s = 0; s < rows; ++s) {
            Complex fit{};
            for (std::size_t j = 0; j < modes; ++j) fit += A[s * modes + j] * X[j * outputs + r];
            residual += std::norm(Y[s * outputs + r] - fit);
            norm += std::norm(Y[s * outputs + r]);
        }
        out.relativeResidual.push_back(std::sqrt(residual / std::max(norm, 1e-300)));
        out.exact = out.exact && out.relativeResidual.back() <= tolerance;
    }
    std::vector<Complex> rest(outputs * K);
    for (std::size_t r = 0; r < outputs; ++r)
        for (std::size_t k = 0; k < K; ++k) {
            Complex v = G[r * K + k];
            for (std::size_t j = 0; j < modes; ++j) v -= out.coefficients[r * modes + j] * C[j * K + k];
            rest[r * K + k] = v;
        }
    out.residual = FunctionExpansion(merged, std::move(rest), g.coordinates(), g.outputNames());
    return out;
}

FunctionExpansion stackExpansions(const std::vector<FunctionExpansion>& parts) {
    if (parts.empty()) throw std::invalid_argument("nothing to stack");
    auto dictionary = parts.front().dictionaryHandle();
    for (std::size_t k = 1; k < parts.size(); ++k) dictionary = mergeDictionaries(*dictionary, parts[k].dictionary());
    std::vector<std::complex<double>> c;
    std::vector<std::string> names;
    for (const auto& part : parts) {
        const auto moved = part.inDictionary(dictionary);
        c.insert(c.end(), moved.coefficients().begin(), moved.coefficients().end());
        for (std::size_t r = 0; r < part.outputCount(); ++r)
            names.push_back(r < part.outputNames().size() ? part.outputNames()[r] : "g" + std::to_string(names.size()));
    }
    return FunctionExpansion(dictionary, std::move(c), parts.front().coordinates(), std::move(names));
}

namespace {
const FunctionExpansion& expansionOf(const Observable& o) {
    const auto* e = dynamic_cast<const FunctionExpansion*>(o.representation().get());
    if (!e) throw std::invalid_argument(o.descriptor().name + " is not represented by an expansion");
    return *e;
}
}  // namespace

ExpansionProjection projectObservable(const Observable& observable, const std::vector<Eigenfunction>& onto,
                                      const std::vector<double>& points, double tolerance) {
    std::vector<FunctionExpansion> modes;
    for (const auto& psi : onto) modes.push_back(expansionOf(psi));
    return projectExpansion(expansionOf(observable), stackExpansions(modes), points, tolerance);
}

}  // namespace dynsyslab::core
