// Interpolant methods for scattered eigenfunction samples. All work on the k nearest training
// samples of the query point (one shared k-d tree search):
//   nearest           value of the nearest sample
//   inverse_distance  Shepard weights 1/d^p (Shepard, Proc. ACM 1968); p = 2 is the default
//   gaussian_kernel   Nadaraya-Watson weights exp(-d^2 / 2h^2)
//   local_linear      weighted affine least squares about the query (moving least squares,
//                     Lancaster & Salkauskas, Math. Comp. 37, 1981)
//   rbf_gaussian, rbf_multiquadric, rbf_thin_plate
//                     local radial-basis interpolation through the k samples; thin-plate splines
//                     carry an affine term (Buhmann, Radial Basis Functions, 2003)
// Width h (or the RBF shape length) defaults to the mean neighbour distance. A local system that
// is singular falls back to inverse-distance weights.
#include <dynsyslab/core/FunctionInterpolation.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<std::string> interpolationMethods() {
    return {"inverse_distance", "nearest", "gaussian_kernel", "local_linear", "rbf_gaussian", "rbf_multiquadric",
            "rbf_thin_plate"};
}

json::Value interpolationMethodJson(const InterpolationMethod& m, std::size_t neighbors) {
    return json::Value::Object{{"name", m.name}, {"neighbors", double(neighbors)}, {"power", m.power}, {"width", m.width}};
}

InterpolationMethod interpolationMethodFromJson(const json::Value& value, std::size_t* neighbors) {
    InterpolationMethod m;
    m.name = value.at("name").string("inverse_distance");
    m.power = value.at("power").number(2.0);
    m.width = value.at("width").number(0.0);
    const auto names = interpolationMethods();
    if (std::find(names.begin(), names.end(), m.name) == names.end())
        throw std::invalid_argument("unknown interpolation method " + m.name);
    if (!(m.power > 0) || !(m.width >= 0)) throw std::invalid_argument("interpolation power must be positive and width nonnegative");
    if (neighbors && !value.at("neighbors").isNull()) *neighbors = std::size_t(value.at("neighbors").number());
    return m;
}

namespace {

using Complex = std::complex<double>;

// Solves A X = B in place (A n x n row-major, B n x m), partial pivoting. False if singular.
bool solve(std::vector<double> A, std::size_t n, std::vector<Complex>& B, std::size_t m) {
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r) if (std::abs(A[r * n + c]) > std::abs(A[p * n + c])) p = r;
        if (!(std::abs(A[p * n + c]) > 1e-14)) return false;
        if (p != c) {
            for (std::size_t j = 0; j < n; ++j) std::swap(A[c * n + j], A[p * n + j]);
            for (std::size_t j = 0; j < m; ++j) std::swap(B[c * m + j], B[p * m + j]);
        }
        for (std::size_t r = c + 1; r < n; ++r) {
            const double f = A[r * n + c] / A[c * n + c];
            if (f == 0.0) continue;
            for (std::size_t j = c; j < n; ++j) A[r * n + j] -= f * A[c * n + j];
            for (std::size_t j = 0; j < m; ++j) B[r * m + j] -= f * B[c * m + j];
        }
    }
    for (std::size_t r = n; r-- > 0;) {
        for (std::size_t j = 0; j < m; ++j) {
            Complex v = B[r * m + j];
            for (std::size_t k = r + 1; k < n; ++k) v -= A[r * n + k] * B[k * m + j];
            B[r * m + j] = v / A[r * n + r];
        }
    }
    return true;
}

std::vector<Complex> weighted(const InterpolationData& d, const std::vector<std::pair<double, std::size_t>>& nb,
                              const std::vector<double>& w) {
    const std::size_t modes = d.outputs();
    std::vector<Complex> out(modes);
    double total = 0;
    for (std::size_t i = 0; i < nb.size(); ++i) {
        total += w[i];
        for (std::size_t j = 0; j < modes; ++j) out[j] += w[i] * d.trainingValues[nb[i].second * modes + j];
    }
    for (auto& v : out) v /= total;
    return out;
}

std::vector<Complex> inverseDistance(const InterpolationData& d, const std::vector<std::pair<double, std::size_t>>& nb,
                                     double power) {
    std::vector<double> w;
    for (const auto& n : nb) w.push_back(1.0 / std::pow(std::max(n.first, 1.0e-26), power / 2.0));
    return weighted(d, nb, w);
}

double meanDistance(const std::vector<std::pair<double, std::size_t>>& nb) {
    double sum = 0;
    for (const auto& n : nb) sum += std::sqrt(n.first);
    return std::max(sum / double(nb.size()), 1e-12);
}

double kernel(const std::string& name, double r, double width) {
    const double z = r / width;
    if (name == "rbf_gaussian") return std::exp(-z * z);
    if (name == "rbf_multiquadric") return std::sqrt(1 + z * z);
    return r > 0 ? r * r * std::log(r) : 0.0;  // thin plate
}

}  // namespace

std::vector<std::complex<double>> interpolateScatteredNeighbours(const InterpolationData& d,
                                           const std::vector<std::pair<double, std::size_t>>& nb, const double* point) {
    const std::size_t modes = d.outputs(), dim = d.stateDimension, k = nb.size();
    const auto& m = d.method;
    if (m.name == "nearest" || nb.front().first <= 1.0e-26)
        return std::vector<Complex>(d.trainingValues.begin() + std::ptrdiff_t(nb.front().second * modes),
                                    d.trainingValues.begin() + std::ptrdiff_t((nb.front().second + 1) * modes));
    if (m.name == "inverse_distance") return inverseDistance(d, nb, m.power);
    const double width = m.width > 0 ? m.width : meanDistance(nb);
    if (m.name == "gaussian_kernel") {
        std::vector<double> w;
        for (const auto& n : nb) w.push_back(std::exp(-n.first / (2 * width * width)));
        if (*std::max_element(w.begin(), w.end()) > 0) return weighted(d, nb, w);
        return inverseDistance(d, nb, 2.0);
    }
    const auto coordinate = [&](std::size_t sample, std::size_t axis) { return d.trainingStates[sample * dim + axis]; };
    if (m.name == "local_linear") {  // value of the weighted affine fit at the query point
        const std::size_t p = dim + 1;
        if (k < p) return inverseDistance(d, nb, 2.0);
        std::vector<double> A(p * p, 0.0);
        std::vector<Complex> B(p * modes);
        for (const auto& n : nb) {
            const double w = std::exp(-n.first / (2 * width * width)) + 1e-12;
            std::vector<double> f{1.0};
            for (std::size_t a = 0; a < dim; ++a) f.push_back(coordinate(n.second, a) - point[a]);
            for (std::size_t r = 0; r < p; ++r) {
                for (std::size_t c = 0; c < p; ++c) A[r * p + c] += w * f[r] * f[c];
                for (std::size_t j = 0; j < modes; ++j) B[r * modes + j] += w * f[r] * d.trainingValues[n.second * modes + j];
            }
        }
        if (!solve(A, p, B, modes)) return inverseDistance(d, nb, 2.0);
        return std::vector<Complex>(B.begin(), B.begin() + std::ptrdiff_t(modes));
    }
    // Local RBF interpolation through the k samples (thin plate with an affine term).
    const bool thin = m.name == "rbf_thin_plate";
    const std::size_t extra = thin ? dim + 1 : 0, n = k + extra;
    std::vector<double> A(n * n, 0.0);
    std::vector<Complex> B(n * modes);
    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) {
            double r = 0;
            for (std::size_t a = 0; a < dim; ++a) r = std::hypot(r, coordinate(nb[i].second, a) - coordinate(nb[j].second, a));
            A[i * n + j] = kernel(m.name, r, width);
        }
        if (thin) {
            A[i * n + k] = A[k * n + i] = 1.0;
            for (std::size_t a = 0; a < dim; ++a) A[i * n + k + 1 + a] = A[(k + 1 + a) * n + i] = coordinate(nb[i].second, a) - point[a];
        }
        for (std::size_t j = 0; j < modes; ++j) B[i * modes + j] = d.trainingValues[nb[i].second * modes + j];
    }
    if (!solve(A, n, B, modes)) return inverseDistance(d, nb, 2.0);
    std::vector<Complex> out(modes);
    for (std::size_t i = 0; i < k; ++i) {
        const double phi = kernel(m.name, std::sqrt(nb[i].first), width);
        for (std::size_t j = 0; j < modes; ++j) out[j] += B[i * modes + j] * phi;
    }
    if (thin) for (std::size_t j = 0; j < modes; ++j) out[j] += B[k * modes + j];  // affine term at the query (offsets 0)
    return out;
}

}  // namespace dynsyslab::core
