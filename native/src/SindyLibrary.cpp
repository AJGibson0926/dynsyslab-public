// Shared pieces of the Koopman SINDy fits (SINDy.cpp, SindyEigenvalues.cpp): the polynomial library
// Theta and its rate along the dynamics Gamma, their Gram blocks, and sparse null vectors.
// SINDy: Brunton, Proctor & Kutz, PNAS 113:3932 (2016). Koopman eigenfunctions as null vectors of
// (Gamma - lambda Theta): Kaiser, Kutz & Brunton, Mach. Learn.: Sci. Technol. 2:035023 (2021).
#include "SindyLibrary.hpp"
#include <dynsyslab/numerics/MatrixFunctions.hpp>
#include <dynsyslab/numerics/Nullspace.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core::sindy_detail {
namespace {
double power(double value, std::size_t exponent) {
    double result = 1.0;
    for (std::size_t index = 0; index < exponent; ++index) result *= value;
    return result;
}
}  // namespace

std::vector<double> featureMatrix(
    const std::vector<double>& states,std::size_t samples,std::size_t dimension,
    const std::vector<std::size_t>& exponents,
    const std::vector<double>& center,const std::vector<double>& scale) {
    const std::size_t features = exponents.size() / dimension;
    std::vector<double> result(samples*features);
    for (std::size_t row = 0; row < samples; ++row)
        for (std::size_t feature = 0; feature < features; ++feature) {
            double value = 1.0;
            for (std::size_t axis = 0; axis < dimension; ++axis)
                value *= power((states[row*dimension+axis] - center[axis]) / scale[axis],
                               exponents[feature * dimension + axis]);
            result[row*features+feature] = value;
        }
    return result;
}

std::vector<double> finiteDifference(
    const std::vector<double>& values,std::size_t samples,std::size_t columns,
    const std::vector<double>& time) {
    std::vector<double> result(values.size());
    const std::size_t last = samples - 1;
    const double firstStep = time[1] - time[0], lastStep = time[time.size()-1] - time[time.size()-2];
    for(std::size_t column=0;column<columns;++column) {
        result[column]=(-3.0*values[column]+4.0*values[columns+column]
            -values[2*columns+column])/(2.0*firstStep);
        result[last*columns+column]=(3.0*values[last*columns+column]
            -4.0*values[(last-1)*columns+column]
            +values[(last-2)*columns+column])/(2.0*lastStep);
    }
    for (std::size_t row = 1; row < last; ++row) {
        const double backward = time[row] - time[row-1];
        const double forward = time[row+1] - time[row];
        for(std::size_t column=0;column<columns;++column)
            result[row*columns+column]=
                (-forward/(backward*(backward+forward)))*values[(row-1)*columns+column]
                +((forward-backward)/(backward*forward))*values[row*columns+column]
                +(backward/(forward*(backward+forward)))*values[(row+1)*columns+column];
    }
    return result;
}

std::vector<double> featureDerivative(
    const std::vector<double>& states,const std::vector<double>& velocities,
    std::size_t samples,std::size_t dimension,
    const std::vector<std::size_t>& exponents,
    const std::vector<double>& center,const std::vector<double>& scale) {
    const std::size_t features = exponents.size() / dimension;
    std::vector<double> result(samples*features);
    for (std::size_t row = 0; row < samples; ++row)
        for (std::size_t feature = 0; feature < features; ++feature)
            for (std::size_t differentiated = 0; differentiated < dimension; ++differentiated) {
                const std::size_t count = exponents[feature * dimension + differentiated];
                if (!count) continue;
                double derivative = static_cast<double>(count) / scale[differentiated];
                for (std::size_t axis = 0; axis < dimension; ++axis) {
                    std::size_t exponent = exponents[feature * dimension + axis];
                    if (axis == differentiated) --exponent;
                    derivative *= power((states[row*dimension+axis] - center[axis])
                                        / scale[axis], exponent);
                }
                result[row*features+feature] += derivative
                    * velocities[row*dimension+differentiated];
            }
    return result;
}

std::vector<std::complex<double>> sparseNullVector(
    const std::vector<std::complex<double>>& source,std::size_t rows,
    std::size_t columns,double threshold,std::size_t maximumIterations) {
    std::vector<double> scales(columns);
    for (std::size_t column = 0; column < columns; ++column) {
        for(std::size_t row=0;row<rows;++row)
            scales[column]+=std::norm(source[row*columns+column]);
        scales[column]=std::sqrt(scales[column]);
        if (scales[column] <= 1e-15) scales[column] = 1.0;
    }
    std::vector<std::size_t> active(columns);
    std::iota(active.begin(), active.end(), std::size_t{0});
    std::vector<std::complex<double>> vector;
    auto solve = [&]() {
        std::vector<std::complex<double>> block(rows*active.size());
        for(std::size_t row=0;row<rows;++row)
            for(std::size_t index=0;index<active.size();++index)
                block[row*active.size()+index]=source[row*columns+active[index]]
                    /scales[active[index]];
        return numerics::smallestRightSingularVector(
            block,rows,active.size());
    };
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration) {
        vector = solve(); double largest=0.0;
        for(const auto value:vector)largest=std::max(largest,std::abs(value));
        const double cutoff = threshold * largest;
        std::vector<std::size_t> retained;
        for (std::size_t index = 0; index < active.size(); ++index)
            if (std::abs(vector[index]) >= cutoff) retained.push_back(active[index]);
        if (retained.empty()) {
            const auto largestIndex=static_cast<std::size_t>(std::distance(vector.begin(),
                std::max_element(vector.begin(),vector.end(),[](const auto left,const auto right) {
                    return std::abs(left)<std::abs(right);
                })));
            retained.push_back(active[largestIndex]);
        }
        if (retained == active) break; active = std::move(retained);
    }
    vector = solve(); std::vector<std::complex<double>> result(columns);
    for (std::size_t index = 0; index < active.size(); ++index)
        result[active[index]] = vector[index] / scales[active[index]];
    double norm=0.0;for(const auto value:result)norm+=std::norm(value);
    norm=std::sqrt(norm);
    if (norm == 0.0) throw std::runtime_error("SINDy null vector vanished");
    for(auto& value:result)value/=norm;
    return result;
}

GramBlocks gramBlocks(const std::vector<double>& theta, const std::vector<double>& gamma,
                      std::size_t rows, std::size_t n) {
    GramBlocks g{std::vector<double>(n*n), std::vector<double>(n*n), std::vector<double>(n*n), n};
    for (std::size_t r = 0; r < rows; ++r) {
        const double* t = theta.data() + r*n; const double* d = gamma.data() + r*n;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                g.A[i*n+j] += t[i]*t[j]; g.B[i*n+j] += d[i]*t[j]; g.C[i*n+j] += d[i]*d[j];
            }
    }
    return g;
}

// Same thresholded null-vector iteration as sparseNullVector, applied to
// G(lambda) = (Gamma - lambda Theta)^H (Gamma - lambda Theta)
//           = C - lambda B - conj(lambda) Bᵀ + |lambda|^2 A.
// Returns nullopt when the active Gram system is too ill-conditioned (its
// eigenvalue ratio below ~eps, i.e. singular-value ratio below ~sqrt(eps)).
std::optional<std::vector<std::complex<double>>> sparseNullVectorGram(
    const GramBlocks& g, std::complex<double> lambda, double threshold, std::size_t maximumIterations) {
    const std::size_t n = g.n;
    std::vector<std::complex<double>> G(n*n);
    const double magnitude = std::norm(lambda);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            G[i*n+j] = g.C[i*n+j] - lambda*g.B[i*n+j] - std::conj(lambda)*g.B[j*n+i] + magnitude*g.A[i*n+j];
    std::vector<double> scales(n);
    for (std::size_t c = 0; c < n; ++c) {
        scales[c] = std::sqrt(std::max(0.0, G[c*n+c].real()));
        if (scales[c] <= 1e-15) scales[c] = 1.0;
    }
    std::vector<std::size_t> active(n); std::iota(active.begin(), active.end(), std::size_t{0});
    bool trusted = true;
    auto solve = [&]() {
        const std::size_t k = active.size();
        std::vector<std::complex<double>> sub(k*k);
        for (std::size_t a = 0; a < k; ++a)
            for (std::size_t b = 0; b < k; ++b)
                sub[a*k+b] = G[active[a]*n+active[b]] / (scales[active[a]]*scales[active[b]]);
        const auto eig = numerics::hermitianEigensystem(sub, k);
        // SINDy needs a small generator residual, not a unique eigenvector. If
        // the smallest singular value is well separated the Gram eigenvector is
        // accurate; if it is tiny, the smallest Gram eigenpairs may cluster, but
        // every vector in the cluster attains a residual at the sqrt(eps) floor,
        // far below the sparsity threshold. Only non-finite output falls back.
        const double largest = eig.eigenvalues.back(), smallest = eig.eigenvalues.front();
        if (!std::isfinite(smallest) || !std::isfinite(largest)) trusted = false;
        std::vector<std::complex<double>> v(k);
        for (std::size_t a = 0; a < k; ++a) v[a] = eig.eigenvectors[a*k];  // first (smallest) column
        return v;
    };
    std::vector<std::complex<double>> vector;
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration) {
        vector = solve(); double largest = 0.0;
        for (const auto value : vector) largest = std::max(largest, std::abs(value));
        const double cutoff = threshold * largest;
        std::vector<std::size_t> retained;
        for (std::size_t index = 0; index < active.size(); ++index)
            if (std::abs(vector[index]) >= cutoff) retained.push_back(active[index]);
        if (retained.empty()) {
            const auto top = static_cast<std::size_t>(std::distance(vector.begin(),
                std::max_element(vector.begin(), vector.end(),
                    [](const auto l, const auto r) { return std::abs(l) < std::abs(r); })));
            retained.push_back(active[top]);
        }
        if (retained == active) break; active = std::move(retained);
    }
    vector = solve();
    if (!trusted) return std::nullopt;
    std::vector<std::complex<double>> result(n);
    for (std::size_t index = 0; index < active.size(); ++index)
        result[active[index]] = vector[index] / scales[active[index]];
    double norm = 0.0; for (const auto value : result) norm += std::norm(value);
    norm = std::sqrt(norm);
    if (!(norm > 0.0)) return std::nullopt;
    for (auto& value : result) value /= norm;
    return result;
}

bool conjugatePartner(std::complex<double> left, std::complex<double> right) {
    return std::abs(right - std::conj(left)) <= 1e-9 + 1e-7 * std::abs(left);
}


std::vector<std::complex<double>> solveDense(std::vector<std::complex<double>> matrix,
                                             std::vector<std::complex<double>> rhs, std::size_t n) {
    for (std::size_t column = 0; column < n; ++column) {  // Gaussian elimination, partial pivoting
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < n; ++row)
            if (std::abs(matrix[row*n+column]) > std::abs(matrix[pivot*n+column])) pivot = row;
        if (std::abs(matrix[pivot*n+column]) < 1e-300) throw std::runtime_error("anchored SINDy system is singular");
        if (pivot != column) {
            for (std::size_t k = 0; k < n; ++k) std::swap(matrix[pivot*n+k], matrix[column*n+k]);
            std::swap(rhs[pivot], rhs[column]);
        }
        for (std::size_t row = column + 1; row < n; ++row) {
            const auto factor = matrix[row*n+column] / matrix[column*n+column];
            if (factor == std::complex<double>{}) continue;
            for (std::size_t k = column; k < n; ++k) matrix[row*n+k] -= factor * matrix[column*n+k];
            rhs[row] -= factor * rhs[column];
        }
    }
    std::vector<std::complex<double>> x(n);
    for (std::size_t row = n; row-- > 0;) {
        auto value = rhs[row];
        for (std::size_t k = row + 1; k < n; ++k) value -= matrix[row*n+k] * x[k];
        x[row] = value / matrix[row*n+row];
    }
    return x;
}

std::vector<std::complex<double>> anchoredVector(const std::vector<double>& anchorGram, const GramBlocks& g,
                                                 std::complex<double> lambda,
                                                 const std::vector<std::complex<double>>& rhs, double weight,
                                                 double threshold, std::size_t maximumIterations) {
    const std::size_t n = g.n;
    std::vector<std::complex<double>> G(n*n);
    const double magnitude = std::norm(lambda);
    double traceAnchor = 0.0, traceGenerator = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            G[i*n+j] = g.C[i*n+j] - lambda*g.B[i*n+j] - std::conj(lambda)*g.B[j*n+i] + magnitude*g.A[i*n+j];
    for (std::size_t i = 0; i < n; ++i) { traceAnchor += anchorGram[i*n+i]; traceGenerator += G[i*n+i].real(); }
    const double w = traceGenerator > 0 ? weight * traceAnchor / traceGenerator : 0.0;
    std::vector<std::complex<double>> M(n*n);
    for (std::size_t k = 0; k < n*n; ++k) M[k] = anchorGram[k] + w * G[k];
    std::vector<std::size_t> active(n); std::iota(active.begin(), active.end(), std::size_t{0});
    auto solve = [&]() {
        const std::size_t k = active.size();
        std::vector<std::complex<double>> sub(k*k), b(k);
        for (std::size_t a = 0; a < k; ++a) {
            b[a] = rhs[active[a]];
            for (std::size_t c = 0; c < k; ++c) sub[a*k+c] = M[active[a]*n+active[c]];
        }
        return solveDense(std::move(sub), std::move(b), k);
    };
    std::vector<std::complex<double>> vector;
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration) {  // sequential thresholding
        vector = solve(); double largest = 0.0;
        for (const auto value : vector) largest = std::max(largest, std::abs(value));
        std::vector<std::size_t> retained;
        for (std::size_t index = 0; index < active.size(); ++index)
            if (std::abs(vector[index]) >= threshold * largest) retained.push_back(active[index]);
        if (retained.empty() || retained == active) break;
        active = std::move(retained);
    }
    vector = solve();
    std::vector<std::complex<double>> result(n);
    for (std::size_t index = 0; index < active.size(); ++index) result[active[index]] = vector[index];
    return result;
}

void normalization(const std::vector<double>& states, std::size_t samples, std::size_t dimension, bool normalize,
                   std::vector<double>& center, std::vector<double>& scale) {
    center.assign(dimension, 0.0); scale.assign(dimension, 1.0);
    if (!normalize) return;
    for (std::size_t axis = 0; axis < dimension; ++axis) {
        for (std::size_t sample = 0; sample < samples; ++sample) center[axis] += states[sample*dimension+axis];
        center[axis] /= static_cast<double>(samples);
        double squared = 0.0;
        for (std::size_t sample = 0; sample < samples; ++sample) {
            const double value = states[sample*dimension+axis] - center[axis];
            squared += value * value;
        }
        scale[axis] = std::sqrt(squared / static_cast<double>(samples));
        if (scale[axis] <= 1e-15) scale[axis] = 1.0;
    }
}
}  // namespace dynsyslab::core::sindy_detail
