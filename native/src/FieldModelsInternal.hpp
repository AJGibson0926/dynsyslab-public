#pragma once
// Shared pieces of the field models (FieldModels.hpp).
#include <dynsyslab/core/FieldModels.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace dynsyslab::core::detail {

json::Value numbersJson(const std::vector<double>& v);
std::vector<double> jsonNumbers(const json::Value& v);

// Indices and distances of the k nearest rows of `points` (rows x d) to q, nearest first.
std::vector<std::pair<double, std::size_t>> nearest(const std::vector<double>& points, std::size_t d, const double* q, std::size_t k);

// Normalized Wendland weights (1 - r)^4 (1 + 4 r), r = d_i / d_k, over the first k of k + 1 nearest
// (equal weights when every ratio reaches 1).
std::vector<double> wendlandWeights(const std::vector<std::pair<double, std::size_t>>& hits);

// Solves the symmetric positive-definite system A x = b (n x n, n x m) by Cholesky factorization.
std::vector<double> choleskySolve(std::vector<double> A, std::size_t n, std::vector<double> b, std::size_t m);

std::vector<double> trapezoidWeights(const std::vector<double>& time);

Observable compactRadialField(const json::Value& definition);
Observable localTangentField(const json::Value& definition);
Observable pairedPhaseField(const json::Value& definition);

}  // namespace dynsyslab::core::detail
