#pragma once
// First-class discrete bilinear reduced model in lifted coordinates z:
//
//   z_plus = A z + B0 u + sum_i u_i B_i z,    z(target) = 0, u(target) = 0.
//
// There is no affine term, so (z, u) = (0, 0) is an equilibrium by construction.
// Each data pair must hold u constant across its sampling interval.
//
// References: R. Strasser, K. Worthmann, et al., "SafEDMD", arXiv:2402.03145
// (constant-input data, lifted coordinates vanishing at the target);
// D. Goswami & D. A. Paley, "Bilinearization, reachability, and optimal control
// of control-affine nonlinear systems: A Koopman spectral approach", IEEE Trans.
// Automat. Control 67 (2022) 2715-2728; D. Bruder, X. Fu & R. Vasudevan,
// "Advantages of bilinear Koopman realizations for the modeling and control of
// systems with unknown dynamics", IEEE Robot. Autom. Lett. 6 (2021) 4369-4376.
#include <dynsyslab/core/Json.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace dynsyslab::core {

struct BilinearKrom {
    std::size_t stateDimension{}, inputDimension{};
    double sampleTime{};                // 0 when not declared
    std::vector<double> A;              // n x n, row-major
    std::vector<double> B0;             // n x m, row-major
    std::vector<std::vector<double>> B; // m matrices, each n x n, row-major
    std::string provenance;
};

// Sample-major pairs (z_k, u_k) -> zNext_k. Pairs may pool separate trajectories.
struct BilinearDataset {
    std::size_t samples{};
    std::vector<double> z, zNext; // samples x n
    std::vector<double> u;        // samples x m
};

struct BilinearFitOptions {
    double ridge{0.0};            // Tikhonov weight relative to the largest squared singular value
    double rankTolerance{1e-10};  // relative to the largest singular value of the scaled design
    bool allowRankDeficient{false};
    double sampleTime{0.0};
    std::string provenance;
};

struct BilinearResidualStats {
    std::vector<double> norms;  // ||residual_k||
    double rootMeanSquare{}, maximum{}, relative{};  // relative = ||R||_F / ||zNext||_F
};

struct BilinearFit {
    BilinearKrom model;
    std::size_t regressors{}, rank{};
    double conditionNumber{};      // of the column-scaled design
    bool identifiable{false};      // full column rank
    std::vector<double> residuals; // samples x n, zNext - prediction
    BilinearResidualStats stats;
};

struct FrozenInputLinearization {
    std::vector<double> A;      // n x n: A + sum_i u_i B_i
    std::vector<double> B;      // n x m: B0 + [B_1 z ... B_m z]
    std::vector<double> offset; // n: f(z, u) - A z - B u
};

void validateBilinearKrom(const BilinearKrom& model);
// Throws when the design is rank deficient unless allowRankDeficient is set.
[[nodiscard]] BilinearFit fitBilinearKrom(const BilinearDataset& data, std::size_t stateDimension,
    std::size_t inputDimension, const BilinearFitOptions& options = {});
// Frozen-A fit: regresses only B0 and B_i on zNext - A z. A is copied bit for bit into
// the result and never estimated, so autonomous (unforced) dynamics stay those of the
// source K-ROM. Diagnostics match fitBilinearKrom; regressors counts B0 and B_i only.
[[nodiscard]] BilinearFit fitBilinearKromFrozenA(const BilinearDataset& data,
    const std::vector<double>& frozenA, std::size_t stateDimension, std::size_t inputDimension,
    const BilinearFitOptions& options = {});
[[nodiscard]] std::vector<double> stepBilinearKrom(const BilinearKrom& model,
    const std::vector<double>& z, const std::vector<double>& u);
// inputs: steps x m. Returns (steps + 1) x n, starting at z0.
[[nodiscard]] std::vector<double> simulateBilinearKrom(const BilinearKrom& model,
    const std::vector<double>& z0, const std::vector<double>& inputs, std::size_t steps);
[[nodiscard]] FrozenInputLinearization linearizeBilinearKrom(const BilinearKrom& model,
    const std::vector<double>& z, const std::vector<double>& u);
// Returns samples x n residuals; fills stats when supplied.
[[nodiscard]] std::vector<double> bilinearResiduals(const BilinearKrom& model,
    const BilinearDataset& data, BilinearResidualStats* stats = nullptr);

// Schema dynsyslab.bilinear-krom/v1. Matrices are arrays of rows.
[[nodiscard]] json::Value bilinearKromToJson(const BilinearKrom& model);
[[nodiscard]] BilinearKrom bilinearKromFromJson(const json::Value& value);
[[nodiscard]] json::Value denseRowsToJson(const std::vector<double>& data, std::size_t rows, std::size_t columns);
// Reads an array of equal-length numeric rows; updates rows and columns.
[[nodiscard]] std::vector<double> denseRowsFromJson(const json::Value& value, std::size_t& rows, std::size_t& columns);

}  // namespace dynsyslab::core
