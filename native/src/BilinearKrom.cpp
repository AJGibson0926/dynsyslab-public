// Least-squares fit and evaluation of z_plus = A z + B0 u + sum_i u_i B_i z.
// The constant-input regression follows SafEDMD (Strasser, Worthmann et al.,
// arXiv:2402.03145); the bilinear lifted form follows Goswami & Paley, IEEE TAC
// 67 (2022) 2715-2728, and Bruder, Fu & Vasudevan, IEEE RA-L 6 (2021) 4369-4376.
#include <dynsyslab/core/BilinearKrom.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ConstRowMap = Eigen::Map<const RowMatrix>;
using Index = Eigen::Index;

Index idx(std::size_t value) { return static_cast<Index>(value); }

void requireFinite(const std::vector<double>& values, const char* name) {
    for (double v : values)
        if (!std::isfinite(v)) throw std::invalid_argument(std::string(name) + " must be finite");
}

void validateData(const BilinearDataset& d, std::size_t n, std::size_t m) {
    if (n == 0 || m == 0 || d.samples == 0) throw std::invalid_argument("Bilinear data needs samples, states and inputs");
    if (d.z.size() != d.samples * n || d.zNext.size() != d.samples * n || d.u.size() != d.samples * m)
        throw std::invalid_argument("Bilinear data arrays have inconsistent sizes");
    requireFinite(d.z, "z"); requireFinite(d.zNext, "zNext"); requireFinite(d.u, "u");
}

std::vector<double> toVector(const RowMatrix& m) {
    return std::vector<double>(m.data(), m.data() + m.size());
}

struct ScaledSolve {
    MatrixXd coefficients;  // regressors x n
    std::size_t rank{};
    double conditionNumber{};
    bool identifiable{false};
};

// Column-scaled truncated-SVD least squares shared by the free and frozen-A fits.
ScaledSolve solveColumnScaled(MatrixXd design, const MatrixXd& target, const BilinearFitOptions& options) {
    const Index p = design.cols();
    VectorXd columnScale(p);
    for (Index j = 0; j < p; ++j) {
        const double norm = design.col(j).norm();
        columnScale(j) = norm > 0.0 ? norm : 1.0;
        design.col(j) /= columnScale(j);
    }
    const Eigen::BDCSVD<MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(design);
    if (svd.info() != Eigen::Success) throw std::runtime_error("Bilinear least-squares SVD failed");
    const VectorXd& sigma = svd.singularValues();
    const double largest = sigma.size() ? sigma(0) : 0.0;
    std::size_t rank = 0;
    for (Index i = 0; i < sigma.size(); ++i)
        if (sigma(i) > options.rankTolerance * largest) ++rank;
    const bool identifiable = rank == static_cast<std::size_t>(p);
    if (!identifiable && !options.allowRankDeficient)
        throw std::invalid_argument("Bilinear design is rank deficient: excite every input and state direction");

    VectorXd filter = VectorXd::Zero(sigma.size());
    const double penalty = options.ridge * largest * largest;
    for (Index i = 0; i < sigma.size(); ++i)
        if (sigma(i) > options.rankTolerance * largest) filter(i) = sigma(i) / (sigma(i) * sigma(i) + penalty);
    MatrixXd coefficients = svd.matrixV() * filter.asDiagonal() * (svd.matrixU().transpose() * target);
    for (Index j = 0; j < p; ++j) coefficients.row(j) /= columnScale(j);
    const double smallest = sigma.size() ? sigma(sigma.size() - 1) : 0.0;
    return {std::move(coefficients), rank, identifiable && smallest > 0.0 ? largest / smallest : HUGE_VAL, identifiable};
}

}  // namespace

void validateBilinearKrom(const BilinearKrom& model) {
    const auto n = model.stateDimension, m = model.inputDimension;
    if (n == 0 || m == 0) throw std::invalid_argument("Bilinear model needs positive dimensions");
    if (model.A.size() != n * n || model.B0.size() != n * m || model.B.size() != m)
        throw std::invalid_argument("Bilinear model matrices have the wrong size");
    requireFinite(model.A, "A"); requireFinite(model.B0, "B0");
    for (const auto& matrix : model.B) {
        if (matrix.size() != n * n) throw std::invalid_argument("Bilinear operator has the wrong size");
        requireFinite(matrix, "B_i");
    }
    if (!std::isfinite(model.sampleTime) || model.sampleTime < 0.0)
        throw std::invalid_argument("Bilinear sample time must be finite and non-negative");
}

BilinearFit fitBilinearKrom(const BilinearDataset& data, std::size_t n, std::size_t m,
    const BilinearFitOptions& options) {
    validateData(data, n, m);
    if (!(options.ridge >= 0.0) || !(options.rankTolerance > 0.0))
        throw std::invalid_argument("Invalid bilinear fit options");
    const std::size_t p = n + m + m * n;
    if (data.samples < p && !options.allowRankDeficient)
        throw std::invalid_argument("Bilinear fit needs at least as many pairs as regressors");

    const ConstRowMap Z(data.z.data(), idx(data.samples), idx(n));
    const ConstRowMap U(data.u.data(), idx(data.samples), idx(m));
    const ConstRowMap Y(data.zNext.data(), idx(data.samples), idx(n));
    MatrixXd design(idx(data.samples), idx(p));
    design.leftCols(idx(n)) = Z;
    design.middleCols(idx(n), idx(m)) = U;
    for (std::size_t i = 0; i < m; ++i)
        design.middleCols(idx(n + m + i * n), idx(n)) = U.col(idx(i)).asDiagonal() * Z;
    const auto solved = solveColumnScaled(std::move(design), Y, options);
    const MatrixXd& coefficients = solved.coefficients;

    BilinearFit fit;
    auto& model = fit.model;
    model.stateDimension = n; model.inputDimension = m;
    model.sampleTime = options.sampleTime;
    model.provenance = options.provenance.empty() ? "bilinear_fit constant-input least squares" : options.provenance;
    model.A = toVector(coefficients.topRows(idx(n)).transpose());
    model.B0 = toVector(coefficients.middleRows(idx(n), idx(m)).transpose());
    for (std::size_t i = 0; i < m; ++i)
        model.B.push_back(toVector(coefficients.middleRows(idx(n + m + i * n), idx(n)).transpose()));
    validateBilinearKrom(model);
    fit.regressors = p;
    fit.rank = solved.rank;
    fit.conditionNumber = solved.conditionNumber;
    fit.identifiable = solved.identifiable;
    fit.residuals = bilinearResiduals(model, data, &fit.stats);
    return fit;
}

BilinearFit fitBilinearKromFrozenA(const BilinearDataset& data, const std::vector<double>& frozenA,
    std::size_t n, std::size_t m, const BilinearFitOptions& options) {
    validateData(data, n, m);
    if (frozenA.size() != n * n) throw std::invalid_argument("Frozen A has the wrong size");
    requireFinite(frozenA, "frozen A");
    if (!(options.ridge >= 0.0) || !(options.rankTolerance > 0.0))
        throw std::invalid_argument("Invalid bilinear fit options");
    const std::size_t p = m + m * n;
    if (data.samples < p && !options.allowRankDeficient)
        throw std::invalid_argument("Frozen-A bilinear fit needs at least as many pairs as regressors");

    const ConstRowMap Z(data.z.data(), idx(data.samples), idx(n));
    const ConstRowMap U(data.u.data(), idx(data.samples), idx(m));
    const ConstRowMap Y(data.zNext.data(), idx(data.samples), idx(n));
    const ConstRowMap A(frozenA.data(), idx(n), idx(n));
    const MatrixXd target = Y - Z * A.transpose();
    MatrixXd design(idx(data.samples), idx(p));
    design.leftCols(idx(m)) = U;
    for (std::size_t i = 0; i < m; ++i)
        design.middleCols(idx(m + i * n), idx(n)) = U.col(idx(i)).asDiagonal() * Z;
    const auto solved = solveColumnScaled(std::move(design), target, options);
    const MatrixXd& coefficients = solved.coefficients;

    BilinearFit fit;
    auto& model = fit.model;
    model.stateDimension = n; model.inputDimension = m;
    model.sampleTime = options.sampleTime;
    model.provenance = options.provenance.empty() ? "bilinear_fit frozen-A constant-input least squares" : options.provenance;
    model.A = frozenA;  // copied, never recomputed
    model.B0 = toVector(coefficients.topRows(idx(m)).transpose());
    for (std::size_t i = 0; i < m; ++i)
        model.B.push_back(toVector(coefficients.middleRows(idx(m + i * n), idx(n)).transpose()));
    validateBilinearKrom(model);
    fit.regressors = p;
    fit.rank = solved.rank;
    fit.conditionNumber = solved.conditionNumber;
    fit.identifiable = solved.identifiable;
    fit.residuals = bilinearResiduals(model, data, &fit.stats);
    return fit;
}

std::vector<double> stepBilinearKrom(const BilinearKrom& model, const std::vector<double>& z,
    const std::vector<double>& u) {
    validateBilinearKrom(model);
    const auto n = model.stateDimension, m = model.inputDimension;
    if (z.size() != n || u.size() != m) throw std::invalid_argument("Bilinear step has wrong z or u size");
    const Eigen::Map<const VectorXd> zv(z.data(), idx(n)), uv(u.data(), idx(m));
    VectorXd next = ConstRowMap(model.A.data(), idx(n), idx(n)) * zv + ConstRowMap(model.B0.data(), idx(n), idx(m)) * uv;
    for (std::size_t i = 0; i < m; ++i) next += u[i] * (ConstRowMap(model.B[i].data(), idx(n), idx(n)) * zv);
    return std::vector<double>(next.data(), next.data() + next.size());
}

std::vector<double> simulateBilinearKrom(const BilinearKrom& model, const std::vector<double>& z0,
    const std::vector<double>& inputs, std::size_t steps) {
    validateBilinearKrom(model);
    const auto n = model.stateDimension, m = model.inputDimension;
    if (z0.size() != n || inputs.size() != steps * m) throw std::invalid_argument("Bilinear simulation sizes are inconsistent");
    std::vector<double> out(z0);
    out.reserve((steps + 1) * n);
    std::vector<double> z(z0), u(m);
    for (std::size_t k = 0; k < steps; ++k) {
        std::copy_n(inputs.begin() + static_cast<std::ptrdiff_t>(k * m), m, u.begin());
        z = stepBilinearKrom(model, z, u);
        out.insert(out.end(), z.begin(), z.end());
    }
    return out;
}

FrozenInputLinearization linearizeBilinearKrom(const BilinearKrom& model, const std::vector<double>& z,
    const std::vector<double>& u) {
    validateBilinearKrom(model);
    const auto n = model.stateDimension, m = model.inputDimension;
    if (z.size() != n || u.size() != m) throw std::invalid_argument("Linearization point has wrong size");
    const Eigen::Map<const VectorXd> zv(z.data(), idx(n)), uv(u.data(), idx(m));
    RowMatrix frozen = ConstRowMap(model.A.data(), idx(n), idx(n));
    RowMatrix input = ConstRowMap(model.B0.data(), idx(n), idx(m));
    for (std::size_t i = 0; i < m; ++i) {
        const ConstRowMap Bi(model.B[i].data(), idx(n), idx(n));
        frozen += u[i] * Bi;
        input.col(idx(i)) += Bi * zv;
    }
    const auto next = stepBilinearKrom(model, z, u);
    const VectorXd offset = Eigen::Map<const VectorXd>(next.data(), idx(n)) - frozen * zv - input * uv;
    return {toVector(frozen), toVector(input), std::vector<double>(offset.data(), offset.data() + offset.size())};
}

std::vector<double> bilinearResiduals(const BilinearKrom& model, const BilinearDataset& data,
    BilinearResidualStats* stats) {
    validateBilinearKrom(model);
    const auto n = model.stateDimension, m = model.inputDimension;
    validateData(data, n, m);
    std::vector<double> residual(data.samples * n), z(n), u(m);
    double total = 0.0, reference = 0.0, maximum = 0.0;
    std::vector<double> norms(data.samples);
    for (std::size_t k = 0; k < data.samples; ++k) {
        std::copy_n(data.z.begin() + static_cast<std::ptrdiff_t>(k * n), n, z.begin());
        std::copy_n(data.u.begin() + static_cast<std::ptrdiff_t>(k * m), m, u.begin());
        const auto prediction = stepBilinearKrom(model, z, u);
        double squared = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double target = data.zNext[k * n + i];
            const double r = target - prediction[i];
            residual[k * n + i] = r;
            squared += r * r;
            reference += target * target;
        }
        total += squared;
        norms[k] = std::sqrt(squared);
        maximum = std::max(maximum, norms[k]);
    }
    if (stats) {
        stats->norms = std::move(norms);
        stats->rootMeanSquare = std::sqrt(total / double(data.samples));
        stats->maximum = maximum;
        stats->relative = reference > 0.0 ? std::sqrt(total / reference) : std::sqrt(total);
    }
    return residual;
}

json::Value denseRowsToJson(const std::vector<double>& data, std::size_t rows, std::size_t columns) {
    if (data.size() != rows * columns) throw std::invalid_argument("Matrix size mismatch");
    json::Value::Array out;
    out.reserve(rows);
    for (std::size_t i = 0; i < rows; ++i) {
        json::Value::Array row;
        row.reserve(columns);
        for (std::size_t j = 0; j < columns; ++j) row.emplace_back(data[i * columns + j]);
        out.emplace_back(std::move(row));
    }
    return out;
}

std::vector<double> denseRowsFromJson(const json::Value& value, std::size_t& rows, std::size_t& columns) {
    const auto& array = value.array();
    rows = array.size();
    columns = rows ? array.front().array().size() : 0;
    std::vector<double> out;
    out.reserve(rows * columns);
    for (const auto& row : array) {
        const auto& entries = row.array();
        if (entries.size() != columns) throw std::invalid_argument("Matrix rows have unequal length");
        for (const auto& entry : entries) {
            if (!std::holds_alternative<double>(entry.storage())) throw std::invalid_argument("Matrix entry must be numeric");
            out.push_back(entry.number());
        }
    }
    requireFinite(out, "matrix");
    return out;
}

json::Value bilinearKromToJson(const BilinearKrom& model) {
    validateBilinearKrom(model);
    const auto n = model.stateDimension, m = model.inputDimension;
    json::Value::Array operators;
    for (const auto& matrix : model.B) operators.push_back(denseRowsToJson(matrix, n, n));
    return json::Value::Object{
        {"schema", "dynsyslab.bilinear-krom/v1"},
        {"convention", "z_plus = A z + B0 u + sum_i u_i B_i z; z(target) = 0, u(target) = 0"},
        {"state_dimension", double(n)}, {"input_dimension", double(m)},
        {"sample_time", model.sampleTime},
        {"A", denseRowsToJson(model.A, n, n)}, {"B0", denseRowsToJson(model.B0, n, m)},
        {"B", std::move(operators)}, {"provenance", model.provenance}};
}

BilinearKrom bilinearKromFromJson(const json::Value& value) {
    if (value.at("schema").string() != "dynsyslab.bilinear-krom/v1")
        throw std::invalid_argument("Expected dynsyslab.bilinear-krom/v1");
    const double nd = value.at("state_dimension").number(-1.0), md = value.at("input_dimension").number(-1.0);
    if (nd < 1.0 || md < 1.0 || nd != std::floor(nd) || md != std::floor(md) || nd > 1e6 || md > 1e6)
        throw std::invalid_argument("Bilinear model dimensions must be positive integers");
    BilinearKrom model;
    model.stateDimension = static_cast<std::size_t>(nd);
    model.inputDimension = static_cast<std::size_t>(md);
    model.sampleTime = value.at("sample_time").number(0.0);
    model.provenance = value.at("provenance").string();
    std::size_t rows = 0, columns = 0;
    model.A = denseRowsFromJson(value.at("A"), rows, columns);
    if (rows != model.stateDimension || columns != model.stateDimension) throw std::invalid_argument("A has the wrong shape");
    model.B0 = denseRowsFromJson(value.at("B0"), rows, columns);
    if (rows != model.stateDimension || columns != model.inputDimension) throw std::invalid_argument("B0 has the wrong shape");
    for (const auto& matrix : value.at("B").array()) {
        model.B.push_back(denseRowsFromJson(matrix, rows, columns));
        if (rows != model.stateDimension || columns != model.stateDimension) throw std::invalid_argument("B_i has the wrong shape");
    }
    validateBilinearKrom(model);
    return model;
}

}  // namespace dynsyslab::core
