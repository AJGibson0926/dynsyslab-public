#include "EigenNumericsLinearSystems.hpp"

#include <unsupported/Eigen/MatrixFunctions>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace dynsyslab::core {
namespace linear_detail {

Matrix mapped(const std::vector<double>& values, std::size_t rows,
              std::size_t columns, const char* name) {
    if (values.size() != rows * columns)
        throw std::invalid_argument(std::string(name) + " dimensions are inconsistent");
    Matrix result(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(columns));
    if (!values.empty()) result = Eigen::Map<const Matrix>(values.data(),
        static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(columns));
    if (!result.allFinite()) throw std::invalid_argument(std::string(name) + " must be finite");
    return result;
}

std::vector<double> flat(const Matrix& matrix) {
    return {matrix.data(), matrix.data() + matrix.size()};
}

Realization realization(const StateSpaceModel& model) {
    validateStateSpace(model);
    return {mapped(model.A, model.states, model.states, "state-space A"),
            mapped(model.B, model.states, model.inputs, "state-space B"),
            mapped(model.C, model.outputs, model.states, "state-space C"),
            mapped(model.D, model.outputs, model.inputs, "state-space D")};
}

StateSpaceModel fromRealization(const Realization& r, TimeDomain domain, double sampleTime) {
    StateSpaceModel model;
    model.states = static_cast<std::size_t>(r.A.rows());
    model.inputs = static_cast<std::size_t>(r.B.cols());
    model.outputs = static_cast<std::size_t>(r.C.rows());
    model.A = flat(r.A); model.B = flat(r.B); model.C = flat(r.C); model.D = flat(r.D);
    model.domain = domain; model.sampleTime = domain == TimeDomain::discrete ? sampleTime : 0.0;
    return model;
}

std::vector<std::complex<double>> eigenvalues(const Matrix& matrix) {
    if (matrix.rows() == 0) return {};
    Eigen::EigenSolver<Matrix> solver(matrix, false);
    if (solver.info() != Eigen::Success) throw std::runtime_error("eigenvalue solve failed");
    return {solver.eigenvalues().data(), solver.eigenvalues().data() + solver.eigenvalues().size()};
}

bool stable(const Matrix& A, bool discrete, double margin) {
    for (const auto value : eigenvalues(A))
        if (discrete ? std::abs(value) >= 1.0 - margin : value.real() >= -margin) return false;
    return true;
}

Matrix lyapunov(const Matrix& A, const Matrix& Q, bool discrete) {
    const Eigen::Index n = A.rows();
    if (A.cols() != n || Q.rows() != n || Q.cols() != n)
        throw std::invalid_argument("Lyapunov dimensions mismatch");
    if (n == 0) return Matrix(0, 0);
    Eigen::ComplexSchur<CMatrix> schur(A.cast<std::complex<double>>());
    if (schur.info() != Eigen::Success) throw std::runtime_error("Lyapunov Schur decomposition failed");
    const CMatrix& T = schur.matrixT();
    const CMatrix& U = schur.matrixU();
    const CMatrix F = U.adjoint() * Q.cast<std::complex<double>>() * U;
    CMatrix Y = CMatrix::Zero(n, n);
    const double scale = std::max(1.0, T.cwiseAbs().maxCoeff());
    for (Eigen::Index j = n - 1; j >= 0; --j) {
        Eigen::VectorXcd accumulated = Eigen::VectorXcd::Zero(n);
        for (Eigen::Index k = j + 1; k < n; ++k) accumulated += std::conj(T(j, k)) * Y.col(k);
        CMatrix op; Eigen::VectorXcd right;
        if (discrete) {
            op = std::conj(T(j, j)) * T - CMatrix::Identity(n, n);
            right = -F.col(j) - T * accumulated;
        } else {
            op = T + std::conj(T(j, j)) * CMatrix::Identity(n, n);
            right = -F.col(j) - accumulated;
        }
        const double pivot = op.diagonal().cwiseAbs().minCoeff();
        if (!(pivot > 1.0e-14 * scale * (discrete ? scale : 1.0)))
            throw std::invalid_argument(discrete
                ? "Stein equation is singular (eigenvalue pair with mu_i conj(mu_j) = 1)"
                : "Lyapunov equation is singular (eigenvalue pair with lambda_i + conj(lambda_j) = 0)");
        Y.col(j) = op.triangularView<Eigen::Upper>().solve(right);
    }
    Matrix X = (U * Y * U.adjoint()).real();
    X = 0.5 * (X + X.transpose()).eval();
    const Matrix residual = discrete ? Matrix(A * X * A.transpose() - X + Q)
                                     : Matrix(A * X + X * A.transpose() + Q);
    const double size = std::max({1.0e-300, Q.norm(), X.norm() * std::max(1.0, A.norm() * A.norm())});
    if (!X.allFinite() || residual.norm() > 1.0e-8 * size)
        throw std::runtime_error("Lyapunov residual is too large (ill-conditioned equation)");
    return X;
}

Matrix finiteContinuousGramian(const Matrix& A, const Matrix& Q, double horizon) {
    const Eigen::Index n = A.rows();
    if (!(horizon > 0.0) || !std::isfinite(horizon))
        throw std::invalid_argument("continuous Gramian horizon must be positive");
    // exp([[-A, Q], [0, A^T]] h) = [[F11, F12], [0, F22]] and W(h) = F22^T F12.
    const double size = A.cwiseAbs().rowwise().sum().maxCoeff() * horizon;
    if (!std::isfinite(size)) throw std::invalid_argument("Gramian time scale overflow");
    const int squarings = size > 0.5 ? static_cast<int>(std::ceil(std::log2(size / 0.5))) : 0;
    const double step = std::ldexp(horizon, -squarings);
    Matrix block = Matrix::Zero(2 * n, 2 * n);
    block.topLeftCorner(n, n) = -A; block.topRightCorner(n, n) = Q;
    block.bottomRightCorner(n, n) = A.transpose();
    const Matrix E = (block * step).exp();
    Matrix flow = (A * step).exp();
    Matrix W = flow * E.topRightCorner(n, n);
    for (int i = 0; i < squarings; ++i) {
        W = (W + flow * W * flow.transpose()).eval();
        flow = (flow * flow).eval();
    }
    if (!W.allFinite()) throw std::runtime_error("finite continuous Gramian is not finite");
    return 0.5 * (W + W.transpose());
}

Matrix psdFactor(const Matrix& W) {
    if (W.rows() == 0) return Matrix(0, 0);
    Eigen::SelfAdjointEigenSolver<Matrix> solver(0.5 * (W + W.transpose()));
    if (solver.info() != Eigen::Success) throw std::runtime_error("gramian factorization failed");
    const Vector root = solver.eigenvalues().cwiseMax(0.0).cwiseSqrt();
    return solver.eigenvectors() * root.asDiagonal();
}

Matrix orthonormalBasis(const Matrix& M, double relativeTolerance) {
    if (M.rows() == 0 || M.cols() == 0) return Matrix(M.rows(), 0);
    Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeThinU);
    const auto& s = svd.singularValues();
    const double threshold = relativeTolerance * std::max(s.size() ? s[0] : 0.0,
                                                          std::numeric_limits<double>::min());
    Eigen::Index rank = 0;
    while (rank < s.size() && s[rank] > threshold) ++rank;
    if (s.size() && s[0] <= std::numeric_limits<double>::min()) rank = 0;
    return svd.matrixU().leftCols(rank);
}

Matrix orthogonalComplement(const Matrix& Q, Eigen::Index dimension) {
    if (Q.cols() == 0) return Matrix::Identity(dimension, dimension);
    if (Q.cols() >= dimension) return Matrix(dimension, 0);
    Eigen::JacobiSVD<Matrix> svd(Q, Eigen::ComputeFullU);
    return svd.matrixU().rightCols(dimension - Q.cols());
}

}  // namespace linear_detail

using namespace linear_detail;

const char* timeDomainName(TimeDomain domain) noexcept {
    return domain == TimeDomain::discrete ? "discrete" : "continuous";
}

void validateStateSpace(const StateSpaceModel& m) {
    if (!m.states) throw std::invalid_argument("state-space model needs at least one state");
    if (m.A.size() != m.states * m.states || m.B.size() != m.states * m.inputs
        || m.C.size() != m.outputs * m.states || m.D.size() != m.outputs * m.inputs)
        throw std::invalid_argument("state-space matrix dimensions are inconsistent");
    if (m.domain == TimeDomain::discrete && !(m.sampleTime > 0.0 && std::isfinite(m.sampleTime)))
        throw std::invalid_argument("discrete state-space model needs a positive sample time");
    for (const auto* values : {&m.A, &m.B, &m.C, &m.D})
        for (double value : *values)
            if (!std::isfinite(value)) throw std::invalid_argument("state-space matrices must be finite");
}

StateSpaceModel makeStateSpace(std::vector<double> A, std::vector<double> B,
    std::vector<double> C, std::vector<double> D, std::size_t states,
    std::size_t inputs, std::size_t outputs, TimeDomain domain, double sampleTime) {
    StateSpaceModel m{std::move(A), std::move(B), std::move(C), std::move(D),
                      states, inputs, outputs, domain, sampleTime};
    if (m.D.empty()) m.D.assign(outputs * inputs, 0.0);
    if (domain == TimeDomain::continuous) m.sampleTime = 0.0;
    validateStateSpace(m);
    return m;
}

std::vector<std::complex<double>> systemPoles(const StateSpaceModel& model) {
    return eigenvalues(realization(model).A);
}

bool isAsymptoticallyStable(const StateSpaceModel& model, double margin) {
    return stable(realization(model).A, model.domain == TimeDomain::discrete, margin);
}

std::vector<double> solveContinuousLyapunov(const std::vector<double>& A,
    const std::vector<double>& Q, std::size_t n) {
    return flat(lyapunov(mapped(A, n, n, "Lyapunov A"), mapped(Q, n, n, "Lyapunov Q"), false));
}

std::vector<double> solveDiscreteLyapunov(const std::vector<double>& A,
    const std::vector<double>& Q, std::size_t n) {
    return flat(lyapunov(mapped(A, n, n, "Stein A"), mapped(Q, n, n, "Stein Q"), true));
}

namespace {
Matrix infiniteGramian(const Matrix& A, const Matrix& BB, bool discrete) {
    if (!stable(A, discrete))
        throw std::invalid_argument("infinite-horizon gramian requires an asymptotically stable model");
    return lyapunov(A, BB, discrete);
}

Matrix finiteGramian(const Matrix& A, const Matrix& BB, bool discrete, double horizon) {
    const Eigen::Index n = A.rows();
    if (!(horizon > 0.0) || !std::isfinite(horizon))
        throw std::invalid_argument("finite gramian horizon must be positive");
    if (discrete) {
        if (std::floor(horizon) != horizon || horizon > 1.0e7)
            throw std::invalid_argument("discrete finite gramian horizon must be a bounded step count");
        Matrix W = Matrix::Zero(n, n), power = Matrix::Identity(n, n);
        for (long k = 0; k < static_cast<long>(horizon); ++k) {
            W += power * BB * power.transpose();
            power = (A * power).eval();
            if (!W.allFinite()) throw std::runtime_error("finite discrete gramian overflowed");
        }
        return 0.5 * (W + W.transpose());
    }
    return finiteContinuousGramian(A, BB, horizon);
}
}  // namespace

std::vector<double> controllabilityGramian(const StateSpaceModel& model) {
    const auto r = realization(model);
    return flat(infiniteGramian(r.A, r.B * r.B.transpose(), model.domain == TimeDomain::discrete));
}

std::vector<double> observabilityGramian(const StateSpaceModel& model) {
    const auto r = realization(model);
    return flat(infiniteGramian(r.A.transpose(), r.C.transpose() * r.C,
                                model.domain == TimeDomain::discrete));
}

std::vector<double> finiteControllabilityGramian(const StateSpaceModel& model, double horizon) {
    const auto r = realization(model);
    return flat(finiteGramian(r.A, r.B * r.B.transpose(), model.domain == TimeDomain::discrete, horizon));
}

std::vector<double> finiteObservabilityGramian(const StateSpaceModel& model, double horizon) {
    const auto r = realization(model);
    return flat(finiteGramian(r.A.transpose(), r.C.transpose() * r.C,
                              model.domain == TimeDomain::discrete, horizon));
}

StateSpaceModel parallelConnection(const StateSpaceModel& a, const StateSpaceModel& b, double sign) {
    if (a.inputs != b.inputs || a.outputs != b.outputs || a.domain != b.domain
        || (a.domain == TimeDomain::discrete && std::abs(a.sampleTime - b.sampleTime) > 1e-12 * a.sampleTime))
        throw std::invalid_argument("parallel connection needs matching inputs, outputs and time domain");
    const std::size_t n1 = a.states, n2 = b.states, n = n1 + n2, m = a.inputs, p = a.outputs;
    if (a.D.size() != p * m || b.D.size() != p * m) throw std::invalid_argument("parallel connection D mismatch");
    StateSpaceModel out; out.states = n; out.inputs = m; out.outputs = p;
    out.domain = a.domain; out.sampleTime = a.sampleTime;
    out.A.assign(n * n, 0.0); out.B.assign(n * m, 0.0); out.C.assign(p * n, 0.0); out.D.assign(p * m, 0.0);
    for (std::size_t i = 0; i < n1; ++i) {
        for (std::size_t j = 0; j < n1; ++j) out.A[i * n + j] = a.A[i * n1 + j];
        for (std::size_t j = 0; j < m; ++j) out.B[i * m + j] = a.B[i * m + j];
    }
    for (std::size_t i = 0; i < n2; ++i) {
        for (std::size_t j = 0; j < n2; ++j) out.A[(n1 + i) * n + n1 + j] = b.A[i * n2 + j];
        for (std::size_t j = 0; j < m; ++j) out.B[(n1 + i) * m + j] = b.B[i * m + j];
    }
    for (std::size_t i = 0; i < p; ++i) {
        for (std::size_t j = 0; j < n1; ++j) out.C[i * n + j] = a.C[i * n1 + j];
        for (std::size_t j = 0; j < n2; ++j) out.C[i * n + n1 + j] = sign * b.C[i * n2 + j];
        for (std::size_t j = 0; j < m; ++j) out.D[i * m + j] = a.D[i * m + j] + sign * b.D[i * m + j];
    }
    if (n) validateStateSpace(out);
    return out;
}

StateSpaceModel similarityTransform(const StateSpaceModel& model,
    const std::vector<double>& transform, const std::vector<double>& inverse) {
    const auto r = realization(model);
    const Matrix T = mapped(transform, model.states, model.states, "similarity transform");
    const Matrix Ti = mapped(inverse, model.states, model.states, "similarity inverse");
    if ((Ti * T - Matrix::Identity(T.rows(), T.cols())).norm() > 1.0e-8 * std::sqrt(double(model.states)))
        throw std::invalid_argument("similarity inverse does not invert the transform");
    return fromRealization({Ti * r.A * T, Ti * r.B, r.C * T, r.D}, model.domain, model.sampleTime);
}

}  // namespace dynsyslab::core
