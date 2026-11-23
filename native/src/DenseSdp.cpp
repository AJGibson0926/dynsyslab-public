// Dense homogeneous self-dual splitting for LMIs. Reference: O'Donoghue, Chu,
// Parikh & Boyd, J. Optim. Theory Appl. 169 (2016) 1042-1068, Algorithm (SCS):
//   u_tilde = (I + Q)^-1 (u + v)
//   u       = Pi_C(alpha u_tilde + (1 - alpha) u - v)
//   v       = v - alpha u_tilde - (1 - alpha) u_old + u
// with Q = [0 A' c; -A 0 b; -c' -b' 0] and C = R^n x K x R_+.
// PSD blocks use the isometric lower-triangular vectorization (off-diagonals
// scaled by sqrt 2), so the cone stays self-dual in the Euclidean product.
#include <dynsyslab/core/DenseSdp.hpp>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using RowMap = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

const double kRootTwo = std::sqrt(2.0);

std::size_t triangle(std::size_t d) { return d * (d + 1) / 2; }

void validate(const SdpProblem& p) {
    if (!p.objective.empty() && p.objective.size() != p.variables)
        throw std::invalid_argument("SDP objective size differs from the variable count");
    for (double v : p.objective)
        if (!std::isfinite(v)) throw std::invalid_argument("SDP objective must be finite");
    if (p.blocks.empty()) throw std::invalid_argument("SDP needs at least one LMI block");
    const auto checkMatrix = [](const std::vector<double>& m, std::size_t d) {
        if (m.size() != d * d) throw std::invalid_argument("LMI matrix has the wrong size");
        double scale = 0.0;
        for (double v : m) {
            if (!std::isfinite(v)) throw std::invalid_argument("LMI matrix must be finite");
            scale = std::max(scale, std::abs(v));
        }
        for (std::size_t i = 0; i < d; ++i)
            for (std::size_t j = 0; j < i; ++j)
                if (std::abs(m[i * d + j] - m[j * d + i]) > 1e-10 * (1.0 + scale))
                    throw std::invalid_argument("LMI matrix must be symmetric");
    };
    for (const auto& block : p.blocks) {
        if (block.dimension == 0) throw std::invalid_argument("LMI block dimension must be positive");
        if (block.coefficients.size() != p.variables)
            throw std::invalid_argument("LMI block needs one coefficient per variable");
        checkMatrix(block.constant, block.dimension);
        for (const auto& f : block.coefficients)
            if (!f.empty()) checkMatrix(f, block.dimension);
    }
}

void packInto(const std::vector<double>& m, std::size_t d, double sign, Eigen::Ref<VectorXd> out) {
    Eigen::Index k = 0;
    for (std::size_t j = 0; j < d; ++j)
        for (std::size_t i = j; i < d; ++i, ++k)
            out(k) = sign * m[i * d + j] * (i == j ? 1.0 : kRootTwo);
}

void projectPsd(Eigen::Ref<VectorXd> v, std::size_t d, MatrixXd& work) {
    if (d == 1) { v(0) = std::max(v(0), 0.0); return; }
    const auto n = static_cast<Eigen::Index>(d);
    work.resize(n, n);
    Eigen::Index k = 0;
    for (Eigen::Index j = 0; j < n; ++j)
        for (Eigen::Index i = j; i < n; ++i, ++k) {
            const double value = i == j ? v(k) : v(k) / kRootTwo;
            work(i, j) = value; work(j, i) = value;
        }
    Eigen::SelfAdjointEigenSolver<MatrixXd> eig(work);
    if (eig.info() != Eigen::Success) throw std::runtime_error("PSD projection failed");
    const auto& values = eig.eigenvalues();
    if (values(0) >= 0.0) return;
    if (values(n - 1) <= 0.0) { v.setZero(); return; }
    Eigen::Index first = 0;
    while (first < n && values(first) < 0.0) ++first;
    const auto count = n - first;
    const MatrixXd vectors = eig.eigenvectors().rightCols(count);
    const MatrixXd projected = vectors * values.tail(count).asDiagonal() * vectors.transpose();
    k = 0;
    for (Eigen::Index j = 0; j < n; ++j)
        for (Eigen::Index i = j; i < n; ++i, ++k)
            v(k) = i == j ? projected(i, j) : projected(i, j) * kRootTwo;
}

struct ConicForm {
    MatrixXd A;
    VectorXd b, c;
    std::vector<Eigen::Index> offsets;
    std::vector<std::size_t> dims;
};

ConicForm conicForm(const SdpProblem& p) {
    ConicForm f;
    Eigen::Index rows = 0;
    for (const auto& block : p.blocks) {
        f.offsets.push_back(rows);
        f.dims.push_back(block.dimension);
        rows += static_cast<Eigen::Index>(triangle(block.dimension));
    }
    const auto n = static_cast<Eigen::Index>(p.variables);
    f.A = MatrixXd::Zero(rows, n);
    f.b.resize(rows);
    f.c = p.objective.empty() ? VectorXd(VectorXd::Zero(n))
                              : VectorXd(Eigen::Map<const VectorXd>(p.objective.data(), n));
    for (std::size_t b = 0; b < p.blocks.size(); ++b) {
        const auto& block = p.blocks[b];
        const auto size = static_cast<Eigen::Index>(triangle(block.dimension));
        packInto(block.constant, block.dimension, 1.0, f.b.segment(f.offsets[b], size));
        for (Eigen::Index i = 0; i < n; ++i) {
            const auto& coefficient = block.coefficients[static_cast<std::size_t>(i)];
            if (!coefficient.empty())
                packInto(coefficient, block.dimension, -1.0, f.A.col(i).segment(f.offsets[b], size));
        }
    }
    return f;
}

// Ruiz-style equilibration. Rows share one positive scalar per PSD block, so
// the scaled slack stays in the same cone.
struct Scaling {
    VectorXd column, row;
    double b{1.0}, c{1.0};
};

Scaling equilibrate(const ConicForm& f, MatrixXd& A, VectorXd& b, VectorXd& c, bool enabled) {
    Scaling s;
    A = f.A;
    s.column = VectorXd::Ones(A.cols());
    s.row = VectorXd::Ones(A.rows());
    if (enabled) {
        for (int pass = 0; pass < 10; ++pass) {
            for (Eigen::Index j = 0; j < A.cols(); ++j) {
                const double norm = A.col(j).norm();
                if (norm > 1e-12) {
                    const double factor = 1.0 / std::sqrt(norm);
                    A.col(j) *= factor; s.column(j) *= factor;
                }
            }
            for (std::size_t block = 0; block < f.dims.size(); ++block) {
                const auto size = static_cast<Eigen::Index>(triangle(f.dims[block]));
                const double norm = A.middleRows(f.offsets[block], size).norm() / std::sqrt(double(size));
                if (norm > 1e-12) {
                    const double factor = 1.0 / std::sqrt(norm);
                    A.middleRows(f.offsets[block], size) *= factor;
                    s.row.segment(f.offsets[block], size).array() *= factor;
                }
            }
        }
    }
    b = s.row.cwiseProduct(f.b);
    c = s.column.cwiseProduct(f.c);
    s.b = 1.0 / std::max(1.0, b.norm());
    s.c = 1.0 / std::max(1.0, c.norm());
    b *= s.b; c *= s.c;
    return s;
}

std::vector<double> unpackBlock(const VectorXd& y, Eigen::Index offset, std::size_t d) {
    std::vector<double> m(d * d);
    Eigen::Index k = offset;
    for (std::size_t j = 0; j < d; ++j)
        for (std::size_t i = j; i < d; ++i, ++k) {
            const double value = i == j ? y(k) : y(k) / kRootTwo;
            m[i * d + j] = value; m[j * d + i] = value;
        }
    return m;
}

}  // namespace

const char* sdpStatusName(SdpStatus status) {
    switch (status) {
        case SdpStatus::converged: return "converged";
        case SdpStatus::primalInfeasible: return "primal_infeasible";
        case SdpStatus::dualInfeasible: return "dual_infeasible";
        default: return "iteration_limit";
    }
}

SdpResult solveDenseSdp(const SdpProblem& problem, const SdpOptions& options) {
    validate(problem);
    if (problem.variables == 0) throw std::invalid_argument("SDP needs at least one variable");
    if (!(options.tolerance > 0.0) || !(options.infeasibilityTolerance > 0.0) ||
        !(options.relaxation > 0.0 && options.relaxation < 2.0) || options.checkInterval == 0)
        throw std::invalid_argument("Invalid SDP options");

    const ConicForm form = conicForm(problem);
    MatrixXd A; VectorXd b, c;
    const Scaling scale = equilibrate(form, A, b, c, options.equilibrate);
    const auto n = A.cols();
    const auto m = A.rows();
    const Eigen::LLT<MatrixXd> factor(MatrixXd::Identity(n, n) + A.transpose() * A);
    if (factor.info() != Eigen::Success) throw std::runtime_error("SDP KKT factorization failed");

    // (x, y) = M^-1 (p, q), with M = [I A'; -A I].
    const auto solveM = [&](const VectorXd& p, const VectorXd& q, VectorXd& x, VectorXd& y) {
        x = factor.solve(p - A.transpose() * q);
        y = q + A * x;
    };
    VectorXd xh, yh;
    solveM(c, b, xh, yh);
    const double denominator = 1.0 + c.dot(xh) + b.dot(yh);

    VectorXd ux = VectorXd::Zero(n), uy = VectorXd::Zero(m), vy = VectorXd::Zero(m);
    double ut = 1.0, vt = 1.0;
    VectorXd xp, yp, wy(m);
    MatrixXd work;
    const double alpha = options.relaxation;
    const double bNorm = form.b.norm(), cNorm = form.c.norm();

    SdpResult result;
    const auto recordIterate = [&](const VectorXd& x, const VectorXd& y) {
        result.x.assign(x.data(), x.data() + x.size());
        result.dual.clear();
        for (std::size_t block = 0; block < form.dims.size(); ++block)
            result.dual.push_back(unpackBlock(y, form.offsets[block], form.dims[block]));
    };

    for (std::size_t it = 1; it <= options.maxIterations; ++it) {
        solveM(ux, uy + vy, xp, yp);
        const double tt = (ut + vt + c.dot(xp) + b.dot(yp)) / denominator;
        const VectorXd xt = xp - tt * xh;
        const VectorXd yt = yp - tt * yh;

        ux = alpha * xt + (1.0 - alpha) * ux;
        wy = alpha * yt + (1.0 - alpha) * uy - vy;
        const double wt = alpha * tt + (1.0 - alpha) * ut - vt;
        uy = wy;
        for (std::size_t block = 0; block < form.dims.size(); ++block) {
            const auto size = static_cast<Eigen::Index>(triangle(form.dims[block]));
            projectPsd(uy.segment(form.offsets[block], size), form.dims[block], work);
        }
        ut = std::max(wt, 0.0);
        vy = uy - wy;
        vt = ut - wt;

        if (it % options.checkInterval != 0 && it != options.maxIterations) continue;
        if (!ux.allFinite() || !uy.allFinite() || !std::isfinite(ut))
            throw std::runtime_error("SDP iteration diverged to a non-finite value");
        result.iterations = it;

        if (ut > 1e-12) {
            const VectorXd x = scale.column.cwiseProduct(ux) / (scale.b * ut);
            const VectorXd y = scale.row.cwiseProduct(uy) / (scale.c * ut);
            const VectorXd s = vy.cwiseQuotient(scale.row) / (scale.b * ut);
            const double cx = form.c.dot(x), by = form.b.dot(y);
            result.primalResidual = (form.A * x + s - form.b).norm() / (1.0 + bNorm);
            result.dualResidual = (form.A.transpose() * y + form.c).norm() / (1.0 + cNorm);
            result.gap = std::abs(cx + by) / (1.0 + std::abs(cx) + std::abs(by));
            result.objective = cx;
            if (result.primalResidual <= options.tolerance && result.dualResidual <= options.tolerance &&
                result.gap <= options.tolerance) {
                result.status = SdpStatus::converged;
                recordIterate(x, y);
                return result;
            }
            if (it == options.maxIterations) recordIterate(x, y);
        }
        if (ut < vt) {  // The embedding is drifting to a certificate of infeasibility.
            const VectorXd y = scale.row.cwiseProduct(uy);
            const double by = form.b.dot(y);
            if (by < 0.0 && (form.A.transpose() * y).norm() / -by <= options.infeasibilityTolerance) {
                result.status = SdpStatus::primalInfeasible;
                result.x.clear();
                return result;
            }
            const VectorXd x = scale.column.cwiseProduct(ux);
            const VectorXd s = vy.cwiseQuotient(scale.row);
            const double cx = form.c.dot(x);
            if (cx < 0.0 && (form.A * x + s).norm() / -cx <= options.infeasibilityTolerance) {
                result.status = SdpStatus::dualInfeasible;
                result.x.clear();
                return result;
            }
        }
    }
    result.status = SdpStatus::iterationLimit;
    return result;
}

LmiCertificate verifyLmiCertificate(const SdpProblem& problem, const std::vector<double>& x, double margin) {
    LmiCertificate out;
    out.margin = margin;
    out.worstMinEigenvalue = std::numeric_limits<double>::infinity();
    out.scope = "floating-point eigenvalue check with a rounding allowance; not interval arithmetic";
    if (!std::isfinite(margin)) throw std::invalid_argument("Certificate margin must be finite");
    if (x.size() != problem.variables) { out.reason = "solution size differs from the variable count"; return out; }
    if (problem.blocks.empty()) { out.reason = "no LMI blocks to verify"; return out; }
    for (double v : x)
        if (!std::isfinite(v)) { out.reason = "solution is not finite"; return out; }

    const double unit = std::numeric_limits<double>::epsilon();
    bool passed = true;
    for (const auto& block : problem.blocks) {
        const auto d = static_cast<Eigen::Index>(block.dimension);
        if (d == 0 || block.constant.size() != block.dimension * block.dimension ||
            block.coefficients.size() != problem.variables) {
            out.reason = "malformed LMI block"; return out;
        }
        MatrixXd F = RowMap(block.constant.data(), d, d);
        double accumulated = F.norm();
        for (std::size_t i = 0; i < problem.variables; ++i) {
            if (block.coefficients[i].empty()) continue;
            if (block.coefficients[i].size() != block.dimension * block.dimension) {
                out.reason = "malformed LMI coefficient"; return out;
            }
            const MatrixXd term = x[i] * RowMap(block.coefficients[i].data(), d, d);
            accumulated += term.norm();
            F += term;
        }
        if (!F.allFinite()) { out.reason = "reconstructed LMI is not finite"; return out; }
        if ((F - F.transpose()).norm() > 1e-9 * (1.0 + accumulated)) {
            out.reason = "reconstructed LMI is not symmetric"; return out;
        }
        const MatrixXd symmetric = 0.5 * (F + F.transpose());
        Eigen::SelfAdjointEigenSolver<MatrixXd> eig(symmetric, Eigen::EigenvaluesOnly);
        if (eig.info() != Eigen::Success) { out.reason = "eigenvalue check failed"; return out; }
        const double minimum = eig.eigenvalues()(0);
        // Covers assembly and backward-stable eigensolver error, to a modest constant.
        const double allowance = 16.0 * double(d) * unit * accumulated;
        out.blockMinEigenvalues.push_back(minimum);
        out.worstMinEigenvalue = std::min(out.worstMinEigenvalue, minimum);
        out.roundingAllowance = std::max(out.roundingAllowance, allowance);
        if (minimum - allowance < margin) passed = false;
    }
    out.certified = passed;
    if (!passed) out.reason = "minimum eigenvalue is below the required margin";
    return out;
}

}  // namespace dynsyslab::core
