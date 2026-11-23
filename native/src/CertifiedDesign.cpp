// Robust LMI design for a bilinear K-ROM; see CertifiedDesign.hpp for the
// inequality and citations (Bernussou, Peres & Geromel 1989; Boyd et al. 1994;
// Strasser, Worthmann et al., arXiv:2402.03145).
//
// Derivation. With r'r <= z'Qz, Q = g_x^2 I + g_u^2 K'K (Young's inequality on the
// sector, g = sqrt(2) c when both coefficients are nonzero), the S-procedure with
// multiplier tau = 1/mu asks that (A_cl z + r)'P(A_cl z + r) <= lambda z'Pz. Two
// Schur complements and the congruence diag(X, I) give the linear block below.
// Vertices suffice because the block is affine in theta and u = K z stays in the box on E.
//
// The state is normalized by R (X_n = X/R^2, L_n = L/R, B0_n = B0/R, g_u,n = g_u/R),
// which multiplies the main block by 1/R^2, so margins are scale-free.
#include <dynsyslab/core/CertifiedDesign.hpp>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Index = Eigen::Index;
using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ConstRowMap = Eigen::Map<const RowMatrix>;

constexpr std::size_t kMaximumInputsForVertices = 12;

struct Normalized {
    Index n{}, m{};
    MatrixXd A, B0;
    std::vector<MatrixXd> B;
    double gx{}, gu{}, lambda{};
    VectorXd bounds;
    std::vector<VectorXd> vertices;
};

void validate(const BilinearKrom& model, const ResidualSectorBound& sector, const CertifiedDesignSpec& s) {
    validateBilinearKrom(model);
    if (!sector.feasible) throw std::invalid_argument("Certified design needs a feasible residual sector");
    if (!(sector.cState >= 0.0) || !(sector.cInput >= 0.0) || !std::isfinite(sector.cState) || !std::isfinite(sector.cInput))
        throw std::invalid_argument("Sector coefficients must be finite and non-negative");
    if (s.inputBounds.size() != model.inputDimension) throw std::invalid_argument("Need one input bound per input");
    if (model.inputDimension > kMaximumInputsForVertices) throw std::invalid_argument("Too many inputs for box-vertex enumeration");
    for (double b : s.inputBounds)
        if (!(b > 0.0) || !std::isfinite(b)) throw std::invalid_argument("Input bounds must be positive and finite");
    if (!(s.stateRadius > 0.0) || !std::isfinite(s.stateRadius)) throw std::invalid_argument("State radius must be positive and finite");
    if (!(s.contraction > 0.0 && s.contraction < 1.0)) throw std::invalid_argument("Contraction must lie in (0, 1)");
    if (!(s.certificateMargin > 0.0) || !(s.solverMargin >= s.certificateMargin))
        throw std::invalid_argument("Need solverMargin >= certificateMargin > 0");
}

Normalized normalize(const BilinearKrom& model, const ResidualSectorBound& sector, const CertifiedDesignSpec& s) {
    Normalized o;
    o.n = static_cast<Index>(model.stateDimension);
    o.m = static_cast<Index>(model.inputDimension);
    o.A = ConstRowMap(model.A.data(), o.n, o.n);
    o.B0 = ConstRowMap(model.B0.data(), o.n, o.m) / s.stateRadius;
    for (const auto& b : model.B) o.B.emplace_back(ConstRowMap(b.data(), o.n, o.n));
    const bool both = sector.cState > 0.0 && sector.cInput > 0.0;
    const double young = both ? std::sqrt(2.0) : 1.0;
    o.gx = young * sector.cState;
    o.gu = young * sector.cInput / s.stateRadius;
    o.lambda = s.contraction;
    o.bounds = Eigen::Map<const VectorXd>(s.inputBounds.data(), o.m);
    const std::size_t count = std::size_t{1} << model.inputDimension;
    for (std::size_t mask = 0; mask < count; ++mask) {
        VectorXd theta(o.m);
        for (Index i = 0; i < o.m; ++i) theta(i) = (mask >> i) & 1U ? o.bounds(i) : -o.bounds(i);
        o.vertices.push_back(std::move(theta));
    }
    return o;
}

MatrixXd vertexBlock(const Normalized& o, const VectorXd& theta, const MatrixXd& X, const MatrixXd& L, double mu) {
    const bool useX = o.gx > 0.0, useU = o.gu > 0.0;
    const Index size = 2 * o.n + (useX ? o.n : 0) + (useU ? o.m : 0);
    MatrixXd F = MatrixXd::Zero(size, size);
    MatrixXd Phi = o.A * X + o.B0 * L;
    for (Index i = 0; i < o.m; ++i) Phi += theta(i) * o.B[static_cast<std::size_t>(i)] * X;
    F.topLeftCorner(o.n, o.n) = o.lambda * X;
    F.block(o.n, 0, o.n, o.n) = Phi;
    F.block(0, o.n, o.n, o.n) = Phi.transpose();
    F.block(o.n, o.n, o.n, o.n) = X - mu * MatrixXd::Identity(o.n, o.n);
    Index next = 2 * o.n;
    if (useX) {
        F.block(next, 0, o.n, o.n) = o.gx * X;
        F.block(0, next, o.n, o.n) = o.gx * X.transpose();
        F.block(next, next, o.n, o.n) = mu * MatrixXd::Identity(o.n, o.n);
        next += o.n;
    }
    if (useU) {
        F.block(next, 0, o.m, o.n) = o.gu * L;
        F.block(0, next, o.n, o.m) = o.gu * L.transpose();
        F.block(next, next, o.m, o.m) = mu * MatrixXd::Identity(o.m, o.m);
    }
    return F;
}

// All inequalities that define the certificate (the objective's t block is separate).
std::vector<MatrixXd> certificateBlocks(const Normalized& o, const MatrixXd& X, const MatrixXd& L, double mu) {
    std::vector<MatrixXd> blocks;
    for (const auto& theta : o.vertices) blocks.push_back(vertexBlock(o, theta, X, L, mu));
    for (Index i = 0; i < o.m; ++i) {  // Congruence by diag(1/ubar_i, I) keeps entries O(1).
        MatrixXd F(o.n + 1, o.n + 1);
        F(0, 0) = 1.0;
        F.block(0, 1, 1, o.n) = L.row(i) / o.bounds(i);
        F.block(1, 0, o.n, 1) = L.row(i).transpose() / o.bounds(i);
        F.block(1, 1, o.n, o.n) = X;
        blocks.push_back(std::move(F));
    }
    blocks.emplace_back(MatrixXd::Identity(o.n, o.n) - X);
    blocks.push_back(X);
    blocks.emplace_back(MatrixXd::Constant(1, 1, mu));
    return blocks;
}

struct Layout {
    Index n, m, triangle;
    Index variables() const { return triangle + n * m + 2; }
    Index muIndex() const { return triangle + n * m; }
    Index tIndex() const { return triangle + n * m + 1; }
    void unpack(const VectorXd& x, MatrixXd& X, MatrixXd& L, double& mu, double& t) const {
        X.resize(n, n); L.resize(m, n);
        Index k = 0;
        for (Index j = 0; j < n; ++j)
            for (Index i = 0; i <= j; ++i, ++k) { X(i, j) = x(k); X(j, i) = x(k); }
        for (Index i = 0; i < m; ++i)
            for (Index j = 0; j < n; ++j, ++k) L(i, j) = x(k);
        mu = x(muIndex());
        t = x(tIndex());
    }
};

std::vector<double> rowMajor(const MatrixXd& m) {
    std::vector<double> out(static_cast<std::size_t>(m.size()));
    Eigen::Map<RowMatrix>(out.data(), m.rows(), m.cols()) = m;
    return out;
}

SdpProblem buildProblem(const Normalized& o, const Layout& layout, double margin) {
    const auto assemble = [&](const VectorXd& x) {
        MatrixXd X, L; double mu, t;
        layout.unpack(x, X, L, mu, t);
        auto blocks = certificateBlocks(o, X, L, mu);
        blocks.emplace_back(X - t * MatrixXd::Identity(o.n, o.n));  // X >= t I
        return blocks;
    };
    const Index count = layout.variables();
    const auto base = assemble(VectorXd::Zero(count));
    SdpProblem p;
    p.variables = static_cast<std::size_t>(count);
    p.objective.assign(p.variables, 0.0);
    p.objective[static_cast<std::size_t>(layout.tIndex())] = -1.0;  // maximize the inscribed ball
    p.blocks.resize(base.size());
    for (std::size_t b = 0; b < base.size(); ++b) {
        const auto d = base[b].rows();
        p.blocks[b].dimension = static_cast<std::size_t>(d);
        p.blocks[b].constant = rowMajor(base[b] - margin * MatrixXd::Identity(d, d));
        p.blocks[b].coefficients.resize(p.variables);
    }
    for (Index v = 0; v < count; ++v) {
        VectorXd e = VectorXd::Zero(count);
        e(v) = 1.0;
        const auto value = assemble(e);
        for (std::size_t b = 0; b < base.size(); ++b) {
            const MatrixXd difference = value[b] - base[b];
            if (difference.cwiseAbs().maxCoeff() > 0.0)
                p.blocks[b].coefficients[static_cast<std::size_t>(v)] = rowMajor(difference);
        }
    }
    return p;
}

LmiCertificate checkBlocks(const std::vector<MatrixXd>& blocks, double margin) {
    SdpProblem p;
    p.variables = 0;
    for (const auto& F : blocks) {
        LmiBlock block;
        block.dimension = static_cast<std::size_t>(F.rows());
        block.constant = rowMajor(F);
        p.blocks.push_back(std::move(block));
    }
    return verifyLmiCertificate(p, {}, margin);
}

}  // namespace

const char* certifiedDesignStatusName(CertifiedDesignStatus status) {
    switch (status) {
        case CertifiedDesignStatus::certified: return "certified";
        case CertifiedDesignStatus::infeasible: return "infeasible";
        case CertifiedDesignStatus::verificationFailed: return "verification_failed";
        default: return "inconclusive";
    }
}

LmiCertificate verifyCertifiedDesign(const BilinearKrom& model, const ResidualSectorBound& sector,
    const CertifiedDesignSpec& spec, const std::vector<double>& gain, const std::vector<double>& X, double mu) {
    validate(model, sector, spec);
    const auto n = model.stateDimension, m = model.inputDimension;
    LmiCertificate rejected;
    rejected.margin = spec.certificateMargin;
    rejected.scope = "floating-point eigenvalue check with a rounding allowance; not interval arithmetic";
    if (gain.size() != m * n || X.size() != n * n) { rejected.reason = "gain or X has the wrong size"; return rejected; }
    if (!std::isfinite(mu) || !(mu > 0.0)) { rejected.reason = "mu must be positive and finite"; return rejected; }
    for (const auto* values : {&gain, &X})
        for (double v : *values)
            if (!std::isfinite(v)) { rejected.reason = "gain or X is not finite"; return rejected; }
    const auto o = normalize(model, sector, spec);
    const double R = spec.stateRadius;
    const MatrixXd Xraw = ConstRowMap(X.data(), o.n, o.n);
    const MatrixXd Kraw = ConstRowMap(gain.data(), o.m, o.n);
    if ((Xraw - Xraw.transpose()).norm() > 1e-9 * (1.0 + Xraw.norm())) {
        rejected.reason = "X is not symmetric"; return rejected;
    }
    const MatrixXd Xn = Xraw / (R * R);
    const MatrixXd Ln = Kraw * Xraw / R;
    return checkBlocks(certificateBlocks(o, Xn, Ln, mu / (R * R)), spec.certificateMargin);
}

CertifiedDesignResult designCertifiedBilinearController(const BilinearKrom& model,
    const ResidualSectorBound& sector, const CertifiedDesignSpec& spec) {
    validate(model, sector, spec);
    const auto o = normalize(model, sector, spec);
    const Layout layout{o.n, o.m, o.n * (o.n + 1) / 2};
    const auto problem = buildProblem(o, layout, spec.solverMargin);
    const auto solved = solveDenseSdp(problem, spec.solver);

    CertifiedDesignResult out;
    out.solverStatus = solved.status;
    out.primalResidual = solved.primalResidual;
    out.dualResidual = solved.dualResidual;
    out.gap = solved.gap;
    out.iterations = solved.iterations;
    if (solved.status == SdpStatus::primalInfeasible) {
        out.status = CertifiedDesignStatus::infeasible;
        out.message = "this sufficient LMI has no solution at the requested sector, box, radius, contraction and margin; this does not prove that no stabilizing controller exists";
        return out;
    }
    if (solved.status != SdpStatus::converged) {
        out.status = CertifiedDesignStatus::inconclusive;
        out.message = std::string("solver did not converge: ") + sdpStatusName(solved.status);
        return out;
    }
    MatrixXd Xn, Ln; double muN = 0.0, t = 0.0;
    layout.unpack(Eigen::Map<const VectorXd>(solved.x.data(), layout.variables()), Xn, Ln, muN, t);
    const Eigen::LLT<MatrixXd> factor(Xn);
    if (factor.info() != Eigen::Success) {
        out.status = CertifiedDesignStatus::verificationFailed;
        out.message = "solver returned X that is not positive definite";
        return out;
    }
    const double R = spec.stateRadius;
    // K = L X^-1 in original units: K = (L_n / X_n) / R.
    const MatrixXd K = factor.solve(Ln.transpose()).transpose() / R;
    const MatrixXd X = Xn * R * R;
    const double mu = muN * R * R;
    const auto gain = rowMajor(K);
    const auto xValues = rowMajor(X);
    out.certificate = verifyCertifiedDesign(model, sector, spec, gain, xValues, mu);
    if (!out.certificate.certified) {
        out.status = CertifiedDesignStatus::verificationFailed;
        out.message = "independent LMI verification rejected the solver output: " + out.certificate.reason;
        return out;
    }
    Eigen::SelfAdjointEigenSolver<MatrixXd> eig(X, Eigen::EigenvaluesOnly);
    out.status = CertifiedDesignStatus::certified;
    out.certified = true;
    out.gain = gain;
    out.X = xValues;
    out.mu = mu;
    out.inscribedRadius = std::sqrt(eig.eigenvalues()(0));
    out.message = "certified for the model with any residual in the stated sector";
    return out;
}

}  // namespace dynsyslab::core
