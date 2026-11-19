#include <dynsyslab/core/BalancedTruncation.hpp>

#include "EigenNumericsLinearSystems.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace dynsyslab::core {
using namespace linear_detail;
namespace {

bool isDiscrete(const StateSpaceModel& m) { return m.domain == TimeDomain::discrete; }

struct SquareRoot { Matrix T, Ti; Vector sigma; Eigen::Index order{}; };

SquareRoot squareRoot(const Realization& r, bool discrete, double tolerance) {
    if (!stable(r.A, discrete))
        throw std::invalid_argument("balancing requires an asymptotically stable model");
    const Matrix Wc = lyapunov(r.A, r.B * r.B.transpose(), discrete);
    const Matrix Wo = lyapunov(r.A.transpose(), r.C.transpose() * r.C, discrete);
    const Matrix Lc = psdFactor(Wc), Lo = psdFactor(Wo);
    Eigen::JacobiSVD<Matrix> svd(Lo.transpose() * Lc, Eigen::ComputeFullU | Eigen::ComputeFullV);
    SquareRoot out; out.sigma = svd.singularValues();
    const double largest = out.sigma.size() ? out.sigma[0] : 0.0;
    while (out.order < out.sigma.size() && out.sigma[out.order] > tolerance * largest
           && out.sigma[out.order] > 0.0) ++out.order;
    const Eigen::Index k = out.order;
    const Vector inverseRoot = out.sigma.head(k).cwiseSqrt().cwiseInverse();
    out.T = Lc * svd.matrixV().leftCols(k) * inverseRoot.asDiagonal();
    out.Ti = inverseRoot.asDiagonal() * svd.matrixU().leftCols(k).transpose() * Lo.transpose();
    return out;
}

Realization reduceBalanced(const Realization& b, Eigen::Index r, BalancedReductionMethod method,
                           bool discrete) {
    const Eigen::Index k = b.A.rows();
    if (r >= k || method == BalancedReductionMethod::truncation)
        return {b.A.topLeftCorner(r, r), b.B.topRows(r), b.C.leftCols(r), b.D};
    const Matrix A12 = b.A.topRightCorner(r, k - r), A21 = b.A.bottomLeftCorner(k - r, r);
    const Matrix A22 = b.A.bottomRightCorner(k - r, k - r);
    const Matrix B2 = b.B.bottomRows(k - r), C2 = b.C.rightCols(k - r);
    // Singular perturbation: continuous sets x2' = 0, discrete sets x2[k+1] = x2[k].
    const Matrix pivot = discrete ? Matrix(Matrix::Identity(k - r, k - r) - A22) : Matrix(-A22);
    Eigen::FullPivLU<Matrix> lu(pivot);
    if (!lu.isInvertible()) throw std::runtime_error("residualization pivot block is singular");
    const Matrix solveA21 = lu.solve(A21), solveB2 = lu.solve(B2);
    return {b.A.topLeftCorner(r, r) + A12 * solveA21, b.B.topRows(r) + A12 * solveB2,
            b.C.leftCols(r) + C2 * solveA21, b.D + C2 * solveB2};
}

}  // namespace

const char* balancedReductionMethodName(BalancedReductionMethod method) noexcept {
    return method == BalancedReductionMethod::residualization ? "residualization" : "truncation";
}

std::vector<double> hankelSingularValues(const StateSpaceModel& model) {
    const auto s = squareRoot(realization(model), isDiscrete(model), 0.0);
    return {s.sigma.data(), s.sigma.data() + s.sigma.size()};
}

BalancedRealization balancedRealization(const StateSpaceModel& model, double tolerance) {
    const auto r = realization(model);
    const auto s = squareRoot(r, isDiscrete(model), tolerance);
    BalancedRealization out;
    out.balanced = fromRealization({s.Ti * r.A * s.T, s.Ti * r.B, r.C * s.T, r.D},
                                   model.domain, model.sampleTime);
    out.hankelSingularValues.assign(s.sigma.data(), s.sigma.data() + s.sigma.size());
    out.transform = flat(s.T); out.inverse = flat(s.Ti);
    out.order = static_cast<std::size_t>(s.order);
    return out;
}

BalancedReduction balancedReduction(const StateSpaceModel& model, const BalancedReductionOptions& o) {
    const auto full = realization(model);
    const bool discrete = isDiscrete(model);
    const Eigen::Index n = full.A.rows();
    const auto split = splitStableUnstable(model, o.stabilityMargin);
    const auto ku = static_cast<Eigen::Index>(split.unstableStates);
    if (ku && !o.retainUnstable)
        throw std::invalid_argument("model has unstable or marginal modes and retainUnstable is false");
    BalancedReduction out; out.method = o.method; out.fullOrder = static_cast<std::size_t>(n);
    out.unstableOrder = static_cast<std::size_t>(ku);
    const Eigen::Index ks = n - ku;
    Realization reducedStable{Matrix(0, 0), Matrix(0, full.B.cols()), Matrix(full.C.rows(), 0), full.D};
    Matrix Ps(0, ks), Ls(ks, 0);
    if (ks > 0) {
        const auto stablePart = realization(split.stable);
        const auto s = squareRoot(stablePart, discrete, o.hankelRelativeTolerance);
        out.hankelSingularValues.assign(s.sigma.data(), s.sigma.data() + s.sigma.size());
        const Vector& sigma = s.sigma;
        auto tail = [&](Eigen::Index r) { return 2.0 * sigma.tail(sigma.size() - r).sum(); };
        Eigen::Index r = s.order;
        if (o.order) {
            if (static_cast<Eigen::Index>(o.order) < ku)
                throw std::invalid_argument("reduced order is smaller than the retained unstable part");
            r = std::min<Eigen::Index>(static_cast<Eigen::Index>(o.order) - ku, s.order);
        } else if (o.errorTolerance > 0.0) {
            r = 0; while (r < s.order && tail(r) > o.errorTolerance) ++r;
        }
        const Realization balanced{s.Ti * stablePart.A * s.T, s.Ti * stablePart.B,
                                   stablePart.C * s.T, stablePart.D};
        reducedStable = reduceBalanced(balanced, r, o.method, discrete);
        out.stablePart = split.stable;
        if (r) out.reducedStablePart = fromRealization(reducedStable, model.domain, model.sampleTime);
        out.errorBound = tail(r);
        out.lowerBound = r < sigma.size() ? sigma[r] : 0.0;
        Ps = s.Ti.topRows(r); Ls = s.T.leftCols(r);
        out.stableOrder = static_cast<std::size_t>(r);
    }
    const Eigen::Index rs = static_cast<Eigen::Index>(out.stableOrder), rt = rs + ku;
    out.order = static_cast<std::size_t>(rt);
    const auto unstable = ku ? realization(split.unstable)
        : Realization{Matrix(0, 0), Matrix(0, full.B.cols()), Matrix(full.C.rows(), 0), full.D};
    Realization reduced{Matrix::Zero(rt, rt), Matrix(rt, full.B.cols()), Matrix(full.C.rows(), rt), reducedStable.D};
    reduced.A.topLeftCorner(rs, rs) = reducedStable.A;
    reduced.A.bottomRightCorner(ku, ku) = unstable.A;
    reduced.B.topRows(rs) = reducedStable.B; reduced.B.bottomRows(ku) = unstable.B;
    reduced.C.leftCols(rs) = reducedStable.C; reduced.C.rightCols(ku) = unstable.C;
    // Coordinates: x = T_split [x_s; x_u]; z_s = P_s x_s; x_u kept.
    const Matrix T = mapped(split.transform, n, n, "split transform");
    const Matrix Ti = mapped(split.inverse, n, n, "split inverse");
    Matrix blockP = Matrix::Zero(rt, n), blockL = Matrix::Zero(n, rt);
    blockP.topLeftCorner(rs, ks) = Ps; blockL.topLeftCorner(ks, rs) = Ls;
    blockP.bottomRightCorner(ku, ku) = Matrix::Identity(ku, ku);
    blockL.bottomRightCorner(ku, ku) = Matrix::Identity(ku, ku);
    out.projection = flat(blockP * Ti); out.lift = flat(T * blockL);
    out.reduced = rt ? fromRealization(reduced, model.domain, model.sampleTime) : StateSpaceModel{};
    if (!rt) { out.reduced.D = flat(full.D); out.reduced.inputs = model.inputs; out.reduced.outputs = model.outputs;
               out.reduced.domain = model.domain; out.reduced.sampleTime = model.sampleTime; }
    return out;
}

EraResult eigensystemRealization(const std::vector<double>& markov, std::size_t samples,
    std::size_t p, std::size_t m, const EraOptions& o) {
    if (!p || !m || markov.size() != samples * p * m)
        throw std::invalid_argument("ERA Markov parameters must be samples x outputs x inputs");
    if (samples < 3) throw std::invalid_argument("ERA needs at least three Markov parameters");
    const std::size_t rows = o.blockRows ? o.blockRows : (samples - 1) / 2;
    const std::size_t cols = o.blockColumns ? o.blockColumns : samples - 1 - rows;
    if (rows + cols + 1 > samples || !rows || !cols)
        throw std::invalid_argument("ERA Hankel blocks exceed the available Markov parameters");
    auto Y = [&](std::size_t k, std::size_t i, std::size_t j) { return markov[(k * p + i) * m + j]; };
    Matrix H0(rows * p, cols * m), H1(rows * p, cols * m);
    for (std::size_t br = 0; br < rows; ++br)
        for (std::size_t bc = 0; bc < cols; ++bc)
            for (std::size_t i = 0; i < p; ++i)
                for (std::size_t j = 0; j < m; ++j) {
                    H0(br * p + i, bc * m + j) = Y(br + bc + 1, i, j);
                    H1(br * p + i, bc * m + j) = Y(br + bc + 2, i, j);
                }
    if (!H0.allFinite() || !H1.allFinite()) throw std::invalid_argument("ERA Markov parameters must be finite");
    Eigen::JacobiSVD<Matrix> svd(H0, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Vector s = svd.singularValues();
    Eigen::Index r = o.order ? static_cast<Eigen::Index>(o.order) : 0;
    if (!r) while (r < s.size() && s[r] > o.relativeTolerance * s[0]) ++r;
    if (r > s.size() || r == 0 || s[r - 1] <= 0.0) throw std::invalid_argument("ERA order exceeds the Hankel rank");
    const Vector root = s.head(r).cwiseSqrt(), inverseRoot = root.cwiseInverse();
    const Matrix U = svd.matrixU().leftCols(r), V = svd.matrixV().leftCols(r);
    Realization era;
    era.A = inverseRoot.asDiagonal() * U.transpose() * H1 * V * inverseRoot.asDiagonal();
    era.B = (root.asDiagonal() * V.transpose()).leftCols(static_cast<Eigen::Index>(m));
    era.C = (U * root.asDiagonal()).topRows(static_cast<Eigen::Index>(p));
    era.D = Matrix(static_cast<Eigen::Index>(p), static_cast<Eigen::Index>(m));
    for (std::size_t i = 0; i < p; ++i) for (std::size_t j = 0; j < m; ++j) era.D(i, j) = Y(0, i, j);
    EraResult out; out.model = fromRealization(era, TimeDomain::discrete, o.sampleTime);
    out.hankelSingularValues.assign(s.data(), s.data() + s.size()); out.order = static_cast<std::size_t>(r);
    return out;
}

std::vector<double> observerKalmanMarkovParameters(const std::vector<double>& u,
    const std::vector<double>& y, std::size_t N, std::size_t m, std::size_t p,
    std::size_t l, std::size_t count) {
    if (!m || !p || !l || !count) throw std::invalid_argument("OKID dimensions must be positive");
    if (u.size() != N * m || y.size() != N * p) throw std::invalid_argument("OKID records must be samples x channels");
    const std::size_t width = m + l * (m + p);
    if (N <= l || N - l < width) throw std::invalid_argument("OKID record is too short for the observer lag");
    const std::size_t rows = N - l;
    Matrix Phi(rows, width), Yt(rows, p);
    for (std::size_t k = l; k < N; ++k) {
        const std::size_t row = k - l;
        for (std::size_t j = 0; j < m; ++j) Phi(row, j) = u[k * m + j];
        for (std::size_t i = 1; i <= l; ++i) {
            const std::size_t base = m + (i - 1) * (m + p);
            for (std::size_t j = 0; j < m; ++j) Phi(row, base + j) = u[(k - i) * m + j];
            for (std::size_t j = 0; j < p; ++j) Phi(row, base + m + j) = y[(k - i) * p + j];
        }
        for (std::size_t j = 0; j < p; ++j) Yt(row, j) = y[k * p + j];
    }
    if (!Phi.allFinite() || !Yt.allFinite()) throw std::invalid_argument("OKID records must be finite");
    const Matrix theta = Phi.completeOrthogonalDecomposition().solve(Yt).transpose();  // p x width
    std::vector<Matrix> Ym(count);
    Ym[0] = theta.leftCols(m);
    for (std::size_t k = 1; k < count; ++k) {
        Matrix value = Matrix::Zero(p, m);
        if (k <= l) value = theta.middleCols(m + (k - 1) * (m + p), m);
        for (std::size_t i = 1; i <= std::min(k, l); ++i)
            value += theta.middleCols(m + (i - 1) * (m + p) + m, p) * Ym[k - i];
        Ym[k] = value;
    }
    std::vector<double> out(count * p * m);
    for (std::size_t k = 0; k < count; ++k)
        for (std::size_t i = 0; i < p; ++i) for (std::size_t j = 0; j < m; ++j)
            out[(k * p + i) * m + j] = Ym[k](i, j);
    return out;
}

BalancedPod balancedPod(const std::vector<double>& primal, const std::vector<double>& adjoint,
    std::size_t n, std::size_t a, std::size_t b, std::size_t order, double tolerance) {
    const Matrix X = mapped(primal, n, a, "balanced POD primal snapshots");
    const Matrix Z = mapped(adjoint, n, b, "balanced POD adjoint snapshots");
    Eigen::JacobiSVD<Matrix> svd(Z.transpose() * X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Vector s = svd.singularValues();
    Eigen::Index r = order ? static_cast<Eigen::Index>(order) : 0;
    if (!r) while (r < s.size() && s[r] > tolerance * s[0]) ++r;
    if (r == 0 || r > s.size() || s[r - 1] <= 0.0) throw std::invalid_argument("balanced POD order exceeds the snapshot rank");
    const Vector inverseRoot = s.head(r).cwiseSqrt().cwiseInverse();
    BalancedPod out; out.states = n; out.order = static_cast<std::size_t>(r);
    out.modes = flat(X * svd.matrixV().leftCols(r) * inverseRoot.asDiagonal());
    out.adjointModes = flat(Z * svd.matrixU().leftCols(r) * inverseRoot.asDiagonal());
    out.hankelSingularValues.assign(s.data(), s.data() + s.size());
    return out;
}

StateSpaceModel projectStateSpace(const StateSpaceModel& model, const std::vector<double>& lift,
    const std::vector<double>& adjointModes, std::size_t order) {
    const auto r = realization(model);
    const Matrix L = mapped(lift, model.states, order, "projection lift");
    const Matrix Pt = mapped(adjointModes, model.states, order, "projection adjoint modes").transpose();
    return fromRealization({Pt * r.A * L, Pt * r.B, r.C * L, r.D}, model.domain, model.sampleTime);
}

}  // namespace dynsyslab::core
