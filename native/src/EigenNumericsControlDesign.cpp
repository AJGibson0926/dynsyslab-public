#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/numerics/Control.hpp>

#include "EigenNumericsLinearSystems.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
using namespace linear_detail;
namespace {

using CVector = Eigen::VectorXcd;
bool isDiscrete(const StateSpaceModel& m) { return m.domain == TimeDomain::discrete; }

void checkPoles(const std::vector<std::complex<double>>& poles, std::size_t n) {
    if (poles.size() != n) throw std::invalid_argument("pole placement needs exactly one pole per state");
    std::vector<bool> used(n, false);
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(poles[i].real()) || !std::isfinite(poles[i].imag()))
            throw std::invalid_argument("poles must be finite");
        const double tolerance = 1.0e-10 * std::max(1.0, std::abs(poles[i]));
        if (std::abs(poles[i].imag()) <= tolerance || used[i]) continue;
        bool paired = false;
        for (std::size_t j = 0; j < n && !paired; ++j)
            if (j != i && !used[j] && std::abs(poles[j] - std::conj(poles[i])) <= tolerance)
                used[i] = used[j] = paired = true;
        if (!paired) throw std::invalid_argument("complex poles must appear in conjugate pairs");
    }
}

void verifyPlacement(const Matrix& closed, std::vector<std::complex<double>> desired) {
    auto actual = eigenvalues(closed);
    const double scale = std::max(1.0, closed.norm());
    for (const auto pole : desired) {  // greedy matching
        auto best = std::min_element(actual.begin(), actual.end(),
            [&](auto a, auto b) { return std::abs(a - pole) < std::abs(b - pole); });
        if (best == actual.end() || std::abs(*best - pole) > 1.0e-5 * scale)
            throw std::runtime_error("pole placement failed: closed-loop spectrum differs from the request "
                                     "(uncontrollable mode or ill-conditioned assignment)");
        actual.erase(best);
    }
}

Matrix ackermann(const Matrix& A, const Matrix& B, const std::vector<std::complex<double>>& poles) {
    const Eigen::Index n = A.rows();
    std::vector<std::complex<double>> coefficients{1.0};  // monic, highest power first
    for (const auto pole : poles) {
        std::vector<std::complex<double>> next(coefficients.size() + 1, 0.0);
        for (std::size_t k = 0; k < coefficients.size(); ++k) {
            next[k] += coefficients[k]; next[k + 1] -= pole * coefficients[k];
        }
        coefficients = std::move(next);
    }
    Matrix phi = Matrix::Zero(n, n);
    for (const auto c : coefficients) phi = (phi * A + c.real() * Matrix::Identity(n, n)).eval();
    Matrix ctrb(n, n); Matrix column = B;
    for (Eigen::Index k = 0; k < n; ++k) { ctrb.col(k) = column; column = A * column; }
    Eigen::FullPivLU<Matrix> lu(ctrb.transpose());
    if (!lu.isInvertible()) throw std::invalid_argument("pole placement requires a controllable pair");
    Vector last = Vector::Zero(n); last[n - 1] = 1.0;
    return lu.solve(last).transpose() * phi;
}

Matrix eigenvectorAssignment(const Matrix& A, const Matrix& B, const std::vector<std::complex<double>>& poles) {
    const Eigen::Index n = A.rows(), m = B.cols();
    std::vector<CMatrix> Nx(n), Nw(n);
    std::vector<int> partner(n, -1);
    std::vector<int> occurrence(n, 0);
    for (Eigen::Index i = 0; i < n; ++i) {
        CMatrix M(n, n + m);
        M.leftCols(n) = A.cast<std::complex<double>>() - poles[i] * CMatrix::Identity(n, n);
        M.rightCols(m) = B.cast<std::complex<double>>();
        Eigen::JacobiSVD<CMatrix> svd(M, Eigen::ComputeFullV);
        const auto& s = svd.singularValues();
        if (s[n - 1] <= 1.0e-12 * std::max(1.0, s[0]))
            throw std::invalid_argument("pole placement: [A - pI, B] loses rank (uncontrollable mode)");
        const CMatrix null = svd.matrixV().rightCols(m);
        Nx[i] = null.topRows(n); Nw[i] = null.bottomRows(m);
        for (Eigen::Index j = 0; j < i; ++j) {
            if (std::abs(poles[j] - poles[i]) <= 1.0e-10 * std::max(1.0, std::abs(poles[i]))) ++occurrence[i];
            if (partner[j] < 0 && std::abs(poles[i].imag()) > 0
                && std::abs(poles[j] - std::conj(poles[i])) <= 1.0e-10 * std::max(1.0, std::abs(poles[i])))
                { partner[j] = static_cast<int>(i); partner[i] = static_cast<int>(j); }
        }
        if (occurrence[i] >= m) throw std::invalid_argument("a pole repeats more often than there are inputs");
    }
    CMatrix X(n, n), W(m, n);
    auto setColumn = [&](Eigen::Index i, const CVector& c) {
        CVector x = Nx[i] * c, w = Nw[i] * c;
        const double size = x.norm();
        if (!(size > 0.0)) throw std::runtime_error("eigenvector assignment produced a zero vector");
        x /= size; w /= size;
        if (std::abs(poles[i].imag()) == 0.0) {  // keep real eigenvectors for real poles
            const double phase = std::arg(x[0] != 0.0 ? x[0] : std::complex<double>(1.0));
            const auto rotation = std::polar(1.0, -phase);
            x = (x * rotation).real().cast<std::complex<double>>(); w = (w * rotation).real().cast<std::complex<double>>();
            const double norm = x.norm(); x /= norm; w /= norm;
        }
        X.col(i) = x; W.col(i) = w;
        if (partner[i] >= 0) { X.col(partner[i]) = x.conjugate(); W.col(partner[i]) = w.conjugate(); }
    };
    for (Eigen::Index i = 0; i < n; ++i)
        if (partner[i] < 0 || partner[i] > i) {
            CVector c = CVector::Zero(m); c[occurrence[i] % m] = 1.0; setColumn(i, c);
        }
    for (int sweep = 0; sweep < 8; ++sweep)  // KNV method 0: make each x_i orthogonal to the others
        for (Eigen::Index i = 0; i < n; ++i) {
            if (partner[i] >= 0 && partner[i] < i) continue;
            CMatrix others(n, n - 1);
            for (Eigen::Index j = 0, k = 0; j < n; ++j) if (j != i) others.col(k++) = X.col(j);
            Eigen::HouseholderQR<CMatrix> qr(others);
            const CVector q = CMatrix(qr.householderQ()).col(n - 1);
            setColumn(i, Nx[i].completeOrthogonalDecomposition().solve(q));
        }
    Eigen::FullPivLU<CMatrix> lu(X);
    if (!lu.isInvertible()) throw std::runtime_error("eigenvector assignment produced a singular eigenbasis");
    const CMatrix K = -W * lu.inverse();
    if (K.imag().norm() > 1.0e-8 * std::max(1.0, K.real().norm()))
        throw std::runtime_error("eigenvector assignment produced a complex gain");
    return K.real();
}

Matrix placeGain(const Matrix& A, const Matrix& B, const std::vector<std::complex<double>>& poles) {
    checkPoles(poles, static_cast<std::size_t>(A.rows()));
    if (B.cols() == 0) throw std::invalid_argument("pole placement needs at least one input");
    const Matrix K = B.cols() == 1 ? ackermann(A, B, poles) : eigenvectorAssignment(A, B, poles);
    verifyPlacement(A - B * K, poles);
    return K;
}

struct Riccati { Matrix P, K; double residual{}; };

// Structure-preserving doubling for P = A^T P A - A^T P B (R + B^T P B)^-1 B^T P A + Q.
Riccati discreteRiccati(const Matrix& A, const Matrix& B, const Matrix& Q, const Matrix& R) {
    const Eigen::Index n = A.rows();
    Eigen::LDLT<Matrix> rFactor(R);
    if (rFactor.info() != Eigen::Success || !rFactor.isPositive() || R.rows() == 0)
        throw std::invalid_argument("DARE input weight must be positive definite");
    Matrix Ak = A, Gk = B * rFactor.solve(B.transpose()), Hk = Q;
    const Matrix I = Matrix::Identity(n, n);
    for (int iteration = 0; iteration < 200; ++iteration) {
        Eigen::PartialPivLU<Matrix> W(I + Gk * Hk);
        const Matrix WA = W.solve(Ak), WG = W.solve(Gk);
        const Matrix nextH = Hk + Ak.transpose() * Hk * WA;
        Gk = (Gk + Ak * WG * Ak.transpose()).eval();
        Ak = (Ak * WA).eval();
        const double change = (nextH - Hk).norm();
        Hk = 0.5 * (nextH + nextH.transpose());
        if (!Hk.allFinite()) throw std::runtime_error("DARE doubling diverged (unstabilizable or undetectable)");
        if (change <= 1.0e-14 * std::max(1.0, Hk.norm())) {
            Riccati out; out.P = Hk;
            const Matrix normal = R + B.transpose() * Hk * B;
            out.K = normal.ldlt().solve(B.transpose() * Hk * A);
            const Matrix residual = A.transpose() * Hk * A - A.transpose() * Hk * B * out.K + Q - Hk;
            out.residual = residual.norm() / std::max({Hk.norm(), Q.norm(), 1.0e-300});
            if (!stable(A - B * out.K, true)) throw std::runtime_error("DARE solution is not stabilizing");
            if (out.residual > 1.0e-8) throw std::runtime_error("DARE residual exceeds tolerance");
            return out;
        }
    }
    throw std::runtime_error("DARE doubling did not converge");
}

}  // namespace

std::vector<double> placePoles(const StateSpaceModel& model, const std::vector<std::complex<double>>& poles) {
    const auto r = realization(model);
    return flat(placeGain(r.A, r.B, poles));
}

std::vector<double> placeObserverPoles(const StateSpaceModel& model, const std::vector<std::complex<double>>& poles) {
    const auto r = realization(model);
    if (r.C.rows() == 0) throw std::invalid_argument("observer placement needs at least one output");
    return flat(placeGain(r.A.transpose(), r.C.transpose(), poles).transpose());
}

numerics::ContinuousLqrResult designLQR(const StateSpaceModel& model, const std::vector<double>& q, const std::vector<double>& rw) {
    const auto r = realization(model);
    const Matrix Q = mapped(q, model.states, model.states, "LQR state penalty");
    const Matrix R = mapped(rw, model.inputs, model.inputs, "LQR input penalty");
    numerics::ContinuousLqrResult out;
    if (isDiscrete(model)) {
        const auto s = discreteRiccati(r.A, r.B, Q, R);
        out.riccati = flat(s.P); out.gain = flat(s.K); out.relativeResidual = s.residual;
        out.closedLoopPoles = eigenvalues(r.A - r.B * s.K);
    } else {
        const auto s = numerics::solveContinuousLqr(model.A, model.B, q, rw, model.states, model.inputs);
        out.riccati = s.riccati; out.gain = s.gain; out.relativeResidual = s.relativeResidual;
        out.closedLoopPoles = s.closedLoopPoles;
    }
    return out;
}

EstimatorDesign designLQE(const StateSpaceModel& model, const std::vector<double>& w, const std::vector<double>& v) {
    const auto r = realization(model);
    const Matrix W = mapped(w, model.states, model.states, "LQE process covariance");
    const Matrix V = mapped(v, model.outputs, model.outputs, "LQE measurement covariance");
    EstimatorDesign out;
    Matrix L;
    if (isDiscrete(model)) {
        const auto s = discreteRiccati(r.A.transpose(), r.C.transpose(), W, V);
        L = s.K.transpose();
        const Matrix innovation = r.C * s.P * r.C.transpose() + V;
        out.filterGain = flat(innovation.ldlt().solve(r.C * s.P).transpose());
        out.covariance = flat(s.P); out.relativeResidual = s.residual;
    } else {
        const auto s = numerics::solveContinuousLqr(flat(r.A.transpose()), flat(r.C.transpose()), w, v,
                                                    model.states, model.outputs);
        L = mapped(s.gain, model.outputs, model.states, "LQE dual gain").transpose();
        out.covariance = s.riccati; out.relativeResidual = s.relativeResidual;
    }
    out.gain = flat(L);
    out.observerPoles = eigenvalues(r.A - L * r.C);
    return out;
}

StateSpaceModel stateFeedbackClosedLoop(const StateSpaceModel& model, const std::vector<double>& gain) {
    const auto r = realization(model);
    const Matrix K = mapped(gain, model.inputs, model.states, "state-feedback gain");
    return fromRealization({r.A - r.B * K, r.B, r.C - r.D * K, r.D}, model.domain, model.sampleTime);
}

LqgDesign designLQG(const StateSpaceModel& model, const std::vector<double>& q, const std::vector<double>& rw,
                    const std::vector<double>& w, const std::vector<double>& v) {
    const auto r = realization(model);
    LqgDesign out;
    out.regulator = designLQR(model, q, rw);
    out.estimator = designLQE(model, w, v);
    const Matrix K = mapped(out.regulator.gain, model.inputs, model.states, "LQG regulator gain");
    const Matrix L = mapped(out.estimator.gain, model.states, model.outputs, "LQG estimator gain");
    const Eigen::Index n = r.A.rows();
    out.compensator = fromRealization({r.A - r.B * K - L * r.C + L * r.D * K, L, -K,
        Matrix::Zero(r.B.cols(), r.C.rows())}, model.domain, model.sampleTime);
    Realization loop{Matrix(2 * n, 2 * n), Matrix(2 * n, r.B.cols()), Matrix(r.C.rows(), 2 * n), r.D};
    loop.A << r.A, -r.B * K, L * r.C, r.A - r.B * K - L * r.C;
    loop.B << r.B, r.B;
    loop.C << r.C, -r.D * K;
    out.closedLoop = fromRealization(loop, model.domain, model.sampleTime);
    out.closedLoopPoles = eigenvalues(loop.A);
    return out;
}

}  // namespace dynsyslab::core
