#include <dynsyslab/core/SystemStructure.hpp>

#include "EigenNumericsLinearSystems.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace dynsyslab::core {
using namespace linear_detail;
namespace {

bool isDiscrete(const StateSpaceModel& m) { return m.domain == TimeDomain::discrete; }

// Horizontal concatenation that tolerates empty blocks.
Matrix columns(std::initializer_list<const Matrix*> blocks, Eigen::Index rows) {
    Eigen::Index total = 0;
    for (const auto* block : blocks) total += block->cols();
    Matrix out(rows, total);
    Eigen::Index at = 0;
    for (const auto* block : blocks) { if (block->cols()) out.middleCols(at, block->cols()) = *block; at += block->cols(); }
    return out;
}

// Orthonormal basis of the smallest A-invariant subspace containing range(B),
// built block by block so each new direction is orthogonal to the previous ones.
StaircaseForm staircase(const Matrix& A, const Matrix& B, double tolerance) {
    const Eigen::Index n = A.rows();
    const double threshold = tolerance * std::max({A.norm(), B.norm(), 1.0e-300});
    StaircaseForm out;
    Matrix basis(n, 0);
    Matrix frontier = B;
    while (basis.cols() < n && frontier.cols() > 0) {
        if (basis.cols()) frontier -= basis * (basis.transpose() * frontier);  // twice: re-orthogonalize
        if (basis.cols()) frontier -= basis * (basis.transpose() * frontier);
        Eigen::JacobiSVD<Matrix> svd(frontier, Eigen::ComputeThinU);
        Eigen::Index added = 0;
        const auto& s = svd.singularValues();
        while (added < s.size() && s[added] > threshold) ++added;
        if (!added) break;
        const Matrix block = svd.matrixU().leftCols(added);
        basis = columns({&basis, &block}, n);
        out.blockSizes.push_back(static_cast<std::size_t>(added));
        frontier = A * block;
    }
    out.reachable = static_cast<std::size_t>(basis.cols());
    Matrix T(n, n); T.leftCols(basis.cols()) = basis;
    T.rightCols(n - basis.cols()) = orthogonalComplement(basis, n);
    out.transform = flat(T);
    return out;
}

Matrix basisOf(const StaircaseForm& form, Eigen::Index n) {
    return mapped(form.transform, n, n, "staircase").leftCols(static_cast<Eigen::Index>(form.reachable));
}

// Orthonormal basis of span(U) intersect span(V), both orthonormal.
Matrix intersection(const Matrix& U, const Matrix& V, double tolerance) {
    const Eigen::Index n = U.rows();
    if (!U.cols() || !V.cols()) return Matrix(n, 0);
    Eigen::JacobiSVD<Matrix> svd(U.transpose() * V, Eigen::ComputeThinU);
    Eigen::Index k = 0;  // principal angles with cos ~ 1
    while (k < svd.singularValues().size() && svd.singularValues()[k] > 1.0 - tolerance) ++k;
    return orthonormalBasis(U * svd.matrixU().leftCols(k), 1.0e-12);
}

// Orthonormal basis of the part of span(Super) orthogonal to span(Sub).
Matrix relativeComplement(const Matrix& Super, const Matrix& Sub) {
    if (!Super.cols()) return Super;
    Matrix residue = Super;
    if (Sub.cols()) { residue -= Sub * (Sub.transpose() * residue); residue -= Sub * (Sub.transpose() * residue); }
    const Eigen::Index expected = Super.cols() - Sub.cols();
    Eigen::JacobiSVD<Matrix> svd(residue, Eigen::ComputeThinU);
    return svd.matrixU().leftCols(std::max<Eigen::Index>(expected, 0));
}

}  // namespace

StaircaseForm controllabilityStaircase(const StateSpaceModel& model, double tolerance) {
    const auto r = realization(model);
    return staircase(r.A, r.B, tolerance);
}

StaircaseForm observabilityStaircase(const StateSpaceModel& model, double tolerance) {
    const auto r = realization(model);
    return staircase(r.A.transpose(), r.C.transpose(), tolerance);
}

KalmanDecomposition kalmanDecomposition(const StateSpaceModel& model, double tolerance) {
    const auto r = realization(model);
    const Eigen::Index n = r.A.rows();
    const Matrix Vc = basisOf(staircase(r.A, r.B, tolerance), n);
    const Matrix Vo = basisOf(staircase(r.A.transpose(), r.C.transpose(), tolerance), n);
    const Matrix Vno = orthogonalComplement(Vo, n);  // unobservable subspace = range(O^T)^perp
    const Matrix Vcno = intersection(Vc, Vno, 1.0e-8);
    const Matrix Vco = relativeComplement(Vc, Vcno);
    const Matrix Vncno = relativeComplement(Vno, Vcno);
    const Matrix used = columns({&Vco, &Vcno, &Vncno}, n);
    const Matrix Vnco = orthogonalComplement(orthonormalBasis(used, 1.0e-12), n);
    const Matrix T = columns({&Vco, &Vcno, &Vnco, &Vncno}, n);
    if (T.cols() != n) throw std::runtime_error("Kalman decomposition subspaces do not span the state space");
    Eigen::FullPivLU<Matrix> lu(T);
    if (!lu.isInvertible()) throw std::runtime_error("Kalman decomposition basis is singular");
    const Matrix Ti = lu.inverse();
    const Matrix At = Ti * r.A * T, Bt = Ti * r.B, Ct = r.C * T;
    KalmanDecomposition out;
    const Eigen::Index a = Vco.cols(), b = Vcno.cols(), c = Vnco.cols(), d = Vncno.cols();
    out.controllableObservable = static_cast<std::size_t>(a);
    out.controllableUnobservable = static_cast<std::size_t>(b);
    out.uncontrollableObservable = static_cast<std::size_t>(c);
    out.uncontrollableUnobservable = static_cast<std::size_t>(d);
    // Required zero blocks: uncontrollable rows of B and of A's controllable columns,
    // unobservable columns of C and of A's observable-to-unobservable coupling.
    double residual = Bt.bottomRows(c + d).norm() + Ct.middleCols(a, b).norm() + Ct.rightCols(d).norm();
    residual += At.block(a + b, 0, c + d, a + b).norm();          // controllable subspace invariant
    residual += At.block(0, a, a, b).norm() + At.block(a + b, a, c, b).norm(); // unobservable invariant
    residual += At.block(0, a + b + c, a, d).norm() + At.block(a + b, a + b + c, c, d).norm();
    out.structureResidual = residual / std::max({r.A.norm(), r.B.norm(), r.C.norm(), 1.0e-300});
    out.transform = flat(T); out.inverse = flat(Ti);
    out.transformed = fromRealization({At, Bt, Ct, r.D}, model.domain, model.sampleTime);
    return out;
}

MinimalRealization minimalRealization(const StateSpaceModel& model, double tolerance) {
    const auto k = kalmanDecomposition(model, tolerance);
    const Eigen::Index n = static_cast<Eigen::Index>(model.states);
    const auto a = static_cast<Eigen::Index>(k.controllableObservable);
    if (!a) throw std::invalid_argument("model has no controllable and observable part (zero transfer function)");
    const auto t = realization(k.transformed);
    MinimalRealization out;
    out.minimal = fromRealization({t.A.topLeftCorner(a, a), t.B.topRows(a), t.C.leftCols(a), t.D},
                                  model.domain, model.sampleTime);
    out.projection = flat(mapped(k.inverse, n, n, "Kalman inverse").topRows(a));
    out.lift = flat(mapped(k.transform, n, n, "Kalman transform").leftCols(a));
    out.removedStates = static_cast<std::size_t>(n - a);
    return out;
}

StructureReport kalmanStructureReport(const StateSpaceModel& model, double tolerance) {
    const auto k = kalmanDecomposition(model, tolerance);
    const auto t = realization(k.transformed);
    const bool discrete = isDiscrete(model);
    const auto a = static_cast<Eigen::Index>(k.controllableObservable), b = static_cast<Eigen::Index>(k.controllableUnobservable);
    const auto c = static_cast<Eigen::Index>(k.uncontrollableObservable), d = static_cast<Eigen::Index>(k.uncontrollableUnobservable);
    StructureReport out;
    out.controllableDimension = static_cast<std::size_t>(a + b);
    out.observableDimension = static_cast<std::size_t>(a + c);
    out.minimalOrder = static_cast<std::size_t>(a);
    out.controllable = c + d == 0; out.observable = b + d == 0;
    // Uncontrollable modes: spectrum of the quotient on [z_nco; z_ncno]; unobservable: [z_cno; z_ncno].
    out.uncontrollableModes = eigenvalues(t.A.bottomRightCorner(c + d, c + d));
    std::vector<Eigen::Index> hidden;
    for (Eigen::Index i = a; i < a + b; ++i) hidden.push_back(i);
    for (Eigen::Index i = a + b + c; i < a + b + c + d; ++i) hidden.push_back(i);
    Matrix unobservable(b + d, b + d);
    for (std::size_t i = 0; i < hidden.size(); ++i)
        for (std::size_t j = 0; j < hidden.size(); ++j)
            unobservable(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = t.A(hidden[i], hidden[j]);
    out.unobservableModes = eigenvalues(unobservable);
    auto stableMode = [&](std::complex<double> v) { return discrete ? std::abs(v) < 1.0 : v.real() < 0.0; };
    out.stabilizable = std::all_of(out.uncontrollableModes.begin(), out.uncontrollableModes.end(), stableMode);
    out.detectable = std::all_of(out.unobservableModes.begin(), out.unobservableModes.end(), stableMode);
    return out;
}

}  // namespace dynsyslab::core
