#include <dynsyslab/core/SystemNorms.hpp>

#include "EigenNumericsLinearSystems.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace dynsyslab::core {
using namespace linear_detail;
namespace {

constexpr double kPi = 3.14159265358979323846;
const double kInfinity = std::numeric_limits<double>::infinity();
bool isDiscrete(const StateSpaceModel& m) { return m.domain == TimeDomain::discrete; }

// Transfer evaluator on the complex Schur form: G(s) = C U (sI - T)^-1 U^H B + D.
class TransferEvaluator {
public:
    explicit TransferEvaluator(const StateSpaceModel& model) : model_(model) {
        const auto r = realization(model);
        Eigen::ComplexSchur<CMatrix> schur(r.A.cast<std::complex<double>>());
        if (schur.info() != Eigen::Success) throw std::runtime_error("transfer evaluation Schur form failed");
        T_ = schur.matrixT();
        B_ = schur.matrixU().adjoint() * r.B.cast<std::complex<double>>();
        C_ = r.C.cast<std::complex<double>>() * schur.matrixU();
        D_ = r.D.cast<std::complex<double>>();
    }
    CMatrix at(std::complex<double> s) const {
        const CMatrix shifted = s * CMatrix::Identity(T_.rows(), T_.cols()) - T_;
        const double pivot = shifted.diagonal().cwiseAbs().minCoeff();
        if (!(pivot > 1.0e-300)) throw std::invalid_argument("frequency coincides with a pole");
        return C_ * shifted.triangularView<Eigen::Upper>().solve(B_) + D_;
    }
    CMatrix atFrequency(double omega) const {
        return at(isDiscrete(model_) ? std::polar(1.0, omega * model_.sampleTime)
                                     : std::complex<double>(0.0, omega));
    }
    double largestSingularValue(double omega) const {
        const CMatrix g = atFrequency(omega);
        if (g.size() == 0) return 0.0;
        Eigen::JacobiSVD<CMatrix> svd(g);
        return svd.singularValues()[0];
    }
private:
    const StateSpaceModel& model_;
    CMatrix T_, B_, C_, D_;
};

double largestSingularValue(const Matrix& M) {
    if (M.size() == 0) return 0.0;
    Eigen::JacobiSVD<Matrix> svd(M); return svd.singularValues()[0];
}

std::vector<double> imaginaryAxisFrequencies(const Matrix& H) {
    Eigen::EigenSolver<Matrix> solver(H, false);
    if (solver.info() != Eigen::Success) throw std::runtime_error("H-infinity Hamiltonian eigensolve failed");
    const double scale = std::max(1.0, H.cwiseAbs().maxCoeff());
    std::vector<double> out;
    for (Eigen::Index i = 0; i < solver.eigenvalues().size(); ++i) {
        const auto v = solver.eigenvalues()[i];
        if (std::abs(v.real()) <= 1.0e-8 * std::max(scale, std::abs(v)) && v.imag() >= 0.0) out.push_back(v.imag());
    }
    std::sort(out.begin(), out.end());
    return out;
}

HInfinityNorm continuousHInfinity(const StateSpaceModel& model, double tolerance) {
    const auto r = realization(model);
    const TransferEvaluator G(model);
    const Eigen::Index m = r.B.cols(), p = r.C.rows();
    HInfinityNorm out;
    out.value = largestSingularValue(r.D); out.peakFrequency = kInfinity;
    std::vector<double> candidates{0.0};
    for (const auto pole : eigenvalues(r.A)) { candidates.push_back(std::abs(pole)); candidates.push_back(std::abs(pole.imag())); }
    for (double w : candidates) {
        const double value = G.largestSingularValue(w);
        if (value > out.value) { out.value = value; out.peakFrequency = w; }
    }
    if (!(out.value > 0.0)) { out.peakFrequency = 0.0; return out; }
    for (out.iterations = 1; out.iterations <= 60; ++out.iterations) {
        const double gamma = (1.0 + 2.0 * tolerance) * out.value;
        const Matrix R = r.D.transpose() * r.D - gamma * gamma * Matrix::Identity(m, m);
        const Matrix S = r.D * r.D.transpose() - gamma * gamma * Matrix::Identity(p, p);
        const Matrix Ri = R.inverse(), Si = S.inverse();
        const Eigen::Index n = r.A.rows();
        Matrix H(2 * n, 2 * n);
        H.topLeftCorner(n, n) = r.A - r.B * Ri * r.D.transpose() * r.C;
        H.topRightCorner(n, n) = -gamma * r.B * Ri * r.B.transpose();
        H.bottomLeftCorner(n, n) = gamma * r.C.transpose() * Si * r.C;
        H.bottomRightCorner(n, n) = -r.A.transpose() + r.C.transpose() * r.D * Ri * r.B.transpose();
        const auto omega = imaginaryAxisFrequencies(H);
        if (omega.empty()) return out;
        std::vector<double> probes;
        if (omega.size() == 1) probes.push_back(omega[0]);
        for (std::size_t i = 0; i + 1 < omega.size(); ++i) probes.push_back(0.5 * (omega[i] + omega[i + 1]));
        double best = out.value;
        for (double w : probes) {
            const double value = G.largestSingularValue(w);
            if (value > best) { best = value; out.peakFrequency = w; }
        }
        if (best <= out.value * (1.0 + tolerance)) { out.value = std::max(best, out.value); return out; }
        out.value = best;
    }
    return out;
}

}  // namespace

FrequencyResponse frequencyResponse(const StateSpaceModel& model, const std::vector<double>& frequencies) {
    validateStateSpace(model);
    const TransferEvaluator G(model);
    FrequencyResponse out; out.frequencies = frequencies;
    out.outputs = model.outputs; out.inputs = model.inputs;
    const std::size_t k = std::min(model.outputs, model.inputs);
    for (double w : frequencies) {
        if (!std::isfinite(w)) throw std::invalid_argument("frequencies must be finite");
        const CMatrix g = G.atFrequency(w);
        for (Eigen::Index i = 0; i < g.rows(); ++i)
            for (Eigen::Index j = 0; j < g.cols(); ++j) out.response.push_back(g(i, j));
        if (k) {
            Eigen::JacobiSVD<CMatrix> svd(g);
            for (std::size_t i = 0; i < k; ++i) out.singularValues.push_back(svd.singularValues()[static_cast<Eigen::Index>(i)]);
        }
    }
    return out;
}

std::vector<double> logarithmicFrequencyGrid(const StateSpaceModel& model, std::size_t count) {
    if (count < 2) throw std::invalid_argument("frequency grid needs at least two points");
    std::vector<double> magnitudes;
    for (auto pole : systemPoles(model)) {
        if (isDiscrete(model)) pole = std::log(pole) / model.sampleTime;
        if (std::abs(pole) > 1.0e-12) magnitudes.push_back(std::abs(pole));
    }
    double low = 1.0e-3, high = 1.0e3;
    if (!magnitudes.empty()) {
        low = *std::min_element(magnitudes.begin(), magnitudes.end()) / 100.0;
        high = *std::max_element(magnitudes.begin(), magnitudes.end()) * 100.0;
    }
    if (isDiscrete(model)) high = std::min(high, kPi / model.sampleTime);
    if (!(high > low)) low = high / 1.0e4;
    std::vector<double> grid(count);
    for (std::size_t i = 0; i < count; ++i)
        grid[i] = low * std::pow(high / low, static_cast<double>(i) / static_cast<double>(count - 1));
    return grid;
}

double h2Norm(const StateSpaceModel& model) {
    const auto r = realization(model);
    const bool discrete = isDiscrete(model);
    if (!stable(r.A, discrete)) return kInfinity;
    if (!discrete && r.D.norm() > 0.0) return kInfinity;
    const Matrix Wc = lyapunov(r.A, r.B * r.B.transpose(), discrete);
    double trace = (r.C * Wc * r.C.transpose()).trace();
    if (discrete) trace += (r.D * r.D.transpose()).trace();
    return std::sqrt(std::max(trace, 0.0));
}

HInfinityNorm hInfinityNorm(const StateSpaceModel& model, double tolerance) {
    validateStateSpace(model);
    if (!isAsymptoticallyStable(model)) return {kInfinity, 0.0, 0};
    if (!model.inputs || !model.outputs) return {0.0, 0.0, 0};
    if (!isDiscrete(model)) return continuousHInfinity(model, tolerance);
    auto bilinear = model; bilinear.sampleTime = 2.0;  // s = (z - 1)/(z + 1) maps |z| = 1 to the imaginary axis
    auto out = continuousHInfinity(continuousFromTustin(bilinear), tolerance);
    out.peakFrequency = std::isinf(out.peakFrequency) ? kPi / model.sampleTime
                                                      : 2.0 * std::atan(out.peakFrequency) / model.sampleTime;
    return out;
}

std::vector<std::complex<double>> transmissionZeros(const StateSpaceModel& model) {
    const auto r = realization(model);
    if (model.inputs != model.outputs || !model.inputs)
        throw std::invalid_argument("transmission zeros are implemented for square systems (inputs == outputs)");
    const Eigen::Index n = r.A.rows(), m = r.B.cols();
    Eigen::MatrixXd M(n + m, n + m), N = Eigen::MatrixXd::Zero(n + m, n + m);
    M.topLeftCorner(n, n) = r.A; M.topRightCorner(n, m) = r.B;
    M.bottomLeftCorner(m, n) = r.C; M.bottomRightCorner(m, m) = r.D;
    N.topLeftCorner(n, n) = Eigen::MatrixXd::Identity(n, n);
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver(M, N);
    if (solver.info() != Eigen::Success) throw std::runtime_error("Rosenbrock pencil QZ failed");
    const double scale = std::max(1.0, M.cwiseAbs().maxCoeff());
    std::vector<std::complex<double>> zeros;
    for (Eigen::Index i = 0; i < solver.alphas().size(); ++i) {
        const auto alpha = solver.alphas()[i];
        const double beta = solver.betas()[i];
        if (std::abs(beta) * 1.0e8 * scale > std::abs(alpha) && std::abs(beta) > 1.0e-13 * scale)
            zeros.push_back(alpha / beta);
    }
    return zeros;
}

std::vector<double> dcGain(const StateSpaceModel& model) {
    const auto r = realization(model);
    const Eigen::Index n = r.A.rows();
    const Matrix pivot = isDiscrete(model) ? Matrix(Matrix::Identity(n, n) - r.A) : Matrix(-r.A);
    Eigen::FullPivLU<Matrix> lu(pivot);
    if (!lu.isInvertible()) throw std::invalid_argument("DC gain is infinite (pole at s = 0 or z = 1)");
    return flat(r.C * lu.solve(r.B) + r.D);
}

std::vector<ModalDamping> modalDamping(const StateSpaceModel& model) {
    std::vector<ModalDamping> out;
    for (auto pole : systemPoles(model)) {
        if (isDiscrete(model)) pole = std::log(pole) / model.sampleTime;
        ModalDamping d; d.pole = pole; d.naturalFrequency = std::abs(pole);
        d.dampingRatio = d.naturalFrequency > 0.0 ? -pole.real() / d.naturalFrequency : -1.0;
        d.timeConstant = pole.real() < 0.0 ? -1.0 / pole.real() : kInfinity;
        out.push_back(d);
    }
    std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.naturalFrequency < b.naturalFrequency; });
    return out;
}

StabilityMargins stabilityMargins(const StateSpaceModel& model, std::size_t points) {
    if (model.inputs != 1 || model.outputs != 1) throw std::invalid_argument("stability margins require a SISO open loop");
    const TransferEvaluator L(model);
    auto value = [&](double w) { return L.atFrequency(w)(0, 0); };
    auto grid = logarithmicFrequencyGrid(model, std::max<std::size_t>(points, 16));
    StabilityMargins out{kInfinity, kInfinity, kInfinity, kInfinity};
    auto refine = [&](double a, double b, auto&& f) {  // bisection on log frequency
        for (int i = 0; i < 80; ++i) {
            const double c = std::sqrt(a * b);
            if ((f(a) < 0.0) == (f(c) < 0.0)) a = c; else b = c;
        }
        return std::sqrt(a * b);
    };
    auto magnitude = [&](double w) { return std::abs(value(w)) - 1.0; };
    auto imaginary = [&](double w) { return value(w).imag(); };
    for (std::size_t i = 0; i + 1 < grid.size(); ++i) {
        const double a = grid[i], b = grid[i + 1];
        if ((magnitude(a) < 0.0) != (magnitude(b) < 0.0)) {
            const double w = refine(a, b, magnitude);
            double phase = std::arg(value(w)) * 180.0 / kPi + 180.0;
            if (phase > 180.0) phase -= 360.0;
            if (std::abs(phase) < std::abs(out.phaseMargin)) { out.phaseMargin = phase; out.gainCrossover = w; }
        }
        if ((imaginary(a) < 0.0) != (imaginary(b) < 0.0)) {
            const double w = refine(a, b, imaginary);
            const auto l = value(w);
            if (l.real() < 0.0 && 1.0 / std::abs(l) < out.gainMargin) { out.gainMargin = 1.0 / std::abs(l); out.phaseCrossover = w; }
        }
    }
    return out;
}

}  // namespace dynsyslab::core
