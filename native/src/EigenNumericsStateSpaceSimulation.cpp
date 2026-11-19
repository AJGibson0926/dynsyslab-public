#include "EigenNumericsLinearSystems.hpp"
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <unsupported/Eigen/MatrixFunctions>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
using namespace linear_detail;
namespace {

bool isDiscrete(const StateSpaceModel& m) { return m.domain == TimeDomain::discrete; }

struct Quadrature { Matrix flow; std::vector<double> weights; std::size_t samples{}; };

// Exact one-step flow and trapezoid (continuous) or unit (discrete) weights.
Quadrature quadrature(const StateSpaceModel& model, const Matrix& A,
                      const EmpiricalGramianOptions& options) {
    if (!(options.horizon > 0.0) || !std::isfinite(options.horizon))
        throw std::invalid_argument("empirical gramian horizon must be positive");
    if (options.scales.empty()) throw std::invalid_argument("empirical gramian needs perturbation scales");
    for (double c : options.scales)
        if (!(c > 0.0) || !std::isfinite(c)) throw std::invalid_argument("empirical gramian scales must be positive");
    Quadrature q;
    if (isDiscrete(model)) {
        if (std::floor(options.horizon) != options.horizon || options.horizon > 1.0e7)
            throw std::invalid_argument("discrete empirical horizon must be a bounded step count");
        q.samples = static_cast<std::size_t>(options.horizon);
        q.flow = A; q.weights.assign(q.samples, 1.0);
        return q;
    }
    if (options.samples < 2) throw std::invalid_argument("empirical gramian needs at least two samples");
    q.samples = options.samples;
    const double dt = options.horizon / static_cast<double>(q.samples - 1);
    q.flow = (A * dt).exp();
    q.weights.assign(q.samples, dt);
    q.weights.front() *= 0.5; q.weights.back() *= 0.5;
    return q;
}

}  // namespace

std::vector<double> empiricalControllabilityGramian(const StateSpaceModel& model,
    const EmpiricalGramianOptions& options) {
    const auto r = realization(model);
    const auto q = quadrature(model, r.A, options);
    const Eigen::Index n = r.A.rows();
    Matrix W = Matrix::Zero(n, n);
    for (double c : options.scales)
        for (Eigen::Index input = 0; input < r.B.cols(); ++input) {
            Vector x = c * r.B.col(input);  // impulse of magnitude c on one input
            for (std::size_t k = 0; k < q.samples; ++k) {
                W += (q.weights[k] / (c * c)) * x * x.transpose();
                x = q.flow * x;
            }
        }
    W /= static_cast<double>(options.scales.size());
    return flat(0.5 * (W + W.transpose()));
}

std::vector<double> empiricalObservabilityGramian(const StateSpaceModel& model,
    const EmpiricalGramianOptions& options) {
    const auto r = realization(model);
    const auto q = quadrature(model, r.A, options);
    const Eigen::Index n = r.A.rows();
    Matrix W = Matrix::Zero(n, n);
    for (double c : options.scales) {
        // Column j holds y(+c e_j) - y(-c e_j) at the current sample.
        Matrix plus = c * Matrix::Identity(n, n), minus = -plus;
        for (std::size_t k = 0; k < q.samples; ++k) {
            const Matrix difference = r.C * plus - r.C * minus;  // p x n
            W += (q.weights[k] / (4.0 * c * c)) * difference.transpose() * difference;
            plus = q.flow * plus; minus = q.flow * minus;
        }
    }
    W /= static_cast<double>(options.scales.size());
    return flat(0.5 * (W + W.transpose()));
}

std::vector<double> snapshotGramian(const std::vector<double>& snapshots,
    std::size_t dimension, std::size_t samples, const std::vector<double>& weights) {
    const Matrix X = mapped(snapshots, dimension, samples, "snapshot matrix");
    if (!weights.empty() && weights.size() != samples)
        throw std::invalid_argument("snapshot weights must match the sample count");
    Vector w = Vector::Ones(static_cast<Eigen::Index>(samples));
    for (std::size_t k = 0; k < weights.size(); ++k) {
        if (!(weights[k] >= 0.0) || !std::isfinite(weights[k]))
            throw std::invalid_argument("snapshot weights must be nonnegative");
        w[static_cast<Eigen::Index>(k)] = weights[k];
    }
    const Matrix W = X * w.asDiagonal() * X.transpose();
    return flat(0.5 * (W + W.transpose()));
}

StateSpaceModel discretize(const StateSpaceModel& model, double h, DiscretizationMethod method) {
    if (isDiscrete(model)) throw std::invalid_argument("model is already discrete");
    if (!(h > 0.0) || !std::isfinite(h)) throw std::invalid_argument("sample time must be positive");
    const auto r = realization(model);
    const Eigen::Index n = r.A.rows(), m = r.B.cols();
    const Matrix I = Matrix::Identity(n, n);
    Realization d{r.A, r.B, r.C, r.D};
    if (method == DiscretizationMethod::zeroOrderHold) {
        // The numerics zero-order hold (one block exponential); no inputs: e^{A h}.
        const auto states = static_cast<std::size_t>(n), inputs = static_cast<std::size_t>(m);
        if (inputs) {
            const auto held = numerics::zeroOrderHold(model.A, model.B, states, inputs, h);
            d.A = mapped(held.transition, states, states, "ZOH transition");
            d.B = mapped(held.response, states, inputs, "ZOH response");
        } else {
            d.A = mapped(numerics::matrixExponential(model.A, states, h), states, states, "ZOH transition");
        }
    } else if (method == DiscretizationMethod::tustin) {
        Eigen::FullPivLU<Matrix> lu(I - 0.5 * h * r.A);
        if (!lu.isInvertible()) throw std::invalid_argument("Tustin map is singular (eigenvalue 2/h)");
        const Matrix M = lu.inverse();
        d.A = M * (I + 0.5 * h * r.A); d.B = h * M * r.B;
        d.C = r.C * M; d.D = r.D + 0.5 * h * r.C * M * r.B;
    } else {
        d.A = I + h * r.A; d.B = h * r.B;
    }
    return fromRealization(d, TimeDomain::discrete, h);
}

StateSpaceModel continuousFromTustin(const StateSpaceModel& model) {
    if (!isDiscrete(model)) throw std::invalid_argument("model is not discrete");
    const auto r = realization(model);
    const double h = model.sampleTime;
    const Eigen::Index n = r.A.rows();
    Eigen::FullPivLU<Matrix> lu(r.A + Matrix::Identity(n, n));
    if (!lu.isInvertible()) throw std::invalid_argument("inverse Tustin map is singular (eigenvalue -1)");
    const Matrix inverse = lu.inverse();
    Realization c;
    c.A = (2.0 / h) * (r.A - Matrix::Identity(n, n)) * inverse;
    c.B = (2.0 / h) * inverse * r.B;
    c.C = 2.0 * r.C * inverse;
    c.D = r.D - r.C * inverse * r.B;
    return fromRealization(c, TimeDomain::continuous, 0.0);
}

LinearResponse simulateLinear(const StateSpaceModel& model, const std::vector<double>& time,
    const std::vector<double>& inputs, const std::vector<double>& initialState) {
    const auto r = realization(model);
    const std::size_t N = time.size(), n = model.states, m = model.inputs, p = model.outputs;
    if (N == 0) throw std::invalid_argument("simulation needs at least one sample");
    if (inputs.size() != N * m) throw std::invalid_argument("simulation inputs must be samples x inputs");
    if (!initialState.empty() && initialState.size() != n)
        throw std::invalid_argument("initial state dimension mismatch");
    Matrix Ad = r.A, Bd = r.B;
    if (!isDiscrete(model) && N > 1) {
        const double dt = time[1] - time[0];
        for (std::size_t k = 1; k < N; ++k)
            if (std::abs((time[k] - time[k - 1]) - dt) > 1.0e-9 * std::max(std::abs(dt), 1.0))
                throw std::invalid_argument("continuous simulation requires a uniform time grid");
        const auto d = realization(discretize(model, dt));
        Ad = d.A; Bd = d.B;
    }
    Vector x = initialState.empty() ? Vector::Zero(static_cast<Eigen::Index>(n))
        : Vector(Eigen::Map<const Vector>(initialState.data(), static_cast<Eigen::Index>(n)));
    LinearResponse out; out.time = time; out.samples = N;
    out.outputs.resize(N * p); out.states.resize(N * n);
    for (std::size_t k = 0; k < N; ++k) {
        const Vector u = m ? Vector(Eigen::Map<const Vector>(inputs.data() + k * m, static_cast<Eigen::Index>(m)))
                           : Vector(0);
        const Vector y = r.C * x + r.D * u;
        for (std::size_t j = 0; j < p; ++j) out.outputs[k * p + j] = y[static_cast<Eigen::Index>(j)];
        for (std::size_t j = 0; j < n; ++j) out.states[k * n + j] = x[static_cast<Eigen::Index>(j)];
        x = Ad * x + Bd * u;
    }
    return out;
}

namespace {
std::vector<double> grid(const StateSpaceModel& model, double horizon, std::size_t samples) {
    if (samples < 2 || (!isDiscrete(model) && !(horizon > 0.0)))
        throw std::invalid_argument("response needs two samples and (continuous) a positive horizon");
    std::vector<double> time(samples);
    const double dt = isDiscrete(model) ? model.sampleTime : horizon / static_cast<double>(samples - 1);
    for (std::size_t k = 0; k < samples; ++k) time[k] = dt * static_cast<double>(k);
    return time;
}
}  // namespace

LinearResponse impulseResponse(const StateSpaceModel& model, std::size_t input,
                               double horizon, std::size_t samples) {
    if (input >= model.inputs) throw std::invalid_argument("impulse input index out of range");
    const auto time = grid(model, horizon, samples);
    std::vector<double> u(samples * model.inputs, 0.0);
    if (isDiscrete(model)) { u[input] = 1.0; return simulateLinear(model, time, u); }
    // Continuous unit impulse: x(0+) = B e_i; the D*delta term is omitted.
    std::vector<double> x0(model.states);
    for (std::size_t i = 0; i < model.states; ++i) x0[i] = model.B[i * model.inputs + input];
    return simulateLinear(model, time, u, x0);
}

LinearResponse stepResponse(const StateSpaceModel& model, std::size_t input,
                            double horizon, std::size_t samples) {
    if (input >= model.inputs) throw std::invalid_argument("step input index out of range");
    const auto time = grid(model, horizon, samples);
    std::vector<double> u(samples * model.inputs, 0.0);
    for (std::size_t k = 0; k < samples; ++k) u[k * model.inputs + input] = 1.0;
    return simulateLinear(model, time, u);
}

namespace {
Matrix matrixSign(Matrix X) {
    const Eigen::Index n = X.rows();
    for (int iteration = 0; iteration < 100; ++iteration) {
        Eigen::PartialPivLU<Matrix> lu(X);
        const double determinant = std::abs(lu.determinant());
        if (!(determinant > 0.0) || !std::isfinite(determinant))
            throw std::invalid_argument("stable/unstable split: eigenvalue on the separating boundary");
        const double c = iteration < 20 ? std::pow(determinant, -1.0 / static_cast<double>(n)) : 1.0;
        const Matrix next = 0.5 * (c * X + lu.inverse() / c);
        const double change = (next - X).norm();
        X = next;
        if (!X.allFinite()) throw std::runtime_error("matrix sign iteration diverged");
        if (change <= 1.0e-13 * X.norm()) return X;
    }
    throw std::runtime_error("matrix sign iteration did not converge");
}
}  // namespace

StableUnstableSplit splitStableUnstable(const StateSpaceModel& model, double margin) {
    const auto r = realization(model);
    const Eigen::Index n = r.A.rows();
    const Matrix I = Matrix::Identity(n, n);
    if (!(margin >= 0.0) || margin >= 1.0) throw std::invalid_argument("split margin must lie in [0, 1)");
    Matrix shifted;
    if (isDiscrete(model)) {
        const double radius = 1.0 - margin;
        Eigen::FullPivLU<Matrix> lu(r.A + radius * I);
        if (!lu.isInvertible()) throw std::invalid_argument("stable/unstable split: eigenvalue at -radius");
        shifted = (r.A - radius * I) * lu.inverse();  // |mu| < radius <=> Re(w) < 0
    } else shifted = r.A + margin * I;
    const Matrix S = matrixSign(shifted);
    const Matrix Ps = 0.5 * (I - S), Pu = 0.5 * (I + S);
    const auto ks = static_cast<Eigen::Index>(std::llround(Ps.trace()));
    if (ks < 0 || ks > n) throw std::runtime_error("stable projector has invalid rank");
    Matrix T(n, n);
    if (ks > 0) { Eigen::JacobiSVD<Matrix> svd(Ps, Eigen::ComputeThinU); T.leftCols(ks) = svd.matrixU().leftCols(ks); }
    if (ks < n) { Eigen::JacobiSVD<Matrix> svd(Pu, Eigen::ComputeThinU); T.rightCols(n - ks) = svd.matrixU().leftCols(n - ks); }
    Eigen::FullPivLU<Matrix> lu(T);
    if (!lu.isInvertible()) throw std::runtime_error("stable and unstable subspaces are not complementary");
    const Matrix Ti = lu.inverse();
    Matrix At = Ti * r.A * T;
    const double coupling = ks > 0 && ks < n
        ? std::max(At.topRightCorner(ks, n - ks).norm(), At.bottomLeftCorner(n - ks, ks).norm()) : 0.0;
    if (coupling > 1.0e-7 * std::max(1.0, r.A.norm()) * std::max(1.0, Ti.norm() * T.norm()))
        throw std::runtime_error("stable/unstable split is ill-conditioned (residual coupling)");
    const Matrix Bt = Ti * r.B, Ct = r.C * T;
    StableUnstableSplit out;
    out.stableStates = static_cast<std::size_t>(ks); out.unstableStates = static_cast<std::size_t>(n - ks);
    out.margin = margin; out.transform = flat(T); out.inverse = flat(Ti);
    const Eigen::Index p = r.C.rows(), m = r.B.cols(), ku = n - ks;
    auto part = [&](Eigen::Index start, Eigen::Index size, bool withD) {
        StateSpaceModel s;
        s.states = static_cast<std::size_t>(size); s.inputs = r.B.cols(); s.outputs = r.C.rows();
        s.A = flat(At.block(start, start, size, size)); s.B = flat(Bt.middleRows(start, size));
        s.C = flat(Ct.middleCols(start, size));
        s.D = withD ? flat(r.D) : std::vector<double>(static_cast<std::size_t>(p * m), 0.0);
        s.domain = model.domain; s.sampleTime = model.sampleTime;
        return s;
    };
    out.stable = part(0, ks, true);
    out.unstable = part(ks, ku, false);
    return out;
}

}  // namespace dynsyslab::core
