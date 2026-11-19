#include <dynsyslab/core/EigenfunctionObservability.hpp>
#include <dynsyslab/core/ControlDesign.hpp>
#include <dynsyslab/core/KoopmanStateSpace.hpp>
#include <dynsyslab/core/ModalObservability.hpp>
#include <dynsyslab/numerics/MatrixFunctions.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
using Complex = std::complex<double>;

std::vector<std::size_t> conjugatePartners(const std::vector<Complex>& lambda) {
    const std::size_t n = lambda.size();
    std::vector<std::size_t> partner(n);
    for (std::size_t i = 0; i < n; ++i) {
        partner[i] = i;
        if (std::abs(lambda[i].imag()) <= 1e-12 * std::max(1.0, std::abs(lambda[i]))) continue;
        double best = 1e300;
        for (std::size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            const double d = std::abs(lambda[j] - std::conj(lambda[i]));
            if (d < best && d <= 1e-8 * std::max(1.0, std::abs(lambda[i]))) { best = d; partner[i] = j; }
        }
    }
    return partner;
}
}  // namespace

EigenfunctionObservability eigenfunctionObservability(const KoopmanModel& model, const std::string& measurement,
                                                      const EigenfunctionObservabilityOptions& o) {
    if (!(o.horizon > 0) || !(o.noiseIntensity > 0) || !(o.regularization >= 0))
        throw std::invalid_argument("eigenfunction observability needs horizon > 0, noise intensity > 0, regularization >= 0");
    ObservabilityOptions g;
    g.method = ObservabilityMethod::continuousFiniteGramian; g.horizon = o.horizon;
    const auto modal = analyzeModalObservability(model, measurement, g);  // G over [0, T]
    const std::size_t n = modal.modes;
    EigenfunctionObservability r;
    r.modes = n; r.horizon = o.horizon; r.noiseIntensity = o.noiseIntensity;
    r.eigenvalues = model.continuousEigenvalues();
    r.partner = conjugatePartners(r.eigenvalues);
    r.gramian = modal.gramian;
    // regularized inverse through the Hermitian eigensystem: G^-1 = U diag(1 / (l + d)) U*
    const auto eig = numerics::hermitianEigensystem(modal.gramian, n);
    const double top = *std::max_element(eig.eigenvalues.begin(), eig.eigenvalues.end());
    const double floor = o.regularization * std::max(top, 0.0);
    std::vector<Complex> inverse(n * n, Complex{});
    for (std::size_t k = 0; k < n; ++k) {
        const double l = std::max(eig.eigenvalues[k], 0.0) + floor;
        if (!(l > 0)) continue;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                inverse[i * n + j] += eig.eigenvectors[i * n + k] * std::conj(eig.eigenvectors[j * n + k]) / l;
    }
    for (std::size_t i = 0; i < n; ++i) {
        const double gii = std::max(modal.gramian[i * n + i].real(), 0.0);
        const double hii = std::max(inverse[i * n + i].real(), 0.0);
        r.signal.push_back(std::sqrt(gii));
        r.distinguishable.push_back(hii > 0 ? 1.0 / std::sqrt(hii) : 0.0);
        r.isolation.push_back(gii > 0 ? std::min(1.0, r.distinguishable.back() / std::sqrt(gii)) : 0.0);
        r.boundStd.push_back(std::sqrt(o.noiseIntensity * hii));
        const std::size_t p = r.partner[i];
        if (p == i) { r.pairBoundStd.push_back(r.boundStd.back()); continue; }
        const double a = inverse[i * n + i].real(), d = inverse[p * n + p].real();
        const double b = std::abs(inverse[i * n + p]);
        const double top2 = 0.5 * (a + d) + std::sqrt(0.25 * (a - d) * (a - d) + b * b);
        r.pairBoundStd.push_back(std::sqrt(o.noiseIntensity * std::max(top2, 0.0)));
    }
    if (o.processCovariance > 0) {  // the KLQE's own steady-state uncertainty per eigenfunction
        KoopmanStateSpaceOptions so; so.measurements = {measurement};
        const auto plant = koopmanStateSpace(model, so);
        const std::size_t ns = plant.states;
        std::vector<double> W(ns * ns, 0.0), V(plant.outputs * plant.outputs, 0.0);
        for (std::size_t k = 0; k < ns; ++k) W[k * ns + k] = o.processCovariance;
        for (std::size_t k = 0; k < plant.outputs; ++k) V[k * plant.outputs + k] = o.noiseIntensity;
        const auto lqe = designLQE(plant, W, V);
        // Split coordinates hold psi_m and its partner independently; their conjugate-inconsistent
        // combinations are unobservable by construction. Use the conjugate-consistent estimate
        // (psi_m + conj(psi_p)) / 2, the part a physical state determines.
        const auto splitOf = [&](std::size_t m, Complex value) {
            std::vector<Complex> unit(n, Complex{}); unit[m] = value;
            return model.splitValues(unit);
        };
        const auto variance = [&](const std::vector<double>& c) {
            double v = 0;
            for (std::size_t i = 0; i < ns; ++i) for (std::size_t j = 0; j < ns; ++j) v += c[i] * lqe.covariance[i * ns + j] * c[j];
            return std::max(v, 0.0);
        };
        for (std::size_t m = 0; m < n; ++m) {
            const std::size_t p = r.partner[m];
            // real part: (Re psi_m + Re psi_p) / 2; imaginary part: (Im psi_m - Im psi_p) / 2
            auto re = splitOf(m, Complex(1.0, 0.0)), im = splitOf(m, Complex(0.0, 1.0));
            if (p != m) {
                const auto rp = splitOf(p, Complex(1.0, 0.0)), ip = splitOf(p, Complex(0.0, 1.0));
                for (std::size_t s = 0; s < ns; ++s) { re[s] = 0.5 * (re[s] + rp[s]); im[s] = 0.5 * (im[s] - ip[s]); }
            }
            r.kalmanStd.push_back(std::sqrt(variance(re) + variance(im)));
        }
    }
    return r;
}

}  // namespace dynsyslab::core
