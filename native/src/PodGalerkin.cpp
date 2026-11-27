#include <dynsyslab/core/PodGalerkin.hpp>

#include <dynsyslab/numerics/Svd.hpp>

#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

std::vector<double> PodGalerkinModel::project(const double* state) const {
    std::vector<double> z(rank, 0.0);
    for (std::size_t i = 0; i < states; ++i)
        for (std::size_t k = 0; k < rank; ++k) z[k] += modes[i * rank + k] * (state[i] - reference[i]);
    return z;
}

std::vector<double> PodGalerkinModel::reconstruct(const double* z) const {
    auto x = reference;
    for (std::size_t i = 0; i < states; ++i)
        for (std::size_t k = 0; k < rank; ++k) x[i] += modes[i * rank + k] * z[k];
    return x;
}

std::vector<double> PodGalerkinModel::reducedRhs(const ODESystem& system, double time, const double* z,
                                                 double input) const {
    const auto x = reconstruct(z);
    std::vector<double> dx(states);
    if (input != 0.0) system.rhsWithInput(time, x.data(), input, dx.data());
    else system.rhs(time, x.data(), dx.data());
    std::vector<double> dz(rank, 0.0);
    for (std::size_t i = 0; i < states; ++i)
        for (std::size_t k = 0; k < rank; ++k) dz[k] += modes[i * rank + k] * dx[i];
    return dz;
}

PodGalerkinModel fitPodGalerkin(const std::vector<double>& snapshots, std::size_t samples, std::size_t n,
                                std::size_t rank, double energy, const std::string& reference,
                                const ODESystem& system) {
    if (snapshots.size() != samples * n || samples < 2) throw std::invalid_argument("POD-Galerkin needs at least two snapshots");
    if (reference != "mean" && reference != "zero") throw std::invalid_argument("POD-Galerkin reference must be mean or zero");
    PodGalerkinModel model;
    model.states = n;
    model.reference.assign(n, 0.0);
    if (reference == "mean")
        for (std::size_t s = 0; s < samples; ++s)
            for (std::size_t i = 0; i < n; ++i) model.reference[i] += snapshots[s * n + i] / double(samples);
    // Snapshot matrix n x samples (columns are states minus the reference).
    std::vector<double> matrix(n * samples);
    for (std::size_t s = 0; s < samples; ++s)
        for (std::size_t i = 0; i < n; ++i) matrix[i * samples + s] = snapshots[s * n + i] - model.reference[i];
    const auto svd = numerics::singularValueDecomposition(matrix, n, samples);
    model.singularValues = svd.singularValues;
    double total = 0.0;
    for (double s : svd.singularValues) total += s * s;
    if (!(total > 0.0)) throw std::invalid_argument("POD-Galerkin snapshots do not vary about the reference");
    if (!rank) {
        if (!(energy > 0.0 && energy <= 1.0)) throw std::invalid_argument("POD-Galerkin energy must lie in (0, 1]");
        double kept = 0.0;
        while (rank < svd.rank && kept < energy * total) kept += svd.singularValues[rank] * svd.singularValues[rank], ++rank;
    }
    if (rank > svd.rank) throw std::invalid_argument("POD-Galerkin rank exceeds the snapshot rank");
    model.rank = rank;
    model.modes.assign(n * rank, 0.0);
    double kept = 0.0;
    for (std::size_t k = 0; k < rank; ++k) kept += svd.singularValues[k] * svd.singularValues[k];
    model.energyFraction = kept / total;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < rank; ++k) model.modes[i * rank + k] = svd.left[i * svd.rank + k];
    const auto J = system.jacobian(0.0, model.reference.data());
    model.A.assign(rank * rank, 0.0);
    for (std::size_t a = 0; a < rank; ++a)
        for (std::size_t b = 0; b < rank; ++b)
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    model.A[a * rank + b] += model.modes[i * rank + a] * J[i * n + j] * model.modes[j * rank + b];
    const auto input = system.inputMatrix(model.reference.data());
    if (!input.empty() && input.size() % n == 0) {
        model.inputs = input.size() / n;
        model.B.assign(rank * model.inputs, 0.0);
        for (std::size_t a = 0; a < rank; ++a)
            for (std::size_t m = 0; m < model.inputs; ++m)
                for (std::size_t i = 0; i < n; ++i) model.B[a * model.inputs + m] += model.modes[i * rank + a] * input[i * model.inputs + m];
    }
    return model;
}

}  // namespace dynsyslab::core
