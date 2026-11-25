// Smooth polynomial representation of scattered (H-DMD) eigenfunctions.
//
// H-DMD eigenfunctions are known only at training states; between them a K-ROM
// interpolates by inverse distance over the nearest neighbours, which is continuous but
// rough (its slope jumps whenever the neighbour set changes). This fit keeps the same
// eigenvalues and the same eigenfunctions -- least squares on the training values in the
// normalized physical state -- and writes a polynomial K-ROM, so estimators and
// controllers evaluate a smooth surface. Each conjugate partner keeps its own stored
// normalization; mode_selection's K-ROM writer later makes it exactly conj(upper).
#include "PostprocessInternal.hpp"

#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/FunctionInterpolation.hpp>
#include <dynsyslab/core/Text.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>

#include <cmath>
#include <functional>
#include <stdexcept>

namespace dynsyslab::workflow::postprocess_detail {
namespace c = dynsyslab::core;

std::filesystem::path fitPolynomialEigenfunctions(const std::filesystem::path& kromPath, std::size_t degree,
                                                  const std::filesystem::path& root, json::Value::Object& report) {
    if (degree < 1 || degree > 12) throw std::invalid_argument("eigenfunction_polynomial_degree must be 1..12");
    const auto model = c::loadKoopmanModel(kromPath);
    const auto* scattered = model.scatteredEigenfunctions();
    if (!scattered) throw std::invalid_argument("eigenfunction_polynomial_degree needs a K-ROM with scattered (H-DMD) eigenfunctions");
    const std::size_t m = scattered->samples, d = scattered->stateDimension, modes = scattered->outputs();
    // Normalized coordinates (x - center) / scale keep high powers well conditioned.
    std::vector<double> center(d, 0.0), scale(d, 0.0);
    for (std::size_t i = 0; i < m; ++i) for (std::size_t j = 0; j < d; ++j) center[j] += scattered->trainingStates[i * d + j] / double(m);
    for (std::size_t i = 0; i < m; ++i) for (std::size_t j = 0; j < d; ++j)
        scale[j] += std::pow(scattered->trainingStates[i * d + j] - center[j], 2) / double(m);
    for (auto& s : scale) s = std::sqrt(s) > 1e-14 ? std::sqrt(s) : 1.0;
    // All monomials of total degree <= degree, feature-major exponents.
    std::vector<std::size_t> exponents, current(d, 0);
    std::function<void(std::size_t, std::size_t)> enumerate = [&](std::size_t axis, std::size_t left) {
        if (axis + 1 == d) { current[axis] = left; exponents.insert(exponents.end(), current.begin(), current.end()); return; }
        for (std::size_t p = 0; p <= left; ++p) { current[axis] = p; enumerate(axis + 1, left - p); }
    };
    for (std::size_t total = 0; total <= degree; ++total) enumerate(0, total);
    const std::size_t features = exponents.size() / d;
    if (features >= m) throw std::invalid_argument("eigenfunction_polynomial_degree has more terms than training samples");
    std::vector<std::complex<double>> design(m * features);
    for (std::size_t i = 0; i < m; ++i) for (std::size_t f = 0; f < features; ++f) {
        double v = 1.0;
        for (std::size_t j = 0; j < d; ++j)
            v *= std::pow((scattered->trainingStates[i * d + j] - center[j]) / scale[j], double(exponents[f * d + j]));
        design[i * features + f] = v;
    }
    const auto solved = numerics::solveLeastSquares(design, m, features, scattered->trainingValues, modes);  // features x modes
    std::vector<std::complex<double>> coefficients(modes * features);
    json::Value::Array errors;
    for (std::size_t k = 0; k < modes; ++k) {
        double residual = 0.0, norm = 0.0;
        for (std::size_t f = 0; f < features; ++f) coefficients[k * features + f] = solved[f * modes + k];
        for (std::size_t i = 0; i < m; ++i) {
            std::complex<double> fit{};
            for (std::size_t f = 0; f < features; ++f) fit += design[i * features + f] * solved[f * modes + k];
            residual += std::norm(fit - scattered->trainingValues[i * modes + k]);
            norm += std::norm(scattered->trainingValues[i * modes + k]);
        }
        errors.emplace_back(std::sqrt(residual / std::max(norm, 1e-300)));
    }
    const auto smooth = c::KoopmanModel::fromPolynomialCoefficients(model.systemConfiguration(), d, exponents,
        coefficients, model.continuousEigenvalues(), center, scale);
    std::map<std::string, c::NpyArray> extra;
    const auto parent = c::readNpz(kromPath);
    for (const char* name : {"training_coordinates", "training_time", "training_states"})
        if (parent.contains(name)) extra[name] = parent.at(name);
    const auto directory = objectFolder(root, "koopman-polynomial-fit");
    std::filesystem::create_directories(directory);
    const auto path = directory / c::wide(safeName(kromPath.stem().u8string() + "__polynomial-fit-degree-"
        + std::to_string(degree) + "__" + timestamp()) + ".npz");
    report = {{"source_krom", kromPath.generic_u8string()}, {"degree", double(degree)}, {"terms", double(features)},
              {"training_samples", double(m)}, {"relative_fit_error", std::move(errors)},
              {"representation", "polynomial in (x - center) / scale, least squares on the scattered training values"}};
    json::Value::Object provenance{{"format", "dynsyslab.derived/v1"}, {"method", "eigenfunction-polynomial-fit"},
        {"label", "polynomial fit of " + kromPath.stem().u8string()}, {"source_krom", kromPath.generic_u8string()},
        {"source_archive", kromPath.generic_u8string()}, {"independent_archive", true},
        {"eigenfunction_polynomial_fit", json::Value(report)}};
    c::saveKoopmanModel(path, smooth, json::dump(json::Value(std::move(provenance))), std::move(extra));
    return path;
}

}  // namespace dynsyslab::workflow::postprocess_detail
