#include <algorithm>
#include <cmath>
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>
#include <stdexcept>

namespace dynsyslab::core {
std::vector<std::complex<double>> KoopmanModel::discreteEigenvalues(double dt) const {
    return spectrum().discreteMultipliers(dt);
}
KoopmanModel KoopmanModel::fromPolynomialCoefficients(
    SystemConfiguration system, std::size_t dimension, std::vector<std::size_t> exponents,
    std::vector<std::complex<double>> coefficients, std::vector<std::complex<double>> eigenvalues,
    std::vector<double> center, std::vector<double> scale) {
    return KoopmanModel(polynomialEigenfunctions(dimension, std::move(exponents), std::move(coefficients),
                                                 Spectrum::fromContinuousRates(std::move(eigenvalues)),
                                                 std::move(center), std::move(scale)),
                        std::move(system));
}
namespace {
void requireEigenfunctions(const FunctionExpansion& functions, const Spectrum& spectrum) {
    if (functions.outputCount() != spectrum.size())
        throw std::invalid_argument("Expansion outputs and spectrum size differ");
    const auto n = functions.dictionary().size();
    const auto& c = functions.coefficients();
    for (std::size_t j = 0; j < spectrum.size(); ++j)
        if (std::all_of(c.begin() + std::ptrdiff_t(j * n), c.begin() + std::ptrdiff_t((j + 1) * n),
                        [](auto z) { return z == std::complex<double>{}; }))
            throw std::invalid_argument("Identically zero expansion is not an eigenfunction candidate");
}
}  // namespace
KoopmanModel::KoopmanModel(EigenfunctionFitResult fit, std::size_t inputs)
    : KoopmanModel(std::move(fit.functions), std::move(fit.spectrum), inputs, std::move(fit.diagnostics)) {}
KoopmanModel::KoopmanModel(EigenfunctionFitResult fit, SystemConfiguration system)
    : KoopmanModel(std::move(fit.functions), std::move(fit.spectrum), std::move(system), std::move(fit.diagnostics)) {}
KoopmanModel::KoopmanModel(FunctionExpansion p, Spectrum s, std::size_t inputs,
                           EigenfunctionFitDiagnostics diagnostics)
    : polynomial_(std::move(p)), spectrum_(std::move(s)), fitDiagnostics_(std::move(diagnostics)),
      suppliedInputDimension_(inputs) {
    requireEigenfunctions(*polynomial_, *spectrum_);
    const auto modes = eigenfunctionCount();
    if (!fitDiagnostics_.residualMaximum.empty() && fitDiagnostics_.residualMaximum.size() != modes)
        throw std::invalid_argument("Fit diagnostics must be mode aligned");
    coordinates().validate();  // the state map (if any) matches the coordinate dimension
    samples_.modes = modes;
    representations_.push_back({representationName(*polynomial_), std::make_shared<const FunctionExpansion>(*polynomial_)});
}
KoopmanModel::KoopmanModel(FunctionExpansion p, Spectrum s, SystemConfiguration system,
                           EigenfunctionFitDiagnostics diagnostics)
    : KoopmanModel(std::move(p), std::move(s), 0, std::move(diagnostics)) {
    attachPhysicalSystem(std::move(system));
}
void KoopmanModel::attachPhysicalSystem(SystemConfiguration configuration) {
    system_ = std::shared_ptr<ODESystem>(createSystem(configuration).release());
    systemConfiguration_ = system_->configuration();
    const auto physical = system_->descriptor().stateDimension;
    if (physicalStateDimension() && *physicalStateDimension() != physical)
        throw std::invalid_argument(
            "KROM coordinate contract and physical system differ in dimension");
    auto q = coordinates().bound(system_->observableCatalog());
    if (q.stateMap && (q.stateMap->stateDimension() != physical || q.stateMap->outputDimension() != q.dimension))
        throw std::invalid_argument("KROM observable and polynomial dimensions differ");
    if (!q.stateMap) q.physicalStateDimension = physical;
    polynomial_ = polynomial_->withCoordinates(std::move(q));
    if (!representations_.empty()) representations_.front().function = std::make_shared<const FunctionExpansion>(*polynomial_);
}
const FunctionExpansion &KoopmanModel::expansion() const {
    if (!polynomial_)
        throw std::logic_error("Scattered eigenfunctions have no polynomial dictionary");
    return *polynomial_;
}
const Spectrum &KoopmanModel::spectrum() const noexcept {
    return *spectrum_;
}
const FunctionDomain &KoopmanModel::coordinates() const noexcept {
    return polynomial_ ? polynomial_->coordinates() : scattered_->coordinates();
}
const EigenfunctionFitDiagnostics &KoopmanModel::fitDiagnostics() const noexcept {
    return fitDiagnostics_;
}
std::optional<std::size_t> KoopmanModel::physicalStateDimension() const noexcept {
    const auto& q = coordinates();
    return q.stateMap ? std::optional<std::size_t>(q.stateMap->stateDimension()) : q.physicalStateDimension;
}
bool KoopmanModel::isRealMode(std::size_t j) const noexcept {
    const auto rate = spectrum().continuousRate(j);
    return std::abs((rate ? *rate : spectrum().values()[j]).imag()) <= 1e-12;
}
void KoopmanModel::setDiscreteInputMap(std::string name, DiscreteInputMap map) {
    if (name.empty() || map.eigenfunctionCount != eigenfunctionCount() ||
        map.inputDimension != inputDimension() || !map.inputDimension ||
        map.matrix.size() != map.eigenfunctionCount * map.inputDimension ||
        !std::isfinite(map.sampleInterval) || map.sampleInterval <= 0)
        throw std::invalid_argument("Discrete input map dimensions/interval do not match KROM");
    if (spectrum().sampleInterval() > 0 &&
        std::abs(map.sampleInterval - spectrum().sampleInterval()) >
            1e-12 * std::max(map.sampleInterval, spectrum().sampleInterval()))
        throw std::invalid_argument("Discrete input map interval differs from spectrum");
    for (auto c : map.matrix)
        if (!std::isfinite(c.real()) || !std::isfinite(c.imag()))
            throw std::invalid_argument("Discrete input map must be finite");
    const auto identity=koopmanModelIdentity(*this);
    if(!map.modelIdentity.empty()&&map.modelIdentity!=identity)
        throw std::invalid_argument("Discrete input map belongs to a different KROM");
    map.modelIdentity=identity;
    inherited_.discreteMaps.erase(name);
    discreteInputMaps_[std::move(name)] = std::move(map);
}
} // namespace dynsyslab::core
