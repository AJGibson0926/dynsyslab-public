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
    return KoopmanModel(KoopmanSINDyResult::fromCoefficients(
                            dimension, std::move(exponents), std::move(coefficients),
                            std::move(eigenvalues), std::move(center), std::move(scale)),
                        std::move(system));
}
KoopmanModel::KoopmanModel(PolynomialEigenfunctionSet p, std::size_t inputs,
                           EigenfunctionFitDiagnostics diagnostics)
    : polynomial_(std::move(p)), fitDiagnostics_(std::move(diagnostics)),
      suppliedInputDimension_(inputs) {
    const auto modes = eigenfunctionCount();
    if ((!fitDiagnostics_.residualRms.empty() && fitDiagnostics_.residualRms.size() != modes) ||
        (!fitDiagnostics_.residualMaximum.empty() &&
         fitDiagnostics_.residualMaximum.size() != modes))
        throw std::invalid_argument("Fit diagnostics must be mode aligned");
    if (spectrum().hasContinuousRates())
        eigenfunctions_ = legacyPolynomialEigenfunctions(*polynomial_, fitDiagnostics_);
    const auto &q = coordinates();
    if (q.kind == EigenfunctionCoordinateKind::analyticObservable &&
        !q.expressionDefinition.empty())
        coordinateObservable_ = std::make_shared<Observable>(
            Observable::fromExpressionDefinition(json::parse(q.expressionDefinition)));
    if (coordinateObservable_ &&
        (coordinateObservable_->outputDimension() != coordinateDimension() ||
         (physicalStateDimension() &&
          coordinateObservable_->stateDimension() != *physicalStateDimension())))
        throw std::invalid_argument("Analytic coordinate definition dimensions differ from KROM");
    if(coordinateObservable_&&!physicalStateDimension()) {
        auto resolved=coordinates();resolved.physicalStateDimension=coordinateObservable_->stateDimension();
        polynomial_=polynomial_->withCoordinates(std::move(resolved));
        if(spectrum().hasContinuousRates())eigenfunctions_=legacyPolynomialEigenfunctions(*polynomial_,fitDiagnostics_);
    }
}
KoopmanModel::KoopmanModel(PolynomialEigenfunctionSet p, SystemConfiguration system,
                           EigenfunctionFitDiagnostics diagnostics)
    : KoopmanModel(std::move(p), 0, std::move(diagnostics)) {
    attachPhysicalSystem(std::move(system));
}
void KoopmanModel::attachPhysicalSystem(SystemConfiguration configuration) {
    system_ = std::shared_ptr<ODESystem>(createSystem(configuration).release());
    systemConfiguration_ = system_->configuration();
    auto q = coordinates();
    const auto physical = system_->descriptor().stateDimension;
    if (q.physicalStateDimension && *q.physicalStateDimension != physical)
        throw std::invalid_argument(
            "KROM coordinate contract and physical system differ in dimension");
    q.physicalStateDimension = physical;
    polynomial_ = polynomial_->withCoordinates(q);
    if (!q.requiresExternalValues()) {
        const auto catalog = system_->observableCatalog();
        auto observable =
            q.kind == EigenfunctionCoordinateKind::physicalState ? catalog.at("state")
            : !q.expressionDefinition.empty()
                ? Observable::fromExpressionDefinition(json::parse(q.expressionDefinition))
                : catalog.coordinates(q.names, "krom_coordinates");
        if (observable.stateDimension() != physical || observable.outputDimension() != q.dimension)
            throw std::invalid_argument("KROM observable and polynomial dimensions differ");
        coordinateObservable_ = std::make_shared<Observable>(std::move(observable));
    }
    if (spectrum().hasContinuousRates())
        eigenfunctions_ = legacyPolynomialEigenfunctions(*polynomial_, fitDiagnostics_);
}
const PolynomialEigenfunctionSet &KoopmanModel::polynomialEigenfunctions() const {
    if (!polynomial_)
        throw std::logic_error("Scattered eigenfunctions have no polynomial dictionary");
    return *polynomial_;
}
const EigenfunctionSpectrum &KoopmanModel::spectrum() const noexcept {
    return polynomial_ ? polynomial_->spectrum() : *scatteredSpectrum_;
}
const EigenfunctionCoordinates &KoopmanModel::coordinates() const noexcept {
    return polynomial_ ? polynomial_->coordinates() : scatteredCoordinates_;
}
const EigenfunctionFitDiagnostics &KoopmanModel::fitDiagnostics() const noexcept {
    return fitDiagnostics_;
}
std::optional<std::size_t> KoopmanModel::physicalStateDimension() const noexcept {
    return coordinates().kind == EigenfunctionCoordinateKind::physicalState
               ? std::optional<std::size_t>(coordinateDimension())
               : coordinates().physicalStateDimension;
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
    discreteInputMaps_[std::move(name)] = std::move(map);
}
} // namespace dynsyslab::core
