#include <dynsyslab/core/EigenfunctionObjects.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {

EigenfunctionSet::EigenfunctionSet(Spectrum spectrum, EigenfunctionSamples samples, DataOrigin origin,
                                   std::vector<std::shared_ptr<const EigenfunctionRepresentation>> representations)
    : spectrum_(std::move(spectrum)), samples_(std::move(samples)), origin_(std::move(origin)),
      representations_(std::move(representations)) {
    if (samples_.modes != spectrum_.size() || samples_.values.size() != samples_.samples() * samples_.modes
        || samples_.sampleIndices.size() != samples_.samples())
        throw std::invalid_argument("eigenfunction samples do not match the spectrum");
    for (const auto& rep : representations_)
        if (!rep) throw std::invalid_argument("eigenfunction representation is empty");
}

const EigenfunctionRepresentation& EigenfunctionSet::representation(const std::string& kind) const {
    if (representations_.empty()) throw std::logic_error("eigenfunction set has no representation");
    if (kind.empty()) return *representations_.front();
    for (const auto& rep : representations_)
        if (rep->kind() == kind) return *rep;
    throw std::invalid_argument("eigenfunction set has no representation " + kind);
}

std::string EigenfunctionSet::recognizeDomain(const std::vector<std::string>& observables,
                                              const std::string& kind) const {
    const auto& rep = representation(kind);
    const auto domains = rep.domains();
    if (observables.empty()) return domains.front();
    for (const auto& domain : domains)
        if (rep.domainObservables(domain) == observables) return domain;
    std::string known;
    for (const auto& domain : domains) known += (known.empty() ? "" : ", ") + domain;
    throw std::invalid_argument("no " + rep.kind() + " domain matches these observables; domains: " + known);
}

std::vector<std::complex<double>> EigenfunctionSet::evaluate(const double* point, const std::string& domain,
                                                             const std::string& kind) const {
    const auto& rep = representation(kind);
    return rep.evaluate(point, domain.empty() ? rep.domains().front() : domain);
}

EigenfunctionErrors EigenfunctionSet::errors(const std::string& kind, const std::string& domainName,
                                             const TimeSeriesRecord* input) const {
    const auto& rep = representation(kind);
    const auto domain = domainName.empty() ? rep.domains().front() : domainName;
    const auto predicted = rep.evaluateTraining(samples_, domain, input);
    const std::size_t n = samples_.samples(), m = samples_.modes;
    if (predicted.size() != n * m) throw std::logic_error("representation returned the wrong number of values");
    const auto& spectrum = spectrum_;
    EigenfunctionErrors out{rep.kind(), domain, n, std::vector<double>(m), std::vector<double>(m)};
    for (std::size_t mode = 0; mode < m; ++mode) {
        double error = 0, norm = 0, residual = 0, base = 0;
        for (std::size_t k = 0; k < n; ++k) {
            const auto value = predicted[k * m + mode], truth = samples_.values[k * m + mode];
            error += std::norm(value - truth);
            norm += std::norm(truth);
            if (k + 1 == n) continue;
            const double elapsed = samples_.time[k + 1] - samples_.time[k];
            residual += std::norm(predicted[(k + 1) * m + mode] - spectrum.transition(mode, elapsed) * value);
            base += std::norm(value);
        }
        out.valueRelative[mode] = std::sqrt(error / std::max(norm, 1e-300));
        out.koopmanResidual[mode] = std::sqrt(residual / std::max(base, 1e-300));
    }
    return out;
}

void EigenfunctionSet::replaceRepresentation(std::shared_ptr<const EigenfunctionRepresentation> representation) {
    if (!representation) throw std::invalid_argument("eigenfunction representation is empty");
    for (auto& rep : representations_)
        if (rep->kind() == representation->kind()) { rep = std::move(representation); return; }
    throw std::invalid_argument("eigenfunction set has no representation " + representation->kind());
}

EigenfunctionSet EigenfunctionSet::extract(const std::vector<std::size_t>& modes) const {
    EigenfunctionSamples subset{samples_.time, samples_.sampleIndices, {}, modes.size()};
    subset.values.reserve(samples_.samples() * modes.size());
    for (std::size_t k = 0; k < samples_.samples(); ++k)
        for (const auto mode : modes) subset.values.push_back(samples_.values[k * samples_.modes + mode]);
    std::vector<std::shared_ptr<const EigenfunctionRepresentation>> reps;
    for (const auto& rep : representations_) reps.push_back(rep->extract(modes));
    return EigenfunctionSet(spectrum_.extract(modes), std::move(subset), origin_, std::move(reps));
}

Eigenfunction EigenfunctionSet::eigenfunction(std::size_t mode, const std::string& domain, const std::string& kind,
                                              const ObservableCatalog* catalog) const {
    if (mode >= size()) throw std::out_of_range("eigenfunction index");
    const auto& rep = representation(kind);
    const auto name = domain.empty() ? rep.domains().front() : domain;
    const auto label = "psi_" + std::to_string(mode);
    return Eigenfunction::fromRepresentation({label, label}, rep.function(name), mode, spectrum_[mode], catalog);
}

SharedRepresentation::SharedRepresentation(std::shared_ptr<const FunctionRepresentation> function)
    : function_(std::move(function)) {
    if (!function_) throw std::invalid_argument("shared representation needs a function");
    domain_ = domainName(function_->coordinates());
}

std::string domainName(const FunctionDomain& domain) {
    if (domain.isPhysicalState()) return "state";
    std::string joined;
    for (const auto& name : domain.names) {
        if (name.empty()) return "coordinates";
        joined += (joined.empty() ? "" : ",") + name;
    }
    return joined.empty() ? "coordinates" : joined;
}

void SharedRepresentation::requireDomain(const std::string& domain) const {
    if (!domain.empty() && domain != domain_)
        throw std::invalid_argument("this " + kind() + " is on domain " + domain_ + ", not " + domain);
}

std::size_t SharedRepresentation::domainDimension(const std::string& domain) const {
    requireDomain(domain);
    return function_->coordinateDimension();
}

std::vector<std::string> SharedRepresentation::domainObservables(const std::string& domain) const {
    requireDomain(domain);
    return function_->coordinates().names;
}

std::vector<std::complex<double>> SharedRepresentation::evaluate(const double* point, const std::string& domain) const {
    requireDomain(domain);
    return function_->evaluate(point);
}

std::vector<std::complex<double>> SharedRepresentation::evaluateTraining(
    const EigenfunctionSamples& samples, const std::string& domain, const TimeSeriesRecord*) const {
    requireDomain(domain);
    const auto* interpolation = dynamic_cast<const FunctionInterpolation*>(function_.get());
    if (!interpolation) throw std::invalid_argument("this " + kind() + " stores no training points");
    const auto& points = interpolation->points();
    const std::size_t width = function_->coordinateDimension();
    if (points.size() != samples.samples() * width) throw std::invalid_argument("interpolation points do not match the samples");
    std::vector<std::complex<double>> out;
    out.reserve(samples.samples() * function_->outputCount());
    for (std::size_t k = 0; k < samples.samples(); ++k) {
        const auto value = function_->evaluate(points.data() + k * width);
        out.insert(out.end(), value.begin(), value.end());
    }
    return out;
}

std::shared_ptr<const EigenfunctionRepresentation> SharedRepresentation::extract(const std::vector<std::size_t>& modes) const {
    return std::make_shared<SharedRepresentation>(function_->selectOutputs(modes));
}

json::Value SharedRepresentation::describe() const {
    json::Value::Object out{{"kind", kind()}, {"domain_name", domain_}, {"domain", functionDomainJson(function_->coordinates())}};
    if (const auto* interpolation = dynamic_cast<const FunctionInterpolation*>(function_.get())) {
        out["method"] = interpolationMethodJson(interpolation->data().method, interpolation->data().neighbors);
        out["points"] = "interpolation_points";
    }
    return out;
}

std::shared_ptr<const FunctionRepresentation> SharedRepresentation::function(const std::string& domain) const {
    requireDomain(domain);
    return function_;
}

const FunctionInterpolation& interpolationOf(const EigenfunctionSet& set) {
    for (const auto& rep : set.representations())
        if (const auto* shared = dynamic_cast<const SharedRepresentation*>(rep.get()))
            if (const auto* interpolation = dynamic_cast<const FunctionInterpolation*>(shared->function().get()))
                return *interpolation;
    throw std::invalid_argument("this eigenfunction set has no interpolation representation");
}

void setInterpolationMethod(EigenfunctionSet& set, InterpolationMethod method, std::size_t neighbors) {
    set.replaceRepresentation(std::make_shared<SharedRepresentation>(
        interpolationOf(set).withMethod(std::move(method), neighbors)));
}

}  // namespace dynsyslab::core
