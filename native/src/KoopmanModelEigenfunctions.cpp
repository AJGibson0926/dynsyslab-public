// A K-ROM's eigenfunctions as data: the spectrum, the values at the training samples, the origin and
// the named representations, of which the first is evaluated (KoopmanModel.hpp). Also the parts a
// model inherits from the krom it was loaded from, and mode extraction.
#include <dynsyslab/core/KoopmanModel.hpp>
#include <dynsyslab/core/ScientificCompatibility.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace dynsyslab::core {
namespace {

// Catalog name for an interpolation domain: one observable, a vector observable split into "name[k]"
// channels, or a comma-joined list of scalar observables; empty when the domain is not catalog
// observables of the state.
std::string domainObservableName(const FunctionDomain& domain) {
    if (!domain.stateMapFromCatalog || domain.names.empty()) return {};
    for (const auto& name : domain.names)
        if (name.empty()) return {};
    if (domain.names.size() == 1) return domain.names.front();
    const auto bracket = domain.names.front().find('[');
    const auto base = domain.names.front().substr(0, bracket);
    bool vector = bracket != std::string::npos;
    for (std::size_t k = 0; vector && k < domain.names.size(); ++k) vector = domain.names[k] == base + "[" + std::to_string(k) + "]";
    if (vector) return base;
    std::string joined;
    for (const auto& name : domain.names) joined += (joined.empty() ? "" : ",") + name;
    return joined;
}

}  // namespace

std::string domainName(const FunctionDomain& domain) {
    if (domain.window) return "window";
    if (domain.isPhysicalState()) return "state";
    std::string joined;
    for (const auto& name : domain.names) {
        if (name.empty()) return "coordinates";
        joined += (joined.empty() ? "" : ",") + name;
    }
    return joined.empty() ? "coordinates" : joined;
}

std::vector<double> empiricalSpectralWeight(const EigenfunctionSamples& samples) {
    const std::size_t n = samples.samples(), m = samples.modes;
    if (!n) throw std::invalid_argument("spectral weights need eigenfunction samples");
    std::vector<double> weight(n * m), total(m, 0.0);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < m; ++j) total[j] += weight[k * m + j] = std::norm(samples.values[k * m + j]);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < m; ++j) weight[k * m + j] = total[j] > 0 ? weight[k * m + j] / total[j] : 1.0 / double(n);
    return weight;
}

KoopmanModel::KoopmanModel(Spectrum spectrum, EigenfunctionSamples samples, DataOrigin origin,
                           std::vector<std::shared_ptr<const FunctionRepresentation>> representations,
                           std::size_t inputDimension)
    : spectrum_(std::move(spectrum)), suppliedInputDimension_(inputDimension),
      samples_(std::move(samples)), origin_(std::move(origin)) {
    if (!samples_.modes) samples_.modes = spectrum_->size();
    if (samples_.modes != spectrum_->size() || samples_.values.size() != samples_.samples() * samples_.modes
        || samples_.sampleIndices.size() != samples_.samples())
        throw std::invalid_argument("eigenfunction samples do not match the spectrum");
    if (representations.empty()) throw std::invalid_argument("a K-ROM needs at least one eigenfunction representation");
    for (auto& function : representations) addRepresentation(std::move(function), false);
}

void KoopmanModel::setSamples(EigenfunctionSamples samples) {
    if (samples.modes != liftedDimension() || samples.values.size() != samples.samples() * samples.modes
        || samples.sampleIndices.size() != samples.samples())
        throw std::invalid_argument("eigenfunction samples do not match the spectrum");
    for (const auto& rep : representations_)
        if (dynamic_cast<const FunctionInterpolation*>(rep.function.get()))
            throw std::invalid_argument("an interpolated K-ROM's samples are its interpolation's");
    samples_ = std::move(samples);
}

std::vector<std::string> KoopmanModel::representationNames() const {
    std::vector<std::string> names;
    for (const auto& rep : representations_) names.push_back(rep.name);
    return names;
}

const KoopmanModel::Representation& KoopmanModel::representation(const std::string& key) const {
    if (representations_.empty()) throw std::logic_error("this K-ROM has no eigenfunction representation");
    if (key.empty()) return representations_.front();
    for (const auto& rep : representations_)
        if (rep.name == key) return rep;
    const Representation* match = nullptr;
    std::size_t matches = 0;
    for (const auto& rep : representations_)
        if (rep.function->kind() == key || (key == "window_map" && rep.name == "time-delay")) { match = &rep; ++matches; }
    std::string names;
    for (const auto& name : representationNames()) names += (names.empty() ? "" : ", ") + name;
    if (matches == 1) return *match;
    if (matches > 1) throw std::invalid_argument("these eigenfunctions have several " + key + " representations; name one of: " + names);
    throw std::invalid_argument("these eigenfunctions have no representation " + key + "; they have: " + names);
}

void KoopmanModel::addRepresentation(std::shared_ptr<const FunctionRepresentation> function, bool first) {
    if (!function) throw std::invalid_argument("eigenfunction representation is empty");
    if (function->outputCount() != liftedDimension())
        throw std::invalid_argument("a representation of these eigenfunctions needs one output per mode");
    Representation rep{representationName(*function), std::move(function)};
    for (const auto& existing : representations_)
        if (existing.name == rep.name) throw std::invalid_argument("these eigenfunctions already have a " + rep.name + " representation");
    const bool activates = first || representations_.empty();
    if (activates) representations_.insert(representations_.begin(), std::move(rep));
    else representations_.push_back(std::move(rep));
    if (activates) activate();
}

void KoopmanModel::replaceRepresentation(std::shared_ptr<const FunctionRepresentation> function) {
    if (!function) throw std::invalid_argument("eigenfunction representation is empty");
    const auto name = representationName(*function);
    for (std::size_t k = 0; k < representations_.size(); ++k)
        if (representations_[k].name == name) {
            representations_[k].function = std::move(function);
            inherited_.representations.erase(name);
            if (k == 0) activate();
            return;
        }
    throw std::invalid_argument("these eigenfunctions have no representation " + name);
}

void KoopmanModel::useRepresentation(const std::string& key) {
    const auto name = representation(key).name;
    const auto at = std::find_if(representations_.begin(), representations_.end(), [&](const Representation& r) { return r.name == name; });
    if (at == representations_.begin()) return;
    auto chosen = std::move(*at);
    representations_.erase(at);
    representations_.insert(representations_.begin(), std::move(chosen));
    activate();
}

void KoopmanModel::activate() {
    const auto& front = representations_.front();
    const auto& domain = front.function->coordinates();
    if (const auto* fx = dynamic_cast<const FunctionExpansion*>(front.function.get())) {
        scattered_.reset();
        polynomial_ = *fx;
        if (system_ && !domain.window) {
            const auto configuration = systemConfiguration_;
            attachPhysicalSystem(configuration);
        }
    } else if (const auto* fi = dynamic_cast<const FunctionInterpolation*>(front.function.get())) {
        polynomial_.reset();
        auto data = fi->data();
        if (!domain.isPhysicalState() && data.coordinateObservable.empty()) data.coordinateObservable = domainObservableName(domain);
        if (!(data.finiteDifferenceRelativeStep > 0)) data.finiteDifferenceRelativeStep = 1.0e-6;
        std::shared_ptr<const FunctionInterpolation> built = std::make_shared<const FunctionInterpolation>(std::move(data), domain);
        if (system_) {
            const std::size_t physical = system_->descriptor().stateDimension;
            if (domain.isPhysicalState() && built->data().stateDimension != physical)
                throw std::invalid_argument("interpolated eigenfunctions and dynamical system differ in dimension");
            if (domain.needsCatalog() || (!domain.isPhysicalState() && domain.stateMapFromCatalog)) {
                auto q = domain.bound(system_->observableCatalog());
                if (q.stateMap->stateDimension() != physical)
                    throw std::invalid_argument("interpolated K-ROM coordinates and their observable differ in dimension");
                built = std::dynamic_pointer_cast<const FunctionInterpolation>(built->onDomain(std::move(q)));
            }
        }
        scattered_ = std::move(built);
    } else {
        throw std::invalid_argument("the " + front.name + " representation cannot be evaluated by a K-ROM");
    }
    representations_.front().function = scattered_ ? std::shared_ptr<const FunctionRepresentation>(scattered_)
                                                   : std::make_shared<const FunctionExpansion>(*polynomial_);
    restampMaps();
}

void KoopmanModel::restampMaps() {
    if (controlMaps_.empty() && discreteInputMaps_.empty()) return;
    const auto identity = koopmanModelIdentity(*this);
    for (auto& entry : controlMaps_) std::visit([&](auto& map) { map.modelIdentity = identity; }, entry.second);
    for (auto& entry : discreteInputMaps_) entry.second.modelIdentity = identity;
}

void KoopmanModel::bindSystem(SystemConfiguration system) {
    system_ = std::shared_ptr<ODESystem>(createSystem(system).release());
    systemConfiguration_ = system_->configuration();
    activate();
}

const FunctionInterpolation& KoopmanModel::interpolation() const {
    for (const auto& rep : representations_)
        if (const auto* fi = dynamic_cast<const FunctionInterpolation*>(rep.function.get())) return *fi;
    throw std::invalid_argument("these eigenfunctions have no interpolation representation");
}

void KoopmanModel::setInterpolationMethod(InterpolationMethod method, std::size_t neighbors) {
    replaceRepresentation(interpolation().withMethod(std::move(method), neighbors));
}

std::string KoopmanModel::recognizeDomain(const std::vector<std::string>& observables, const std::string& key) const {
    const auto& rep = representation(key);
    const auto& domain = rep.function->coordinates();
    const auto name = domainName(domain);
    if (observables.empty()) return name;
    const auto& own = domain.window ? domain.window->channelObservables : domain.names;
    if (own == observables) return name;
    throw std::invalid_argument("the " + rep.name + " representation is on " + name + ", which these observables do not describe");
}

std::vector<std::complex<double>> KoopmanModel::evaluateRepresentation(const double* point, const std::string& key) const {
    return representation(key).function->evaluate(point);
}

EigenfunctionErrors KoopmanModel::errors(const std::string& key, const TimeSeriesRecord* input) const {
    const auto& rep = representation(key);
    const auto& function = *rep.function;
    const std::size_t n = samples_.samples(), m = samples_.modes;
    std::vector<std::complex<double>> predicted;
    predicted.reserve(n * m);
    if (const auto& window = function.coordinates().window) {  // the window at each sample's label
        if (!input) throw std::invalid_argument("window-map errors need the input time series the eigenfunctions were fitted on");
        for (const auto label : samples_.sampleIndices) {
            const auto h = window->window_at(input->series, label);
            const auto psi = function.evaluate(h.data());
            predicted.insert(predicted.end(), psi.begin(), psi.end());
        }
    } else if (const auto* fi = dynamic_cast<const FunctionInterpolation*>(&function)) {
        const auto& points = fi->points();
        const std::size_t width = function.coordinateDimension();
        if (points.size() != n * width) throw std::invalid_argument("interpolation points do not match the samples");
        for (std::size_t k = 0; k < n; ++k) {
            const auto psi = function.evaluate(points.data() + k * width);
            predicted.insert(predicted.end(), psi.begin(), psi.end());
        }
    } else {
        throw std::invalid_argument("the " + rep.name + " representation stores no sample points; its errors need the input states");
    }
    EigenfunctionErrors out{rep.name, domainName(function.coordinates()), n, std::vector<double>(m), std::vector<double>(m)};
    for (std::size_t mode = 0; mode < m; ++mode) {
        double error = 0, norm = 0, residual = 0, base = 0;
        for (std::size_t k = 0; k < n; ++k) {
            const auto value = predicted[k * m + mode], truth = samples_.values[k * m + mode];
            error += std::norm(value - truth);
            norm += std::norm(truth);
            if (k + 1 == n) continue;
            const double elapsed = samples_.time[k + 1] - samples_.time[k];
            if (!(elapsed > 0)) continue;  // time does not advance: the next sample starts another trajectory
            residual += std::norm(predicted[(k + 1) * m + mode] - spectrum().transition(mode, elapsed) * value);
            base += std::norm(value);
        }
        out.valueRelative[mode] = std::sqrt(error / std::max(norm, 1e-300));
        out.koopmanResidual[mode] = std::sqrt(residual / std::max(base, 1e-300));
    }
    return out;
}

Eigenfunction KoopmanModel::eigenfunction(std::size_t mode, const ObservableCatalog* catalog) const {
    if (mode >= eigenfunctionCount()) throw std::out_of_range("K-ROM eigenfunction index");
    const auto own = system_ && !catalog ? std::optional<ObservableCatalog>(system_->observableCatalog()) : std::nullopt;
    if (own) catalog = &*own;
    const auto label = "psi_" + std::to_string(mode);
    auto psi = Eigenfunction::fromRepresentation({label, label}, representations_.front().function, mode,
                                                 spectrum().eigenvalue(mode), catalog);
    for (std::size_t k = 1; k < representations_.size(); ++k)
        psi = psi.withRepresentation(representations_[k].function->selectOutputs({mode}), catalog, false, representations_[k].name);
    return psi;
}

void KoopmanModel::setParent(KromParent parent) {
    parent_ = std::move(parent);
    inherited_ = {};
    for (const auto& rep : representations_) inherited_.representations.insert(rep.name);
    for (const auto& item : reconstructions_) inherited_.readouts.insert(item.first);
    for (const auto& item : controlMaps_) inherited_.controlMaps.insert(item.first);
    for (const auto& item : discreteInputMaps_) inherited_.discreteMaps.insert(item.first);
    ownSpectralEvidence_ = false;
}

void KoopmanModel::clearParent() {
    parent_.reset();
    inherited_ = {};
    ownSpectralEvidence_ = true;
}

Eigenpair KoopmanModel::eigenpair(std::size_t mode, const ObservableCatalog* catalog) const {
    return {spectrum().eigenvalue(mode), eigenfunction(mode, catalog), spectrum().residual(mode)};
}

void KoopmanModel::setResidual(std::size_t mode, std::optional<EigenpairResidual> residual) {
    spectrum_->setResidual(mode, std::move(residual));
    ownSpectralEvidence_ = true;
}

void KoopmanModel::setSpectralEvidence(const Spectrum& evidence) {
    spectrum_->adoptEvidence(evidence);
    ownSpectralEvidence_ = true;
}

KoopmanModel KoopmanModel::extract(const std::vector<std::size_t>& indices) const {
    if (indices.empty()) throw std::invalid_argument("KROM extraction requires a mode");
    std::set<std::size_t> unique(indices.begin(), indices.end());
    if (unique.size() != indices.size() || *unique.rbegin() >= liftedDimension())
        throw std::out_of_range("KROM extraction indices are invalid");
    EigenfunctionSamples subset{samples_.time, samples_.sampleIndices, {}, indices.size()};
    subset.values.reserve(samples_.samples() * indices.size());
    for (std::size_t k = 0; k < samples_.samples(); ++k)
        for (const auto mode : indices) subset.values.push_back(samples_.values[k * samples_.modes + mode]);
    std::vector<std::shared_ptr<const FunctionRepresentation>> functions;
    for (const auto& rep : representations_) functions.push_back(rep.function->selectOutputs(indices));
    KoopmanModel result(spectrum().extract(indices), std::move(subset), origin_, std::move(functions), suppliedInputDimension_);
    if (system_) result.bindSystem(systemConfiguration_);
    if (!fitDiagnostics_.residualMaximum.empty())
        result.fitDiagnostics_ = fitDiagnostics_.extract(indices);
    // The selected eigenpairs keep their residuals, and the spectrum its continuous part and
    // pseudospectrum (Spectrum::extract). The fit's diagnostics refer to the original trial subspace;
    // extraction intentionally clears them.
    result.datasetReferences_ = datasetReferences_;
    result.bindingEvidenceJson_ = bindingEvidenceJson_;
    const std::size_t inputs = inputDimension(), modes = indices.size();
    for (const auto& entry : discreteInputMaps_) {
        auto map = entry.second; map.matrix.clear(); map.eigenfunctionCount = modes; map.modelIdentity.clear();
        for (auto j : indices)
            map.matrix.insert(map.matrix.end(), entry.second.matrix.begin() + j * inputs, entry.second.matrix.begin() + (j + 1) * inputs);
        result.setDiscreteInputMap(entry.first, std::move(map));
    }
    for (const auto& item : controlMaps_) {
        if (const auto* constantSource = std::get_if<ConstantControlMap>(&item.second)) {
            ConstantControlMap map{{}, modes, inputs, constantSource->provenance + "; independent row subset"};
            map.sampleTime = constantSource->sampleTime;
            for (const auto row : indices)
                map.matrix.insert(map.matrix.end(), constantSource->matrix.begin() + row * inputs,
                                  constantSource->matrix.begin() + (row + 1) * inputs);
            if (!constantSource->discreteMatrix.empty())
                for (const auto row : indices)
                    map.discreteMatrix.insert(map.discreteMatrix.end(), constantSource->discreteMatrix.begin() + row * inputs,
                                              constantSource->discreteMatrix.begin() + (row + 1) * inputs);
            result.setControlMap(item.first, std::move(map));
        } else {
            const auto& bilinearSource = std::get<BilinearControlMap>(item.second);
            BilinearControlMap map{{}, {}, {}, modes, inputs, bilinearSource.provenance + "; independent principal restriction"};
            map.operators.resize(inputs * modes * modes); map.offsets.resize(modes * inputs);
            if (!bilinearSource.conjugateOperators.empty()) map.conjugateOperators.resize(map.operators.size());
            for (std::size_t input = 0; input < inputs; ++input)
                for (std::size_t row = 0; row < modes; ++row) {
                    map.offsets[row * inputs + input] = bilinearSource.offsets[indices[row] * inputs + input];
                    for (std::size_t column = 0; column < modes; ++column) {
                        const auto target = (input * modes + row) * modes + column;
                        const auto origin = (input * bilinearSource.modes + indices[row]) * bilinearSource.modes + indices[column];
                        map.operators[target] = bilinearSource.operators[origin];
                        if (!map.conjugateOperators.empty()) map.conjugateOperators[target] = bilinearSource.conjugateOperators[origin];
                    }
                }
            result.setControlMap(item.first, std::move(map));
        }
    }
    for (const auto& item : reconstructions_) {
        ReconstructionResult projection = item.second;
        projection.modes = modes; projection.matrix.assign(projection.outputs * modes, {});
        for (std::size_t output = 0; output < projection.outputs; ++output)
            for (std::size_t column = 0; column < modes; ++column)
                projection.matrix[output * modes + column] = item.second.matrix[output * liftedDimension() + indices[column]];
        result.setReconstruction(item.first, std::move(projection));
    }
    // The same krom, the chosen modes of the modes this model kept; what was inherited stays inherited.
    if (parent_) {
        auto parent = *parent_;
        parent.modes.clear();
        for (const auto j : indices) parent.modes.push_back(parent_->modes.empty() ? j : parent_->modes.at(j));
        result.parent_ = std::move(parent);
        result.inherited_ = inherited_;
    }
    result.ownSpectralEvidence_ = ownSpectralEvidence_;  // residuals it set itself stay its own
    return result;
}

}  // namespace dynsyslab::core
