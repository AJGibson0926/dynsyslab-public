// The parts of a Spectrum beyond its eigenvalues (Spectrum.hpp, SpectralDecomposition.hpp): the
// eigenfunctions of the atomic part, the spectral bands, the decomposition of observables over
// both (Koopman modes and measures, with atoms derived from the modes), and the functions on spectral
// sets (the pseudospectrum's minimizers, the measures' wave packets).
#include <dynsyslab/core/Spectrum.hpp>
#include <dynsyslab/core/FunctionExpansion.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dynsyslab::core {
namespace {
bool knownContinuity(const std::string& value) {
    return value.empty() || value == continuity::absolute || value == continuity::singular;
}
}  // namespace

// ---- eigenfunctions --------------------------------------------------------------------------------

std::vector<std::string> Spectrum::representationNames() const {
    std::vector<std::string> names;
    for (const auto& rep : representations_) names.push_back(rep.name);
    return names;
}

const EigenfunctionRepresentation& Spectrum::representation(const std::string& key) const {
    if (representations_.empty()) throw std::logic_error("these eigenfunctions have no representation");
    if (key.empty()) return representations_.front();
    for (const auto& rep : representations_)
        if (rep.name == key) return rep;
    const EigenfunctionRepresentation* match = nullptr;
    std::size_t matches = 0;
    for (const auto& rep : representations_)
        if (rep.function->kind() == key || (key == "window_map" && rep.name == "time-delay")) { match = &rep; ++matches; }
    std::string names;
    for (const auto& name : representationNames()) names += (names.empty() ? "" : ", ") + name;
    if (matches == 1) return *match;
    if (matches > 1) throw std::invalid_argument("these eigenfunctions have several " + key + " representations; name one of: " + names);
    throw std::invalid_argument("these eigenfunctions have no representation " + key + "; they have: " + names);
}

Spectrum& Spectrum::addRepresentation(std::shared_ptr<const FunctionRepresentation> function, bool first) {
    if (!function) throw std::invalid_argument("eigenfunction representation is empty");
    if (function->outputCount() != size())
        throw std::invalid_argument("a representation of these eigenfunctions needs one output per mode");
    EigenfunctionRepresentation rep{representationName(*function), std::move(function)};
    for (const auto& existing : representations_)
        if (existing.name == rep.name) throw std::invalid_argument("these eigenfunctions already have a " + rep.name + " representation");
    if (first) representations_.insert(representations_.begin(), std::move(rep));
    else representations_.push_back(std::move(rep));
    return *this;
}

std::size_t Spectrum::replaceRepresentation(std::shared_ptr<const FunctionRepresentation> function) {
    if (!function) throw std::invalid_argument("eigenfunction representation is empty");
    if (function->outputCount() != size())
        throw std::invalid_argument("a representation of these eigenfunctions needs one output per mode");
    const auto name = representationName(*function);
    for (std::size_t k = 0; k < representations_.size(); ++k)
        if (representations_[k].name == name) { representations_[k].function = std::move(function); return k; }
    throw std::invalid_argument("these eigenfunctions have no representation " + name);
}

void Spectrum::setFirstRepresentation(std::shared_ptr<const FunctionRepresentation> function) {
    if (representations_.empty() || !function || function->outputCount() != size())
        throw std::invalid_argument("the first representation can only be replaced by the same eigenfunctions");
    representations_.front().function = std::move(function);
}

bool Spectrum::moveRepresentationToFront(const std::string& name) {
    const auto at = std::find_if(representations_.begin(), representations_.end(),
                                 [&](const EigenfunctionRepresentation& r) { return r.name == name; });
    if (at == representations_.end()) throw std::invalid_argument("these eigenfunctions have no representation " + name);
    if (at == representations_.begin()) return false;
    auto chosen = std::move(*at);
    representations_.erase(at);
    representations_.insert(representations_.begin(), std::move(chosen));
    return true;
}

Spectrum& Spectrum::setSamples(EigenfunctionSamples samples) {
    if (!samples.modes) samples.modes = size();
    if (samples.modes != size() || samples.values.size() != samples.samples() * samples.modes
        || samples.sampleIndices.size() != samples.samples())
        throw std::invalid_argument("eigenfunction samples do not match the spectrum");
    samples_ = std::move(samples);
    return *this;
}

Eigenfunction Spectrum::eigenfunction(std::size_t mode, const ObservableCatalog* catalog) const {
    if (mode >= size()) throw std::out_of_range("eigenfunction index");
    if (representations_.empty()) throw std::logic_error("this spectrum has no eigenfunctions");
    const auto label = "psi_" + std::to_string(mode);
    auto psi = Eigenfunction::fromRepresentation({label, label}, representations_.front().function, mode, eigenvalue(mode), catalog);
    for (std::size_t k = 1; k < representations_.size(); ++k)
        psi = psi.withRepresentation(representations_[k].function->selectOutputs({mode}), catalog, false, representations_[k].name);
    return psi;
}

// ---- continuous part ---------------------------------------------------------------------------------

const SpectralBand* Spectrum::band(const std::string& name) const {
    for (const auto& c : bands_)
        if (c.name == name) return &c;
    return nullptr;
}

Spectrum& Spectrum::addBand(SpectralBand c) {
    if (c.name.empty()) throw std::invalid_argument("a spectral band needs a name");
    (void)bandPoint(c, 0.0);  // a known support
    if (const auto& g = c.eigenfunctions; !g.empty()) {
        if (g.levelSet() && (g.functions->outputCount() != 1 || g.level->outputCount() != 1
                             || g.functions->coordinateDimension() != g.level->coordinateDimension()))
            throw std::invalid_argument("level-set generalized eigenfunctions are one profile and one scalar level map on one domain");
        if (!g.levelSet() && g.functions->outputCount() != c.grid.size())
            throw std::invalid_argument("sampled generalized eigenfunctions have one output per grid point of their band");
    }
    if (!knownContinuity(c.continuity)) throw std::invalid_argument("unknown continuity " + c.continuity);
    for (const double point : c.grid)
        if (!std::isfinite(point)) throw std::invalid_argument("a spectral band's grid must be finite");
    if (const auto* existing = band(c.name)) {
        if (existing->continuity != c.continuity || existing->support != c.support || existing->grid != c.grid || existing->method != c.method
            || existing->bandwidth != c.bandwidth || existing->order != c.order)
            throw std::invalid_argument("this spectrum already has a different spectral band named " + c.name);
        return *this;
    }
    bands_.push_back(std::move(c));
    return *this;
}

Spectrum& Spectrum::setBands(std::vector<SpectralBand> bands) {
    for (const auto& m : measures_)
        if (std::none_of(bands.begin(), bands.end(), [&](const SpectralBand& c) { return c.name == m.band; }))
            throw std::invalid_argument("the measure of " + m.observable + " lies on band " + m.band + ", which would be removed");
    bands_.clear();
    for (auto& c : bands) addBand(std::move(c));
    return *this;
}

// ---- observables -------------------------------------------------------------------------------------

const ReconstructionResult* Spectrum::koopmanModes(const std::string& observable) const {
    const auto found = koopmanModes_.find(observable);
    return found == koopmanModes_.end() ? nullptr : &found->second;
}

Spectrum& Spectrum::setKoopmanModes(std::string observable, ReconstructionResult modes) {
    if (observable.empty() || !modes.outputs || modes.modes != size() || modes.matrix.size() != modes.outputs * size()
        || (!modes.offset.empty() && modes.offset.size() != modes.outputs))
        throw std::invalid_argument("observable reconstruction dimensions do not match KROM");
    koopmanModes_[std::move(observable)] = std::move(modes);
    return *this;
}

const SpectralMeasure* Spectrum::measure(const std::string& band, const std::string& observable, std::size_t output) const {
    for (const auto& m : measures_)
        if (m.band == band && m.observable == observable && m.output == output) return &m;
    return nullptr;
}

Spectrum& Spectrum::addMeasure(SpectralMeasure m) {
    const auto* c = band(m.band);
    if (!c) throw std::invalid_argument("a measure needs a spectral band of this spectrum; there is no " + m.band);
    if (m.observable.empty()) throw std::invalid_argument("a measure names its observable");
    if (c->grid.empty() || m.density.size() != c->grid.size()) throw std::invalid_argument("a measure's density lies on its band's grid");
    if (m.packets && m.packets->outputCount() != c->grid.size())
        throw std::invalid_argument("a measure's wave packets have one output per grid point of its band");
    if (measure(m.band, m.observable, m.output))
        throw std::invalid_argument("this spectrum already has a measure of " + m.observable + " output " + std::to_string(m.output)
                                    + " on " + m.band);
    measures_.push_back(std::move(m));
    return *this;
}

Spectrum& Spectrum::setMeasures(std::vector<SpectralMeasure> measures) {
    measures_.clear();
    for (auto& m : measures) addMeasure(std::move(m));
    return *this;
}

std::vector<const SpectralMeasure*> Spectrum::bandMeasures(const std::string& band) const {
    std::vector<const SpectralMeasure*> out;
    for (const auto& m : measures_)
        if (m.band == band) out.push_back(&m);
    return out;
}

ObservableDecomposition Spectrum::decomposition(const std::string& observable) const {
    ObservableDecomposition out{koopmanModes(observable), {}};
    for (const auto& m : measures_)
        if (m.observable == observable) out.measures.push_back(&m);
    return out;
}

std::vector<SpectralAtom> Spectrum::atoms(const std::string& observable, std::size_t output) const {
    const auto* modes = koopmanModes(observable);
    if (!modes) return {};
    if (output >= modes->outputs) throw std::out_of_range("observable output index");
    std::vector<SpectralAtom> out;
    for (std::size_t j = 0; j < modes->modes; ++j) out.push_back({modeId(j), std::norm(modes->matrix[output * modes->modes + j])});
    return out;
}

// ---- evidence, operator, functions on spectral sets --------------------------------------------------

Spectrum& Spectrum::setCompression(std::shared_ptr<const KoopmanCompression> compression) {
    compression_ = std::move(compression);
    return *this;
}

std::vector<std::string> Spectrum::spectralFunctionKeys() const {
    std::vector<std::string> keys;
    for (const auto& b : bands_) {
        if (b.eigenfunctions.functions) keys.push_back(bandFunctionKey(b.name));
        if (b.eigenfunctions.level) keys.push_back(bandLevelKey(b.name));
        if (b.eigenfunctions.family) keys.push_back(bandFamilyKey(b.name));
    }
    if (pseudospectrum_ && pseudospectrum_->minimizers) keys.emplace_back(pseudospectrumFunctionKey());
    for (const auto& m : measures_)
        if (m.packets) keys.push_back(measureFunctionKey(m));
    return keys;
}

const FunctionExpansion* Spectrum::spectralFunction(const std::string& key) const {
    for (const auto& b : bands_) {
        if (key == bandFunctionKey(b.name)) return b.eigenfunctions.functions.get();
        if (key == bandLevelKey(b.name)) return b.eigenfunctions.level.get();
        if (key == bandFamilyKey(b.name)) return b.eigenfunctions.family.get();
    }
    if (key == pseudospectrumFunctionKey()) return pseudospectrum_ ? pseudospectrum_->minimizers.get() : nullptr;
    for (const auto& m : measures_)
        if (measureFunctionKey(m) == key) return m.packets.get();
    return nullptr;
}

Spectrum& Spectrum::setSpectralFunction(const std::string& key, std::shared_ptr<const FunctionExpansion> function) {
    if (!function) throw std::invalid_argument("a spectral function is empty");
    if (key == pseudospectrumFunctionKey()) {
        if (!pseudospectrum_ || function->outputCount() != pseudospectrum_->grid.size())
            throw std::invalid_argument("the pseudospectrum's minimizers have one output per grid point");
        pseudospectrum_->minimizers = std::move(function);
        return *this;
    }
    for (auto& b : bands_) {  // the band's generalized eigenfunctions, attached as read (checked when complete)
        if (key == bandFunctionKey(b.name)) { b.eigenfunctions.functions = std::move(function); return *this; }
        if (key == bandFamilyKey(b.name)) {
            if (function->outputCount() != 1) throw std::invalid_argument("a band's explicit family is one function of the state and p");
            b.eigenfunctions.family = std::move(function);
            return *this;
        }
        if (key == bandLevelKey(b.name)) {
            if (function->outputCount() != 1) throw std::invalid_argument("a band's level map is scalar");
            b.eigenfunctions.level = std::move(function);
            return *this;
        }
    }
    for (auto& m : measures_)
        if (measureFunctionKey(m) == key) {
            if (function->outputCount() != band(m.band)->grid.size())
                throw std::invalid_argument("a measure's wave packets have one output per grid point of its band");
            m.packets = std::move(function);
            return *this;
        }
    throw std::invalid_argument("this spectrum has no spectral set " + key);
}

}  // namespace dynsyslab::core
